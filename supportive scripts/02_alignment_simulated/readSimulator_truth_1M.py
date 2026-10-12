#!/usr/bin/env python3

import argparse
import csv
import gzip
import random
import re
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path


# ============================================================
# FASTA
# ============================================================

def read_fasta(path):
    records = []
    name = None
    seq_parts = []

    with open(path) as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    records.append(
                        (name, "".join(seq_parts).upper())
                    )

                # Use only the first whitespace-delimited field.
                name = line[1:].split()[0]
                seq_parts = []

            else:
                seq_parts.append(line)

    if name is not None:
        records.append(
            (name, "".join(seq_parts).upper())
        )

    if not records:
        raise RuntimeError(
            f"No FASTA records found in {path}"
        )

    return records


# ============================================================
# FASTQ
# ============================================================

def read_fastq(path):
    records = {}

    with open(path) as f:
        while True:
            header = f.readline()

            if not header:
                break

            seq = f.readline()
            plus = f.readline()
            qual = f.readline()

            if not seq or not plus or not qual:
                raise RuntimeError(
                    f"Truncated FASTQ file: {path}"
                )

            header = header.rstrip()
            seq = seq.rstrip()
            plus = plus.rstrip()
            qual = qual.rstrip()

            if not header.startswith("@"):
                raise RuntimeError(
                    f"Invalid FASTQ header in {path}: "
                    f"{header}"
                )

            name = header[1:].split()[0]

            # ART names:
            #
            # NC_000964.3-20/1
            # NC_000964.3-20/2
            #
            # Strip /1 or /2 to obtain the fragment QNAME used
            # in ART's SAM file.
            if name.endswith("/1") or name.endswith("/2"):
                pair_name = name[:-2]
            else:
                pair_name = name

            if pair_name in records:
                raise RuntimeError(
                    f"Duplicate FASTQ pair identifier "
                    f"{pair_name} in {path}"
                )

            records[pair_name] = (
                seq,
                qual,
            )

    if not records:
        raise RuntimeError(
            f"No FASTQ records found in {path}"
        )

    return records


def write_fastq_record(handle, pair_id, mate, record):
    seq, qual = record

    handle.write(
        f"@{pair_id}/{mate}\n"
    )
    handle.write(
        seq + "\n"
    )
    handle.write(
        "+\n"
    )
    handle.write(
        qual + "\n"
    )


# ============================================================
# CIGAR
# ============================================================

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def cigar_reference_length(cigar):
    """
    Return the number of reference bases consumed by a SAM CIGAR.

    Reference-consuming operations:
        M D N = X
    """

    if cigar == "*":
        return 0

    total = 0
    consumed_chars = 0

    for length, op in CIGAR_RE.findall(cigar):
        n = int(length)

        consumed_chars += len(length) + 1

        if op in {"M", "D", "N", "=", "X"}:
            total += n

    if consumed_chars != len(cigar):
        raise RuntimeError(
            f"Could not completely parse CIGAR: {cigar}"
        )

    return total


# ============================================================
# ART SAM
# ============================================================

def parse_art_sam(path):
    """
    Parse ART's SAM output directly.

    ART 2.5.8 can produce SAM files whose header/reference
    information is not fully accepted by samtools. We therefore
    parse the SAM text ourselves.

    Returns:

        {
            qname: {
                1: {...R1...},
                2: {...R2...}
            }
        }

    Mate identity is determined from SAM flags:
        0x40 = read 1
        0x80 = read 2

    Strand:
        0x10 = reverse
    """

    pairs = {}

    with open(path) as f:
        for line in f:

            if not line.strip():
                continue

            if line.startswith("@"):
                continue

            fields = line.rstrip("\n\r").split("\t")

            if len(fields) < 11:
                raise RuntimeError(
                    "Malformed ART SAM record:\n"
                    + line.rstrip()
                )

            qname = fields[0]

            try:
                flag = int(fields[1])
                pos1 = int(fields[3])
                mapq = int(fields[4])
                pnext1 = int(fields[7])
                tlen = int(fields[8])
            except ValueError as exc:
                raise RuntimeError(
                    "Invalid numeric field in ART SAM:\n"
                    + line.rstrip()
                ) from exc

            cigar = fields[5]

            if flag & 64:
                mate = 1
            elif flag & 128:
                mate = 2
            else:
                raise RuntimeError(
                    f"Could not determine mate from ART SAM "
                    f"flag {flag} for {qname}"
                )

            strand = "-" if (flag & 16) else "+"

            if pos1 <= 0:
                raise RuntimeError(
                    f"ART SAM record has invalid POS for "
                    f"{qname}: {pos1}"
                )

            ref_consumed = cigar_reference_length(cigar)

            if ref_consumed <= 0:
                raise RuntimeError(
                    f"ART SAM record has no reference-consuming "
                    f"alignment for {qname}: {cigar}"
                )

            record = {
                "qname": qname,
                "mate": mate,
                "flag": flag,
                "pos1": pos1,
                "start0": pos1 - 1,
                "strand": strand,
                "cigar": cigar,
                "reference_span": ref_consumed,
                "mapq": mapq,
                "pnext1": pnext1,
                "tlen": tlen,
            }

            if qname not in pairs:
                pairs[qname] = {}

            if mate in pairs[qname]:
                raise RuntimeError(
                    f"Duplicate ART SAM record for "
                    f"{qname} mate {mate}"
                )

            pairs[qname][mate] = record

    if not pairs:
        raise RuntimeError(
            f"No alignment records found in ART SAM: {path}"
        )

    for qname, mates in pairs.items():
        if 1 not in mates or 2 not in mates:
            raise RuntimeError(
                f"Incomplete ART SAM pair for {qname}: "
                f"mates={sorted(mates)}"
            )

    return pairs


# ============================================================
# COORDINATE TRANSFORMATION
# ============================================================

def rotated_to_original(
    rotation_start0,
    rotated_start0,
    reference_length,
):
    """
    Convert a 0-based coordinate on the rotated chromosome
    back to the original circular chromosome.

    rotated sequence:

        seq[rotation_start0:] + seq[:rotation_start0]
    """

    return (
        rotation_start0 + rotated_start0
    ) % reference_length


def alignment_crosses_origin(
    original_start0,
    reference_span,
    reference_length,
):
    """
    True when an alignment beginning at original_start0 and
    consuming reference_span bases crosses the original
    L -> 0 chromosome boundary.
    """

    return (
        original_start0 + reference_span
        > reference_length
    )


# ============================================================
# COMMAND EXECUTION
# ============================================================

def run_command(cmd):
    print(
        " ".join(
            shlex.quote(str(x))
            for x in cmd
        )
    )

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if result.returncode != 0:
        raise RuntimeError(
            "Command failed:\n"
            + " ".join(
                shlex.quote(str(x))
                for x in cmd
            )
            + "\n\nSTDOUT:\n"
            + result.stdout
            + "\nSTDERR:\n"
            + result.stderr
        )


# ============================================================
# ART
# ============================================================

def simulate_art(
    art,
    reference_name,
    sequence,
    pair_count,
    read_length,
    opts,
    seed,
    workdir,
    prefix_name,
):
    """
    Run classic ART Illumina and obtain:

        R1 FASTQ
        R2 FASTQ
        SAM truth

    The .aln files produced by ART are intentionally ignored.
    """

    ref_path = workdir / f"{prefix_name}.fasta"
    prefix = workdir / f"{prefix_name}_"

    with open(ref_path, "w") as handle:
        handle.write(
            f">{reference_name}\n"
        )
        handle.write(
            sequence + "\n"
        )

    # Use ART 2.5.8 native short options.
    cmd = [
        art,
        "-ss", "HS25",
        "-sam",
        "-i", str(ref_path),
        "-c", str(pair_count),
        "-l", str(read_length),
        "-p",
    ]

    cmd.extend(opts)

    cmd.extend([
        "-rs", str(seed),
        "-o", str(prefix),
    ])

    run_command(cmd)

    fq1 = Path(
        str(prefix) + "1.fq"
    )

    fq2 = Path(
        str(prefix) + "2.fq"
    )

    sam = Path(
        str(prefix) + ".sam"
    )

    for path in (fq1, fq2, sam):
        if not path.is_file():
            raise RuntimeError(
                f"Expected ART output not found: {path}"
            )

    reads1 = read_fastq(fq1)
    reads2 = read_fastq(fq2)
    truth = parse_art_sam(sam)

    names1 = set(reads1)
    names2 = set(reads2)
    names_sam = set(truth)

    if names1 != names2:
        raise RuntimeError(
            "ART R1 and R2 FASTQ identifiers disagree.\n"
            f"R1 pairs: {len(names1)}\n"
            f"R2 pairs: {len(names2)}"
        )

    if names1 != names_sam:
        only_fastq = sorted(
            names1 - names_sam
        )[:5]

        only_sam = sorted(
            names_sam - names1
        )[:5]

        raise RuntimeError(
            "ART FASTQ and SAM identifiers disagree.\n"
            f"FASTQ pairs: {len(names1)}\n"
            f"SAM pairs:   {len(names_sam)}\n"
            f"Only FASTQ examples: {only_fastq}\n"
            f"Only SAM examples:   {only_sam}"
        )

    return reads1, reads2, truth


# ============================================================
# TRUTH OUTPUT
# ============================================================

TRUTH_FIELDS = [
    "pair_id",

    "ref1",
    "start1_1",
    "strand1",
    "cigar1",
    "reference_span1",

    "ref2",
    "start1_2",
    "strand2",
    "cigar2",
    "reference_span2",

    "pair_type",
    "iteration",

    "rotation_start0",
    "reference_length",

    "art_start1_0",
    "art_start2_0",

    "art_pos1_1",
    "art_pos1_2",

    "art_tlen1",
    "art_tlen2",

    "read1_crosses_origin",
    "read2_crosses_origin",
    "pair_crosses_origin",
]


ROTATION_FIELDS = [
    "contig",
    "pair_type",
    "iteration",
    "rotation_start0",
    "reference_length",
    "pairs_simulated",
    "pairs_retained",
]


def build_truth_record(
    pair_id,
    ref_name,
    mate1,
    mate2,
    rotation_start0,
    ref_len,
    iteration,
    pair_type,
):
    original1_0 = rotated_to_original(
        rotation_start0,
        mate1["start0"],
        ref_len,
    )

    original2_0 = rotated_to_original(
        rotation_start0,
        mate2["start0"],
        ref_len,
    )

    read1_junction = alignment_crosses_origin(
        original1_0,
        mate1["reference_span"],
        ref_len,
    )

    read2_junction = alignment_crosses_origin(
        original2_0,
        mate2["reference_span"],
        ref_len,
    )

    pair_junction = (
        read1_junction
        or read2_junction
    )

    record = {
        "pair_id": pair_id,

        "ref1": ref_name,
        "start1_1": original1_0 + 1,
        "strand1": mate1["strand"],
        "cigar1": mate1["cigar"],
        "reference_span1": mate1["reference_span"],

        "ref2": ref_name,
        "start1_2": original2_0 + 1,
        "strand2": mate2["strand"],
        "cigar2": mate2["cigar"],
        "reference_span2": mate2["reference_span"],

        "pair_type": pair_type,
        "iteration": iteration,

        "rotation_start0": rotation_start0,
        "reference_length": ref_len,

        "art_start1_0": mate1["start0"],
        "art_start2_0": mate2["start0"],

        "art_pos1_1": mate1["pos1"],
        "art_pos1_2": mate2["pos1"],

        "art_tlen1": mate1["tlen"],
        "art_tlen2": mate2["tlen"],

        "read1_crosses_origin": int(
            read1_junction
        ),

        "read2_crosses_origin": int(
            read2_junction
        ),

        "pair_crosses_origin": int(
            pair_junction
        ),
    }

    return record


# ============================================================
# BACKGROUND SIMULATION
# ============================================================

def art_pair_count(
    depth,
    reference_length,
    read_length,
):
    return int(
        round(
            depth
            * reference_length
            / (2.0 * read_length),
            0,
        )
    )


def simulate_background(
    sample,
    ref_name,
    seq,
    art,
    iterations,
    target_pairs,
    read_length,
    opts,
    rng,
    out1,
    out2,
    truth_writer,
    rotation_writer,
    workdir,
):
    """
    Simulate exactly target_pairs ordinary/background fragments.

    The initial requested total is divided across whole-replicon rotations.
    If ART returns fewer usable pairs than requested, additional rotations
    are generated until the exact target is reached.
    """
    ref_len = len(seq)

    if target_pairs <= 0:
        return 0, 0

    base = target_pairs // iterations
    remainder = target_pairs % iterations

    total_pairs = 0
    junction_pairs = 0
    iteration = 0

    while total_pairs < target_pairs:
        iteration += 1

        if iteration <= iterations:
            requested_pairs = (
                base
                + (1 if iteration <= remainder else 0)
            )
        else:
            # Top up any deficit left by ART.
            requested_pairs = target_pairs - total_pairs

        if requested_pairs <= 0:
            continue

        rotation_start0 = rng.randint(
            0,
            ref_len - 1,
        )

        rotated = (
            seq[rotation_start0:]
            + seq[:rotation_start0]
        )

        art_seed = rng.randint(
            1,
            2_147_483_646,
        )

        if iteration <= iterations:
            label = "background iteration"
        else:
            label = "background top-up"

        print(
            f"    {label} {iteration:02d}: "
            f"rotation={rotation_start0:,}, "
            f"pairs={requested_pairs:,}"
        )

        reads1, reads2, truth = simulate_art(
            art=art,
            reference_name=ref_name,
            sequence=rotated,
            pair_count=requested_pairs,
            read_length=read_length,
            opts=opts,
            seed=art_seed,
            workdir=workdir,
            prefix_name=(
                f"{sample}_{ref_name}_"
                f"background_{iteration:03d}"
            ),
        )

        retained = 0

        for art_name in reads1.keys():

            # Never retain more than the exact genome-level
            # target requested for this replicon.
            if total_pairs >= target_pairs:
                break

            mate1 = truth[art_name][1]
            mate2 = truth[art_name][2]

            pair_number = total_pairs + 1

            pair_id = (
                f"{sample}__{ref_name}__"
                f"bg{iteration:03d}__"
                f"pair{pair_number:09d}"
            )

            write_fastq_record(
                out1,
                pair_id,
                1,
                reads1[art_name],
            )

            write_fastq_record(
                out2,
                pair_id,
                2,
                reads2[art_name],
            )

            record = build_truth_record(
                pair_id=pair_id,
                ref_name=ref_name,
                mate1=mate1,
                mate2=mate2,
                rotation_start0=rotation_start0,
                ref_len=ref_len,
                iteration=iteration,
                pair_type="background",
            )

            truth_writer.writerow(
                record
            )

            total_pairs += 1
            retained += 1

            if record["pair_crosses_origin"]:
                junction_pairs += 1

        rotation_writer.writerow({
            "contig": ref_name,
            "pair_type": "background",
            "iteration": iteration,
            "rotation_start0": rotation_start0,
            "reference_length": ref_len,
            "pairs_simulated": len(reads1),
            "pairs_retained": retained,
        })

        # Prevent an infinite loop if ART unexpectedly
        # produces no usable paired records.
        if retained == 0:
            raise RuntimeError(
                f"ART produced no usable background pairs "
                f"for {ref_name} during iteration "
                f"{iteration}."
            )

    if total_pairs != target_pairs:
        raise RuntimeError(
            f"Background simulation produced "
            f"{total_pairs:,} pairs; expected exactly "
            f"{target_pairs:,}."
        )

    return total_pairs, junction_pairs


# ============================================================
# JUNCTION CHALLENGE
# ============================================================

def simulate_junction_challenge(
    sample,
    ref_name,
    seq,
    art,
    target_pairs,
    read_length,
    opts,
    rng,
    out1,
    out2,
    truth_writer,
    rotation_writer,
    workdir,
):
    """
    Produce exactly target_pairs true origin-crossing fragments.

    Instead of simulating across the entire bacterial chromosome and discarding
    almost everything, ART is run on a local linear window taken from the
    circular chromosome.  The biological origin is placed at the centre of
    that window, safely away from the artificial template ends.

    ART SAM coordinates are transformed back to the original chromosome using
    the same modulo coordinate transform used by the full-length rotations.
    Only pairs for which R1 and/or R2 truly crosses L -> 0 are retained.
    """
    if target_pairs <= 0:
        return 0

    ref_len = len(seq)

    # A 2-kb window is comfortably wider than the default ~300-bp fragment
    # while enriching true origin-crossing fragments enough to make 100k
    # challenge pairs practical.  Increase automatically for long inserts.
    window_len = min(
        ref_len,
        max(2000, read_length * 12),
    )

    # Start half a window before coordinate 0.  Therefore original coordinate
    # 0 lies near the centre of the ART template.
    rotation_start0 = (ref_len - window_len // 2) % ref_len

    doubled = seq + seq
    local_template = doubled[
        rotation_start0:rotation_start0 + window_len
    ]

    if len(local_template) != window_len:
        raise RuntimeError("Could not construct circular junction template.")

    collected = 0
    batch = 0

    while collected < target_pairs:
        batch += 1
        remaining = target_pairs - collected

        # For a ~300-bp fragment in a 2-kb template, roughly 10-15% of
        # fragments are expected to overlap the central origin.  Oversample
        # conservatively, then retain only truth-confirmed crossings.
        batch_pairs = max(50_000, remaining * 10)

        art_seed = rng.randint(1, 2_147_483_646)

        print(
            f"    junction batch {batch:02d}: "
            f"window={window_len:,} bp, "
            f"simulating={batch_pairs:,}, "
            f"collected={collected:,}/{target_pairs:,}"
        )

        reads1, reads2, truth = simulate_art(
            art=art,
            reference_name=ref_name,
            sequence=local_template,
            pair_count=batch_pairs,
            read_length=read_length,
            opts=opts,
            seed=art_seed,
            workdir=workdir,
            prefix_name=(
                f"{sample}_{ref_name}_"
                f"junction_{batch:03d}"
            ),
        )

        retained_this_batch = 0

        for art_name in reads1.keys():
            mate1 = truth[art_name][1]
            mate2 = truth[art_name][2]

            temporary_record = build_truth_record(
                pair_id="temporary",
                ref_name=ref_name,
                mate1=mate1,
                mate2=mate2,
                rotation_start0=rotation_start0,
                ref_len=ref_len,
                iteration=batch,
                pair_type="junction",
            )

            if not temporary_record["pair_crosses_origin"]:
                continue

            collected += 1
            retained_this_batch += 1

            pair_id = (
                f"{sample}__{ref_name}__"
                f"junction__pair{collected:09d}"
            )

            write_fastq_record(out1, pair_id, 1, reads1[art_name])
            write_fastq_record(out2, pair_id, 2, reads2[art_name])

            temporary_record["pair_id"] = pair_id
            truth_writer.writerow(temporary_record)

            if collected >= target_pairs:
                break

        rotation_writer.writerow({
            "contig": ref_name,
            "pair_type": "junction_challenge",
            "iteration": batch,
            "rotation_start0": rotation_start0,
            "reference_length": ref_len,
            "pairs_simulated": len(reads1),
            "pairs_retained": retained_this_batch,
        })

        if retained_this_batch == 0:
            raise RuntimeError(
                "Junction simulation retained zero crossing pairs. "
                "Check ART insert-size options and read length."
            )

        if batch >= 1000 and collected < target_pairs:
            raise RuntimeError(
                "Unable to obtain requested number of junction pairs. "
                f"Collected {collected:,}/{target_pairs:,}."
            )

    return collected


# ============================================================
# COMPRESSION
# ============================================================

def gzip_file(path):
    gz_path = Path(
        str(path) + ".gz"
    )

    with open(path, "rb") as src:
        with gzip.open(
            gz_path,
            "wb",
        ) as dst:
            shutil.copyfileobj(
                src,
                dst,
            )

    path.unlink()

    return gz_path


def allocate_counts(total, weights):
    """
    Allocate an exact integer total according to weights.
    """
    if total < 0:
        raise ValueError("total must be >= 0")

    if not weights:
        raise ValueError("weights must not be empty")

    weight_sum = sum(weights)

    if weight_sum <= 0:
        raise ValueError("weights must sum to > 0")

    raw = [
        total * w / weight_sum
        for w in weights
    ]

    counts = [
        int(x)
        for x in raw
    ]

    remainder = total - sum(counts)

    order = sorted(
        range(len(weights)),
        key=lambda i: raw[i] - counts[i],
        reverse=True,
    )

    for i in order[:remainder]:
        counts[i] += 1

    return counts    


# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate a truth-aware paired-end RNA-seq benchmark for one or more "
            "circular bacterial replicons. The script uses whole-replicon rotations "
            "for background simulation, retains a requested number of true "
            "origin-spanning challenge pairs, and derives mate-level ground-truth "
            "coordinates from classic ART SAM output."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples
--------
Generate the manuscript-style 1,000,000-pair benchmark:

  python3 readSimulator_truth_1M.py \\
      --input reference.fna \\
      --outdir simulation \\
      --iterations 10 \\
      --readlen 150 \\
      --total-pairs 1000000 \\
      --junction-pairs 100000 \\
      --opts "-m 300 -s 30" \\
      --seed 12345

Use an explicitly located ART executable:

  python3 readSimulator_truth_1M.py \\
      --input reference.fna \\
      --simulator_path /path/to/art_illumina \\
      --outdir simulation

Principal outputs
-----------------
For an input named reference.fna, the output directory contains files including:

  reference_1.fastq
  reference_2.fastq
  reference_truth.tsv

The truth TSV records expected reference, 1-based start, strand,
origin-crossing status, and reference length for both mates. It is consumed
directly by compare_aligners.py and evaluate_pe_bam.py.

Notes
-----
* --total-pairs is the exact final number of read pairs.
* --junction-pairs specifies how many truth-confirmed origin-spanning challenge
  pairs are retained in that total.
* --opts is passed to classic ART; quote it as one shell argument.
* Use a fixed --seed when exact reproducibility is required.
"""
    )

    parser.add_argument(
        "--input",
        required=True,
        help=(
            "Input FASTA containing one or more circular replicons. "
            "All sequences in the file are treated as circular."
        ),
    )

    parser.add_argument(
        "--simulator_path",
        default="art",
        help=(
            "Classic ART art_illumina executable "
            "(default: art)"
        ),
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help=(
            "Directory for simulated FASTQ files, truth TSV, rotation metadata, "
            "and temporary simulation products."
        ),
    )

    parser.add_argument(
        "--iterations",
        type=int,
        default=10,
        help=(
            "Number of random whole-replicon "
            "rotations for background simulation "
            "(default: 10)"
        ),
    )

    parser.add_argument(
        "--readlen",
        type=int,
        default=150,
        help=(
            "Read length in bp "
            "(default: 150)"
        ),
    )

    parser.add_argument(
        "--total-pairs",
        type=int,
        default=1_000_000,
        help=(
            "Exact total number of output read pairs "
            "(default: 1000000)"
        ),
    )

    parser.add_argument(
        "--junction-pairs",
        type=int,
        default=100_000,
        help=(
            "Additional true origin-spanning "
            "challenge pairs to retain in the dataset "
            "(default: 100000)"
        ),
    )

    parser.add_argument(
        "--opts",
        default="-m 300 -s 30",
        help=(
            "Additional ART paired-end options "
            '(default: "-m 300 -s 30")'
        ),
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=12345,
        help=(
            "Random seed "
            "(default: 12345)"
        ),
    )

    args = parser.parse_args()

    if args.iterations < 1:
        parser.error(
            "--iterations must be >= 1"
        )

    if args.readlen < 1:
        parser.error(
            "--readlen must be >= 1"
        )

    if args.total_pairs < 1:
        parser.error(
            "--total-pairs must be >= 1"
        )

    if args.junction_pairs < 0:
        parser.error(
            "--junction-pairs must be >= 0"
        )

    if args.junction_pairs > args.total_pairs:
        parser.error(
            "--junction-pairs cannot exceed --total-pairs"
        )

    input_path = Path(
        args.input
    )

    if not input_path.is_file():
        parser.error(
            f"Input FASTA not found: "
            f"{input_path}"
        )

    art = shutil.which(
        args.simulator_path
    )

    if art is None:
        candidate = Path(
            args.simulator_path
        )

        if candidate.is_file():
            art = str(
                candidate.resolve()
            )
        else:
            parser.error(
                f"ART executable not found: "
                f"{args.simulator_path}"
            )

    outdir = Path(
        args.outdir
    )

    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    sample = input_path.name

    for suffix in (
        ".fasta.gz",
        ".fna.gz",
        ".fa.gz",
        ".fasta",
        ".fna",
        ".fa",
        ".gz",
    ):
        if sample.endswith(
            suffix
        ):
            sample = sample[
                :-len(suffix)
            ]
            break

    records = read_fasta(
        input_path
    )

    background_target = (
        args.total_pairs - args.junction_pairs
    )

    replicon_lengths = [
    len(seq)
    for _, seq in records
    ]

    # Background reads are distributed in proportion
    # to replicon length.
    background_allocations = allocate_counts(
        background_target,
        replicon_lengths,
    )

    # Junction challenge reads are distributed equally
    # among circular replicons.
    junction_allocations = allocate_counts(
        args.junction_pairs,
        [1] * len(records),
    )

    rng = random.Random(
        args.seed
    )

    opts = shlex.split(
        args.opts
    )

    # Prevent options controlled by this script from being supplied
    # a second time through --opts.
    forbidden_opts = {
        "-i", "--in",
        "-c", "--rcount",
        "-l", "--len",
        "-p", "--paired",
        "-o", "--out",
        "-rs", "--rndSeed",
        "-sam", "--samout",
        "-ss", "--seqSys",
    }

    supplied_forbidden = [
        x
        for x in opts
        if x in forbidden_opts
    ]

    if supplied_forbidden:
        parser.error(
            "--opts contains options controlled "
            "by readSimulator_truth.py: "
            + ", ".join(
                supplied_forbidden
            )
        )

    r1_path = (
        outdir
        / f"{sample}_1.fastq"
    )

    r2_path = (
        outdir
        / f"{sample}_2.fastq"
    )

    truth_path = (
        outdir
        / f"{sample}_truth.tsv"
    )

    rotation_path = (
        outdir
        / f"{sample}_rotations.tsv"
    )

    print()
    print("=" * 70)
    print(
        "Truth-aware circular read simulation"
    )
    print("=" * 70)
    print(
        f"Sample:             {sample}"
    )
    print(
        f"Circular replicons: "
        f"{len(records)}"
    )
    print(
        f"Iterations:         "
        f"{args.iterations}"
    )
    print(
        f"Read length:        "
        f"{args.readlen}"
    )
    print(
        f"Total pairs:        "
        f"{args.total_pairs:,}"
    )
    print(
        f"Background pairs:   "
        f"{background_target:,}"
    )
    print(
        f"Junction pairs:     "
        f"{args.junction_pairs:,}"
    )
    print(
        f"Random seed:        "
        f"{args.seed}"
    )
    print(
        f"ART:                "
        f"{art}"
    )
    print()

    print("Replicon allocation:")
    for i, (ref_name, seq) in enumerate(records):
        print(
            f"  {ref_name}: "
            f"{len(seq):,} bp | "
            f"background={background_allocations[i]:,} | "
            f"junction={junction_allocations[i]:,}"
        )
    print()

    total_background = 0
    total_background_junction = 0
    total_challenge = 0

    with tempfile.TemporaryDirectory(
        prefix="readSimulator_truth_"
    ) as tmpdir:

        workdir = Path(
            tmpdir
        )

        with (
            open(
                r1_path,
                "w",
            ) as out1,

            open(
                r2_path,
                "w",
            ) as out2,

            open(
                truth_path,
                "w",
                newline="",
            ) as truth_handle,

            open(
                rotation_path,
                "w",
                newline="",
            ) as rotation_handle,
        ):

            truth_writer = (
                csv.DictWriter(
                    truth_handle,
                    fieldnames=TRUTH_FIELDS,
                    delimiter="\t",
                    lineterminator="\n",
                )
            )

            truth_writer.writeheader()

            rotation_writer = (
                csv.DictWriter(
                    rotation_handle,
                    fieldnames=ROTATION_FIELDS,
                    delimiter="\t",
                    lineterminator="\n",
                )
            )

            rotation_writer.writeheader()

            for replicon_index, (ref_name, seq) in enumerate(records):

                replicon_background_target = (
                    background_allocations[replicon_index]
                )

                replicon_junction_target = (
                    junction_allocations[replicon_index]
                )

                ref_workdir = (
                    workdir
                    / ref_name
                )
                ref_workdir.mkdir(
                    parents=True,
                    exist_ok=True,
                )

                print(
                    f"Contig: {ref_name} "
                    f"({len(seq):,} bp)"
                )

                (
                    background_pairs,
                    natural_junctions,
                ) = simulate_background(
                    sample=sample,
                    ref_name=ref_name,
                    seq=seq,
                    art=art,
                    iterations=(
                        args.iterations
                    ),
                    target_pairs=(
                        replicon_background_target
                    ),
                    read_length=(
                        args.readlen
                    ),
                    opts=opts,
                    rng=rng,
                    out1=out1,
                    out2=out2,
                    truth_writer=(
                        truth_writer
                    ),
                    rotation_writer=(
                        rotation_writer
                    ),
                    workdir=ref_workdir,
                )

                print(
                    "  background pairs:          "
                    f"{background_pairs:,}"
                )

                print(
                    "  natural junction pairs:    "
                    f"{natural_junctions:,}"
                )

                challenge_pairs = (
                    simulate_junction_challenge(
                        sample=sample,
                        ref_name=ref_name,
                        seq=seq,
                        art=art,
                        target_pairs=(
                            replicon_junction_target
                        ),
                        read_length=(
                            args.readlen
                        ),
                        opts=opts,
                        rng=rng,
                        out1=out1,
                        out2=out2,
                        truth_writer=(
                            truth_writer
                        ),
                        rotation_writer=(
                            rotation_writer
                        ),
                        workdir=ref_workdir,
                    )
                )

                print(
                    "  challenge junction pairs:  "
                    f"{challenge_pairs:,}"
                )

                print()

                total_background += (
                    background_pairs
                )

                total_background_junction += (
                    natural_junctions
                )

                total_challenge += (
                    challenge_pairs
                )

    print(
        "Compressing FASTQ files..."
    )

    r1_gz = gzip_file(
        r1_path
    )

    r2_gz = gzip_file(
        r2_path
    )

    total_pairs = (
        total_background
        + total_challenge
    )

    if total_pairs != args.total_pairs:
        raise RuntimeError(
            f"Final dataset contains {total_pairs:,} pairs; "
            f"expected exactly {args.total_pairs:,}."
        )

    if total_challenge != args.junction_pairs:
        raise RuntimeError(
            f"Final challenge set contains {total_challenge:,} pairs; "
            f"expected exactly {args.junction_pairs:,}."
        )

    print()
    print("=" * 70)
    print(
        "Simulation completed successfully"
    )
    print("=" * 70)

    print(
        "Background pairs:          "
        f"{total_background:,}"
    )

    print(
        "Natural junction pairs:    "
        f"{total_background_junction:,}"
    )

    print(
        "Challenge junction pairs:  "
        f"{total_challenge:,}"
    )

    print(
        "Total output pairs:        "
        f"{total_pairs:,}"
    )

    print(
        f"R1:                        "
        f"{r1_gz}"
    )

    print(
        f"R2:                        "
        f"{r2_gz}"
    )

    print(
        f"Truth table:               "
        f"{truth_path}"
    )

    print(
        f"Rotation table:            "
        f"{rotation_path}"
    )

    print("=" * 70)


if __name__ == "__main__":
    main()
