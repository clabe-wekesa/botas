#!/usr/bin/env python3

import argparse
import gzip
import random
from pathlib import Path


def read_fasta(path):
    """Read a FASTA file into a dictionary keyed by sequence name."""
    seqs = {}
    name = None
    chunks = []

    with open(path) as fh:
        for line in fh:
            line = line.strip()

            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()

                name = line[1:].split()[0]
                chunks = []

            else:
                chunks.append(line)

        if name is not None:
            seqs[name] = "".join(chunks).upper()

    return seqs


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate single-end reads forced to span the circular "
            "boundaries of R. phaseoli replicons."
        )
    )

    parser.add_argument(
        "--ref",
        default=None,
        help=(
            "Reference FASTA. Default: "
            "../references/rphaseoli/genome.fa relative to the archive root"
        ),
    )

    parser.add_argument(
        "--outdir",
        default=None,
        help=(
            "Output directory. Default: "
            "botas_manuscript_data/simulated/circular_wrap/rphaseoli_se"
        ),
    )

    parser.add_argument(
        "--reads-per-contig",
        type=int,
        default=10000,
        help="Number of reads generated per reference sequence (default: 10000)",
    )

    parser.add_argument(
        "--read-length",
        type=int,
        default=150,
        help="Read length in bp (default: 150)",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed (default: 42)",
    )

    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent

    ref = (
        Path(args.ref)
        if args.ref
        else root / "references" / "rphaseoli" / "genome.fa"
    )

    outdir = (
        Path(args.outdir)
        if args.outdir
        else (
            root
            / "botas_manuscript_data"
            / "simulated"
            / "circular_wrap"
            / "rphaseoli_se"
        )
    )

    if not ref.is_file():
        raise FileNotFoundError(f"Reference FASTA not found: {ref}")

    if args.reads_per_contig <= 0:
        raise ValueError("--reads-per-contig must be > 0")

    if args.read_length <= 1:
        raise ValueError("--read-length must be > 1")

    outdir.mkdir(parents=True, exist_ok=True)

    random.seed(args.seed)

    seqs = read_fasta(ref)

    fq = outdir / "rphaseoli_wrap_SE.fq.gz"
    truth = outdir / "rphaseoli_wrap_SE_truth.tsv"

    read_id = 0

    with gzip.open(fq, "wt") as out, open(truth, "w") as tr:

        tr.write("read_id\tref\tstart\n")

        for contig, seq in seqs.items():

            length = len(seq)

            if length < args.read_length:
                raise ValueError(
                    f"Reference sequence {contig} is shorter than "
                    f"the requested read length ({args.read_length} bp)"
                )

            for _ in range(args.reads_per_contig):

                read_id += 1

                # Choose a start position that guarantees that the read
                # crosses the end/origin boundary of the circular replicon.
                start = random.randint(
                    length - args.read_length + 1,
                    length - 1,
                )

                circular_seq = seq + seq

                read = circular_seq[
                    start : start + args.read_length
                ]

                qname = (
                    f"rphaseoli_se_wrap_{read_id}_{contig}"
                )

                out.write(
                    f"@{qname}\n"
                    f"{read}\n"
                    f"+\n"
                    f"{'I' * args.read_length}\n"
                )

                # Coordinates in this truth table are 0-based,
                # matching pysam.AlignmentSegment.reference_start.
                tr.write(
                    f"{qname}\t{contig}\t{start % length}\n"
                )

    print(f"Reference: {ref}")
    print(f"Random seed: {args.seed}")
    print(f"Read length: {args.read_length}")
    print(f"Reads per replicon: {args.reads_per_contig}")
    print(f"Total reads: {read_id}")
    print(f"Wrote FASTQ: {fq}")
    print(f"Wrote truth: {truth}")


if __name__ == "__main__":
    main()