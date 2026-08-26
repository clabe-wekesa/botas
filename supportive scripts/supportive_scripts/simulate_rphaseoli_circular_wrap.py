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

            if not line:
                continue

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


def revcomp(seq):
    """Return the reverse complement of a DNA sequence."""
    table = str.maketrans(
        "ACGTNacgtn",
        "TGCANtgcan",
    )
    return seq.translate(table)[::-1].upper()


def quality_string(length):
    return "I" * length


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate paired-end reads forced to span the circular "
            "boundaries of R. phaseoli replicons."
        )
    )

    parser.add_argument(
        "--ref",
        default=None,
        help=(
            "Reference FASTA. Default: "
            "references/rphaseoli/genome.fa relative to the archive root"
        ),
    )

    parser.add_argument(
        "--outdir",
        default=None,
        help=(
            "Output directory. Default: "
            "botas_manuscript_data/simulated/circular_wrap/rphaseoli"
        ),
    )

    parser.add_argument(
        "--pairs-per-contig",
        type=int,
        default=5000,
        help="Number of read pairs per replicon (default: 5000)",
    )

    parser.add_argument(
        "--read-length",
        type=int,
        default=150,
        help="Read length in bp (default: 150)",
    )

    parser.add_argument(
        "--insert-mean",
        type=float,
        default=300,
        help="Mean fragment length in bp (default: 300)",
    )

    parser.add_argument(
        "--insert-sd",
        type=float,
        default=50,
        help="Fragment-length standard deviation in bp (default: 50)",
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
            / "rphaseoli"
        )
    )

    if not ref.is_file():
        raise FileNotFoundError(
            f"Reference FASTA not found: {ref}"
        )

    if args.pairs_per_contig <= 0:
        raise ValueError("--pairs-per-contig must be > 0")

    if args.read_length <= 0:
        raise ValueError("--read-length must be > 0")

    if args.insert_mean <= 0:
        raise ValueError("--insert-mean must be > 0")

    if args.insert_sd < 0:
        raise ValueError("--insert-sd must be >= 0")

    outdir.mkdir(parents=True, exist_ok=True)

    random.seed(args.seed)

    seqs = read_fasta(ref)

    r1_path = outdir / "rphaseoli_wrap_R1.fq.gz"
    r2_path = outdir / "rphaseoli_wrap_R2.fq.gz"
    truth_path = outdir / "rphaseoli_wrap_truth.tsv"

    pair_id = 0

    with (
        gzip.open(r1_path, "wt") as r1,
        gzip.open(r2_path, "wt") as r2,
        open(truth_path, "w") as truth,
    ):

        truth.write(
            "pair_id\tref1\tstart1_1\tref2\tstart1_2\n"
        )

        for contig, seq in seqs.items():

            length = len(seq)

            for _ in range(args.pairs_per_contig):

                pair_id += 1

                insert = max(
                    args.read_length * 2,
                    int(
                        random.gauss(
                            args.insert_mean,
                            args.insert_sd,
                        )
                    ),
                )

                if insert >= length:
                    raise ValueError(
                        f"Generated insert length ({insert}) is not "
                        f"appropriate for reference sequence {contig} "
                        f"of length {length}."
                    )

                # Choose a fragment start that guarantees that the
                # fragment crosses the end/origin boundary.
                start = random.randint(
                    length - insert + 1,
                    length - 1,
                )

                circular_fragment = (
                    seq + seq
                )[start : start + insert]

                read1 = circular_fragment[: args.read_length]

                read2 = revcomp(
                    circular_fragment[-args.read_length :]
                )

                # Truth coordinates are stored as 0-based positions.
                pos1 = start % length

                pos2 = (
                    start
                    + insert
                    - args.read_length
                ) % length

                qname = (
                    f"rphaseoli_wrap_{pair_id}_{contig}"
                )

                r1.write(
                    f"@{qname}/1\n"
                    f"{read1}\n"
                    f"+\n"
                    f"{quality_string(args.read_length)}\n"
                )

                r2.write(
                    f"@{qname}/2\n"
                    f"{read2}\n"
                    f"+\n"
                    f"{quality_string(args.read_length)}\n"
                )

                truth.write(
                    f"{qname}\t"
                    f"{contig}\t"
                    f"{pos1}\t"
                    f"{contig}\t"
                    f"{pos2}\n"
                )

    print(f"Reference: {ref}")
    print(f"Random seed: {args.seed}")
    print(f"Read length: {args.read_length}")
    print(f"Mean fragment length: {args.insert_mean}")
    print(f"Fragment-length SD: {args.insert_sd}")
    print(f"Pairs per replicon: {args.pairs_per_contig}")
    print(f"Total pairs: {pair_id}")
    print(f"Wrote R1: {r1_path}")
    print(f"Wrote R2: {r2_path}")
    print(f"Wrote truth: {truth_path}")


if __name__ == "__main__":
    main()