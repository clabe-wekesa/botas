#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

import pysam


def load_truth(path):
    """Load expected single-end alignment positions from a TSV file."""
    truth = {}

    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        required = {"read_id", "ref", "start"}
        missing = required - set(reader.fieldnames or [])

        if missing:
            raise ValueError(
                f"Truth file is missing required columns: "
                f"{', '.join(sorted(missing))}"
            )

        for row in reader:
            truth[row["read_id"]] = (
                row["ref"],
                int(row["start"]),
            )

    return truth


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Evaluate single-end BAM alignments against known genomic "
            "positions from a simulated-read truth table."
        )
    )

    parser.add_argument(
        "--bam",
        required=True,
        help="BAM file to evaluate",
    )

    parser.add_argument(
        "--truth",
        required=True,
        help="TSV truth table containing expected positions",
    )

    parser.add_argument(
        "--out",
        required=True,
        help="Output TSV file for evaluation metrics",
    )

    parser.add_argument(
        "--tolerance",
        type=int,
        default=10,
        help="Maximum positional difference allowed in bp (default: 10)",
    )

    args = parser.parse_args()

    bam_path = Path(args.bam)
    truth_path = Path(args.truth)
    out_path = Path(args.out)

    if not bam_path.is_file():
        raise FileNotFoundError(f"BAM file not found: {bam_path}")

    if not truth_path.is_file():
        raise FileNotFoundError(f"Truth file not found: {truth_path}")

    if args.tolerance < 0:
        raise ValueError("--tolerance must be >= 0")

    truth = load_truth(truth_path)

    total = len(truth)
    mapped = 0
    correct = 0

    # Track reads already evaluated so that secondary/supplementary
    # records or duplicate primary records are not counted repeatedly.
    evaluated = set()

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:

        for aln in bam.fetch(until_eof=True):

            if aln.is_secondary or aln.is_supplementary:
                continue

            qname = aln.query_name

            if qname not in truth:
                continue

            if qname in evaluated:
                continue

            evaluated.add(qname)

            if aln.is_unmapped:
                continue

            mapped += 1

            expected_ref, expected_pos = truth[qname]

            # reference_start is 0-based. The supplied SE truth tables
            # use the same coordinate convention.
            observed_pos = aln.reference_start

            if (
                aln.reference_name == expected_ref
                and abs(observed_pos - expected_pos) <= args.tolerance
            ):
                correct += 1

    unmapped = total - mapped
    wrong = mapped - correct

    precision = correct / mapped if mapped else 0.0
    recall = correct / total if total else 0.0

    f1 = (
        2 * precision * recall / (precision + recall)
        if precision + recall > 0
        else 0.0
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"total_reads\t{total}\n")
        fh.write(f"mapped_reads\t{mapped}\n")
        fh.write(f"correct_reads\t{correct}\n")
        fh.write(f"wrong_reads\t{wrong}\n")
        fh.write(f"unmapped_reads\t{unmapped}\n")
        fh.write(f"precision\t{precision:.6f}\n")
        fh.write(f"recall\t{recall:.6f}\n")
        fh.write(f"f1\t{f1:.6f}\n")
        fh.write(f"tolerance_bp\t{args.tolerance}\n")

    print(f"Wrote metrics to: {out_path}")
    print(f"Precision: {precision:.6f}")
    print(f"Recall:    {recall:.6f}")
    print(f"F1:        {f1:.6f}")


if __name__ == "__main__":
    main()