#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

import pysam


def norm_qname(qname):
    """Remove common /1 and /2 mate suffixes from read names."""
    return qname.replace("/1", "").replace("/2", "")


def load_truth(path):
    """Load paired-end alignment truth positions from a TSV file."""
    truth = {}

    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")

        required = {"pair_id", "ref1", "start1_1", "ref2", "start1_2"}
        missing = required - set(reader.fieldnames or [])

        if missing:
            raise ValueError(
                f"Truth file is missing required columns: {', '.join(sorted(missing))}"
            )

        for row in reader:
            truth[row["pair_id"]] = {
                "ref1": row["ref1"],
                "pos1": int(row["start1_1"]),
                "ref2": row["ref2"],
                "pos2": int(row["start1_2"]),
            }

    return truth


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Evaluate paired-end BAM alignments against known genomic "
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

    if not bam_path.is_file():
        raise FileNotFoundError(f"BAM file not found: {bam_path}")

    if not truth_path.is_file():
        raise FileNotFoundError(f"Truth file not found: {truth_path}")

    if args.tolerance < 0:
        raise ValueError("--tolerance must be >= 0")

    truth = load_truth(truth_path)

    mate_status = {}

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for aln in bam.fetch(until_eof=True):

            # Only primary alignments are evaluated.
            if aln.is_secondary or aln.is_supplementary:
                continue

            pid = norm_qname(aln.query_name)

            if pid not in truth:
                continue

            if aln.is_read1:
                mate = 1
            elif aln.is_read2:
                mate = 2
            else:
                continue

            key = (pid, mate)

            # Evaluate only the first primary record observed for each mate.
            if key in mate_status:
                continue

            if aln.is_unmapped:
                mate_status[key] = "unmapped"
                continue

            t = truth[pid]

            if mate == 1:
                expected_ref = t["ref1"]
                expected_pos = t["pos1"]
            else:
                expected_ref = t["ref2"]
                expected_pos = t["pos2"]

            observed_ref = bam.get_reference_name(aln.reference_id)

            # pysam uses 0-based coordinates; truth tables use 1-based starts.
            observed_pos = aln.reference_start + 1

            if (
                observed_ref == expected_ref
                and abs(observed_pos - expected_pos) <= args.tolerance
            ):
                mate_status[key] = "correct"
            else:
                mate_status[key] = "wrong"

    total_pairs = len(truth)
    total_reads = total_pairs * 2

    correct_reads = sum(
        1 for status in mate_status.values()
        if status == "correct"
    )

    wrong_reads = sum(
        1 for status in mate_status.values()
        if status == "wrong"
    )

    observed_unmapped = sum(
        1 for status in mate_status.values()
        if status == "unmapped"
    )

    missing_reads = total_reads - len(mate_status)
    unmapped_reads = observed_unmapped + missing_reads

    mapped_reads = correct_reads + wrong_reads

    precision = (
        correct_reads / mapped_reads
        if mapped_reads
        else 0.0
    )

    recall = (
        correct_reads / total_reads
        if total_reads
        else 0.0
    )

    f1 = (
        2 * precision * recall / (precision + recall)
        if precision + recall > 0
        else 0.0
    )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    with open(out, "w") as w:
        w.write("metric\tvalue\n")
        w.write(f"total_pairs\t{total_pairs}\n")
        w.write(f"total_reads\t{total_reads}\n")
        w.write(f"mapped_reads\t{mapped_reads}\n")
        w.write(f"correct_reads\t{correct_reads}\n")
        w.write(f"wrong_reads\t{wrong_reads}\n")
        w.write(f"unmapped_reads\t{unmapped_reads}\n")
        w.write(f"precision\t{precision:.6f}\n")
        w.write(f"recall\t{recall:.6f}\n")
        w.write(f"f1\t{f1:.6f}\n")
        w.write(f"tolerance_bp\t{args.tolerance}\n")

    print(f"Wrote metrics to: {out}")
    print(f"Precision: {precision:.6f}")
    print(f"Recall:    {recall:.6f}")
    print(f"F1:        {f1:.6f}")


if __name__ == "__main__":
    main()