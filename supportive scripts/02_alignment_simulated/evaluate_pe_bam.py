#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

import pysam


# ============================================================
# Read-name handling
# ============================================================

def norm_qname(qname):
    """
    Remove a terminal /1 or /2 mate suffix.

    Examples:
        read123/1 -> read123
        read123/2 -> read123
        read123   -> read123
    """
    if qname.endswith("/1") or qname.endswith("/2"):
        return qname[:-2]

    return qname


# ============================================================
# Truth loading
# ============================================================

def load_truth_tsv(path):
    """
    Load paired-end truth from a TSV file.

    Required columns:
        pair_id
        ref1
        start1_1
        ref2
        start1_2

    Truth positions are expected to be 1-based.
    """

    truth = {}

    with open(path, newline="") as f:

        reader = csv.DictReader(f, delimiter="\t")

        required = {
            "pair_id",
            "ref1",
            "start1_1",
            "ref2",
            "start1_2",
        }

        missing = required - set(reader.fieldnames or [])

        if missing:
            raise ValueError(
                "Truth TSV is missing required columns: "
                + ", ".join(sorted(missing))
            )

        for row in reader:

            pid = norm_qname(row["pair_id"])

            if pid in truth:
                raise ValueError(
                    f"Duplicate pair_id in truth TSV: {pid}"
                )

            truth[pid] = {
                "ref1": row["ref1"],
                "pos1": int(row["start1_1"]),
                "ref2": row["ref2"],
                "pos2": int(row["start1_2"]),
            }

    return truth


def load_truth_alignment(path):
    """
    Load paired-end truth directly from SAM or BAM.

    SAM/BAM reference_start is zero-based internally.
    It is converted here to 1-based coordinates so that
    the resulting representation matches the existing
    TSV truth convention.
    """

    truth = {}

    suffix = path.suffix.lower()

    if suffix == ".bam":
        mode = "rb"
    else:
        mode = "r"

    with pysam.AlignmentFile(str(path), mode) as aln_file:

        for aln in aln_file.fetch(until_eof=True):

            # Truth should be defined from primary alignments only.
            if aln.is_secondary or aln.is_supplementary:
                continue

            if aln.is_unmapped:
                continue

            # We require explicit paired-end mate designation.
            if aln.is_read1:
                mate = 1
            elif aln.is_read2:
                mate = 2
            else:
                continue

            pid = norm_qname(aln.query_name)

            ref = aln_file.get_reference_name(
                aln.reference_id
            )

            # Convert BAM/SAM 0-based coordinate to the
            # 1-based convention used by the TSV truth.
            pos = aln.reference_start + 1

            if pid not in truth:
                truth[pid] = {}

            if mate == 1:

                if "ref1" in truth[pid]:
                    raise ValueError(
                        f"Duplicate primary truth alignment "
                        f"for mate 1: {pid}"
                    )

                truth[pid]["ref1"] = ref
                truth[pid]["pos1"] = pos

            else:

                if "ref2" in truth[pid]:
                    raise ValueError(
                        f"Duplicate primary truth alignment "
                        f"for mate 2: {pid}"
                    )

                truth[pid]["ref2"] = ref
                truth[pid]["pos2"] = pos

    # Require complete paired-end truth.
    incomplete = []

    for pid, data in truth.items():

        required = {
            "ref1",
            "pos1",
            "ref2",
            "pos2",
        }

        if not required.issubset(data):
            incomplete.append(pid)

    if incomplete:

        examples = ", ".join(incomplete[:5])

        raise ValueError(
            f"Truth alignment contains "
            f"{len(incomplete)} incomplete read pairs. "
            f"Examples: {examples}"
        )

    return truth


def load_truth(path):
    """
    Automatically detect supported truth format.

    Supported:
        .tsv
        .sam
        .bam
    """

    suffix = path.suffix.lower()

    if suffix == ".tsv":
        return load_truth_tsv(path)

    if suffix in {".sam", ".bam"}:
        return load_truth_alignment(path)

    raise ValueError(
        f"Unsupported truth format: {path}\n"
        "Supported formats are TSV, SAM and BAM."
    )


# ============================================================
# BAM evaluation
# ============================================================

def evaluate_bam(bam_path, truth, tolerance):
    """
    Evaluate primary paired-end BAM alignments against truth.
    """

    mate_status = {}

    with pysam.AlignmentFile(
        str(bam_path), "rb"
    ) as bam:

        for aln in bam.fetch(until_eof=True):

            # Ignore non-primary records.
            if aln.is_secondary or aln.is_supplementary:
                continue

            pid = norm_qname(aln.query_name)

            # Ignore reads not belonging to the truth dataset.
            if pid not in truth:
                continue

            if aln.is_read1:
                mate = 1

            elif aln.is_read2:
                mate = 2

            else:
                continue

            key = (pid, mate)

            # Only one primary alignment should normally exist.
            # Preserve the behaviour of the previous evaluator:
            # use the first primary record encountered.
            if key in mate_status:
                continue

            if aln.is_unmapped:

                mate_status[key] = "unmapped"
                continue

            expected = truth[pid]

            if mate == 1:

                expected_ref = expected["ref1"]
                expected_pos = expected["pos1"]

            else:

                expected_ref = expected["ref2"]
                expected_pos = expected["pos2"]

            observed_ref = bam.get_reference_name(
                aln.reference_id
            )

            # pysam reference_start is zero-based.
            observed_pos = aln.reference_start + 1

            if (
                observed_ref == expected_ref
                and abs(
                    observed_pos - expected_pos
                ) <= tolerance
            ):
                mate_status[key] = "correct"

            else:
                mate_status[key] = "wrong"

    return mate_status


# ============================================================
# Metric calculation
# ============================================================

def calculate_metrics(truth, mate_status):
    """
    Calculate the same read-level metrics used by the
    previous BOTAS paired-end evaluator.
    """

    total_pairs = len(truth)
    total_reads = total_pairs * 2

    correct_reads = sum(
        status == "correct"
        for status in mate_status.values()
    )

    wrong_reads = sum(
        status == "wrong"
        for status in mate_status.values()
    )

    observed_unmapped = sum(
        status == "unmapped"
        for status in mate_status.values()
    )

    missing_reads = (
        total_reads - len(mate_status)
    )

    if missing_reads < 0:
        raise RuntimeError(
            "Observed more evaluated mates than expected "
            "from the truth dataset."
        )

    unmapped_reads = (
        observed_unmapped + missing_reads
    )

    mapped_reads = (
        correct_reads + wrong_reads
    )

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
        2 * precision * recall
        / (precision + recall)
        if precision + recall > 0
        else 0.0
    )

    return {
        "total_pairs": total_pairs,
        "total_reads": total_reads,
        "mapped_reads": mapped_reads,
        "correct_reads": correct_reads,
        "wrong_reads": wrong_reads,
        "unmapped_reads": unmapped_reads,
        "missing_reads": missing_reads,
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }


# ============================================================
# Output
# ============================================================

def write_metrics(path, metrics, tolerance, truth_format):

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with open(path, "w") as w:

        w.write("metric\tvalue\n")

        for key in [
            "total_pairs",
            "total_reads",
            "mapped_reads",
            "correct_reads",
            "wrong_reads",
            "unmapped_reads",
            "missing_reads",
        ]:
            w.write(
                f"{key}\t{metrics[key]}\n"
            )

        w.write(
            f"precision\t"
            f"{metrics['precision']:.6f}\n"
        )

        w.write(
            f"recall\t"
            f"{metrics['recall']:.6f}\n"
        )

        w.write(
            f"f1\t"
            f"{metrics['f1']:.6f}\n"
        )

        w.write(
            f"tolerance_bp\t{tolerance}\n"
        )

        w.write(
            f"truth_format\t{truth_format}\n"
        )


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Evaluate one paired-end BAM against known mate-level genomic truth. "
            "Truth may be supplied as the BOTAS simulator TSV, SAM, or BAM."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples
--------
Evaluate a BAM using the BOTAS simulator truth TSV:

  python3 evaluate_pe_bam.py \\
      --bam botas.sorted.bam \\
      --truth simulation/reference_truth.tsv \\
      --out botas_evaluation.tsv \\
      --tolerance 10

Evaluate against SAM truth:

  python3 evaluate_pe_bam.py \\
      --bam alignments.sorted.bam \\
      --truth truth.sam \\
      --out evaluation.tsv

Purpose
-------
Use this script for single-BAM validation and a compact machine-readable metric
file. For direct comparison of several aligners with separate JUNCTION and
NON_JUNCTION summaries, use compare_aligners.py.

Supported truth formats
-----------------------
.tsv   BOTAS truth TSV generated by readSimulator_truth_1M.py
.sam   paired-end SAM containing known alignments
.bam   paired-end BAM containing known alignments

--out is a tab-delimited metrics file containing counts plus precision, recall
and F1. A concise summary is also printed to standard output.
"""
    )

    parser.add_argument(
        "--bam",
        required=True,
        help="Primary paired-end BAM file to evaluate",
    )

    parser.add_argument(
        "--truth",
        required=True,
        help=(
            "Truth file. Supported formats: "
            ".tsv, .sam, .bam"
        ),
    )

    parser.add_argument(
        "--out",
        required=True,
        help="Output tab-delimited metrics file",
    )

    parser.add_argument(
        "--tolerance",
        type=int,
        default=10,
        help=(
            "Maximum allowed difference between observed "
            "and true alignment start in bp (default: 10)"
        ),
    )

    args = parser.parse_args()

    bam_path = Path(args.bam)
    truth_path = Path(args.truth)
    out_path = Path(args.out)

    if not bam_path.is_file():
        raise FileNotFoundError(
            f"BAM file not found: {bam_path}"
        )

    if not truth_path.is_file():
        raise FileNotFoundError(
            f"Truth file not found: {truth_path}"
        )

    if args.tolerance < 0:
        raise ValueError(
            "--tolerance must be >= 0"
        )

    print(
        f"Loading truth: {truth_path}"
    )

    truth = load_truth(truth_path)

    print(
        f"Truth pairs: {len(truth):,}"
    )

    print(
        f"Evaluating: {bam_path}"
    )

    mate_status = evaluate_bam(
        bam_path=bam_path,
        truth=truth,
        tolerance=args.tolerance,
    )

    metrics = calculate_metrics(
        truth=truth,
        mate_status=mate_status,
    )

    truth_format = (
        truth_path.suffix.lower()
        .lstrip(".")
    )

    write_metrics(
        path=out_path,
        metrics=metrics,
        tolerance=args.tolerance,
        truth_format=truth_format,
    )

    print()
    print(f"Wrote metrics to: {out_path}")
    print(
        f"Total pairs: {metrics['total_pairs']:,}"
    )
    print(
        f"Correct:     {metrics['correct_reads']:,}"
    )
    print(
        f"Wrong:       {metrics['wrong_reads']:,}"
    )
    print(
        f"Unmapped:    {metrics['unmapped_reads']:,}"
    )
    print(
        f"Precision:   {metrics['precision']:.6f}"
    )
    print(
        f"Recall:      {metrics['recall']:.6f}"
    )
    print(
        f"F1:          {metrics['f1']:.6f}"
    )


if __name__ == "__main__":
    main()
