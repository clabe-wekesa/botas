#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Convert the raw 0-based R. phaseoli circular-wrap truth table "
            "to the 1-based format used by evaluate_pe_bam.py."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Raw circular-wrap truth TSV",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output 1-based evaluation truth TSV",
    )

    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.is_file():
        raise FileNotFoundError(f"Input truth file not found: {input_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(input_path, newline="") as infile, open(
        output_path, "w", newline=""
    ) as outfile:

        reader = csv.DictReader(infile, delimiter="\t")

        required = {
            "qname",
            "contig",
            "pos1",
            "pos2",
        }

        missing = required - set(reader.fieldnames or [])

        if missing:
            raise ValueError(
                "Input truth file is missing required columns: "
                + ", ".join(sorted(missing))
            )

        writer = csv.writer(outfile, delimiter="\t")

        writer.writerow(
            [
                "pair_id",
                "ref1",
                "start1_1",
                "ref2",
                "start1_2",
            ]
        )

        for row in reader:
            writer.writerow(
                [
                    row["qname"],
                    row["contig"],
                    int(row["pos1"]) + 1,
                    row["contig"],
                    int(row["pos2"]) + 1,
                ]
            )

    print(f"Wrote evaluation truth table: {output_path}")


if __name__ == "__main__":
    main()