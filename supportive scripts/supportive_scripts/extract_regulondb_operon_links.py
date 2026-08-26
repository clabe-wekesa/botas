#!/usr/bin/env python3

import argparse
from pathlib import Path

from bson import decode_file_iter


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract transcription units and adjacent-gene operon links "
            "from the RegulonDB operonDatamart BSON file."
        )
    )

    root = Path(__file__).resolve().parent.parent
    default_outdir = root / "results" / "operons" / "regulondb"

    parser.add_argument(
        "--bson",
        required=True,
        help="RegulonDB operonDatamart BSON file",
    )

    parser.add_argument(
        "--out-links",
        default=str(default_outdir / "regulondb_adjacent_gene_links.tsv"),
        help=(
            "Output TSV containing adjacent gene links "
            "(default: results/operons/regulondb/"
            "regulondb_adjacent_gene_links.tsv)"
        ),
    )

    parser.add_argument(
        "--out-operons",
        default=str(default_outdir / "regulondb_transcription_units.tsv"),
        help=(
            "Output TSV containing transcription units "
            "(default: results/operons/regulondb/"
            "regulondb_transcription_units.tsv)"
        ),
    )

    args = parser.parse_args()

    bson_path = Path(args.bson)
    links_path = Path(args.out_links)
    operons_path = Path(args.out_operons)

    if not bson_path.is_file():
        raise FileNotFoundError(
            f"RegulonDB BSON file not found: {bson_path}"
        )

    links_path.parent.mkdir(parents=True, exist_ok=True)
    operons_path.parent.mkdir(parents=True, exist_ok=True)

    links = set()

    multi_gene_tus = 0
    tus_written = 0

    with open(operons_path, "w") as out_operons:

        out_operons.write(
            "operon_name\ttu_name\tn_genes\tgenes\n"
        )

        with open(bson_path, "rb") as fh:

            for doc in decode_file_iter(fh):

                operon = doc.get("operon", {})
                operon_name = operon.get("name", "")

                for tu in doc.get("transcriptionUnits", []):

                    tu_name = tu.get("name", "")

                    genes = [
                        gene.get("name")
                        for gene in tu.get("genes", [])
                        if gene.get("name")
                    ]

                    if not genes:
                        continue

                    out_operons.write(
                        f"{operon_name}\t"
                        f"{tu_name}\t"
                        f"{len(genes)}\t"
                        f"{','.join(genes)}\n"
                    )

                    tus_written += 1

                    if len(genes) > 1:
                        multi_gene_tus += 1

                    # Consecutive genes within each RegulonDB
                    # transcription unit define reference links.
                    for gene_a, gene_b in zip(
                        genes,
                        genes[1:],
                    ):
                        links.add(
                            (
                                gene_a,
                                gene_b,
                                operon_name,
                                tu_name,
                            )
                        )

    with open(links_path, "w") as out:

        out.write(
            "geneA\tgeneB\toperon_name\ttu_name\n"
        )

        for (
            gene_a,
            gene_b,
            operon_name,
            tu_name,
        ) in sorted(links):

            out.write(
                f"{gene_a}\t"
                f"{gene_b}\t"
                f"{operon_name}\t"
                f"{tu_name}\n"
            )

    print(f"Transcription units written: {tus_written}")
    print(f"Multi-gene transcription units: {multi_gene_tus}")
    print(f"Reference adjacent-gene links: {len(links)}")
    print(f"Transcription units: {operons_path}")
    print(f"Adjacent-gene links: {links_path}")


if __name__ == "__main__":
    main()