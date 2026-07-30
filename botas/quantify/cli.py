# botas/quantify/cli.py
from __future__ import annotations

import argparse
import logging

from .core import quantify_bam
from .io import compute_rpkm, compute_tpm, write_gene_counts_tsv, write_operon_counts_tsv

QUANTIFY_EPILOG = """
Examples
--------
Gene-level quantification:

  botas quantify --bam sample.bam --gff genome.gff --feature-type gene --id-attribute locus_tag --out gene_counts.tsv

Operon-level quantification:

  botas quantify --bam sample.bam --gff operons.gff --feature-type operon --id-attribute ID --out operon_counts.tsv

Notes
-----
- BAM input must be coordinate-sorted paired-end data.
- The default feature type is gene.
- The default identifier is locus_tag for genes and ID for operons.
"""


def add_quantify_args(parser: argparse.ArgumentParser) -> None:
    """
    Add arguments for the `botas quantify` command.
    """

    parser.add_argument(
        "-b",
        "--bam",
        required=True,
        metavar="BAM",
        help="Coordinate-sorted paired-end BAM file.",
    )

    parser.add_argument(
        "-g",
        "--gff",
        required=True,
        metavar="GFF",
        help="GFF3 annotation file.",
    )

    parser.add_argument(
        "-o",
        "--out",
        required=True,
        metavar="TSV",
        help="Output quantification table.",
    )

    parser.add_argument(
        "--feature-type",
        choices=["gene", "operon"],
        default="gene",
        help="GFF feature type to quantify (default: gene).",
    )

    parser.add_argument(
        "--id-attribute",
        default=None,
        metavar="ATTR",
        help=(
            "GFF attribute used as the feature identifier. "
            "Defaults to locus_tag for genes and ID for operons."
        ),
    )

    parser.add_argument(
        "--min-overlap",
        type=int,
        default=1,
        metavar="N",
        help="Minimum number of overlapping bases required (default: 1).",
    )

    parser.add_argument(
        "--largest-overlap",
        action="store_true",
        help=(
            "Assign an ambiguous fragment to the feature with the largest "
            "overlap. Equal largest overlaps remain ambiguous."
        ),
    )

    parser.add_argument(
        "--log",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging level (default: INFO).",
    )

    parser.set_defaults(func=_run_quantify)


def _run_quantify(args: argparse.Namespace) -> int:
    logging.basicConfig(
        level=getattr(logging, args.log),
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    logger = logging.getLogger("botas.quantify")

    id_attribute = args.id_attribute

    if id_attribute is None:
        if args.feature_type == "operon":
            id_attribute = "ID"
        else:
            id_attribute = "locus_tag"

    logger.info("BAM: %s", args.bam)
    logger.info("GFF: %s", args.gff)
    logger.info("Feature type: %s", args.feature_type)
    logger.info("Identifier attribute: %s", id_attribute)

    result = quantify_bam(
        bam_path=args.bam,
        gff_path=args.gff,
        feature_types=args.feature_type,
        id_attribute=id_attribute,
        min_overlap_bases=args.min_overlap,
        largest_overlap=args.largest_overlap,
    )

    counts = result.gene_counts
    features = result.genes

    lengths = {
        feature_id: feature.length
        for feature_id, feature in features.items()
    }

    total_assigned = float(result.total_assigned_counts)

    rpkm = compute_rpkm(
        counts=counts,
        lengths_bp=lengths,
        total_mapped=total_assigned,
    )

    tpm = compute_tpm(
        counts=counts,
        lengths_bp=lengths,
    )

    if args.feature_type == "gene":
        write_gene_counts_tsv(
            out_tsv=args.out,
            counts=counts,
            lengths_bp=lengths,
            rpkm=rpkm,
            tpm=tpm,
        )

    else:
        operon_metadata = {
            operon_id: {
                "contig": feature.contig,
                "start": feature.start + 1,
                "end": feature.end,
                "strand": feature.strand,
                "genes": "",
            }
            for operon_id, feature in features.items()
        }

        write_operon_counts_tsv(
            out_tsv=args.out,
            op_meta=operon_metadata,
            counts=counts,
            lengths_bp=lengths,
            rpkm=rpkm,
            tpm=tpm,
        )

    summary = result.summary

    logger.info("Output written to: %s", args.out)
    logger.info("Total fragments: %d", summary.total_fragments)
    logger.info("Assigned: %d", summary.assigned)
    logger.info("No feature: %d", summary.no_feature)
    logger.info("Ambiguous: %d", summary.ambiguous)
    logger.info(
        "Insufficient overlap: %d",
        summary.insufficient_overlap,
    )
    logger.info("Unmapped: %d", summary.unmapped)
    logger.info(
        "Assignment rate: %.2f%%",
        summary.assignment_rate * 100.0,
    )

    return 0