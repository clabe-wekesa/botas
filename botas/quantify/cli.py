# botas/quantify/cli.py
from __future__ import annotations

import argparse
import logging
from pathlib import Path

from .core import quantify_bam
from .io import (compute_rpkm, compute_tpm, write_gene_counts_tsv, write_gene_expression_matrix_tsv,
    write_operon_counts_tsv, write_operon_expression_matrix_tsv)

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

    parser.add_argument("-b", "--bam", nargs="+", required=True, metavar="BAM",
        help="One or more coordinate-sorted paired-end BAM files.")

    parser.add_argument("-g", "--gff", required=True, metavar="GFF", help="GFF3 annotation file.")

    parser.add_argument("-o", "--out", required=False, metavar="TSV", help="Output quantification table.")

    parser.add_argument("--feature-type", choices=["gene", "operon"], default="gene",
        help="GFF feature type to quantify (default: gene).")

    parser.add_argument("--id-attribute", default=None, metavar="ATTR",
        help=(
            "GFF attribute used as the feature identifier. "
            "Defaults to locus_tag for genes and ID for operons."
        ),
    )

    parser.add_argument("--min-overlap", type=int, default=1, metavar="N",
        help="Minimum number of overlapping bases required (default: 1).",
    )

    parser.add_argument("--largest-overlap", action="store_true",
        help=(
            "Assign an ambiguous fragment to the feature with the largest "
            "overlap. Equal largest overlaps remain ambiguous."
        ),
    )

    parser.add_argument("--log", choices=["DEBUG", "INFO", "WARNING", "ERROR"], default="INFO",
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

    logger.info("BAM files: %d", len(args.bam))
    logger.info("GFF: %s", args.gff)
    logger.info("Feature type: %s", args.feature_type)
    logger.info("Identifier attribute: %s", id_attribute)

    counts_by_sample = {}
    rpkm_by_sample = {}
    tpm_by_sample = {}
    lengths = {}
    features = {}

    for bam_path in args.bam:
        sample_name = Path(bam_path).name

        if sample_name.endswith(".sorted.bam"):
            sample_name = sample_name[:-len(".sorted.bam")]
        elif sample_name.endswith(".bam"):
            sample_name = sample_name[:-len(".bam")]

        if sample_name in counts_by_sample:
            raise ValueError(
                f"Duplicate sample name derived from BAM files: {sample_name}"
            )

        logger.info("Quantifying sample: %s", sample_name)
        logger.info("BAM: %s", bam_path)

        result = quantify_bam(
            bam_path=bam_path,
            gff_path=args.gff,
            feature_types=args.feature_type,
            id_attribute=id_attribute,
            min_overlap_bases=args.min_overlap,
            largest_overlap=args.largest_overlap,
        )

        sample_counts = result.gene_counts
        sample_features = result.genes

        sample_lengths = {
            feature_id: feature.length
            for feature_id, feature in sample_features.items()
        }

        if not lengths:
            lengths = sample_lengths
            features = sample_features
        elif sample_lengths != lengths:
            raise ValueError(
                "Feature identifiers or lengths differ between BAM files"
            )

        total_assigned = float(result.total_assigned_counts)

        counts_by_sample[sample_name] = sample_counts
        rpkm_by_sample[sample_name] = compute_rpkm(
            counts=sample_counts,
            lengths_bp=lengths,
            total_mapped=total_assigned,
        )
        tpm_by_sample[sample_name] = compute_tpm(
            counts=sample_counts,
            lengths_bp=lengths,
        )

    # Single BAM: preserve the existing output format
    if len(args.bam) == 1:
        sample_name = next(iter(counts_by_sample))

        if args.feature_type == "gene":
            write_gene_counts_tsv(
                out_tsv=args.out,
                counts=counts_by_sample[sample_name],
                lengths_bp=lengths,
                rpkm=rpkm_by_sample[sample_name],
                tpm=tpm_by_sample[sample_name],
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
                counts=counts_by_sample[sample_name],
                lengths_bp=lengths,
                rpkm=rpkm_by_sample[sample_name],
                tpm=tpm_by_sample[sample_name],
            )

    # Multiple BAMs: write one expression matrix
    elif args.feature_type == "gene":
        write_gene_expression_matrix_tsv(
            out_tsv=args.out,
            counts_by_sample=counts_by_sample,
            lengths_bp=lengths,
            rpkm_by_sample=rpkm_by_sample,
            tpm_by_sample=tpm_by_sample,
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

        write_operon_expression_matrix_tsv(
            out_tsv=args.out,
            op_meta=operon_metadata,
            counts_by_sample=counts_by_sample,
            lengths_bp=lengths,
            rpkm_by_sample=rpkm_by_sample,
            tpm_by_sample=tpm_by_sample,
        )

    logger.info("Output: %s", args.out)

    return 0