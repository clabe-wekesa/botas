# botas/quantify/cli.py
from __future__ import annotations

import argparse
import logging
from dataclasses import replace
from .genes import load_gene_features_from_gff
from .core import QuantParams, quantify_genes_single_bam
from .io import compute_rpkm, compute_tpm, write_gene_counts_tsv, write_operon_counts_tsv
from .infer_operons import infer_operons
from .operons import load_operons_tsv, aggregate_operon_counts, aggregate_operon_lengths, compute_operon_metadata
from .core import QuantParams, quantify_genes_single_bam, quantify_genes_bams


QUANTIFY_EPILOG = """
Examples
--------
1) Gene-level quantification (single sample)

  botas quantify --bam sample.bam --gff genome.gff --gff-gene-attribute locus_tag --strand 1 --level gene -o gene_counts.tsv

2) Operon-level quantification using inferred operons

  botas quantify --bam sample.bam --gff genome.gff --gff-gene-attribute ID --level operon --operons sample.operons.tsv --operon-gene-column genes -o operon_counts.tsv

3) Multi-condition operon expression matrix

  botas quantify --bam ctrl.bam stress.bam --gff genome.gff --gff-gene-attribute locus_tag --level operon --operons ctrl.operons.tsv --operon-gene-column genes -o operon_matrix.tsv

Notes
-----
- Gene identifiers are taken from the specified GFF attribute.
- Multi-mapping reads are handled using the NH tag.
- Operon-level counts are computed by aggregating gene counts.
"""

def add_quantify_args(p: argparse.ArgumentParser) -> None:  
    # =========================================================
    # Input / output
    # =========================================================
    io = p.add_argument_group("Input / Output")
    io.add_argument("-b","--bam", required=True, nargs="+", metavar="BAM",
        help=(
            "One or more coordinate-sorted BAM files produced by `botas align`. "
            "If multiple BAMs are provided, a count matrix is produced."))
    io.add_argument("-g","--gff", required=True, metavar="GFF", help="Genome annotation in GFF3 format (gene features required).")
    io.add_argument("-o","--out", required=True, metavar="TSV",
        help=(
            "Output TSV file. For multiple BAMs, this will be a count matrix "
            "(features × samples)."))
    # =========================================================
    # Feature definition
    # =========================================================
    feat = p.add_argument_group("Feature definition")
    feat.add_argument("--level", choices=["gene", "operon"], default="gene", 
        help="Quantification level (default: gene).")
    feat.add_argument("--gff-gene-attribute", required=True, metavar="ATTR",
        help=(
            "GFF attribute used as gene identifier "
            "(e.g. ID, locus_tag, gene_id)."))
    # =========================================================
    # Operon handling
    # =========================================================
    oper = p.add_argument_group("Operon handling")
    oper.add_argument("--operons", metavar="TSV", help="Operon definition TSV produced by `botas getOperons`.")

    oper.add_argument("--infer-operons", action="store_true",
        help=(
            "Infer operons on-the-fly if no operon TSV is provided. "
            "Uses default `getOperons` parameters."))

    oper.add_argument("--operon-gene-column", default=None, metavar="COL",
        help=(
            "Column name in operon TSV containing gene IDs "
            "(required when --operons is used)."))

    # =========================================================
    # Read assignment rules
    # =========================================================
    assign = p.add_argument_group("Read assignment")

    assign.add_argument("--strand", type=int, choices=[0, 1, 2], default=0, metavar="MODE",
        help=(
            "Strand specificity:\n"
            "  0 = unstranded\n"
            "  1 = stranded (same strand)\n"
            "  2 = reverse-stranded\n"
            "(default: 0)"))

    assign.add_argument("--mapq", type=int, default=0, metavar="Q", help="Minimum MAPQ required for read assignment (default: 0).")

    assign.add_argument("--multi", choices=["ignore", "unique", "fractional"], default="ignore",
        help=(
            "Handling of multi-mapping reads based on NH tag:\n"
            "  ignore     = discard multi-mappers\n"
            "  unique     = keep reads with NH=1 only\n"
            "  fractional = distribute 1/N across targets\n"
            "(default: ignore)"))

    # =========================================================
    # Logging
    # =========================================================
    p.add_argument("--log", default="INFO", choices=["DEBUG", "INFO", "WARN", "ERROR"], help="Log verbosity level.")

    p.set_defaults(func=_run_quantify)


def _run_quantify(args: argparse.Namespace) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    logger = logging.getLogger("botas.quantify")
    logger.info("Loading genes from GFF: %s", args.gff)
    genes = load_gene_features_from_gff(args.gff, feature_type="gene", id_from=args.gff_gene_attribute)

    # ---- HARD VALIDATION: gene features must exist ----
    if not genes:
        raise ValueError(
            "No gene features were loaded from the GFF.\n"
            "BOTAS quantify requires GFF 'gene' features.\n"
            "Check that:\n"
            "  1) The GFF contains feature type 'gene'\n"
            "  2) --gff-gene-attribute matches the gene IDs used in operons\n"
        )

    params = QuantParams(strand_mode=args.strand, mapq_min=args.mapq, multi_mode=args.multi, use_nh_tag=True)
    logger.info("Quantifying genes from BAM: %s", args.bam)
    bam_paths = args.bam if isinstance(args.bam, list) else [args.bam]

    if len(bam_paths) == 1:
        counts = quantify_genes_single_bam(
            bam_paths[0],
            genes,
            params,
        )
    else:
        counts_by_sample = quantify_genes_bams(
            bam_paths,
            genes,
            params,
        )
    lengths = {gid: gf.length for gid, gf in genes.items()}
    total_mapped = sum(counts.values())
    rpkm = compute_rpkm(counts, lengths, total_mapped=total_mapped)
    tpm = compute_tpm(counts, lengths)

    # -------------------------------------------------
    # GENE-LEVEL OUTPUT (default)
    # -------------------------------------------------
    if args.level == "gene":
        logger.info("Writing gene-level counts: %s", args.out)
        write_gene_counts_tsv(args.out, counts, lengths, rpkm, tpm)
        logger.info(
            "Done (gene-level). Total assigned (weighted) counts: %.3f",
            total_mapped,
        )
        return 0


    # -------------------------------------------------
    # OPERON-LEVEL OUTPUT
    # -------------------------------------------------
    logger.info("Operon-level quantification requested")

    if args.operons and args.infer_operons:
        raise ValueError(
            "Use either --operons <file> OR --infer-operons, not both."
        )

    if args.operons:
        logger.info("Loading operons from file: %s", args.operons)
        operon_tsv = args.operons

        if not args.operon_gene_column:
            raise ValueError("--operon-gene-column is required when using --operons <file>")
        gene_col = args.operon_gene_column

    elif args.infer_operons:
        logger.info("Inferring operons from BAM + GFF")
        operon_tsv = infer_operons(bam=args.bam, gff=args.gff, min_score=0.5)

        # BOTAS controls the schema of getOperons output
        gene_col = "gene_id"

    else:
        raise ValueError(
            "Operon-level quantification requested but no operons provided.\n"
            "Use either:\n"
            "  --operons <operons.tsv>\n"
            "  --infer-operons"
        )

    # 2) Load operons and aggregate
    operons = load_operons_tsv(operon_tsv, gene_col=gene_col)

    # ---- HARD VALIDATION: operon genes must exist in GFF gene set ----
    missing = set()

    for op_id, gids in operons.items():
        for g in gids:
            if g not in genes:
                missing.add(g)

    if missing:
        raise ValueError(
            f"{len(missing)} operon genes were not found in the GFF gene features.\n"
            f"Example missing genes: {list(missing)[:5]}\n"
            "Ensure --gff-gene-attribute matches the operon gene IDs."
        )


    # compute operon genomic metadata from GFF-derived genes
    op_meta = compute_operon_metadata(operons, genes)

    op_counts = aggregate_operon_counts(counts, operons)
    op_lengths = aggregate_operon_lengths(lengths, operons)    
    op_total = total_mapped
    op_rpkm = compute_rpkm(op_counts, op_lengths, total_mapped=op_total)
    op_tpm = compute_tpm(op_counts, op_lengths)

    # 3) Write operon output
    logger.info("Writing operon-level counts: %s", args.out)
    write_operon_counts_tsv(args.out, op_meta, op_counts, op_lengths, op_rpkm, op_tpm)
    logger.info("Done (operon-level). Total assigned (weighted) counts: %.3f", op_total)
    return 0

