import argparse
import pysam
from pathlib import Path
from botas.operons.io import load_genes, compute_gene_coverage, operon_stats, write_operons_gff
from botas.operons.merge import merge_pairs
from botas.operons.features import operon_igds
from botas.operons.classifier import same_operon, operon_score, operon_confidence
from botas.operons.consensus import build_pair_support, merge_by_support


GETOPERONS_EPILOG = """
Examples
--------
1) Infer operons from a single RNA-seq experiment

2) Infer operons across multiple conditions (consensus)

3) Write operons as GFF features (for genome browsers)



Notes
-----
- Genes are linked into operons if they are on the same strand,
  separated by at most --max-igd base pairs, and show compatible
  RNA-seq coverage profiles.
- When multiple BAM files are provided, consensus operons are
  constructed by retaining gene pairs supported in at least
  --min-support fraction of samples.
"""

def _directons(genes):
    by = {}
    for g in genes:
        by.setdefault((g["chrom"], g["strand"]), []).append(g)

    directons = []
    for (chrom, strand), gs in by.items():
        if strand == "+":
            gs.sort(key=lambda x: x["start"])
        else:
            gs.sort(key=lambda x: x["end"], reverse=True)
        directons.append(gs)

    directons.sort(key=lambda d: (d[0]["chrom"], d[0]["start"]))
    return directons


def run_get_operons(args):
    if args.threads < 1:
        raise ValueError("--threads must be at least 1")

    genes = load_genes(args.gff, feature_types=("gene",))

    # -------------------------------------------------
    # Resolve output paths inside the BOTAS results dir
    # -------------------------------------------------
    results_dir = Path(args._resultsdir)

    if args.prefix:
        prefix_name = Path(args.prefix).name
        out_tsv = results_dir / f"{prefix_name}_operons.tsv"
        out_gff = results_dir / f"{prefix_name}_operons.gff"
    else:
        out_tsv = Path(args.out)
        out_gff = out_tsv.with_suffix(out_tsv.suffix + ".gff")

    for bam in args.bam:

        try:
            with pysam.AlignmentFile(bam, "rb") as bf:
                # minimal index check via pileup
                next(
                    bf.pileup(
                        stepper="all",
                        min_base_quality=0,
                        truncate=True,
                        max_depth=1,
                    ),
                    None,
                )
        except Exception as e:
            raise SystemExit(
                f"[ERROR] BAM file is not indexed or unreadable: {bam}\n"
                f"Underlying error: {e}\n"
                f"Please run: samtools index {bam}"
            )
    
    if not genes:
        raise ValueError("No genes loaded from GFF")

    cov_by_bam = [compute_gene_coverage(bam, genes, max_workers=args.threads) for bam in args.bam]

    if args.consensus:
        cov = {gene["id"]: sum(bam_cov.get(gene["id"], 0.0) for bam_cov in cov_by_bam) / len(cov_by_bam)
            for gene in genes
        }
    else:
        cov = cov_by_bam[0]

    directons = _directons(genes)
    operons_all = []

    if not args.consensus:
        # single-BAM behavior (unchanged)
        for ds in directons:
            if len(ds) == 1:
                operons_all.append(ds)
                continue

            labels = []
            for g1, g2 in zip(ds[:-1], ds[1:]):
                labels.append(
                    same_operon(
                        g1, g2, cov,
                        max_igd=args.max_igd,
                        min_coverage=args.min_coverage,
                        min_cov_ratio=args.min_cov_ratio,
                    )
                )

            operons_all.extend(merge_pairs(ds, labels))

    else:
        # multi-BAM consensus
        pair_support = build_pair_support(directons, cov_by_bam, args)

        for ds in directons:
            operons_all.extend(
                merge_by_support(ds, pair_support, args.min_support)
            )

    retained_operons = []

    with open(out_tsv, "w", encoding="utf-8") as out:
        out.write(
            "operon_id\tchrom\tstrand\tstart\tend\tn_genes\tgene_ids\tigds\t"
            "mean_coverage\tmin_coverage\tcoverage_cv\tscore\tconfidence\n"
        )

        for op in operons_all:
            chrom = op[0]["chrom"]
            strand = op[0]["strand"]
            start = min(g["start"] for g in op)
            end = max(g["end"] for g in op)
            ids = ",".join(g["id"] for g in op)

            igds = operon_igds(op)
            mean_cov, min_cov, cv = operon_stats(op, cov)
            score = operon_score(igds, mean_cov, min_cov, cv, args.max_igd)

            if score < args.min_score:
                continue

            retained_operons.append(op)
            operon_id = len(retained_operons)
            conf = operon_confidence(score)

            out.write(
                f"operon_{operon_id}\t{chrom}\t{strand}\t{start}\t{end}\t"
                f"{len(op)}\t{ids}\t{','.join(map(str, igds))}\t"
                f"{mean_cov:.3f}\t{min_cov:.3f}\t{cv:.3f}\t"
                f"{score:.3f}\t{conf}\n"
            )

    if args.write_gff:
        write_operons_gff(
            out_gff,
            retained_operons,
            cov,
            args.max_igd,
        )


def add_operon_args(p: argparse.ArgumentParser) -> None:
    # =========================================================
    # Input / output
    # =========================================================
    io = p.add_argument_group("Input / Output")

    io.add_argument("--bam","-b", required=True, nargs="+", metavar="BAM",
        help=(
            "One or more coordinate-sorted BAM files. "
            "If multiple BAMs are provided, operons can be inferred per-condition "
            "or as a consensus."))

    io.add_argument("--gff","-g", required=True, metavar="GFF",
        help="Genome annotation in GFF3 format (gene features required).")

    io.add_argument("--out","-o", required=False, metavar="TSV",
        help="Output TSV file containing predicted operons.")

    io.add_argument("--threads", "-t", type=int, default=1, metavar="INT",
        help=(
            "Number of worker processes used for gene coverage calculation "
            "(default: 1)."))

    # =========================================================
    # Operon geometry (genomic constraints)
    # =========================================================
    geom = p.add_argument_group("Operon geometry")

    geom.add_argument("--max-igd", type=int, default=70,
        help=(
            "Maximum intergenic distance (bp) between adjacent genes "
            "to allow operon linkage (default: 70)."))


    # =========================================================
    # Expression / coverage constraints
    # =========================================================
    cov = p.add_argument_group("Coverage constraints")

    cov.add_argument("--min-coverage", type=float, default=0.5, metavar="X",
        help=(
            "Minimum mean per-base coverage required for each gene "
            "to be considered expressed (default: 0.5)."))

    cov.add_argument("--min-cov-ratio", type=float, default=0.25, metavar="R",
        help=(
            "Minimum coverage ratio between adjacent genes "
            "(min/ max) to support co-transcription (default: 0.25)."))

    # =========================================================
    # Multi-BAM consensus inference
    # =========================================================
    consensus = p.add_argument_group("Consensus inference")

    consensus.add_argument("--consensus", action="store_true",
        help=(
            "Infer consensus operons across multiple BAM files "
            "instead of reporting per-BAM operons."))

    consensus.add_argument("--min-support", type=float, default=0.6, metavar="F",
        help=(
            "Minimum fraction of BAM files supporting a gene pair "
            "for inclusion in the consensus operon (default: 0.6)."))

    # =========================================================
    # Scoring and filtering
    # =========================================================
    score = p.add_argument_group("Scoring and filtering")

    score.add_argument("--min-score", type=float, default=0.0, metavar="S",
        help=(
            "Minimum operon confidence score required to report "
            "(default: 0.0)."))

    # =========================================================
    # Output options
    # =========================================================
    out = p.add_argument_group("Output options")

    out.add_argument("--write-gff", action="store_true", help=(
            "Write predicted operons as GFF features "
            "in addition to the TSV output."))

    out.add_argument("--prefix", default=None,
        help=(
            "Optional output prefix. "
            "If provided, multiple output files may be generated "
            "(e.g. prefix.operons.tsv, prefix.operons.gff)."))
