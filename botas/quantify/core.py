# botas/quantify/core.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple
import bisect
import logging
import time
import pysam

try:
    from botas.quantify._core_fast import quantify_bam_fast

    CYTHON_BACKEND_AVAILABLE = True

except ImportError:
    quantify_bam_fast = None
    CYTHON_BACKEND_AVAILABLE = False

from botas.quantify.genes import GeneFeature


logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------
# Type aliases
# ---------------------------------------------------------------------

GeneCounts = Dict[str, float]
MultiBamCounts = Dict[str, GeneCounts]

# start, end, gene_id, strand
IndexedInterval = Tuple[int, int, str, str]


# ---------------------------------------------------------------------
# Quantification parameters and statistics
# ---------------------------------------------------------------------

@dataclass(frozen=True)
class QuantParams:
    """
    Parameters controlling alignment filtering and gene assignment.

    strand_mode
        0 = unstranded
        1 = alignment strand must match gene strand
        2 = alignment strand must be opposite to gene strand

    multi_mode
        ignore:
            Discard alignments explicitly marked as multimapping by NH > 1.
            If NH is absent, treat the alignment as countable.

        unique:
            Count alignments with NH == 1. If NH is absent, treat the
            primary alignment as unique.

        fractional:
            Count alignments with weight 1 / NH. If NH is absent, use 1.
    """

    strand_mode: int = 0
    mapq_min: int = 0
    multi_mode: str = "ignore"
    use_nh_tag: bool = True
    log_every: int = 1_000_000

    def __post_init__(self) -> None:
        if self.strand_mode not in {0, 1, 2}:
            raise ValueError("strand_mode must be 0, 1, or 2")

        if self.mapq_min < 0:
            raise ValueError("mapq_min must be >= 0")

        if self.multi_mode not in {"ignore", "unique", "fractional"}:
            raise ValueError(
                "multi_mode must be one of: ignore, unique, fractional"
            )

        if self.log_every < 0:
            raise ValueError("log_every must be >= 0")


@dataclass
class QuantStats:
    """
    Summary statistics for one quantified BAM file.
    """

    bam_path: str
    total_records: int = 0
    primary_mapped_records: int = 0
    below_mapq: int = 0
    multimapping_discarded: int = 0
    no_annotated_contig: int = 0
    no_feature: int = 0
    ambiguous: int = 0
    assigned_records: int = 0
    assigned_weight: float = 0.0
    elapsed_seconds: float = 0.0


# ---------------------------------------------------------------------
# Gene interval index
# ---------------------------------------------------------------------

@dataclass(frozen=True)
class ContigIntervalIndex:
    """
    Interval index for one reference contig.

    intervals
        Gene intervals sorted by start coordinate.

    starts
        Start coordinates corresponding to intervals.

    prefix_max_ends
        prefix_max_ends[i] is the greatest interval end among
        intervals[0:i+1]. This allows the backward overlap search to stop
        as soon as no earlier interval can overlap the query.
    """

    intervals: Tuple[IndexedInterval, ...]
    starts: Tuple[int, ...]
    prefix_max_ends: Tuple[int, ...]


GeneIntervalIndex = Dict[str, ContigIntervalIndex]


def build_gene_interval_index(
    genes: Mapping[str, GeneFeature],
) -> GeneIntervalIndex:
    """
    Build an interval index from gene features.

    Each block belonging to a gene is inserted separately. Intervals are
    sorted by genomic start coordinate. Prefix maximum ends permit
    correct overlap lookup without relying on an arbitrary search window.
    """

    intervals_by_contig: Dict[str, List[IndexedInterval]] = {}

    for gene_id, feature in genes.items():
        for start, end in feature.blocks:
            if end <= start:
                logger.warning(
                    "Skipping invalid block for gene %s: [%d, %d)",
                    gene_id,
                    start,
                    end,
                )
                continue

            intervals_by_contig.setdefault(feature.contig, []).append(
                (start, end, gene_id, feature.strand)
            )

    index: GeneIntervalIndex = {}

    for contig, intervals in intervals_by_contig.items():
        intervals.sort(key=lambda item: (item[0], item[1], item[2]))

        starts: List[int] = []
        prefix_max_ends: List[int] = []
        maximum_end = -1

        for start, end, _gene_id, _strand in intervals:
            starts.append(start)
            maximum_end = max(maximum_end, end)
            prefix_max_ends.append(maximum_end)

        index[contig] = ContigIntervalIndex(
            intervals=tuple(intervals),
            starts=tuple(starts),
            prefix_max_ends=tuple(prefix_max_ends),
        )

    logger.debug(
        "Built gene interval index for %d contigs and %d genes",
        len(index),
        len(genes),
    )

    return index


# ---------------------------------------------------------------------
# Alignment filtering and weighting
# ---------------------------------------------------------------------

def _read_passes_filters(
    aln: pysam.AlignedSegment,
    mapq_min: int,
) -> bool:
    """
    Return True for mapped primary alignments meeting the MAPQ threshold.
    """

    if aln.is_unmapped:
        return False

    if aln.is_secondary or aln.is_supplementary:
        return False

    if aln.mapping_quality < mapq_min:
        return False

    return True


def _nh(aln: pysam.AlignedSegment) -> Optional[int]:
    """
    Return the NH tag without using exceptions for normally absent tags.
    """

    if not aln.has_tag("NH"):
        return None

    value = aln.get_tag("NH")

    try:
        return int(value)
    except (TypeError, ValueError):
        logger.debug(
            "Invalid NH tag for alignment %s: %r",
            aln.query_name,
            value,
        )
        return None


def _read_weight(
    aln: pysam.AlignedSegment,
    params: QuantParams,
) -> float:
    """
    Determine the weight contributed by one primary alignment.
    """

    if not params.use_nh_tag:
        return 1.0

    nh = _nh(aln)

    if params.multi_mode == "ignore":
        if nh is not None and nh > 1:
            return 0.0
        return 1.0

    if params.multi_mode == "unique":
        # Some aligners omit NH for primary unique alignments.
        if nh is None:
            return 1.0
        return 1.0 if nh == 1 else 0.0

    if params.multi_mode == "fractional":
        if nh is None or nh <= 1:
            return 1.0
        return 1.0 / float(nh)

    # QuantParams validates this before execution.
    raise ValueError(f"Unknown multi_mode: {params.multi_mode}")


def _strand_compatible(
    aln: pysam.AlignedSegment,
    gene_strand: str,
    strand_mode: int,
) -> bool:
    """
    Test alignment/gene strand compatibility.
    """

    if strand_mode == 0:
        return True

    read_strand = "-" if aln.is_reverse else "+"

    if strand_mode == 1:
        return read_strand == gene_strand

    if strand_mode == 2:
        return read_strand != gene_strand

    raise ValueError("strand_mode must be 0, 1, or 2")


# ---------------------------------------------------------------------
# Interval querying
# ---------------------------------------------------------------------

def _genes_overlapping_block(
    block_start: int,
    block_end: int,
    contig_index: ContigIntervalIndex,
    aln: pysam.AlignedSegment,
    strand_mode: int,
    hits: set[str],
) -> bool:
    """
    Add genes overlapping one alignment block to ``hits``.

    Returns True as soon as more than one distinct gene has been found,
    allowing ambiguous alignments to terminate early.
    """

    if block_end <= block_start:
        return False

    intervals = contig_index.intervals
    starts = contig_index.starts
    prefix_max_ends = contig_index.prefix_max_ends

    # Only intervals beginning before block_end can overlap.
    position = bisect.bisect_left(starts, block_end) - 1

    while position >= 0:
        # If the maximum end among this and all previous intervals is not
        # beyond block_start, no earlier interval can overlap.
        if prefix_max_ends[position] <= block_start:
            break

        start, end, gene_id, gene_strand = intervals[position]

        if end > block_start and start < block_end:
            if _strand_compatible(aln, gene_strand, strand_mode):
                hits.add(gene_id)

                if len(hits) > 1:
                    return True

        position -= 1

    return False


def _assign_alignment_to_gene(
    aln: pysam.AlignedSegment,
    contig_index: ContigIntervalIndex,
    strand_mode: int,
) -> Tuple[Optional[str], bool]:
    """
    Assign one alignment to a gene.

    The CIGAR-derived aligned blocks are checked independently. This
    avoids treating skipped reference regions as covered by the read.

    Returns
    -------
    gene_id
        Gene ID when exactly one gene is overlapped, otherwise None.

    ambiguous
        True when more than one distinct gene is overlapped.
    """

    hits: set[str] = set()

    blocks = aln.get_blocks()

    if not blocks:
        return None, False

    for block_start, block_end in blocks:
        ambiguous = _genes_overlapping_block(
            block_start=block_start,
            block_end=block_end,
            contig_index=contig_index,
            aln=aln,
            strand_mode=strand_mode,
            hits=hits,
        )

        if ambiguous:
            return None, True

    if len(hits) == 1:
        return next(iter(hits)), False

    return None, False


# ---------------------------------------------------------------------
# Single-BAM quantification
# ---------------------------------------------------------------------

def _build_tid_index(
    bam: pysam.AlignmentFile,
    gene_index: GeneIntervalIndex,
) -> Dict[int, ContigIntervalIndex]:
    """
    Convert GFF contig names to BAM numeric reference IDs once.
    """

    index_by_tid: Dict[int, ContigIntervalIndex] = {}

    for contig, contig_index in gene_index.items():
        tid = bam.get_tid(contig)

        if tid < 0:
            logger.warning(
                "GFF contig %s is absent from BAM header",
                contig,
            )
            continue

        index_by_tid[tid] = contig_index

    return index_by_tid


def _quantify_one_bam(
    bam_path: str,
    genes: Mapping[str, GeneFeature],
    params: QuantParams,
    gene_index: GeneIntervalIndex,
) -> Tuple[GeneCounts, QuantStats]:
    """
    Internal implementation for one BAM using a pre-built gene index.

    The compiled Cython backend is used when available. Otherwise, the
    original Python implementation is used as a fallback.
    """

    counts: GeneCounts = dict.fromkeys(genes, 0.0)
    stats = QuantStats(bam_path=str(bam_path))
    start_time = time.perf_counter()

    # ---------------------------------------------------------------
    # Compiled Cython implementation
    # ---------------------------------------------------------------

    if CYTHON_BACKEND_AVAILABLE:
        multi_mode_code = {
            "ignore": 0,
            "unique": 1,
            "fractional": 2,
        }[params.multi_mode]

        counts, raw_stats = quantify_bam_fast(
            str(bam_path),
            tuple(genes),
            gene_index,
            params.strand_mode,
            params.mapq_min,
            multi_mode_code,
            params.use_nh_tag,
        )

        (
            stats.total_records,
            stats.primary_mapped_records,
            stats.below_mapq,
            stats.multimapping_discarded,
            stats.no_annotated_contig,
            stats.no_feature,
            stats.ambiguous,
            stats.assigned_records,
            stats.assigned_weight,
        ) = raw_stats

        stats.elapsed_seconds = time.perf_counter() - start_time

        logger.info(
            "Quantification finished for %s: "
            "records=%d, primary_mapped=%d, assigned=%d, "
            "assigned_weight=%.3f, ambiguous=%d, no_feature=%d, "
            "elapsed=%.3f s",
            bam_path,
            stats.total_records,
            stats.primary_mapped_records,
            stats.assigned_records,
            stats.assigned_weight,
            stats.ambiguous,
            stats.no_feature,
            stats.elapsed_seconds,
        )

        return counts, stats

    # ---------------------------------------------------------------
    # Pure-Python fallback
    # ---------------------------------------------------------------

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        gene_index_by_tid = _build_tid_index(
            bam=bam,
            gene_index=gene_index,
        )

        for aln in bam.fetch(until_eof=True):
            stats.total_records += 1

            if aln.is_unmapped:
                continue

            if aln.is_secondary:
                continue

            if aln.is_supplementary:
                continue

            stats.primary_mapped_records += 1

            if aln.mapping_quality < params.mapq_min:
                stats.below_mapq += 1
                continue

            contig_index = gene_index_by_tid.get(
                aln.reference_id
            )

            if contig_index is None:
                stats.no_annotated_contig += 1
                continue

            weight = _read_weight(
                aln=aln,
                params=params,
            )

            if weight == 0.0:
                stats.multimapping_discarded += 1
                continue

            gene_id, ambiguous = _assign_alignment_to_gene(
                aln=aln,
                contig_index=contig_index,
                strand_mode=params.strand_mode,
            )

            if ambiguous:
                stats.ambiguous += 1
                continue

            if gene_id is None:
                stats.no_feature += 1
                continue

            counts[gene_id] += weight
            stats.assigned_records += 1
            stats.assigned_weight += weight

            if (
                params.log_every > 0
                and logger.isEnabledFor(logging.INFO)
                and stats.total_records % params.log_every == 0
            ):
                elapsed = time.perf_counter() - start_time

                logger.info(
                    "Quantify %s: records=%d, assigned=%d, "
                    "elapsed=%.2f s",
                    Path(bam_path).name,
                    stats.total_records,
                    stats.assigned_records,
                    elapsed,
                )

    stats.elapsed_seconds = time.perf_counter() - start_time

    logger.info(
        "Quantification finished for %s: "
        "records=%d, primary_mapped=%d, assigned=%d, "
        "assigned_weight=%.3f, ambiguous=%d, no_feature=%d, "
        "elapsed=%.3f s",
        bam_path,
        stats.total_records,
        stats.primary_mapped_records,
        stats.assigned_records,
        stats.assigned_weight,
        stats.ambiguous,
        stats.no_feature,
        stats.elapsed_seconds,
    )

    return counts, stats


def quantify_genes_single_bam(
    bam_path: str,
    genes: Dict[str, GeneFeature],
    p: QuantParams,
) -> GeneCounts:
    """
    Quantify genes from one BAM file.

    This function preserves the original BOTAS public interface and
    returns only the count dictionary.
    """

    gene_index = build_gene_interval_index(genes)

    counts, _stats = _quantify_one_bam(
        bam_path=bam_path,
        genes=genes,
        params=p,
        gene_index=gene_index,
    )

    return counts


def quantify_genes_single_bam_with_stats(
    bam_path: str,
    genes: Dict[str, GeneFeature],
    p: QuantParams,
) -> Tuple[GeneCounts, QuantStats]:
    """
    Quantify one BAM and also return assignment statistics.
    """

    gene_index = build_gene_interval_index(genes)

    return _quantify_one_bam(
        bam_path=bam_path,
        genes=genes,
        params=p,
        gene_index=gene_index,
    )


# ---------------------------------------------------------------------
# Multiple-BAM quantification
# ---------------------------------------------------------------------

def _sample_name_from_bam(
    bam_path: str,
    existing_names: set[str],
) -> str:
    """
    Generate a unique sample name from a BAM filename.
    """

    path = Path(bam_path)

    name = path.name

    if name.lower().endswith(".bam"):
        name = name[:-4]

    base_name = name
    suffix = 2

    while name in existing_names:
        name = f"{base_name}_{suffix}"
        suffix += 1

    existing_names.add(name)
    return name


def quantify_genes_bams(
    bam_paths: Sequence[str],
    genes: Dict[str, GeneFeature],
    p: QuantParams,
) -> MultiBamCounts:
    """
    Quantify multiple BAM files.

    Returns
    -------
    dict
        Mapping:

            sample_name -> {gene_id: count}

    The gene interval index is built once and reused for every BAM.
    """

    counts_by_sample, _stats_by_sample = quantify_genes_bams_with_stats(
        bam_paths=bam_paths,
        genes=genes,
        p=p,
    )

    return counts_by_sample


def quantify_genes_bams_with_stats(
    bam_paths: Sequence[str],
    genes: Dict[str, GeneFeature],
    p: QuantParams,
) -> Tuple[MultiBamCounts, Dict[str, QuantStats]]:
    """
    Quantify multiple BAM files and return per-sample statistics.
    """

    if not bam_paths:
        raise ValueError("At least one BAM file is required")

    gene_index = build_gene_interval_index(genes)

    counts_by_sample: MultiBamCounts = {}
    stats_by_sample: Dict[str, QuantStats] = {}
    used_sample_names: set[str] = set()

    for bam_path in bam_paths:
        sample_name = _sample_name_from_bam(
            bam_path,
            used_sample_names,
        )

        logger.info(
            "Quantifying sample %s from BAM: %s",
            sample_name,
            bam_path,
        )

        counts, stats = _quantify_one_bam(
            bam_path=bam_path,
            genes=genes,
            params=p,
            gene_index=gene_index,
        )

        counts_by_sample[sample_name] = counts
        stats_by_sample[sample_name] = stats

    return counts_by_sample, stats_by_sample