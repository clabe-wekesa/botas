from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Mapping, Tuple
import bisect

from botas.quantify.fragments import Fragment
from botas.quantify.genes import GeneFeature


Interval = Tuple[int, int]


@dataclass(frozen=True, slots=True)
class ContigIntervalIndex:
    """
    Gene intervals for one contig, sorted by genomic start.
    """

    intervals: Tuple[Tuple[int, int, str], ...]
    starts: Tuple[int, ...]
    prefix_max_ends: Tuple[int, ...]


@dataclass(frozen=True, slots=True)
class GeneIntervalIndex:
    """
    Searchable genomic index together with the original gene features.
    """

    contigs: Mapping[str, ContigIntervalIndex]
    genes: Mapping[str, GeneFeature]

@dataclass(frozen=True, slots=True)
class FragmentGeneOverlap:
    """
    Per-gene overlap evidence from the two mates of a fragment.
    """

    gene_id: str
    read1_overlap: int
    read2_overlap: int

    @property
    def total_overlap(self) -> int:
        """
        featureCounts-style overlap total.

        Overlap bases from read 1 and read 2 are summed. Bases covered by
        both mates may therefore contribute once for each read.
        """

        return self.read1_overlap + self.read2_overlap

    @property
    def supporting_reads(self) -> int:
        return int(self.read1_overlap > 0) + int(self.read2_overlap > 0)

def normalized_blocks_gene_overlaps(
    blocks: Iterable[Interval],
    contig_index: ContigIntervalIndex,
) -> Dict[str, int]:
    """
    Calculate overlap lengths for sorted, non-overlapping BAM blocks.
    """

    intervals = contig_index.intervals
    starts = contig_index.starts
    prefix_max_ends = contig_index.prefix_max_ends
    bisect_left = bisect.bisect_left

    overlaps: Dict[str, int] = {}
    overlaps_get = overlaps.get

    for block_start, block_end in blocks:
        position = bisect_left(starts, block_end) - 1

        while position >= 0:
            if prefix_max_ends[position] <= block_start:
                break

            gene_start, gene_end, gene_id = intervals[position]

            # bisect_left(starts, block_end) guarantees gene_start < block_end.
            if gene_end > block_start:
                if block_start > gene_start:
                    overlap_start = block_start
                else:
                    overlap_start = gene_start

                if block_end < gene_end:
                    overlap_end = block_end
                else:
                    overlap_end = gene_end

                overlap_length = overlap_end - overlap_start

                if overlap_length > 0:
                    overlaps[gene_id] = (
                        overlaps_get(gene_id, 0)
                        + overlap_length
                    )

            position -= 1

    return overlaps

def blocks_gene_overlaps(
    *,
    contig: str,
    blocks: Iterable[Interval],
    index: GeneIntervalIndex,
    blocks_are_normalized: bool = False,
) -> Dict[str, int]:
    """
    Calculate gene overlaps for one alignment.

    Overlap lengths are accumulated directly from the interval index,
    avoiding a separate candidate-discovery and exact-overlap pass.

    Parameters
    ----------
    contig
        Reference contig containing the alignment.

    blocks
        Half-open alignment blocks.

    index
        Gene interval index.

    blocks_are_normalized
        Set to True when blocks are already sorted and non-overlapping,
        such as blocks returned by ``pysam.AlignedSegment.get_blocks()``.
    """

    if blocks_are_normalized:
        normalized_blocks = blocks
    else:
        normalized_blocks = merge_blocks(blocks)

    if not normalized_blocks:
        return {}

    contig_index = index.contigs.get(contig)

    if contig_index is None:
        return {}

    intervals = contig_index.intervals
    starts = contig_index.starts
    prefix_max_ends = contig_index.prefix_max_ends

    overlaps: Dict[str, int] = {}

    for block_start, block_end in normalized_blocks:
        if block_end <= block_start:
            continue

        position = bisect.bisect_left(
            starts,
            block_end,
        ) - 1

        while position >= 0:
            if prefix_max_ends[position] <= block_start:
                break

            gene_start, gene_end, gene_id = intervals[position]

            if gene_start < block_end and gene_end > block_start:
                if block_start > gene_start:
                    overlap_start = block_start
                else:
                    overlap_start = gene_start

                if block_end < gene_end:
                    overlap_end = block_end
                else:
                    overlap_end = gene_end

                overlap_length = overlap_end - overlap_start

                if overlap_length > 0:
                    overlaps[gene_id] = (
                        overlaps.get(gene_id, 0)
                        + overlap_length
                    )

            position -= 1

    return overlaps


def fragment_gene_overlap_details(
    fragment: Fragment,
    index: GeneIntervalIndex,
) -> Dict[str, FragmentGeneOverlap]:
    """
    Return read-specific gene-overlap evidence for one fragment.

    The two mates are evaluated independently. This preserves whether a gene
    overlaps one read or both reads.
    """

    read1_overlaps: Dict[str, int] = {}
    read2_overlaps: Dict[str, int] = {}

    if (
        fragment.read1 is not None
        and not fragment.read1.is_unmapped
        and fragment.read1.contig is not None
    ):
        read1_overlaps = blocks_gene_overlaps(
            contig=fragment.read1.contig,
            blocks=fragment.read1.blocks,
            index=index,
        )

    if (
        fragment.read2 is not None
        and not fragment.read2.is_unmapped
        and fragment.read2.contig is not None
    ):
        read2_overlaps = blocks_gene_overlaps(
            contig=fragment.read2.contig,
            blocks=fragment.read2.blocks,
            index=index,
        )

    gene_ids = set(read1_overlaps) | set(read2_overlaps)

    return {
        gene_id: FragmentGeneOverlap(
            gene_id=gene_id,
            read1_overlap=read1_overlaps.get(gene_id, 0),
            read2_overlap=read2_overlaps.get(gene_id, 0),
        )
        for gene_id in gene_ids
    }


def merge_blocks(blocks: Iterable[Interval]) -> Tuple[Interval, ...]:
    """
    Return sorted, merged, non-overlapping half-open intervals.
    """

    ordered = sorted(
        (int(start), int(end))
        for start, end in blocks
        if end > start
    )

    if not ordered:
        return ()

    merged = [ordered[0]]

    for start, end in ordered[1:]:
        previous_start, previous_end = merged[-1]

        if start <= previous_end:
            merged[-1] = (
                previous_start,
                max(previous_end, end),
            )
        else:
            merged.append((start, end))

    return tuple(merged)


def fragment_blocks_by_contig(
    fragment: Fragment,
) -> Dict[str, Tuple[Interval, ...]]:
    """
    Return merged fragment blocks grouped by reference contig.
    """

    raw: Dict[str, list[Interval]] = {}

    for part in fragment.mapped_parts:
        if part.contig is None:
            continue

        raw.setdefault(part.contig, []).extend(part.blocks)

    return {
        contig: merge_blocks(blocks)
        for contig, blocks in raw.items()
    }


def build_gene_interval_index(
    genes: Mapping[str, GeneFeature],
) -> GeneIntervalIndex:
    """
    Build an interval index from gene blocks.

    Each gene block is indexed separately. The original GeneFeature objects
    remain available for exact overlap calculation after candidate lookup.
    """

    if not genes:
        raise ValueError("Cannot build an interval index without genes")

    by_contig: Dict[str, list[Tuple[int, int, str]]] = {}

    for gene_id, gene in genes.items():
        if gene_id != gene.gene_id:
            raise ValueError(
                f"Gene mapping key {gene_id!r} does not match "
                f"GeneFeature.gene_id {gene.gene_id!r}"
            )

        for start, end in gene.blocks:
            if end <= start:
                raise ValueError(
                    f"Gene {gene_id!r} has invalid block [{start}, {end})"
                )

            by_contig.setdefault(gene.contig, []).append(
                (start, end, gene_id)
            )

    contig_indexes: Dict[str, ContigIntervalIndex] = {}

    for contig, intervals in by_contig.items():
        intervals.sort(key=lambda item: (item[0], item[1], item[2]))

        starts: list[int] = []
        prefix_max_ends: list[int] = []
        maximum_end = -1

        for start, end, _gene_id in intervals:
            starts.append(start)
            maximum_end = max(maximum_end, end)
            prefix_max_ends.append(maximum_end)

        contig_indexes[contig] = ContigIntervalIndex(
            intervals=tuple(intervals),
            starts=tuple(starts),
            prefix_max_ends=tuple(prefix_max_ends),
        )

    return GeneIntervalIndex(
        contigs=contig_indexes,
        genes=genes,
    )


def _candidate_gene_ids_for_block(
    block_start: int,
    block_end: int,
    contig_index: ContigIntervalIndex,
) -> set[str]:
    """
    Find genes having at least one indexed block overlapping the query block.
    """

    if block_end <= block_start:
        return set()

    intervals = contig_index.intervals
    starts = contig_index.starts
    prefix_max_ends = contig_index.prefix_max_ends

    position = bisect.bisect_left(starts, block_end) - 1
    candidates: set[str] = set()

    while position >= 0:
        if prefix_max_ends[position] <= block_start:
            break

        start, end, gene_id = intervals[position]

        if start < block_end and end > block_start:
            candidates.add(gene_id)

        position -= 1

    return candidates


def interval_overlap_length(
    first: Tuple[Interval, ...],
    second: Tuple[Interval, ...],
) -> int:
    """
    Return the number of bases shared by two sorted, non-overlapping
    collections of half-open intervals.
    """

    i = 0
    j = 0
    overlap = 0

    while i < len(first) and j < len(second):
        first_start, first_end = first[i]
        second_start, second_end = second[j]

        overlap += max(
            0,
            min(first_end, second_end)
            - max(first_start, second_start),
        )

        if first_end <= second_end:
            i += 1
        else:
            j += 1

    return overlap


def fragment_gene_overlaps(
    fragment: Fragment,
    index: GeneIntervalIndex,
) -> Dict[str, int]:
    """
    Calculate exact fragment overlap length for every touched gene.

    Returns
    -------
    dict
        ``gene_id -> number of overlapping genomic bases``

    Notes
    -----
    This function makes no assignment decision. Multiple genes may be
    returned. Strand, MAPQ and multimapping are handled elsewhere.
    """

    fragment_by_contig = fragment_blocks_by_contig(fragment)
    overlaps: Dict[str, int] = {}

    for contig, fragment_blocks in fragment_by_contig.items():
        contig_index = index.contigs.get(contig)

        if contig_index is None:
            continue

        candidate_ids: set[str] = set()

        for block_start, block_end in fragment_blocks:
            candidate_ids.update(
                _candidate_gene_ids_for_block(
                    block_start,
                    block_end,
                    contig_index,
                )
            )

        for gene_id in candidate_ids:
            gene = index.genes[gene_id]

            overlap = interval_overlap_length(
                fragment_blocks,
                gene.blocks,
            )

            if overlap > 0:
                overlaps[gene_id] = overlap

    return overlaps
