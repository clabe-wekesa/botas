from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, Optional, Tuple
import logging

import pysam


logger = logging.getLogger(__name__)

Interval = Tuple[int, int]


@dataclass(frozen=True, slots=True)
class AlignmentPart:
    """
    One primary BAM alignment belonging to a fragment.

    No feature-assignment decisions are made here.
    """

    query_name: str
    mate_number: int
    contig: Optional[str]
    blocks: Tuple[Interval, ...]
    mapping_quality: int
    is_reverse: bool
    is_unmapped: bool
    is_proper_pair: bool
    mate_is_unmapped: bool
    nh: Optional[int]

    @property
    def start(self) -> Optional[int]:
        if not self.blocks:
            return None
        return self.blocks[0][0]

    @property
    def end(self) -> Optional[int]:
        if not self.blocks:
            return None
        return self.blocks[-1][1]


@dataclass(frozen=True, slots=True)
class Fragment:
    """
    A sequencing fragment represented by one or two primary alignments.
    """

    query_name: str
    read1: Optional[AlignmentPart]
    read2: Optional[AlignmentPart]

    @property
    def parts(self) -> Tuple[AlignmentPart, ...]:
        return tuple(
            part for part in (self.read1, self.read2) if part is not None
        )

    @property
    def mapped_parts(self) -> Tuple[AlignmentPart, ...]:
        return tuple(part for part in self.parts if not part.is_unmapped)

    @property
    def is_paired(self) -> bool:
        return self.read1 is not None and self.read2 is not None

    @property
    def is_orphan(self) -> bool:
        return not self.is_paired

    @property
    def has_mapped_alignment(self) -> bool:
        return bool(self.mapped_parts)

    @property
    def mapped_contigs(self) -> frozenset[str]:
        return frozenset(
            part.contig
            for part in self.mapped_parts
            if part.contig is not None
        )

    @property
    def is_discordant_contig(self) -> bool:
        return len(self.mapped_contigs) > 1


@dataclass(slots=True)
class FragmentStats:
    """
    Statistics describing BAM records and emitted fragments.
    """

    total_records: int = 0
    secondary_records: int = 0
    supplementary_records: int = 0
    primary_records: int = 0
    paired_fragments: int = 0
    orphan_fragments: int = 0
    fragments_with_no_mapped_part: int = 0
    discordant_contig_fragments: int = 0
    fully_mapped_fragments: int = 0
    partially_mapped_fragments: int = 0
    fully_unmapped_fragments: int = 0

    @property
    def total_fragments(self) -> int:
        return self.paired_fragments + self.orphan_fragments


@dataclass(slots=True)
class _PendingFragment:
    read1: Optional[AlignmentPart] = None
    read2: Optional[AlignmentPart] = None


def _read_nh(alignment: pysam.AlignedSegment) -> Optional[int]:
    if not alignment.has_tag("NH"):
        return None

    value = alignment.get_tag("NH")

    try:
        return int(value)
    except (TypeError, ValueError):
        logger.warning(
            "Ignoring invalid NH tag for %s: %r",
            alignment.query_name,
            value,
        )
        return None


def _alignment_to_part(
    alignment: pysam.AlignedSegment,
) -> AlignmentPart:
    if alignment.is_read1:
        mate_number = 1
    elif alignment.is_read2:
        mate_number = 2
    else:
        # A non-paired or unlabelled alignment is represented as read1.
        mate_number = 1

    if alignment.is_unmapped:
        contig = None
        blocks: Tuple[Interval, ...] = ()
    else:
        contig = alignment.reference_name
        blocks = tuple(
            (int(start), int(end))
            for start, end in alignment.get_blocks()
            if end > start
        )

    return AlignmentPart(
        query_name=alignment.query_name,
        mate_number=mate_number,
        contig=contig,
        blocks=blocks,
        mapping_quality=int(alignment.mapping_quality),
        is_reverse=bool(alignment.is_reverse),
        is_unmapped=bool(alignment.is_unmapped),
        is_proper_pair=bool(alignment.is_proper_pair),
        mate_is_unmapped=bool(alignment.mate_is_unmapped),
        nh=_read_nh(alignment),
    )


def _make_fragment(
    query_name: str,
    pending: _PendingFragment,
    stats: FragmentStats,
) -> Fragment:
    fragment = Fragment(
        query_name=query_name,
        read1=pending.read1,
        read2=pending.read2,
    )

    if fragment.is_paired:
        stats.paired_fragments += 1
    else:
        stats.orphan_fragments += 1

    mapped_part_count = len(fragment.mapped_parts)

    if mapped_part_count == 0:
        stats.fragments_with_no_mapped_part += 1
        stats.fully_unmapped_fragments += 1
    elif mapped_part_count == 1:
        stats.partially_mapped_fragments += 1
    else:
        stats.fully_mapped_fragments += 1

    if fragment.is_discordant_contig:
        stats.discordant_contig_fragments += 1

    return fragment


def iter_fragments(
    bam_path: str | Path,
    *,
    stats: Optional[FragmentStats] = None,
    require_coordinate_sorted: bool = True,
) -> Iterator[Fragment]:
    """
    Yield one Fragment per query name from a BAM file.

    Secondary and supplementary alignments are excluded. Primary mates are
    paired by query name, so coordinate-sorted BAM files are supported.

    Orphan primary alignments are emitted at the end of the BAM.
    """

    path = Path(bam_path)

    if not path.is_file():
        raise FileNotFoundError(f"BAM file does not exist: {path}")

    if stats is None:
        stats = FragmentStats()

    pending: Dict[str, _PendingFragment] = {}

    with pysam.AlignmentFile(str(path), "rb") as bam:
        sort_order = bam.header.to_dict().get("HD", {}).get("SO")

        if require_coordinate_sorted and sort_order != "coordinate":
            raise ValueError(
                f"BAM must be coordinate-sorted, but header SO={sort_order!r}: "
                f"{path}"
            )

        for alignment in bam.fetch(until_eof=True):
            stats.total_records += 1

            if alignment.is_secondary:
                stats.secondary_records += 1
                continue

            if alignment.is_supplementary:
                stats.supplementary_records += 1
                continue

            stats.primary_records += 1

            query_name = alignment.query_name

            if not query_name:
                raise ValueError(
                    "Encountered a primary BAM alignment without a query name"
                )

            part = _alignment_to_part(alignment)
            current = pending.setdefault(query_name, _PendingFragment())

            if part.mate_number == 1:
                if current.read1 is not None:
                    raise ValueError(
                        f"Multiple primary read1 alignments found for "
                        f"fragment {query_name!r}"
                    )
                current.read1 = part
            else:
                if current.read2 is not None:
                    raise ValueError(
                        f"Multiple primary read2 alignments found for "
                        f"fragment {query_name!r}"
                    )
                current.read2 = part

            if current.read1 is not None and current.read2 is not None:
                yield _make_fragment(query_name, current, stats)
                del pending[query_name]

    # Remaining entries have only one observed primary mate.
    for query_name, current in pending.items():
        yield _make_fragment(query_name, current, stats)
