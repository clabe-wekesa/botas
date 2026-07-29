from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Mapping, Optional

import pysam

from botas.quantify.fragments import FragmentStats
from botas.quantify.genes import (
    GeneFeature,
    load_gene_features_from_gff,
)
from botas.quantify.intervals import (
    ContigIntervalIndex,
    build_gene_interval_index,
    normalized_blocks_gene_overlaps,
)


@dataclass(frozen=True, slots=True)
class QuantificationConfig:
    """
    Configuration used for one quantification run.
    """

    bam_path: str
    gff_path: str
    feature_types: str = "gene"
    id_attribute: str = "locus_tag"
    min_overlap_bases: int = 1
    largest_overlap: bool = False

    def __post_init__(self) -> None:
        if not self.bam_path:
            raise ValueError("bam_path cannot be empty")

        if not self.gff_path:
            raise ValueError("gff_path cannot be empty")

        if not self.feature_types:
            raise ValueError("feature_types cannot be empty")

        if not self.id_attribute:
            raise ValueError("id_attribute cannot be empty")

        if self.min_overlap_bases < 1:
            raise ValueError("min_overlap_bases must be at least 1")


@dataclass(frozen=True, slots=True)
class AssignmentSummary:
    """
    Fragment-level assignment categories.
    """

    assigned: int = 0
    unmapped: int = 0
    no_feature: int = 0
    ambiguous: int = 0
    insufficient_overlap: int = 0

    @property
    def total_fragments(self) -> int:
        return (
            self.assigned
            + self.unmapped
            + self.no_feature
            + self.ambiguous
            + self.insufficient_overlap
        )

    @property
    def assignment_rate(self) -> float:
        if self.total_fragments == 0:
            return 0.0

        return self.assigned / self.total_fragments

    def as_dict(self) -> Dict[str, int]:
        return {
            "assigned": self.assigned,
            "unmapped": self.unmapped,
            "no_feature": self.no_feature,
            "ambiguous": self.ambiguous,
            "insufficient_overlap": self.insufficient_overlap,
        }


@dataclass(frozen=True, slots=True)
class QuantificationResult:
    """
    Complete output of one quantification run.
    """

    gene_counts: Mapping[str, int]
    genes: Mapping[str, GeneFeature]
    summary: AssignmentSummary
    fragment_stats: FragmentStats
    config: QuantificationConfig

    @property
    def total_assigned_counts(self) -> int:
        return sum(self.gene_counts.values())

    def validate(self) -> None:
        """
        Check internal consistency of the completed run.
        """

        if set(self.gene_counts) != set(self.genes):
            missing_counts = set(self.genes) - set(self.gene_counts)
            unknown_counts = set(self.gene_counts) - set(self.genes)

            raise ValueError(
                "Gene count identifiers do not match loaded genes. "
                f"Missing counts: {sorted(missing_counts)[:5]}; "
                f"unknown counts: {sorted(unknown_counts)[:5]}"
            )

        if any(count < 0 for count in self.gene_counts.values()):
            raise ValueError("Gene counts cannot be negative")

        if self.total_assigned_counts != self.summary.assigned:
            raise ValueError(
                "Sum of gene counts does not match assigned fragments: "
                f"{self.total_assigned_counts} != "
                f"{self.summary.assigned}"
            )

        if (
            self.summary.total_fragments
            != self.fragment_stats.total_fragments
        ):
            raise ValueError(
                "Assignment summary total does not match fragment total: "
                f"{self.summary.total_fragments} != "
                f"{self.fragment_stats.total_fragments}"
            )

        if (
            self.summary.unmapped
            != self.fragment_stats.fully_unmapped_fragments
        ):
            raise ValueError(
                "Unmapped assignment count does not match completely "
                "unmapped fragments"
            )


def _validate_input_path(
    path: str,
    *,
    description: str,
) -> None:
    input_path = Path(path)

    if not input_path.exists():
        raise FileNotFoundError(
            f"{description} does not exist: {path}"
        )

    if not input_path.is_file():
        raise ValueError(
            f"{description} is not a regular file: {path}"
        )

_NO_FEATURE = object()
_AMBIGUOUS = object()


def _assign_default_overlaps(
    read1_overlaps: Dict[str, int],
    read2_overlaps: Dict[str, int],
) -> object:
    """
    Fast assignment for:

        min_overlap_bases=1
        largest_overlap=False
    """

    if not read1_overlaps:
        if not read2_overlaps:
            return _NO_FEATURE

        if len(read2_overlaps) == 1:
            return next(iter(read2_overlaps))

        return _AMBIGUOUS

    if not read2_overlaps:
        if len(read1_overlaps) == 1:
            return next(iter(read1_overlaps))

        return _AMBIGUOUS

    # Genes supported by both mates take priority.
    if len(read1_overlaps) <= len(read2_overlaps):
        smaller = read1_overlaps
        larger = read2_overlaps
    else:
        smaller = read2_overlaps
        larger = read1_overlaps

    shared_gene: Optional[str] = None

    for gene_id in smaller:
        if gene_id not in larger:
            continue

        if shared_gene is not None:
            return _AMBIGUOUS

        shared_gene = gene_id

    if shared_gene is not None:
        return shared_gene

    # Both dictionaries are non-empty and have no shared genes.
    # Therefore, their union contains at least two genes.
    return _AMBIGUOUS


def _assign_overlaps(
    read1_overlaps: Dict[str, int],
    read2_overlaps: Dict[str, int],
    *,
    min_overlap_bases: int,
    largest_overlap: bool,
) -> tuple[str, Optional[str]]:
    """
    Apply the validated paired-read assignment rule directly.

    Genes overlapped by both mates are preferred. If no gene is supported
    by both mates, the union of genes supported by either mate is used.
    """

    if not read1_overlaps:
        if not read2_overlaps:
            return "no_feature", None

        if not largest_overlap:
            eligible_gene: Optional[str] = None

            for gene_id, overlap in read2_overlaps.items():
                if overlap < min_overlap_bases:
                    continue

                if eligible_gene is not None:
                    return "ambiguous", None

                eligible_gene = gene_id

            if eligible_gene is None:
                return "insufficient_overlap", None

            return "assigned", eligible_gene

        largest_gene: Optional[str] = None
        largest_value = -1
        largest_tied = False

        for gene_id, overlap in read2_overlaps.items():
            if overlap < min_overlap_bases:
                continue

            if overlap > largest_value:
                largest_value = overlap
                largest_gene = gene_id
                largest_tied = False
            elif overlap == largest_value:
                largest_tied = True

        if largest_gene is None:
            return "insufficient_overlap", None

        if largest_tied:
            return "ambiguous", None

        return "assigned", largest_gene

    if not read2_overlaps:
        if not largest_overlap:
            eligible_gene: Optional[str] = None

            for gene_id, overlap in read1_overlaps.items():
                if overlap < min_overlap_bases:
                    continue

                if eligible_gene is not None:
                    return "ambiguous", None

                eligible_gene = gene_id

            if eligible_gene is None:
                return "insufficient_overlap", None

            return "assigned", eligible_gene

        largest_gene: Optional[str] = None
        largest_value = -1
        largest_tied = False

        for gene_id, overlap in read1_overlaps.items():
            if overlap < min_overlap_bases:
                continue

            if overlap > largest_value:
                largest_value = overlap
                largest_gene = gene_id
                largest_tied = False
            elif overlap == largest_value:
                largest_tied = True

        if largest_gene is None:
            return "insufficient_overlap", None

        if largest_tied:
            return "ambiguous", None

        return "assigned", largest_gene

    # Search the smaller dictionary for genes present in both mates.
    if len(read1_overlaps) <= len(read2_overlaps):
        smaller = read1_overlaps
        larger = read2_overlaps
    else:
        smaller = read2_overlaps
        larger = read1_overlaps

    shared_found = False

    if not largest_overlap:
        eligible_gene: Optional[str] = None

        for gene_id, overlap in smaller.items():
            mate_overlap = larger.get(gene_id)

            if mate_overlap is None:
                continue

            shared_found = True
            overlap += mate_overlap

            if overlap < min_overlap_bases:
                continue

            if eligible_gene is not None:
                return "ambiguous", None

            eligible_gene = gene_id

        if shared_found:
            if eligible_gene is None:
                return "insufficient_overlap", None

            return "assigned", eligible_gene

        # No shared genes: evaluate the union. Because no shared gene exists,
        # the dictionaries can be traversed independently.
        for gene_id, overlap in read1_overlaps.items():
            if overlap < min_overlap_bases:
                continue

            if eligible_gene is not None:
                return "ambiguous", None

            eligible_gene = gene_id

        for gene_id, overlap in read2_overlaps.items():
            if overlap < min_overlap_bases:
                continue

            if eligible_gene is not None:
                return "ambiguous", None

            eligible_gene = gene_id

        if eligible_gene is None:
            return "insufficient_overlap", None

        return "assigned", eligible_gene

    largest_gene: Optional[str] = None
    largest_value = -1
    largest_tied = False

    for gene_id, overlap in smaller.items():
        mate_overlap = larger.get(gene_id)

        if mate_overlap is None:
            continue

        shared_found = True
        overlap += mate_overlap

        if overlap < min_overlap_bases:
            continue

        if overlap > largest_value:
            largest_value = overlap
            largest_gene = gene_id
            largest_tied = False
        elif overlap == largest_value:
            largest_tied = True

    if shared_found:
        if largest_gene is None:
            return "insufficient_overlap", None

        if largest_tied:
            return "ambiguous", None

        return "assigned", largest_gene

    # No shared genes: evaluate the union.
    for gene_id, overlap in read1_overlaps.items():
        if overlap < min_overlap_bases:
            continue

        if overlap > largest_value:
            largest_value = overlap
            largest_gene = gene_id
            largest_tied = False
        elif overlap == largest_value:
            largest_tied = True

    for gene_id, overlap in read2_overlaps.items():
        if overlap < min_overlap_bases:
            continue

        if overlap > largest_value:
            largest_value = overlap
            largest_gene = gene_id
            largest_tied = False
        elif overlap == largest_value:
            largest_tied = True

    if largest_gene is None:
        return "insufficient_overlap", None

    if largest_tied:
        return "ambiguous", None

    return "assigned", largest_gene


def quantify_bam(
    bam_path: str,
    gff_path: str,
    *,
    feature_types: str = "gene",
    id_attribute: str = "locus_tag",
    min_overlap_bases: int = 1,
    largest_overlap: bool = False,
) -> QuantificationResult:
    """
    Quantify paired-end fragments against GFF gene annotations.

    The validated default behaviour corresponds to featureCounts with:

        -p --countReadPairs -s 0 --maxMOp 1000
    """

    _validate_input_path(
        bam_path,
        description="BAM file",
    )
    _validate_input_path(
        gff_path,
        description="GFF file",
    )

    config = QuantificationConfig(
        bam_path=str(bam_path),
        gff_path=str(gff_path),
        feature_types=feature_types,
        id_attribute=id_attribute,
        min_overlap_bases=min_overlap_bases,
        largest_overlap=largest_overlap,
    )

    genes = load_gene_features_from_gff(
        gff_path,
        feature_types=feature_types,
        id_attribute=id_attribute,
    )

    index = build_gene_interval_index(genes)

    gene_counts = {
        gene_id: 0
        for gene_id in genes
    }

    assigned = 0
    unmapped = 0
    no_feature = 0
    ambiguous = 0
    insufficient_overlap = 0

    fragment_stats = FragmentStats()

    # Coordinate-sorted BAM records may have many alignments between mates.
    pending: Dict[str, pysam.AlignedSegment] = {}

    reference_contig_indexes: tuple[Optional[ContigIntervalIndex], ...] = ()
    overlap_function = normalized_blocks_gene_overlaps

    use_default_assignment = (
        min_overlap_bases == 1
        and not largest_overlap
    )

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        reference_contig_indexes = tuple(
            index.contigs.get(reference_name)
            for reference_name in bam.references
        )

        for read in bam.fetch(until_eof=True):
            fragment_stats.total_records += 1

            if read.is_secondary:
                fragment_stats.secondary_records += 1
                continue

            if read.is_supplementary:
                fragment_stats.supplementary_records += 1
                continue

            fragment_stats.primary_records += 1

            query_name = read.query_name
            mate = pending.pop(query_name, None)

            if mate is None:
                pending[query_name] = read
                continue

            fragment_stats.paired_fragments += 1

            if mate.is_read1:
                read1 = mate
                read2 = read
            elif read.is_read1:
                read1 = read
                read2 = mate
            else:
                read1 = mate
                read2 = read

            read1_mapped = not read1.is_unmapped
            read2_mapped = not read2.is_unmapped

            if read1_mapped and read2_mapped:
                fragment_stats.fully_mapped_fragments += 1

                if read1.reference_id != read2.reference_id:
                    fragment_stats.discordant_contig_fragments += 1

            elif read1_mapped or read2_mapped:
                fragment_stats.partially_mapped_fragments += 1

            else:
                fragment_stats.fully_unmapped_fragments += 1
                fragment_stats.fragments_with_no_mapped_part += 1
                unmapped += 1
                continue

            if read1_mapped:
                read1_blocks = read1.get_blocks()
                read1_contig_index = reference_contig_indexes[
                    read1.reference_id
                ]

                if read1_blocks and read1_contig_index is not None:
                    read1_overlaps = overlap_function(
                        read1_blocks,
                        read1_contig_index,
                    )
                else:
                    read1_overlaps = {}
            else:
                read1_overlaps = {}

            if read2_mapped:
                read2_blocks = read2.get_blocks()
                read2_contig_index = reference_contig_indexes[
                    read2.reference_id
                ]

                if read2_blocks and read2_contig_index is not None:
                    read2_overlaps = overlap_function(
                        read2_blocks,
                        read2_contig_index,
                    )
                else:
                    read2_overlaps = {}
            else:
                read2_overlaps = {}

            if use_default_assignment:
                assignment = _assign_default_overlaps(
                    read1_overlaps,
                    read2_overlaps,
                )

                if assignment is _NO_FEATURE:
                    no_feature += 1

                elif assignment is _AMBIGUOUS:
                    ambiguous += 1

                else:
                    gene_counts[assignment] += 1
                    assigned += 1

            else:
                status, gene_id = _assign_overlaps(
                    read1_overlaps,
                    read2_overlaps,
                    min_overlap_bases=min_overlap_bases,
                    largest_overlap=largest_overlap,
                )

                if status == "assigned":
                    if gene_id is None:
                        raise RuntimeError(
                            "Assigned fragment has no gene identifier"
                        )

                    gene_counts[gene_id] += 1
                    assigned += 1

                elif status == "no_feature":
                    no_feature += 1

                elif status == "ambiguous":
                    ambiguous += 1

                elif status == "insufficient_overlap":
                    insufficient_overlap += 1

                else:
                    raise RuntimeError(
                        f"Unsupported assignment status: {status!r}"
                    )

    # Preserve fragment accounting for incomplete pairs.
    for read in pending.values():
        fragment_stats.orphan_fragments += 1

        if read.is_unmapped:
            fragment_stats.fully_unmapped_fragments += 1
            fragment_stats.fragments_with_no_mapped_part += 1
            unmapped += 1
            continue

        fragment_stats.partially_mapped_fragments += 1

        blocks = read.get_blocks()
        contig_index = reference_contig_indexes[read.reference_id]

        if blocks and contig_index is not None:
            overlaps = overlap_function(
                blocks,
                contig_index,
            )
        else:
            overlaps = {}

        if read.is_read2:
            read1_overlaps = {}
            read2_overlaps = overlaps
        else:
            read1_overlaps = overlaps
            read2_overlaps = {}

        status, gene_id = _assign_overlaps(
            read1_overlaps,
            read2_overlaps,
            min_overlap_bases=min_overlap_bases,
            largest_overlap=largest_overlap,
        )

        if status == "assigned":
            if gene_id is None:
                raise RuntimeError(
                    "Assigned orphan has no gene identifier"
                )

            gene_counts[gene_id] += 1
            assigned += 1

        elif status == "no_feature":
            no_feature += 1

        elif status == "ambiguous":
            ambiguous += 1

        elif status == "insufficient_overlap":
            insufficient_overlap += 1

        else:
            raise RuntimeError(
                f"Unsupported assignment status: {status!r}"
            )

    summary = AssignmentSummary(
        assigned=assigned,
        unmapped=unmapped,
        no_feature=no_feature,
        ambiguous=ambiguous,
        insufficient_overlap=insufficient_overlap,
    )

    quantification = QuantificationResult(
        gene_counts=gene_counts,
        genes=genes,
        summary=summary,
        fragment_stats=fragment_stats,
        config=config,
    )

    quantification.validate()

    return quantification