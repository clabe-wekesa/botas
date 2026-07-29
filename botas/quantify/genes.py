from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Mapping, Sequence, Tuple
from urllib.parse import unquote
import logging


logger = logging.getLogger(__name__)

Interval = Tuple[int, int]


@dataclass(frozen=True, slots=True)
class GeneFeature:
    """
    Genomic feature used for quantification.

    Coordinates are stored as 0-based, half-open intervals:
    [start, end)
    """

    gene_id: str
    contig: str
    strand: str
    blocks: Tuple[Interval, ...]
    length: int

    @property
    def start(self) -> int:
        return self.blocks[0][0]

    @property
    def end(self) -> int:
        return self.blocks[-1][1]


def parse_gff3_attributes(text: str) -> Dict[str, str]:
    """
    Parse the ninth column of a GFF3 record.

    Percent-encoded values are decoded according to the GFF3 format.
    Empty attributes and malformed entries without '=' are ignored.
    """

    attributes: Dict[str, str] = {}

    for raw_field in text.strip().strip(";").split(";"):
        field = raw_field.strip()

        if not field or "=" not in field:
            continue

        key, value = field.split("=", 1)
        key = unquote(key.strip())
        value = unquote(value.strip())

        if key:
            attributes[key] = value

    return attributes


def merge_intervals(intervals: Iterable[Interval]) -> Tuple[Interval, ...]:
    """
    Merge overlapping or directly adjacent half-open intervals.
    """

    ordered = sorted(intervals)

    if not ordered:
        return ()

    merged = [ordered[0]]

    for start, end in ordered[1:]:
        previous_start, previous_end = merged[-1]

        if start <= previous_end:
            merged[-1] = (previous_start, max(previous_end, end))
        else:
            merged.append((start, end))

    return tuple(merged)


def _normalise_feature_types(
    feature_types: str | Sequence[str],
) -> frozenset[str]:
    if isinstance(feature_types, str):
        values = [feature_types]
    else:
        values = list(feature_types)

    cleaned = frozenset(value.strip() for value in values if value.strip())

    if not cleaned:
        raise ValueError("At least one GFF feature type must be provided")

    return cleaned


def load_gene_features_from_gff(
    gff_path: str | Path,
    *,
    feature_types: str | Sequence[str] = "gene",
    id_attribute: str = "locus_tag",
) -> Dict[str, GeneFeature]:
    """
    Load quantifiable features from a GFF3 annotation.

    Parameters
    ----------
    gff_path
        Input GFF3 file.
    feature_types
        GFF feature type or types to load, normally ``gene``.
    id_attribute
        Attribute used as the unique feature identifier, for example
        ``locus_tag``, ``ID`` or ``gene_id``.

    Returns
    -------
    dict
        Mapping from feature identifier to GeneFeature.

    Raises
    ------
    ValueError
        If coordinates are invalid, identifiers are duplicated across
        incompatible contigs or strands, or no features can be loaded.
    """

    path = Path(gff_path)

    if not path.is_file():
        raise FileNotFoundError(f"GFF file does not exist: {path}")

    accepted_types = _normalise_feature_types(feature_types)

    # gene_id -> list of (contig, strand, start, end, line_number)
    records: Dict[str, list[tuple[str, str, int, int, int]]] = {}

    matching_records = 0
    missing_identifier = 0

    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")

            if not line or line.startswith("#"):
                continue

            columns = line.split("\t")

            if len(columns) != 9:
                raise ValueError(
                    f"{path}:{line_number}: expected 9 tab-separated "
                    f"GFF columns, found {len(columns)}"
                )

            (
                contig,
                _source,
                feature_type,
                start_text,
                end_text,
                _score,
                strand,
                _phase,
                attribute_text,
            ) = columns

            if feature_type not in accepted_types:
                continue

            matching_records += 1

            if strand not in {"+", "-"}:
                raise ValueError(
                    f"{path}:{line_number}: feature has unsupported strand "
                    f"{strand!r}; expected '+' or '-'"
                )

            try:
                start_1based = int(start_text)
                end_1based = int(end_text)
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: start and end must be integers"
                ) from error

            if start_1based < 1 or end_1based < start_1based:
                raise ValueError(
                    f"{path}:{line_number}: invalid coordinates "
                    f"{start_1based}-{end_1based}"
                )

            attributes = parse_gff3_attributes(attribute_text)
            gene_id = attributes.get(id_attribute)

            if not gene_id:
                missing_identifier += 1
                continue

            # GFF3: 1-based inclusive -> Python: 0-based half-open.
            start_0based = start_1based - 1
            end_0based = end_1based

            records.setdefault(gene_id, []).append(
                (
                    contig,
                    strand,
                    start_0based,
                    end_0based,
                    line_number,
                )
            )

    if matching_records == 0:
        requested = ", ".join(sorted(accepted_types))
        raise ValueError(
            f"No records of GFF feature type(s) {requested} were found in {path}"
        )

    if not records:
        raise ValueError(
            f"No features could be loaded from {path} using attribute "
            f"{id_attribute!r}. Matching records without that attribute: "
            f"{missing_identifier}"
        )

    features: Dict[str, GeneFeature] = {}

    for gene_id, gene_records in records.items():
        first_contig, first_strand, _, _, _ = gene_records[0]
        intervals: list[Interval] = []

        for contig, strand, start, end, line_number in gene_records:
            if contig != first_contig:
                raise ValueError(
                    f"Feature {gene_id!r} occurs on multiple contigs: "
                    f"{first_contig!r} and {contig!r} "
                    f"(detected at line {line_number})"
                )

            if strand != first_strand:
                raise ValueError(
                    f"Feature {gene_id!r} occurs on multiple strands: "
                    f"{first_strand!r} and {strand!r} "
                    f"(detected at line {line_number})"
                )

            intervals.append((start, end))

        blocks = merge_intervals(intervals)
        length = sum(end - start for start, end in blocks)

        if length <= 0:
            raise ValueError(
                f"Feature {gene_id!r} has zero or negative merged length"
            )

        features[gene_id] = GeneFeature(
            gene_id=gene_id,
            contig=first_contig,
            strand=first_strand,
            blocks=blocks,
            length=length,
        )

    logger.info(
        "Loaded %d features from %s using type(s) %s and attribute %s; "
        "%d matching records lacked the identifier",
        len(features),
        path,
        ",".join(sorted(accepted_types)),
        id_attribute,
        missing_identifier,
    )

    return features
