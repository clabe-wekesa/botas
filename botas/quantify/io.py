# botas/quantify/io.py
from __future__ import annotations

import csv
from typing import Dict, Mapping, Optional, Sequence, Tuple


# ---------------------------------------------------------------------
# Type aliases
# ---------------------------------------------------------------------

FeatureCounts = Mapping[str, float]
FeatureLengths = Mapping[str, int]

# sample_name -> feature_id -> value
SampleFeatureValues = Mapping[str, Mapping[str, float]]


# ---------------------------------------------------------------------
# Normalization
# ---------------------------------------------------------------------

def _safe_div(numerator: float, denominator: float) -> float:
    """
    Divide safely, returning 0.0 when the denominator is zero.
    """

    if denominator == 0:
        return 0.0

    return numerator / denominator


def compute_rpkm(
    counts: FeatureCounts,
    lengths_bp: FeatureLengths,
    total_mapped: float,
) -> Dict[str, float]:
    """
    Compute RPKM values.

    RPKM = 1e9 × count / (total mapped or assigned weight × length_bp)

    Parameters
    ----------
    counts
        Feature-level counts.

    lengths_bp
        Feature lengths in base pairs.

    total_mapped
        Normalization denominator. In the current BOTAS workflow, this is
        normally the total assigned weight, calculated as sum(counts.values()).

    Returns
    -------
    dict
        feature_id -> RPKM
    """

    output: Dict[str, float] = {}

    if total_mapped <= 0:
        return {feature_id: 0.0 for feature_id in counts}

    for feature_id, count in counts.items():
        length = float(lengths_bp.get(feature_id, 0))

        if length <= 0:
            output[feature_id] = 0.0
            continue

        output[feature_id] = _safe_div(
            1e9 * float(count),
            total_mapped * length,
        )

    return output


def compute_tpm(
    counts: FeatureCounts,
    lengths_bp: FeatureLengths,
) -> Dict[str, float]:
    """
    Compute TPM values.

    TPM is calculated by first computing reads per kilobase:

        RPK = count / (length_bp / 1000)

    and then scaling all RPK values so that they sum to one million.
    """

    rpk: Dict[str, float] = {}

    for feature_id, count in counts.items():
        length = float(lengths_bp.get(feature_id, 0))

        if length <= 0:
            rpk[feature_id] = 0.0
            continue

        rpk[feature_id] = _safe_div(
            float(count),
            length / 1_000.0,
        )

    denominator = sum(rpk.values())

    if denominator <= 0:
        return {feature_id: 0.0 for feature_id in counts}

    return {
        feature_id: (value / denominator) * 1_000_000.0
        for feature_id, value in rpk.items()
    }


def compute_multi_sample_normalization(
    counts_by_sample: SampleFeatureValues,
    lengths_bp: FeatureLengths,
) -> Tuple[
    Dict[str, Dict[str, float]],
    Dict[str, Dict[str, float]],
]:
    """
    Compute RPKM and TPM independently for every sample.

    Parameters
    ----------
    counts_by_sample
        Mapping:

            sample_name -> {feature_id: count}

    lengths_bp
        Feature lengths shared by all samples.

    Returns
    -------
    tuple
        rpkm_by_sample, tpm_by_sample
    """

    rpkm_by_sample: Dict[str, Dict[str, float]] = {}
    tpm_by_sample: Dict[str, Dict[str, float]] = {}

    for sample_name, counts in counts_by_sample.items():
        total_assigned = float(sum(counts.values()))

        rpkm_by_sample[sample_name] = compute_rpkm(
            counts=counts,
            lengths_bp=lengths_bp,
            total_mapped=total_assigned,
        )

        tpm_by_sample[sample_name] = compute_tpm(
            counts=counts,
            lengths_bp=lengths_bp,
        )

    return rpkm_by_sample, tpm_by_sample


# ---------------------------------------------------------------------
# Validation and formatting helpers
# ---------------------------------------------------------------------

def _format_float(value: float) -> str:
    """
    Format numeric output consistently.
    """

    return f"{float(value):.6f}"


def _validate_sample_tables(
    counts_by_sample: SampleFeatureValues,
    rpkm_by_sample: Optional[SampleFeatureValues] = None,
    tpm_by_sample: Optional[SampleFeatureValues] = None,
) -> Sequence[str]:
    """
    Validate that multi-sample tables contain compatible sample names.

    Returns sample names in input insertion order.
    """

    if not counts_by_sample:
        raise ValueError("counts_by_sample is empty")

    sample_names = list(counts_by_sample.keys())
    expected = set(sample_names)

    if rpkm_by_sample is not None:
        observed = set(rpkm_by_sample.keys())

        if observed != expected:
            missing = expected - observed
            extra = observed - expected

            raise ValueError(
                "RPKM sample names do not match count sample names. "
                f"Missing={sorted(missing)}, extra={sorted(extra)}"
            )

    if tpm_by_sample is not None:
        observed = set(tpm_by_sample.keys())

        if observed != expected:
            missing = expected - observed
            extra = observed - expected

            raise ValueError(
                "TPM sample names do not match count sample names. "
                f"Missing={sorted(missing)}, extra={sorted(extra)}"
            )

    return sample_names


def _collect_feature_ids(
    values_by_sample: SampleFeatureValues,
) -> Sequence[str]:
    """
    Collect the union of all feature IDs across samples.
    """

    feature_ids = {
        feature_id
        for sample_values in values_by_sample.values()
        for feature_id in sample_values
    }

    return sorted(feature_ids)


# ---------------------------------------------------------------------
# Single-sample gene output
# ---------------------------------------------------------------------

def write_gene_counts_tsv(
    out_tsv: str,
    counts: FeatureCounts,
    lengths_bp: FeatureLengths,
    rpkm: FeatureCounts,
    tpm: FeatureCounts,
) -> None:
    """
    Write a single-sample gene quantification table.
    """

    with open(
        out_tsv,
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")

        writer.writerow([
            "gene_id",
            "length_bp",
            "count",
            "RPKM",
            "TPM",
        ])

        for gene_id in sorted(counts):
            writer.writerow([
                gene_id,
                lengths_bp.get(gene_id, 0),
                _format_float(counts.get(gene_id, 0.0)),
                _format_float(rpkm.get(gene_id, 0.0)),
                _format_float(tpm.get(gene_id, 0.0)),
            ])


# ---------------------------------------------------------------------
# Multi-sample gene output
# ---------------------------------------------------------------------

def write_gene_count_matrix_tsv(
    out_tsv: str,
    counts_by_sample: SampleFeatureValues,
    lengths_bp: FeatureLengths,
) -> None:
    """
    Write a gene count matrix for multiple BAM files.

    Output structure:

        gene_id    length_bp    sample_1    sample_2    ...
    """

    sample_names = _validate_sample_tables(counts_by_sample)
    gene_ids = _collect_feature_ids(counts_by_sample)

    with open(
        out_tsv,
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")

        writer.writerow([
            "gene_id",
            "length_bp",
            *sample_names,
        ])

        for gene_id in gene_ids:
            writer.writerow([
                gene_id,
                lengths_bp.get(gene_id, 0),
                *[
                    _format_float(
                        counts_by_sample[sample_name].get(
                            gene_id,
                            0.0,
                        )
                    )
                    for sample_name in sample_names
                ],
            ])


def write_gene_expression_matrix_tsv(
    out_tsv: str,
    counts_by_sample: SampleFeatureValues,
    lengths_bp: FeatureLengths,
    rpkm_by_sample: Optional[SampleFeatureValues] = None,
    tpm_by_sample: Optional[SampleFeatureValues] = None,
) -> None:
    """
    Write counts, RPKM, and TPM for multiple samples in one table.

    For every sample, the output contains three columns:

        sample.count
        sample.RPKM
        sample.TPM

    When RPKM or TPM tables are omitted, they are computed automatically.
    """

    if rpkm_by_sample is None or tpm_by_sample is None:
        computed_rpkm, computed_tpm = compute_multi_sample_normalization(
            counts_by_sample=counts_by_sample,
            lengths_bp=lengths_bp,
        )

        if rpkm_by_sample is None:
            rpkm_by_sample = computed_rpkm

        if tpm_by_sample is None:
            tpm_by_sample = computed_tpm

    sample_names = _validate_sample_tables(
        counts_by_sample=counts_by_sample,
        rpkm_by_sample=rpkm_by_sample,
        tpm_by_sample=tpm_by_sample,
    )

    gene_ids = _collect_feature_ids(counts_by_sample)

    header = [
        "gene_id",
        "length_bp",
    ]

    for sample_name in sample_names:
        header.extend([
            f"{sample_name}.count",
            f"{sample_name}.RPKM",
            f"{sample_name}.TPM",
        ])

    with open(
        out_tsv,
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)

        for gene_id in gene_ids:
            row = [
                gene_id,
                lengths_bp.get(gene_id, 0),
            ]

            for sample_name in sample_names:
                row.extend([
                    _format_float(
                        counts_by_sample[sample_name].get(
                            gene_id,
                            0.0,
                        )
                    ),
                    _format_float(
                        rpkm_by_sample[sample_name].get(
                            gene_id,
                            0.0,
                        )
                    ),
                    _format_float(
                        tpm_by_sample[sample_name].get(
                            gene_id,
                            0.0,
                        )
                    ),
                ])

            writer.writerow(row)


# ---------------------------------------------------------------------
# Single-sample operon output
# ---------------------------------------------------------------------

def write_operon_counts_tsv(
    out_tsv: str,
    op_meta: Mapping[str, Mapping],
    counts: FeatureCounts,
    lengths_bp: FeatureLengths,
    rpkm: FeatureCounts,
    tpm: FeatureCounts,
) -> None:
    """
    Write a single-sample operon quantification table.
    """

    with open(
        out_tsv,
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")

        writer.writerow([
            "operon_id",
            "contig",
            "start",
            "end",
            "strand",
            "genes",
            "length_bp",
            "count",
            "RPKM",
            "TPM",
        ])

        for operon_id in sorted(counts):
            metadata = op_meta.get(operon_id, {})

            writer.writerow([
                operon_id,
                metadata.get("contig", ""),
                metadata.get("start", -1),
                metadata.get("end", -1),
                metadata.get("strand", "."),
                metadata.get("genes", ""),
                lengths_bp.get(operon_id, 0),
                _format_float(counts.get(operon_id, 0.0)),
                _format_float(rpkm.get(operon_id, 0.0)),
                _format_float(tpm.get(operon_id, 0.0)),
            ])


# ---------------------------------------------------------------------
# Multi-sample operon output
# ---------------------------------------------------------------------

def write_operon_count_matrix_tsv(
    out_tsv: str,
    op_meta: Mapping[str, Mapping],
    counts_by_sample: SampleFeatureValues,
    lengths_bp: FeatureLengths,
) -> None:
    """
    Write an operon count matrix for multiple samples.
    """

    sample_names = _validate_sample_tables(counts_by_sample)
    operon_ids = _collect_feature_ids(counts_by_sample)

    with open(
        out_tsv,
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")

        writer.writerow([
            "operon_id",
            "contig",
            "start",
            "end",
            "strand",
            "genes",
            "length_bp",
            *sample_names,
        ])

        for operon_id in operon_ids:
            metadata = op_meta.get(operon_id, {})

            writer.writerow([
                operon_id,
                metadata.get("contig", ""),
                metadata.get("start", -1),
                metadata.get("end", -1),
                metadata.get("strand", "."),
                metadata.get("genes", ""),
                lengths_bp.get(operon_id, 0),
                *[
                    _format_float(
                        counts_by_sample[sample_name].get(
                            operon_id,
                            0.0,
                        )
                    )
                    for sample_name in sample_names
                ],
            ])


def write_operon_expression_matrix_tsv(
    out_tsv: str,
    op_meta: Mapping[str, Mapping],
    counts_by_sample: SampleFeatureValues,
    lengths_bp: FeatureLengths,
    rpkm_by_sample: Optional[SampleFeatureValues] = None,
    tpm_by_sample: Optional[SampleFeatureValues] = None,
) -> None:
    """
    Write counts, RPKM, and TPM for multiple operon samples.
    """

    if rpkm_by_sample is None or tpm_by_sample is None:
        computed_rpkm, computed_tpm = compute_multi_sample_normalization(
            counts_by_sample=counts_by_sample,
            lengths_bp=lengths_bp,
        )

        if rpkm_by_sample is None:
            rpkm_by_sample = computed_rpkm

        if tpm_by_sample is None:
            tpm_by_sample = computed_tpm

    sample_names = _validate_sample_tables(
        counts_by_sample=counts_by_sample,
        rpkm_by_sample=rpkm_by_sample,
        tpm_by_sample=tpm_by_sample,
    )

    operon_ids = _collect_feature_ids(counts_by_sample)

    header = [
        "operon_id",
        "contig",
        "start",
        "end",
        "strand",
        "genes",
        "length_bp",
    ]

    for sample_name in sample_names:
        header.extend([
            f"{sample_name}.count",
            f"{sample_name}.RPKM",
            f"{sample_name}.TPM",
        ])

    with open(
        out_tsv,
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)

        for operon_id in operon_ids:
            metadata = op_meta.get(operon_id, {})

            row = [
                operon_id,
                metadata.get("contig", ""),
                metadata.get("start", -1),
                metadata.get("end", -1),
                metadata.get("strand", "."),
                metadata.get("genes", ""),
                lengths_bp.get(operon_id, 0),
            ]

            for sample_name in sample_names:
                row.extend([
                    _format_float(
                        counts_by_sample[sample_name].get(
                            operon_id,
                            0.0,
                        )
                    ),
                    _format_float(
                        rpkm_by_sample[sample_name].get(
                            operon_id,
                            0.0,
                        )
                    ),
                    _format_float(
                        tpm_by_sample[sample_name].get(
                            operon_id,
                            0.0,
                        )
                    ),
                ])

            writer.writerow(row)