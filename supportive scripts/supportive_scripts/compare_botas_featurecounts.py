#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import (
    mean_absolute_error,
    mean_squared_error,
    r2_score,
)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare gene-level counts produced by BOTAS and featureCounts "
            "and calculate agreement statistics."
        )
    )

    root = Path(__file__).resolve().parent.parent
    result_dir = root / "results" / "quantification"

    parser.add_argument(
        "--botas",
        default=str(
            result_dir
            / "botas_quant.botas"
            / "results"
            / "botas_quant.tsv"
        ),
        help="BOTAS gene-count TSV file",
    )

    parser.add_argument(
        "--featurecounts",
        default=str(result_dir / "featurecounts.tsv"),
        help="featureCounts output file",
    )

    parser.add_argument(
        "--out",
        default=str(result_dir / "botas_featurecounts_comparison.tsv"),
        help="Output TSV containing comparison statistics",
    )

    parser.add_argument(
        "--merged-out",
        default=str(result_dir / "botas_featurecounts_gene_counts.tsv"),
        help="Output TSV containing gene-by-gene counts",
    )

    args = parser.parse_args()

    botas_path = Path(args.botas)
    fc_path = Path(args.featurecounts)
    out_path = Path(args.out)
    merged_out = Path(args.merged_out)

    if not botas_path.is_file():
        raise FileNotFoundError(
            f"BOTAS count file not found: {botas_path}"
        )

    if not fc_path.is_file():
        raise FileNotFoundError(
            f"featureCounts file not found: {fc_path}"
        )

    # --------------------------------------------------------
    # Read BOTAS counts
    # --------------------------------------------------------

    botas = pd.read_csv(
        botas_path,
        sep="\t",
    )

    required_botas = {"gene_id", "count"}

    missing = required_botas - set(botas.columns)

    if missing:
        raise ValueError(
            "BOTAS file is missing required columns: "
            + ", ".join(sorted(missing))
        )

    botas = botas[["gene_id", "count"]].rename(
        columns={"count": "botas"}
    )

    # --------------------------------------------------------
    # Read featureCounts counts
    # --------------------------------------------------------

    fc = pd.read_csv(
        fc_path,
        sep="\t",
        comment="#",
    )

    if "Geneid" not in fc.columns:
        raise ValueError(
            "featureCounts output does not contain a Geneid column"
        )

    # The final column contains counts for the input BAM.
    count_column = fc.columns[-1]

    fc = fc[["Geneid", count_column]].rename(
        columns={
            "Geneid": "gene_id",
            count_column: "featurecounts",
        }
    )

    # --------------------------------------------------------
    # Match genes
    # --------------------------------------------------------

    merged = botas.merge(
        fc,
        on="gene_id",
        how="inner",
        validate="one_to_one",
    )

    if merged.empty:
        raise ValueError(
            "No matching gene IDs were found between BOTAS "
            "and featureCounts."
        )

    merged["botas"] = pd.to_numeric(
        merged["botas"],
        errors="raise",
    )

    merged["featurecounts"] = pd.to_numeric(
        merged["featurecounts"],
        errors="raise",
    )

    # --------------------------------------------------------
    # Agreement statistics
    # --------------------------------------------------------

    pearson_r, _ = pearsonr(
        merged["botas"],
        merged["featurecounts"],
    )

    spearman_r, _ = spearmanr(
        merged["botas"],
        merged["featurecounts"],
    )

    r2 = r2_score(
        merged["featurecounts"],
        merged["botas"],
    )

    mae = mean_absolute_error(
        merged["featurecounts"],
        merged["botas"],
    )

    rmse = np.sqrt(
        mean_squared_error(
            merged["featurecounts"],
            merged["botas"],
        )
    )

    identical = int(
        (merged["botas"] == merged["featurecounts"]).sum()
    )

    max_abs_difference = float(
        (
            merged["botas"] - merged["featurecounts"]
        )
        .abs()
        .max()
    )

    # --------------------------------------------------------
    # Write outputs
    # --------------------------------------------------------

    out_path.parent.mkdir(parents=True, exist_ok=True)
    merged_out.parent.mkdir(parents=True, exist_ok=True)

    merged.to_csv(
        merged_out,
        sep="\t",
        index=False,
    )

    metrics = pd.DataFrame(
        {
            "metric": [
                "genes_compared",
                "identical_gene_counts",
                "pearson_r",
                "spearman_r",
                "r2",
                "mae",
                "rmse",
                "max_absolute_difference",
            ],
            "value": [
                len(merged),
                identical,
                pearson_r,
                spearman_r,
                r2,
                mae,
                rmse,
                max_abs_difference,
            ],
        }
    )

    metrics.to_csv(
        out_path,
        sep="\t",
        index=False,
    )

    # --------------------------------------------------------
    # Console summary
    # --------------------------------------------------------

    print(f"Genes compared:          {len(merged)}")
    print(f"Identical gene counts:   {identical}")
    print(f"Pearson r:               {pearson_r:.6f}")
    print(f"Spearman r:              {spearman_r:.6f}")
    print(f"R²:                      {r2:.6f}")
    print(f"MAE:                     {mae:.3f}")
    print(f"RMSE:                    {rmse:.3f}")
    print(
        f"Maximum absolute diff.:  {max_abs_difference:.3f}"
    )

    print()
    print(f"Statistics written to: {out_path}")
    print(f"Gene-level comparison: {merged_out}")


if __name__ == "__main__":
    main()