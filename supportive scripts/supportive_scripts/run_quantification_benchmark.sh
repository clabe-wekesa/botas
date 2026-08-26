#!/usr/bin/env bash
set -euo pipefail

export LC_ALL=C

THREADS="${THREADS:-6}"

# ------------------------------------------------------------
# Project paths
# ------------------------------------------------------------

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

REF_DIR="$ROOT/references/ecoli"
RESULTS_DIR="$ROOT/results"

# Experimental E. coli BAM used for the quantification comparison.
# This BAM should be generated from the experimental dataset described
# in the manuscript.
BAM="${BAM:-$RESULTS_DIR/experimental/ecoli/ecoli.sorted.bam}"

GFF="${GFF:-$REF_DIR/annotation.gff}"

OUTDIR="$RESULTS_DIR/quantification"
mkdir -p "$OUTDIR"

TIMEFMT="runtime_seconds\t%e\npeak_ram_kb\t%M"


# ------------------------------------------------------------
# Utility functions
# ------------------------------------------------------------

check_file() {
    [[ -f "$1" ]] || {
        echo "ERROR: Missing file: $1"
        exit 1
    }
}

check_command() {
    command -v "$1" >/dev/null 2>&1 || {
        echo "ERROR: Required command not found: $1"
        exit 1
    }
}


# ------------------------------------------------------------
# Check requirements
# ------------------------------------------------------------

check_file "$BAM"
check_file "$GFF"

check_command botas
check_command featureCounts


echo "========================================="
echo "BOTAS quantification"
echo "========================================="

/usr/bin/time \
    -f "$TIMEFMT" \
    -o "$OUTDIR/botas_quant.time.tsv" \
botas quantify \
    -b "$BAM" \
    -g "$GFF" \
    -o "$OUTDIR/botas_quant" \
    --gff-gene-attribute ID

echo
cat "$OUTDIR/botas_quant.time.tsv"


echo
echo "========================================="
echo "featureCounts"
echo "========================================="

/usr/bin/time \
    -f "$TIMEFMT" \
    -o "$OUTDIR/featurecounts.time.tsv" \
featureCounts \
    -T "$THREADS" \
    -F GFF \
    -g ID \
    -p \
    --countReadPairs \
    -a "$GFF" \
    -o "$OUTDIR/featurecounts.tsv" \
    "$BAM"

echo
cat "$OUTDIR/featurecounts.time.tsv"


# ------------------------------------------------------------
# Runtime and memory summary
# ------------------------------------------------------------

SUMMARY="$OUTDIR/quantification_performance.tsv"

printf \
"tool\tthreads\truntime_seconds\tpeak_ram_kb\tpeak_ram_mb\n" \
> "$SUMMARY"

for tool in botas_quant featurecounts
do
    runtime=$(
        awk -F'\t' \
        '$1=="runtime_seconds"{print $2}' \
        "$OUTDIR/${tool}.time.tsv"
    )

    ram=$(
        awk -F'\t' \
        '$1=="peak_ram_kb"{print $2}' \
        "$OUTDIR/${tool}.time.tsv"
    )

    ram_mb=$(
        awk -v kb="$ram" \
        'BEGIN {printf "%.2f", kb / 1024}'
    )

    printf "%s\t%s\t%s\t%s\t%s\n" \
        "$tool" \
        "$THREADS" \
        "$runtime" \
        "$ram" \
        "$ram_mb" \
        >> "$SUMMARY"
done


echo
echo "========================================="
echo "Performance summary"
echo "========================================="

if command -v column >/dev/null 2>&1; then
    column -t -s $'\t' "$SUMMARY"
else
    cat "$SUMMARY"
fi

echo
echo "BOTAS output:"
echo "$OUTDIR/botas_quant"

echo
echo "featureCounts output:"
echo "$OUTDIR/featurecounts.tsv"

echo
echo "Performance summary:"
echo "$SUMMARY"