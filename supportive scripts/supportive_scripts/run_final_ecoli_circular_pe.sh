#!/usr/bin/env bash
set -euo pipefail

export LC_ALL=C

THREADS="${THREADS:-6}"
TOLERANCE="${TOLERANCE:-10}"
DATASET="ecoli_circular_pe"

# ------------------------------------------------------------
# Project paths
# ------------------------------------------------------------

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

DATA_DIR="$ROOT/botas_manuscript_data"
SCRIPT_DIR="$ROOT/botas_manuscript_scripts"
REF_DIR="$ROOT/references/ecoli"
INDEX_DIR="$ROOT/index/ecoli"
RESULTS_DIR="$ROOT/results"

REF="${REF:-$REF_DIR/genome.fa}"

R1="$DATA_DIR/simulated/circular_wrap/ecoli/ecoli_wrap_20k_R1.fq.gz"
R2="$DATA_DIR/simulated/circular_wrap/ecoli/ecoli_wrap_20k_R2.fq.gz"
TRUTH="$DATA_DIR/simulated/circular_wrap/ecoli/ecoli_wrap_20k_truth.tsv"

OUTDIR="$RESULTS_DIR/alignment/$DATASET"
METRIC_DIR="$RESULTS_DIR/metrics/$DATASET"
TABLE_DIR="$RESULTS_DIR/tables"

mkdir -p "$OUTDIR" "$METRIC_DIR" "$TABLE_DIR"

SUMMARY="$TABLE_DIR/${DATASET}_summary.tsv"

printf \
"tool\tdataset\tthreads\truntime_seconds\tpeak_ram_kb\tpeak_ram_mb\tprecision\trecall\tf1\n" \
> "$SUMMARY"


# ------------------------------------------------------------
# Utility functions
# ------------------------------------------------------------

check_file() {
    [[ -f "$1" ]] || {
        echo "ERROR: Missing file: $1"
        exit 1
    }
}

extract_metric() {
    local file="$1"
    local key="$2"

    awk -F'\t' -v k="$key" '$1==k {print $2}' "$file"
}

add_summary() {
    local tool="$1"
    local timefile="$2"
    local metricfile="$3"

    local runtime ram ram_mb precision recall f1

    runtime=$(awk -F'\t' '$1=="runtime_seconds"{print $2}' "$timefile")
    ram=$(awk -F'\t' '$1=="peak_ram_kb"{print $2}' "$timefile")

    ram_mb=$(awk -v kb="$ram" 'BEGIN {printf "%.2f", kb / 1024}')

    precision=$(extract_metric "$metricfile" "precision")
    recall=$(extract_metric "$metricfile" "recall")
    f1=$(extract_metric "$metricfile" "f1")

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$tool" \
        "$DATASET" \
        "$THREADS" \
        "$runtime" \
        "$ram" \
        "$ram_mb" \
        "$precision" \
        "$recall" \
        "$f1" \
        >> "$SUMMARY"
}

run_eval() {
    local bam="$1"
    local out="$2"

    python "$SCRIPT_DIR/evaluate_pe_bam.py" \
        --bam "$bam" \
        --truth "$TRUTH" \
        --out "$out" \
        --tolerance "$TOLERANCE"
}


# ------------------------------------------------------------
# Check input files
# ------------------------------------------------------------

check_file "$REF"
check_file "$R1"
check_file "$R2"
check_file "$TRUTH"

echo "=== Final E. coli circular-wrap PE benchmark ==="


# ============================================================
# BOTAS
# ============================================================

BOTAS_INDEX="$INDEX_DIR/botas/genome.circular5.botas/results/genome.circular5.idx"
check_file "$BOTAS_INDEX"

echo "=== BOTAS ==="

/usr/bin/time \
    -f "runtime_seconds\t%e\npeak_ram_kb\t%M" \
    -o "$METRIC_DIR/botas.time.tsv" \
botas align \
    --ref "$REF" \
    --index "$BOTAS_INDEX" \
    --fq1 "$R1" \
    --fq2 "$R2" \
    --circular \
    --pool \
    -t "$THREADS" \
    -o "$OUTDIR/botas.bam"

BOTAS_BAM="$OUTDIR/botas.botas/results/botas.bam"
check_file "$BOTAS_BAM"

run_eval \
    "$BOTAS_BAM" \
    "$METRIC_DIR/botas.metrics.tsv"

add_summary \
    "BOTAS" \
    "$METRIC_DIR/botas.time.tsv" \
    "$METRIC_DIR/botas.metrics.tsv"


# ============================================================
# Bowtie2
# ============================================================

BOWTIE2_INDEX="$INDEX_DIR/bowtie2/ecoli"
check_file "${BOWTIE2_INDEX}.1.bt2"

echo "=== Bowtie2 ==="

/usr/bin/time \
    -f "runtime_seconds\t%e\npeak_ram_kb\t%M" \
    -o "$METRIC_DIR/bowtie2.time.tsv" \
bash -o pipefail -c "
bowtie2 \
    -x '$BOWTIE2_INDEX' \
    -1 '$R1' \
    -2 '$R2' \
    -p '$THREADS' \
    2> '$METRIC_DIR/bowtie2.log' \
| samtools view \
    -@ '$THREADS' \
    -b \
    -o '$OUTDIR/bowtie2.bam'
"

run_eval \
    "$OUTDIR/bowtie2.bam" \
    "$METRIC_DIR/bowtie2.metrics.tsv"

add_summary \
    "Bowtie2" \
    "$METRIC_DIR/bowtie2.time.tsv" \
    "$METRIC_DIR/bowtie2.metrics.tsv"


# ============================================================
# BWA-MEM2
# ============================================================

BWA_REF="$INDEX_DIR/bwa/genome.fa"
check_file "$BWA_REF"

echo "=== BWA-MEM2 ==="

/usr/bin/time \
    -f "runtime_seconds\t%e\npeak_ram_kb\t%M" \
    -o "$METRIC_DIR/bwa_mem2.time.tsv" \
bash -o pipefail -c "
bwa-mem2 mem \
    -t '$THREADS' \
    '$BWA_REF' \
    '$R1' \
    '$R2' \
    2> '$METRIC_DIR/bwa_mem2.log' \
| samtools view \
    -@ '$THREADS' \
    -b \
    -o '$OUTDIR/bwa_mem2.bam'
"

run_eval \
    "$OUTDIR/bwa_mem2.bam" \
    "$METRIC_DIR/bwa_mem2.metrics.tsv"

add_summary \
    "BWA-MEM2" \
    "$METRIC_DIR/bwa_mem2.time.tsv" \
    "$METRIC_DIR/bwa_mem2.metrics.tsv"


# ============================================================
# Minimap2
# ============================================================

MMI="$INDEX_DIR/minimap2/ecoli.mmi"
check_file "$MMI"

echo "=== Minimap2 ==="

/usr/bin/time \
    -f "runtime_seconds\t%e\npeak_ram_kb\t%M" \
    -o "$METRIC_DIR/minimap2.time.tsv" \
bash -o pipefail -c "
minimap2 \
    -ax sr \
    -t '$THREADS' \
    '$MMI' \
    '$R1' \
    '$R2' \
    2> '$METRIC_DIR/minimap2.log' \
| samtools view \
    -@ '$THREADS' \
    -b \
    -o '$OUTDIR/minimap2.bam'
"

run_eval \
    "$OUTDIR/minimap2.bam" \
    "$METRIC_DIR/minimap2.metrics.tsv"

add_summary \
    "Minimap2" \
    "$METRIC_DIR/minimap2.time.tsv" \
    "$METRIC_DIR/minimap2.metrics.tsv"


# ============================================================
# HISAT2
# ============================================================

HISAT2_INDEX="$INDEX_DIR/hisat2/ecoli"
check_file "${HISAT2_INDEX}.1.ht2"

echo "=== HISAT2 ==="

/usr/bin/time \
    -f "runtime_seconds\t%e\npeak_ram_kb\t%M" \
    -o "$METRIC_DIR/hisat2.time.tsv" \
bash -o pipefail -c "
hisat2 \
    -x '$HISAT2_INDEX' \
    -1 '$R1' \
    -2 '$R2' \
    -p '$THREADS' \
    2> '$METRIC_DIR/hisat2.log' \
| samtools view \
    -@ '$THREADS' \
    -b \
    -o '$OUTDIR/hisat2.bam'
"

run_eval \
    "$OUTDIR/hisat2.bam" \
    "$METRIC_DIR/hisat2.metrics.tsv"

add_summary \
    "HISAT2" \
    "$METRIC_DIR/hisat2.time.tsv" \
    "$METRIC_DIR/hisat2.metrics.tsv"


# ============================================================
# STAR
# ============================================================

STAR_INDEX="$INDEX_DIR/star"
check_file "$STAR_INDEX/Genome"

echo "=== STAR ==="

STAR_OUT="$OUTDIR/star_"
rm -f "${STAR_OUT}"*

/usr/bin/time \
    -f "runtime_seconds\t%e\npeak_ram_kb\t%M" \
    -o "$METRIC_DIR/star.time.tsv" \
STAR \
    --runThreadN "$THREADS" \
    --genomeDir "$STAR_INDEX" \
    --readFilesIn "$R1" "$R2" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$STAR_OUT" \
    --outSAMtype BAM Unsorted \
    --outSAMattributes NH HI AS nM NM \
    > "$METRIC_DIR/star.stdout.log" \
    2> "$METRIC_DIR/star.stderr.log"

STAR_BAM="${STAR_OUT}Aligned.out.bam"
check_file "$STAR_BAM"

run_eval \
    "$STAR_BAM" \
    "$METRIC_DIR/star.metrics.tsv"

add_summary \
    "STAR" \
    "$METRIC_DIR/star.time.tsv" \
    "$METRIC_DIR/star.metrics.tsv"


# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------

echo
echo "=== Final summary ==="

if command -v column >/dev/null 2>&1; then
    column -t -s $'\t' "$SUMMARY"
else
    cat "$SUMMARY"
fi

echo
echo "Summary written to:"
echo "$SUMMARY"