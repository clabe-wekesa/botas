#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Unified multi-genome, multi-aligner benchmark
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
THREADS=6; TOTAL_PAIRS=1000000; JUNCTION_PAIRS=100000; ITERATIONS=10; READLEN=150; SEED=12345; TOLERANCE=10; OVERHANG_PERCENT=5
REF_DIR=""; OUT_ROOT=""; SIMULATOR="$SCRIPT_DIR/readSimulator_truth_1M.py"; EVALUATOR="$SCRIPT_DIR/compare_aligners.py"; CIRCULAR_MODE="all"
usage(){ cat <<'EOF'
Usage: run_all_alignment_benchmarks.sh --reference-dir DIR --outdir DIR [options]
  --threads N
  --total-pairs N
  --junction-pairs N
  --iterations N
  --read-length N
  --seed N
  --tolerance N
  --overhang-percent N
  --circular all|none   BOTAS circular treatment (default: all)
  --simulator FILE
  --evaluator FILE
EOF
}
while [[ $# -gt 0 ]]; do
 case "$1" in
  --reference-dir) REF_DIR="$2";shift 2;;
  --outdir) OUT_ROOT="$2";shift 2;;
  --threads) THREADS="$2";shift 2;;
  --total-pairs) TOTAL_PAIRS="$2";shift 2;;
  --junction-pairs) JUNCTION_PAIRS="$2";shift 2;;
  --iterations) ITERATIONS="$2";shift 2;;
  --read-length) READLEN="$2";shift 2;;
  --seed) SEED="$2";shift 2;;
  --tolerance) TOLERANCE="$2";shift 2;;
  --overhang-percent) OVERHANG_PERCENT="$2";shift 2;;
  --circular) CIRCULAR_MODE="$2";shift 2;;
  --simulator) SIMULATOR="$2";shift 2;;
  --evaluator) EVALUATOR="$2";shift 2;;
  -h|--help) usage;exit 0;;
  *) echo "Unknown option: $1" >&2;usage >&2;exit 2;;
 esac
done
[[ -n "$REF_DIR" && -n "$OUT_ROOT" ]] || { usage >&2; exit 2; }

mkdir -p "$OUT_ROOT"

MASTER_METRICS="$OUT_ROOT/master_alignment_metrics.tsv"
MASTER_RESOURCES="$OUT_ROOT/master_resource_usage.tsv"
MASTER_FLAGSTAT="$OUT_ROOT/master_flagstat.tsv"

printf 'genome\taligner\tcategory\ttotal\tcorrect\twrong\tunmapped\tprecision\trecall\tf1\tjunction_tags\n' > "$MASTER_METRICS"
printf 'genome\taligner\truntime_seconds\tpeak_ram_kb\tpeak_ram_mb\n' > "$MASTER_RESOURCES"
printf 'genome\taligner\ttotal_reads\tmapped_reads\tproperly_paired\tsingletons\n' > "$MASTER_FLAGSTAT"

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || {
        echo "ERROR: required command not found: $1" >&2
        exit 1
    }
}

for cmd in python3 samtools bowtie2 bowtie2-build bwa-mem2 minimap2 botas /usr/bin/time; do
    require_cmd "$cmd"
done

[[ -f "$SIMULATOR" ]] || { echo "ERROR: simulator not found: $SIMULATOR" >&2; exit 1; }
[[ -f "$EVALUATOR" ]] || { echo "ERROR: evaluator not found: $EVALUATOR" >&2; exit 1; }
[[ -d "$REF_DIR" ]] || { echo "ERROR: reference directory not found: $REF_DIR" >&2; exit 1; }

mapfile -t REFS < <(find "$REF_DIR" -type f \( -name '*.fna' -o -name '*.fa' -o -name '*.fasta' -o -name '*.fna.gz' -o -name '*.fa.gz' -o -name '*.fasta.gz' \) | sort)
if (( ${#REFS[@]} == 0 )); then
    echo "ERROR: no reference FASTA files found under $REF_DIR" >&2
    exit 1
fi

echo "Found ${#REFS[@]} reference file(s)."
echo "Output root: $OUT_ROOT"

run_timed() {
    local label="$1"
    local timefile="$2"
    shift 2

    echo "[$label] starting"
    /usr/bin/time \
        -f $'runtime_seconds\t%e\npeak_ram_kb\t%M' \
        -o "$timefile" \
        "$@"
    echo "[$label] done"
}

append_resources() {
    local genome="$1"
    local aligner="$2"
    local timefile="$3"

    local runtime ram ram_mb
    runtime=$(awk -F'\t' '$1=="runtime_seconds"{print $2}' "$timefile")
    ram=$(awk -F'\t' '$1=="peak_ram_kb"{print $2}' "$timefile")
    ram_mb=$(awk -v kb="$ram" 'BEGIN{printf "%.2f", kb/1024}')

    printf '%s\t%s\t%s\t%s\t%s\n' \
        "$genome" "$aligner" "$runtime" "$ram" "$ram_mb" \
        >> "$MASTER_RESOURCES"
}

append_flagstat() {
    local genome="$1"
    local aligner="$2"
    local bam="$3"
    local outfile="$4"

    samtools flagstat "$bam" > "$outfile"

    local total mapped proper singletons
    total=$(awk '$3=="in" && $4=="total" {print $1; exit}' "$outfile")
    mapped=$(awk '$3=="mapped" {print $1; exit}' "$outfile")
    proper=$(awk '$3=="properly" && $4=="paired" {print $1; exit}' "$outfile")
    singletons=$(awk '$3=="singletons" {print $1; exit}' "$outfile")

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$genome" "$aligner" "${total:-NA}" "${mapped:-NA}" "${proper:-NA}" "${singletons:-NA}" \
        >> "$MASTER_FLAGSTAT"
}

for REF_SRC in "${REFS[@]}"; do
    base=$(basename "$REF_SRC")
    base=${base%.gz}
    stem=${base%.*}
    genome="$stem"

    WORK="$OUT_ROOT/$genome"
    REF="$WORK/reference.fna"
    DATA="$WORK/data"
    IDX="$WORK/index"
    ALN="$WORK/alignment"
    MET="$WORK/metrics"

    mkdir -p "$DATA" "$IDX" "$ALN" "$MET"

    echo
    echo "======================================================================"
    echo "GENOME: $genome"
    echo "======================================================================"

    # ------------------------------------------------------------
    # Reference
    # ------------------------------------------------------------
    if [[ "$REF_SRC" == *.gz ]]; then
        gzip -cd "$REF_SRC" > "$REF"
    else
        cp "$REF_SRC" "$REF"
    fi

    # Optionally treat every FASTA record as circular for BOTAS.
    # Use --circular none for linear references.
    CIRCULAR_CONTIGS=""
    if [[ "$CIRCULAR_MODE" == "all" ]]; then
        CIRCULAR_CONTIGS=$(awk '/^>/{sub(/^>/,""); split($0,a,/[^[:graph:]]/); print a[1]}' "$REF" | paste -sd, -)
        [[ -n "$CIRCULAR_CONTIGS" ]] || { echo "ERROR: no FASTA sequence identifiers found in $REF" >&2; exit 1; }
    elif [[ "$CIRCULAR_MODE" != "none" ]]; then
        echo "ERROR: --circular must be all or none" >&2; exit 2
    fi

    # ------------------------------------------------------------
    # Simulation
    # ------------------------------------------------------------
    R1="$DATA/${genome}_1.fastq.gz"
    R2="$DATA/${genome}_2.fastq.gz"
    TRUTH="$DATA/${genome}_truth.tsv"

    if [[ ! -s "$R1" || ! -s "$R2" || ! -s "$TRUTH" ]]; then
        rm -rf "$DATA"
        mkdir -p "$DATA"

        python3 "$SIMULATOR" \
            --input "$REF" \
            --outdir "$DATA" \
            --total-pairs "$TOTAL_PAIRS" \
            --junction-pairs "$JUNCTION_PAIRS" \
            --iterations "$ITERATIONS" \
            --readlen "$READLEN" \
            --seed "$SEED"

        # Simulator names outputs from input stem (reference), not genome stem.
        # Normalize the names once so the rest of this driver stays generic.
        sim_r1=$(find "$DATA" -maxdepth 1 -name '*_1.fastq.gz' | head -1)
        sim_r2=$(find "$DATA" -maxdepth 1 -name '*_2.fastq.gz' | head -1)
        sim_truth=$(find "$DATA" -maxdepth 1 -name '*_truth.tsv' | head -1)
        sim_rot=$(find "$DATA" -maxdepth 1 -name '*_rotations.tsv' | head -1 || true)

        mv "$sim_r1" "$R1"
        mv "$sim_r2" "$R2"
        mv "$sim_truth" "$TRUTH"
        [[ -n "${sim_rot:-}" && -f "$sim_rot" ]] && mv "$sim_rot" "$DATA/${genome}_rotations.tsv"
    else
        echo "[simulation] existing dataset found; reusing it"
    fi

    # ------------------------------------------------------------
    # BOTAS index + alignment
    # ------------------------------------------------------------
    BOTAS_PREFIX="$IDX/botas/index"
    BOTAS_INDEX="${BOTAS_PREFIX}.botas/results/index.botas.idx"
    BOTAS_BAM="$ALN/botas.bam"

    if [[ ! -s "$BOTAS_INDEX" ]]; then
        mkdir -p "$IDX/botas"
        if [[ -n "$CIRCULAR_CONTIGS" ]]; then
            botas index -r "$REF" --circular-contigs "$CIRCULAR_CONTIGS" --circular-overhang-percent "$OVERHANG_PERCENT" -o "$BOTAS_PREFIX"
        else
            botas index -r "$REF" -o "$BOTAS_PREFIX"
        fi
    fi

    if [[ ! -s "$BOTAS_BAM" ]]; then
        rm -rf "$ALN/botas.botas"
        CIRC_ARGS=()
        [[ -n "$CIRCULAR_CONTIGS" ]] && CIRC_ARGS=(--circular-contigs "$CIRCULAR_CONTIGS")
        run_timed "BOTAS $genome" "$MET/botas.time.tsv" \
            botas align \
                --ref "$REF" \
                --index "$BOTAS_INDEX" \
                --fq1 "$R1" \
                --fq2 "$R2" \
                "${CIRC_ARGS[@]}" \
                --pool \
                -t "$THREADS" \
                -o "$ALN/botas.bam"

        produced=$(find "$ALN" -type f -path '*/botas.botas/results/*.bam' | head -1 || true)
        if [[ -n "${produced:-}" && "$produced" != "$BOTAS_BAM" ]]; then
            cp "$produced" "$BOTAS_BAM"
        fi
    else
        echo "[BOTAS] existing BAM found; reusing it"
    fi

    # ------------------------------------------------------------
    # Bowtie2
    # ------------------------------------------------------------
    BT2_PREFIX="$IDX/bowtie2/index"
    BT2_BAM="$ALN/bowtie2.bam"

    if [[ ! -s "${BT2_PREFIX}.1.bt2" && ! -s "${BT2_PREFIX}.1.bt2l" ]]; then
        mkdir -p "$IDX/bowtie2"
        bowtie2-build "$REF" "$BT2_PREFIX" > "$MET/bowtie2-build.log" 2>&1
    fi

    if [[ ! -s "$BT2_BAM" ]]; then
        run_timed "Bowtie2 $genome" "$MET/bowtie2.time.tsv" \
            bash -o pipefail -c \
            "bowtie2 -x '$BT2_PREFIX' -1 '$R1' -2 '$R2' -p '$THREADS' 2> '$MET/bowtie2.log' | samtools view -@ '$THREADS' -b -o '$BT2_BAM'"
    else
        echo "[Bowtie2] existing BAM found; reusing it"
    fi

    # ------------------------------------------------------------
    # BWA-MEM2
    # ------------------------------------------------------------
    BWA_REF="$IDX/bwa/reference.fna"
    BWA_BAM="$ALN/bwa_mem2.bam"

    if [[ ! -s "$BWA_REF" ]]; then
        mkdir -p "$IDX/bwa"
        cp "$REF" "$BWA_REF"
    fi

    if [[ ! -e "${BWA_REF}.0123" && ! -e "${BWA_REF}.bwt.2bit.64" ]]; then
        bwa-mem2 index "$BWA_REF" > "$MET/bwa_mem2-index.log" 2>&1
    fi

    if [[ ! -s "$BWA_BAM" ]]; then
        run_timed "BWA-MEM2 $genome" "$MET/bwa_mem2.time.tsv" \
            bash -o pipefail -c \
            "bwa-mem2 mem -t '$THREADS' '$BWA_REF' '$R1' '$R2' 2> '$MET/bwa_mem2.log' | samtools view -@ '$THREADS' -b -o '$BWA_BAM'"
    else
        echo "[BWA-MEM2] existing BAM found; reusing it"
    fi

    # ------------------------------------------------------------
    # Minimap2
    # ------------------------------------------------------------
    MMI="$IDX/minimap2/reference.mmi"
    MM2_BAM="$ALN/minimap2.bam"

    if [[ ! -s "$MMI" ]]; then
        mkdir -p "$IDX/minimap2"
        minimap2 -d "$MMI" "$REF" > "$MET/minimap2-index.log" 2>&1
    fi

    if [[ ! -s "$MM2_BAM" ]]; then
        run_timed "Minimap2 $genome" "$MET/minimap2.time.tsv" \
            bash -o pipefail -c \
            "minimap2 -ax sr -t '$THREADS' '$MMI' '$R1' '$R2' 2> '$MET/minimap2.log' | samtools view -@ '$THREADS' -b -o '$MM2_BAM'"
    else
        echo "[Minimap2] existing BAM found; reusing it"
    fi

    # ------------------------------------------------------------
    # Flagstat + resource summaries
    # ------------------------------------------------------------
    for spec in \
        "BOTAS|$BOTAS_BAM|$MET/botas.time.tsv" \
        "Bowtie2|$BT2_BAM|$MET/bowtie2.time.tsv" \
        "BWA-MEM2|$BWA_BAM|$MET/bwa_mem2.time.tsv" \
        "Minimap2|$MM2_BAM|$MET/minimap2.time.tsv"
    do
        IFS='|' read -r label bam timefile <<< "$spec"
        append_flagstat "$genome" "$label" "$bam" "$MET/${label// /_}.flagstat.txt"
        [[ -s "$timefile" ]] && append_resources "$genome" "$label" "$timefile"
    done

    # ------------------------------------------------------------
    # Accuracy evaluation
    # ------------------------------------------------------------
    python3 "$EVALUATOR" \
        --truth "$TRUTH" \
        --tolerance "$TOLERANCE" \
        --bam BOTAS "$BOTAS_BAM" \
        --bam Bowtie2 "$BT2_BAM" \
        --bam BWA-MEM2 "$BWA_BAM" \
        --bam Minimap2 "$MM2_BAM" \
        --tsv "$MET/alignment_metrics.tsv" \
        > "$MET/alignment_metrics.txt"

    # Generic evaluator TSV has one row per aligner/category. Prefix genome.
    awk -v genome="$genome" 'BEGIN{FS=OFS="\t"} NR==1{next} {print genome,$0}' \
        "$MET/alignment_metrics.tsv" >> "$MASTER_METRICS"

    echo "[$genome] complete"
done

echo
echo "======================================================================"
echo "ALL BENCHMARKS COMPLETE"
echo "======================================================================"
echo "Accuracy summary:  $MASTER_METRICS"
echo "Resource summary:  $MASTER_RESOURCES"
echo "Flagstat summary:  $MASTER_FLAGSTAT"
