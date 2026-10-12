#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Build simulated rRNA benchmark datasets
#
# For each genome:
#   50,000 true canonical-rRNA PE pairs
#   50,000 true non-rRNA PE pairs
#   100,000 total PE pairs
#
# Positive rRNA source:
#   5S, 16S and 23S ribosomal RNA only
#
# Negative source:
#   whole genome with ALL annotated rRNA features masked
#
# Output:
#   rrna_sim/<genome>/
# ============================================================


# ------------------------------------------------------------
# Command-line configuration
# ------------------------------------------------------------
usage() {
  cat <<'EOF'
Usage: build_rrna_benchmark.sh --manifest samples.tsv --outdir DIR [options]

Required manifest columns: sample, reference, gff
Options:
  --rrna-pairs N        true rRNA pairs per sample (default: 50000)
  --non-rrna-pairs N    true non-rRNA pairs per sample (default: 50000)
  --read-length N       read length (default: 150)
  --insert-mean N       fragment mean (default: 300)
  --insert-sd N         fragment SD (default: 50)
  --threads N           ART threads (default: 6)
EOF
}
MANIFEST=""; OUT_ROOT=""
TARGET_RRNA_PAIRS=50000; TARGET_NON_RRNA_PAIRS=50000
READ_LEN=150; INSERT_MEAN=300; INSERT_SD=50; THREADS=6
RRNA_FCOV_START=700; NON_RRNA_FCOV_START=5; MAX_ATTEMPTS=10
RRNA_ART_SEED=12345; NON_RRNA_ART_SEED=54321; SUBSAMPLE_SEED=12345; MIX_SEED=67890
while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2;;
    --outdir) OUT_ROOT="$2"; shift 2;;
    --rrna-pairs) TARGET_RRNA_PAIRS="$2"; shift 2;;
    --non-rrna-pairs) TARGET_NON_RRNA_PAIRS="$2"; shift 2;;
    --read-length) READ_LEN="$2"; shift 2;;
    --insert-mean) INSERT_MEAN="$2"; shift 2;;
    --insert-sd) INSERT_SD="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2;;
  esac
done
[[ -n "$MANIFEST" && -n "$OUT_ROOT" ]] || { usage >&2; exit 2; }
[[ -f "$MANIFEST" ]] || { echo "ERROR: manifest not found: $MANIFEST" >&2; exit 1; }


# ------------------------------------------------------------
# Check dependencies
# ------------------------------------------------------------

for cmd in art_modern bedtools python3 samtools; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: required command not found: $cmd" >&2
        exit 1
    fi
done


mkdir -p "$OUT_ROOT"


SUMMARY="${OUT_ROOT}/rrna_benchmark_summary.tsv"

printf \
"genome\tall_rRNA_features\tcanonical_rRNA_features\tmasked_rRNA_bp\tcanonical_rRNA_bp\ttrue_rRNA_pairs\ttrue_non_rRNA_pairs\ttotal_pairs\n" \
> "$SUMMARY"


# ============================================================
# Helper: count pairs in an interleaved FASTQ
# ============================================================

count_pairs() {
    local fq="$1"

    if [[ ! -s "$fq" ]]; then
        echo 0
        return
    fi

    local lines
    lines=$(wc -l < "$fq")

    echo $(( lines / 8 ))
}


# ============================================================
# Helper: simulate until at least requested number of pairs
# ============================================================

simulate_until_enough() {

    local source_fasta="$1"
    local output_fastq="$2"
    local start_fcov="$3"
    local target_pairs="$4"
    local art_seed="$5"
    local label="$6"

    local fcov="$start_fcov"
    local attempt=1

    while (( attempt <= MAX_ATTEMPTS )); do

        echo
        echo "  [$label] attempt=$attempt fcov=$fcov"

        rm -f "$output_fastq"

        art_modern \
            --mode wgs \
            --lc pe \
            --i-file "$source_fasta" \
            --i-fcov "$fcov" \
            --i-seed "$art_seed" \
            --builtin_qual_file HiSeq2500_150bp \
            --read_len "$READ_LEN" \
            --pe_frag_dist_mean "$INSERT_MEAN" \
            --pe_frag_dist_std_dev "$INSERT_SD" \
            --o-fastq "$output_fastq" \
            --parallel "$THREADS"

        local observed
        observed=$(count_pairs "$output_fastq")

        echo "  [$label] generated_pairs=$observed"

        if (( observed >= target_pairs )); then
            return 0
        fi

        fcov=$(( fcov * 2 ))
        attempt=$(( attempt + 1 ))
    done

    echo "ERROR: could not generate ${target_pairs} pairs from ${source_fasta}" >&2
    exit 1
}


# ============================================================
# Process genomes
# ============================================================

while IFS=$'\t' read -r GENOME REF GFF; do
    [[ -n "${GENOME:-}" ]] || continue
    [[ -f "$REF" ]] || { echo "ERROR: reference not found for $GENOME: $REF" >&2; exit 1; }
    [[ -f "$GFF" ]] || { echo "ERROR: GFF not found for $GENOME: $GFF" >&2; exit 1; }

    OUT="${OUT_ROOT}/${GENOME}"

    echo
    echo "============================================================"
    echo "Genome: $GENOME"
    echo "============================================================"

    # Start clean so old benchmark files cannot contaminate
    # the regenerated dataset.
    rm -rf "$OUT"
    mkdir -p "$OUT"


    # --------------------------------------------------------
    # File names
    # --------------------------------------------------------

    ALL_RRNA_BED="${OUT}/${GENOME}.all_rRNA.bed"
    CANONICAL_BED="${OUT}/${GENOME}.rRNA.bed"

    RRNA_FASTA="${OUT}/${GENOME}.rRNA.fasta"
    NON_RRNA_FASTA="${OUT}/${GENOME}.non_rRNA_masked.fasta"

    RRNA_RAW="${OUT}/rRNA.raw.fastq"
    NON_RAW="${OUT}/non_rRNA.raw.fastq"

    RRNA_50K="${OUT}/rRNA_50k_pairs.fastq"
    NON_50K="${OUT}/non_rRNA_50k_pairs.fastq"

    MIXED="${OUT}/${GENOME}_rrna_benchmark_100k.fastq"

    FINAL_R1="${OUT}/${GENOME}_rrna_benchmark_100k_R1.fastq"
    FINAL_R2="${OUT}/${GENOME}_rrna_benchmark_100k_R2.fastq"

    TRUTH="${OUT}/${GENOME}_rrna_benchmark_truth.tsv"


    # ========================================================
    # 1. Extract ALL annotated rRNA features
    #
    # These are used for masking the negative genome.
    # ========================================================

    awk -F'\t' '
    BEGIN {OFS="\t"}

    $0 !~ /^#/ && $3=="rRNA" {
        print $1, $4-1, $5, $9, ".", $7
    }
    ' "$GFF" > "$ALL_RRNA_BED"


    ALL_FEATURES=$(wc -l < "$ALL_RRNA_BED")

    if (( ALL_FEATURES == 0 )); then
        echo "ERROR: no rRNA annotations found for $GENOME" >&2
        exit 1
    fi


    # ========================================================
    # 2. Extract canonical bacterial rRNAs only
    #
    # Positive class:
    #   5S
    #   16S
    #   23S
    #
    # This deliberately excludes 4.5S and other unusual
    # features from the positive simulation source.
    # ========================================================

    awk -F'\t' '
    BEGIN {OFS="\t"}

    $0 !~ /^#/ &&
    $3=="rRNA" &&
    (
        $9 ~ /product=5S ribosomal RNA/  ||
        $9 ~ /product=16S ribosomal RNA/ ||
        $9 ~ /product=23S ribosomal RNA/
    )
    {
        print $1, $4-1, $5, $9, ".", $7
    }
    ' "$GFF" > "$CANONICAL_BED"


    CANONICAL_FEATURES=$(wc -l < "$CANONICAL_BED")

    if (( CANONICAL_FEATURES == 0 )); then
        echo "ERROR: no canonical 5S/16S/23S rRNA annotations found for $GENOME" >&2
        exit 1
    fi


    echo "  all rRNA features       : $ALL_FEATURES"
    echo "  canonical rRNA features : $CANONICAL_FEATURES"


    # ========================================================
    # 3. Calculate annotated rRNA lengths
    # ========================================================

    MASKED_BP=$(awk '
        {sum += $3-$2}
        END {print sum+0}
    ' "$ALL_RRNA_BED")


    CANONICAL_BP=$(awk '
        {sum += $3-$2}
        END {print sum+0}
    ' "$CANONICAL_BED")


    echo "  all rRNA bp             : $MASKED_BP"
    echo "  canonical rRNA bp       : $CANONICAL_BP"


    # ========================================================
    # 4. Extract canonical rRNA sequences
    # ========================================================

    bedtools getfasta \
        -fi "$REF" \
        -bed "$CANONICAL_BED" \
        -fo "$RRNA_FASTA" \
        -s


    # ========================================================
    # 5. Mask ALL annotated rRNA from negative source
    #
    # This prevents 4.5S or other annotated rRNA sequences
    # from accidentally becoming "non-rRNA" truth reads.
    # ========================================================

    bedtools maskfasta \
        -fi "$REF" \
        -bed "$ALL_RRNA_BED" \
        -fo "$NON_RRNA_FASTA"


    # ========================================================
    # 6. Simulate canonical-rRNA pairs
    # ========================================================

    simulate_until_enough \
        "$RRNA_FASTA" \
        "$RRNA_RAW" \
        "$RRNA_FCOV_START" \
        "$TARGET_RRNA_PAIRS" \
        "$RRNA_ART_SEED" \
        "rRNA"


    # ========================================================
    # 7. Simulate non-rRNA pairs
    #
    # ART default max_n=0 prevents fragments containing masked
    # N positions from being emitted.
    # ========================================================

    simulate_until_enough \
        "$NON_RRNA_FASTA" \
        "$NON_RAW" \
        "$NON_RRNA_FCOV_START" \
        "$TARGET_NON_RRNA_PAIRS" \
        "$NON_RRNA_ART_SEED" \
        "non-rRNA"


    # ========================================================
    # 8. Pair-aware subsampling
    #
    # ART PE output is interleaved:
    #   4 lines mate 1
    #   4 lines mate 2
    #
    # Therefore one pair = 8 FASTQ lines.
    # ========================================================

    python3 - \
        "$RRNA_RAW" \
        "$RRNA_50K" \
        "$TARGET_RRNA_PAIRS" \
        "$SUBSAMPLE_SEED" <<'PY'
import sys
import random

src = sys.argv[1]
dst = sys.argv[2]
target = int(sys.argv[3])
seed = int(sys.argv[4])

pairs = []

with open(src) as f:
    while True:
        block = [f.readline() for _ in range(8)]

        if not block[0]:
            break

        if any(x == "" for x in block):
            raise RuntimeError(
                f"Incomplete paired FASTQ block in {src}"
            )

        pairs.append(block)

if len(pairs) < target:
    raise RuntimeError(
        f"{src}: only {len(pairs)} pairs available; "
        f"{target} requested"
    )

rng = random.Random(seed)
indices = rng.sample(range(len(pairs)), target)

with open(dst, "w") as out:
    for i in indices:
        out.writelines(pairs[i])

print(f"Subsampled {target} pairs from {len(pairs)}")
PY


    python3 - \
        "$NON_RAW" \
        "$NON_50K" \
        "$TARGET_NON_RRNA_PAIRS" \
        "$((SUBSAMPLE_SEED + 1))" <<'PY'
import sys
import random

src = sys.argv[1]
dst = sys.argv[2]
target = int(sys.argv[3])
seed = int(sys.argv[4])

pairs = []

with open(src) as f:
    while True:
        block = [f.readline() for _ in range(8)]

        if not block[0]:
            break

        if any(x == "" for x in block):
            raise RuntimeError(
                f"Incomplete paired FASTQ block in {src}"
            )

        pairs.append(block)

if len(pairs) < target:
    raise RuntimeError(
        f"{src}: only {len(pairs)} pairs available; "
        f"{target} requested"
    )

rng = random.Random(seed)
indices = rng.sample(range(len(pairs)), target)

with open(dst, "w") as out:
    for i in indices:
        out.writelines(pairs[i])

print(f"Subsampled {target} pairs from {len(pairs)}")
PY


    # ========================================================
    # 9. Mix the two truth classes and rename pair IDs
    #
    # Generates:
    #   interleaved FASTQ
    #   R1 FASTQ
    #   R2 FASTQ
    #   pair-level truth table
    # ========================================================

    python3 - \
        "$GENOME" \
        "$RRNA_50K" \
        "$NON_50K" \
        "$MIXED" \
        "$FINAL_R1" \
        "$FINAL_R2" \
        "$TRUTH" \
        "$MIX_SEED" <<'PY'

import sys
import random

genome = sys.argv[1]

rrna_file = sys.argv[2]
non_file = sys.argv[3]

mixed_file = sys.argv[4]
r1_file = sys.argv[5]
r2_file = sys.argv[6]
truth_file = sys.argv[7]

seed = int(sys.argv[8])


def read_pairs(path, truth_class):

    result = []

    with open(path) as f:

        while True:

            r1 = [f.readline() for _ in range(4)]

            if not r1[0]:
                break

            r2 = [f.readline() for _ in range(4)]

            if any(x == "" for x in r1 + r2):
                raise RuntimeError(
                    f"Incomplete PE FASTQ block: {path}"
                )

            result.append(
                (truth_class, r1, r2)
            )

    return result


pairs = []

pairs.extend(
    read_pairs(rrna_file, "rRNA")
)

pairs.extend(
    read_pairs(non_file, "non_rRNA")
)


rng = random.Random(seed)
rng.shuffle(pairs)


with (
    open(mixed_file, "w") as mixed,
    open(r1_file, "w") as out1,
    open(r2_file, "w") as out2,
    open(truth_file, "w") as truth
):

    truth.write(
        "pair_id\ttruth_class\n"
    )

    for i, (klass, read1, read2) in enumerate(
        pairs,
        start=1
    ):

        pair_id = (
            f"{genome}_rrnaBench_{i:06d}"
        )

        read1[0] = f"@{pair_id}/1\n"
        read2[0] = f"@{pair_id}/2\n"

        mixed.writelines(read1)
        mixed.writelines(read2)

        out1.writelines(read1)
        out2.writelines(read2)

        truth.write(
            f"{pair_id}\t{klass}\n"
        )


print(
    f"Final benchmark: {len(pairs)} pairs"
)
PY


    # ========================================================
    # 10. Validate final counts
    # ========================================================

    FINAL_R1_COUNT=$(( $(wc -l < "$FINAL_R1") / 4 ))
    FINAL_R2_COUNT=$(( $(wc -l < "$FINAL_R2") / 4 ))

    TRUE_RRNA=$(
        awk -F'\t' '
        NR>1 && $2=="rRNA" {n++}
        END {print n+0}
        ' "$TRUTH"
    )

    TRUE_NON=$(
        awk -F'\t' '
        NR>1 && $2=="non_rRNA" {n++}
        END {print n+0}
        ' "$TRUTH"
    )


    if (( FINAL_R1_COUNT != TARGET_RRNA_PAIRS + TARGET_NON_RRNA_PAIRS )); then
        echo "ERROR: $GENOME R1 has $FINAL_R1_COUNT pairs" >&2
        exit 1
    fi

    if (( FINAL_R2_COUNT != TARGET_RRNA_PAIRS + TARGET_NON_RRNA_PAIRS )); then
        echo "ERROR: $GENOME R2 has $FINAL_R2_COUNT pairs" >&2
        exit 1
    fi

    if (( TRUE_RRNA != TARGET_RRNA_PAIRS )); then
        echo "ERROR: $GENOME rRNA truth count=$TRUE_RRNA" >&2
        exit 1
    fi

    if (( TRUE_NON != TARGET_NON_RRNA_PAIRS )); then
        echo "ERROR: $GENOME non-rRNA truth count=$TRUE_NON" >&2
        exit 1
    fi


    TOTAL=$(( TRUE_RRNA + TRUE_NON ))


    # ========================================================
    # 11. Summary table
    # ========================================================

    printf \
    "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" \
    "$GENOME" \
    "$ALL_FEATURES" \
    "$CANONICAL_FEATURES" \
    "$MASKED_BP" \
    "$CANONICAL_BP" \
    "$TRUE_RRNA" \
    "$TRUE_NON" \
    "$TOTAL" \
    >> "$SUMMARY"


    echo
    echo "Completed: $GENOME"
    echo "  true rRNA pairs     : $TRUE_RRNA"
    echo "  true non-rRNA pairs : $TRUE_NON"
    echo "  total pairs         : $TOTAL"

done < <(python3 - "$MANIFEST" <<'PY_MANIFEST'
import csv,sys
with open(sys.argv[1],newline='') as fh:
    r=csv.DictReader(fh,delimiter='\t')
    req={'sample','reference','gff'}
    if not req.issubset(r.fieldnames or []):
        raise SystemExit('manifest must contain columns: sample, reference, gff')
    for row in r:
        if row['sample'].strip():
            print(row['sample'].strip(), row['reference'].strip(), row['gff'].strip(), sep='\t')
PY_MANIFEST
)


echo
echo "============================================================"
echo "All benchmarks completed"
echo "============================================================"
echo
column -t -s $'\t' "$SUMMARY"
