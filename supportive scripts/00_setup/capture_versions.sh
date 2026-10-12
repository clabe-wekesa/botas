#!/usr/bin/env bash
set -u
OUT="${1:-software_versions.txt}"
{
  date -Is
  python3 --version 2>&1
  botas --version 2>&1 || true
  samtools --version 2>&1 | head -n2 || true
  bowtie2 --version 2>&1 | head -n1 || true
  bwa-mem2 version 2>&1 || true
  minimap2 --version 2>&1 || true
  featureCounts -v 2>&1 || true
  htseq-count --version 2>&1 || true
  sortmerna --version 2>&1 || true
  ribodetector_cpu --version 2>&1 || true
  java -version 2>&1 | head -n1 || true
} > "$OUT"
echo "Wrote $OUT"
