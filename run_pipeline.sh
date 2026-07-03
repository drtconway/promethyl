#!/usr/bin/env bash
set -euo pipefail

# ---------------------------
# Promethyl pipeline runner
# Option A: CLI wrapper with sample.yml input
# ---------------------------

# Config file passed as first argument (default: sample.yml)
CONFIG="${1:-sample.yml}"

# Project root (location of this script)
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Paths inside repo
MODKIT_DIR="${ROOT_DIR}/output/modkit"
ANNOTATION="${ROOT_DIR}/reference/CpGs_with_promoters.bed"
INCLUDE_BED="${ROOT_DIR}/reference/CPGIslandsBED3.bed"
OUTPUT="${ROOT_DIR}/output/probands_modkit_cohort_methylation.tsv"

# External reference genome (NOT in repo)
REF="/path/to/reference/hg38.fa"

# Threads (fixed for Option A)
THREADS=10

echo "==================================="
echo "Triomethyl pipeline starting"
echo "Config: $CONFIG"
echo "Repo root: $ROOT_DIR"
echo "Threads: $THREADS"
echo "==================================="

# Ensure output directories exist
mkdir -p "${ROOT_DIR}/output/modkit"

# Run pipeline
python "${ROOT_DIR}/src/main.py" \
    --config "$CONFIG" \
    --modkit-dir "$MODKIT_DIR" \
    --annotation "$ANNOTATION" \
    --ref "$REF" \
    --output "$OUTPUT" \
    --include-bed "$INCLUDE_BED" \
    --threads "$THREADS"

echo "==================================="
echo "Pipeline finished successfully"
echo "==================================="