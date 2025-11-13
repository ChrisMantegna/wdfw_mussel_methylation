#!/usr/bin/env bash
set -euo pipefail

# Configuration - loacation of sequences
RAW_URL="https://owl.fish.washington.edu/nightingales/M_trossulus/"

# Local folders (adjust if needed)
RAW_DIR="../data/raw"
OUT_DIR="../output/qc_reports"

# Create Directories
mkdir -p "$RAW_DIR" "$OUT_DIR"

# Download Files
echo "Downloading"
wget --recursive \
  --no-check-certificate \
  --continue \
  --cut-dirs=3 \
  --no-host-directories \
  --no-parent \
  --quiet \
  --level=1 \
  --accept "*.fastq.gz,*.fastq.gz.md5sum" \
  -P "$RAW_DIR" \
  "${RAW_URL}"

echo "Downloads complete. Files are in $RAW_DIR"

# Run FastQC
echo "Running FastQC..."
/home/shared/ fastqc -o "$OUT_DIR" "$RAW_DIR"/*.fastq.gz
echo "FastQC complete."

# Run MultiQC
echo "Running MultiQC..."
multiqc "$OUT_DIR" -o "$OUT_DIR"
echo "MultiQC report generated at $OUT_DIR"

echo "Done! Check your FastQC and MultiQC reports."
