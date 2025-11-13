#!/usr/bin/env bash
set -euo pipefail

# Configuration - loacation of sequences
RAW_URL="https://owl.fish.washington.edu/nightingales/M_trossulus/"

# Local folders (adjust if needed)
RAW_DIR="../data/raw"
OUT_DIR="../output/qc_reports"

# Find fastqc and multiqc
FASTQC_BIN="${FASTQC_BIN:-$(command -v fastqc || true)}"
MULTIQC_BIN="${MULTIQC_BIN:-$(command -v multiqc || true)}"

# If not in PATH, search shared tool area
if [ -z "$FASTQC_BIN" ]; then
  FASTQC_BIN="$(find /home/shared -type f -iname fastqc -perm -111 2>/dev/null | head -n1 || true)"
fi
if [ -z "$MULTIQC_BIN" ]; then
  for p in \
    /home/shared/MultiQC*/multiqc \
    /home/shared/*/bin/multiqc \
    /home/shared/*/MultiQC*/multiqc
  do
    if [ -x "$p" ]; then
      MULTIQC_BIN="$p"
      break
    fi
  done
fi

# final broad search fallback
if [ -z "$MULTIQC_BIN" ]; then
  MULTIQC_BIN="$(find /home/shared -type f \( -iname 'multiqc' -o -iname 'multiqc*' \) -perm -111 2>/dev/null | head -n1 || true)"
fi

# guard rail
[ -n "$MULTIQC_BIN" ] && [ -x "$MULTIQC_BIN" ] || {
  echo "ERROR: Could not find an executable MultiQC. Set MULTIQC_BIN=/full/path/to/multiqc"
  exit 1
}


# Guard rails
if [ -z "$FASTQC_BIN" ] || [ ! -x "$FASTQC_BIN" ]; then
  echo "ERROR: Could not find an executable FastQC. Try adding it to PATH or set FASTQC_BIN=/full/path/to/fastqc"
  exit 1
fi

# Run FastQC
echo "Running FastQC..."
"$FASTQC_BIN" -o "$OUT_DIR" "$RAW_DIR"/*.fastq.gz
echo "FastQC complete."

# MultiQC
echo "Running MultiQC..."
"$MULTIQC_BIN" "$OUT_DIR" -o "$OUT_DIR"
echo "MultiQC report generated at $OUT_DIR"
