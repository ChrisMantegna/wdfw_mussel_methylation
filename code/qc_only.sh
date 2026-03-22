#!/usr/bin/env bash
set -euo pipefail

# Configuration
RAW_URL="https://owl.fish.washington.edu/nightingales/M_trossulus/"
RAW_DIR="../data/raw"
OUT_DIR="../output/qc_reports"

mkdir -p "$OUT_DIR"

# Locate FastQC
FASTQC_BIN="${FASTQC_BIN:-$(command -v fastqc || true)}"
if [ -z "$FASTQC_BIN" ]; then
  FASTQC_BIN="$(find /home/shared -type f -iname fastqc -perm -111 2>/dev/null | head -n1 || true)"
fi

# Locate MultiQC
MULTIQC_BIN="${MULTIQC_BIN:-$(command -v multiqc || true)}"
if [ -z "$MULTIQC_BIN" ]; then
  for p in \
    /home/shared/MultiQC*/bin/multiqc \
    /home/shared/*/bin/multiqc \
    /home/shared/*/MultiQC*/bin/multiqc \
    /home/shared/*/MultiQC*/multiqc
  do
    if [ -f "$p" ] && [ -x "$p" ]; then
      MULTIQC_BIN="$p"
      break
    fi
  done
fi

# Broad fallback search for MultiQC
if [ -z "$MULTIQC_BIN" ]; then
  MULTIQC_BIN="$(find /home/shared -type f \( -iname 'multiqc' -o -iname 'multiqc*' \) -perm -111 2>/dev/null | head -n1 || true)"
fi

# Fail clearly if tools were not found
if [ -z "$FASTQC_BIN" ]; then
  echo "Error: fastqc not found in PATH or /home/shared" >&2
  exit 1
fi

if [ -z "$MULTIQC_BIN" ]; then
  echo "Error: multiqc not found in PATH or /home/shared" >&2
  exit 1
fi

echo "Using FastQC:  $FASTQC_BIN"
echo "Using MultiQC: $MULTIQC_BIN"

# Gather input files safely
shopt -s nullglob
FASTQ_FILES=("$RAW_DIR"/*.fastq.gz)

if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
  echo "Error: no .fastq.gz files found in $RAW_DIR" >&2
  exit 1
fi

# Run FastQC first
echo "Running FastQC..."
"$FASTQC_BIN" -o "$OUT_DIR" "${FASTQ_FILES[@]}"
echo "FastQC complete."

# Run MultiQC after FastQC
echo "Running MultiQC..."
"$MULTIQC_BIN" "$OUT_DIR" -o "$OUT_DIR"
echo "MultiQC report generated at $OUT_DIR"