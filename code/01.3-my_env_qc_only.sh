#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# -----------------------------
# Configuration
# -----------------------------
RAW_DIR="../data/raw"
OUT_DIR="../output/qc_reports"

# Local venv-installed MultiQC
MULTIQC_BIN="../.multiqc_env/bin/multiqc"

# Optional: force rerun FastQC
FORCE_FASTQC="${FORCE_FASTQC:-0}"

# -----------------------------
# Make output folders
# -----------------------------
mkdir -p "$OUT_DIR"

# -----------------------------
# Find executable helper
# -----------------------------
find_executable() {
  local candidate
  for candidate in "$@"; do
    if [[ -n "${candidate:-}" && -f "$candidate" && -x "$candidate" ]]; then
      echo "$candidate"
      return 0
    fi
  done
  return 1
}

# -----------------------------
# Find FastQC
# -----------------------------
FASTQC_BIN="${FASTQC_BIN:-}"

if [[ -z "$FASTQC_BIN" ]]; then
  FASTQC_BIN="$(command -v fastqc 2>/dev/null || true)"
fi

FASTQC_BIN="$(find_executable \
  "$FASTQC_BIN" \
  "/home/shared/FastQC-0.12.1/fastqc" \
  "/home/shared/FastQC-0.11.9/fastqc" \
  "/home/shared/FastQC/fastqc" \
  "/home/shared/fastqc" \
  "/usr/local/bin/fastqc" \
  "/usr/bin/fastqc" \
)" || {
  echo "ERROR: Could not find a runnable FastQC executable."
  echo "Set it manually, for example:"
  echo "  export FASTQC_BIN=/full/path/to/fastqc"
  exit 1
}

echo "Using FastQC: $FASTQC_BIN"

# -----------------------------
# Check MultiQC
# -----------------------------
if [[ ! -x "$MULTIQC_BIN" ]]; then
  echo "ERROR: MultiQC not found or not executable at $MULTIQC_BIN"
  exit 1
fi

echo "Using MultiQC: $MULTIQC_BIN"

# -----------------------------
# Check FASTQ inputs
# -----------------------------
fastq_files=( "$RAW_DIR"/*.fastq.gz )

if (( ${#fastq_files[@]} == 0 )); then
  echo "ERROR: No FASTQ files found in $RAW_DIR"
  exit 1
fi

echo "Found ${#fastq_files[@]} FASTQ files."

# -----------------------------
# Run FastQC with skip feature
# -----------------------------
echo "Running FastQC with skip check..."

run_count=0
skip_count=0

for fq in "${fastq_files[@]}"; do
  base="$(basename "$fq")"
  sample="${base%.fastq.gz}"

  html_out="$OUT_DIR/${sample}_fastqc.html"
  zip_out="$OUT_DIR/${sample}_fastqc.zip"

  if [[ "$FORCE_FASTQC" -eq 0 && -f "$html_out" && -f "$zip_out" ]]; then
    echo "Skipping FastQC for $base (outputs already exist)"
    ((skip_count+=1))
    continue
  fi

  echo "Running FastQC on $base"
  "$FASTQC_BIN" -o "$OUT_DIR" "$fq"
  ((run_count+=1))
done

echo "FastQC complete."
echo "FastQC run on:  $run_count file(s)"
echo "FastQC skipped: $skip_count file(s)"

# -----------------------------
# Run MultiQC
# -----------------------------
echo "Running MultiQC..."
"$MULTIQC_BIN" "$OUT_DIR" -o "$OUT_DIR" --force

if [[ ! -f "$OUT_DIR/multiqc_report.html" ]]; then
  echo "ERROR: MultiQC finished, but no report was created."
  exit 1
fi

echo "MultiQC complete."
echo "Report created at: $OUT_DIR/multiqc_report.html"

# -----------------------------
# Done
# -----------------------------
echo "QC workflow complete."