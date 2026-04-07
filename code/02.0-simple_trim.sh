#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# -----------------------------
# Paths
# -----------------------------
RAW_DIR="../data/raw"
TRIM_ROOT="../data/trimmed"
FASTP_DIR="$TRIM_ROOT/fastp"
FASTQC_DIR="$TRIM_ROOT/fastqc"
MULTIQC_DIR="$TRIM_ROOT/multiqc"

mkdir -p "$TRIM_ROOT" "$FASTP_DIR" "$FASTQC_DIR"

# Set to 1 only if/when you want MultiQC to run
RUN_MULTIQC="${RUN_MULTIQC:-1}"

# Local venv-installed MultiQC
# Assumes script is run from project code/ directory
MULTIQC_BIN_DEFAULT="../.multiqc_env/bin/multiqc"

# -----------------------------
# Helper: find executable
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
# Find tools
# -----------------------------
FASTP_BIN="${FASTP_BIN:-$(command -v fastp || true)}"
FASTQC_BIN="${FASTQC_BIN:-$(command -v fastqc || true)}"
MULTIQC_BIN="${MULTIQC_BIN:-$MULTIQC_BIN_DEFAULT}"

FASTP_BIN="$(find_executable \
  "$FASTP_BIN" \
  "/home/shared/fastp" \
  "/home/shared/bin/fastp" \
  "/usr/local/bin/fastp" \
  "/usr/bin/fastp" \
)" || {
  echo "ERROR: fastp executable not found."
  echo "Set it manually, for example:"
  echo "  export FASTP_BIN=/full/path/to/fastp"
  exit 1
}

FASTQC_BIN="$(find_executable \
  "$FASTQC_BIN" \
  "/home/shared/FastQC-0.12.1/fastqc" \
  "/home/shared/FastQC-0.11.9/fastqc" \
  "/home/shared/FastQC/fastqc" \
  "/home/shared/fastqc" \
  "/usr/local/bin/fastqc" \
  "/usr/bin/fastqc" \
)" || {
  echo "ERROR: fastqc executable not found."
  echo "Set it manually, for example:"
  echo "  export FASTQC_BIN=/full/path/to/fastqc"
  exit 1
}

if [[ "$RUN_MULTIQC" == "1" ]]; then
  if [[ ! -f "$MULTIQC_BIN" || ! -x "$MULTIQC_BIN" ]]; then
    echo "Warning: MultiQC requested, but not found at $MULTIQC_BIN"
    echo "Skipping MultiQC for this run."
    RUN_MULTIQC=0
  fi
fi

echo "Using fastp:   $FASTP_BIN"
echo "Using fastqc:  $FASTQC_BIN"
if [[ "$RUN_MULTIQC" == "1" ]]; then
  echo "Using multiqc: $MULTIQC_BIN"
else
  echo "MultiQC disabled for this run."
fi

# -----------------------------
# Gather input files
# -----------------------------
R1_FILES=( "$RAW_DIR"/*_1.fastq.gz )

if (( ${#R1_FILES[@]} == 0 )); then
  echo "ERROR: no *_1.fastq.gz files found in $RAW_DIR"
  exit 1
fi

: > "$FASTP_DIR/fastp.stderr"

# -----------------------------
# Trim paired reads
# -----------------------------
echo "=== Trimming paired reads with fastp ==="

for R1 in "${R1_FILES[@]}"; do
  base="$(basename "$R1")"
  sample="${base%_1.fastq.gz}"
  R2="$RAW_DIR/${sample}_2.fastq.gz"

  if [[ ! -f "$R2" ]]; then
    echo "Skipping $sample: missing pair $R2" >&2
    continue
  fi

  echo "--- Sample: $sample ---"

  "$FASTP_BIN" \
    --in1 "$R1" \
    --in2 "$R2" \
    --detect_adapter_for_pe \
    --trim_front1 20 \
    --trim_front2 20 \
    --trim_tail1 20 \
    --trim_tail2 20 \
    --html "$FASTP_DIR/${sample}.fastp-trim.report.html" \
    --json "$FASTP_DIR/${sample}.fastp-trim.report.json" \
    --out1 "$FASTP_DIR/${sample}_R1.fastp-trim.fq.gz" \
    --out2 "$FASTP_DIR/${sample}_R2.fastp-trim.fq.gz" \
    2>> "$FASTP_DIR/fastp.stderr"

  md5sum "$FASTP_DIR/${sample}_R1.fastp-trim.fq.gz" > "$FASTP_DIR/${sample}_R1.fastp-trim.fq.gz.md5"
  md5sum "$FASTP_DIR/${sample}_R2.fastp-trim.fq.gz" > "$FASTP_DIR/${sample}_R2.fastp-trim.fq.gz.md5"
done

# -----------------------------
# FastQC on trimmed reads
# -----------------------------
TRIMMED_FASTQ=( "$FASTP_DIR"/*.fastp-trim.fq.gz )

if (( ${#TRIMMED_FASTQ[@]} == 0 )); then
  echo "ERROR: no trimmed fastq files were produced."
  exit 1
fi

echo "=== Running FastQC on trimmed reads ==="
"$FASTQC_BIN" -o "$FASTQC_DIR" "${TRIMMED_FASTQ[@]}"

# -----------------------------
# Optional MultiQC
# -----------------------------
if [[ "$RUN_MULTIQC" == "1" ]]; then
  mkdir -p "$MULTIQC_DIR"
  echo "=== Running MultiQC ==="
  "$MULTIQC_BIN" "$FASTQC_DIR" -m fastqc -o "$MULTIQC_DIR" --force

  if [[ ! -f "$MULTIQC_DIR/multiqc_report.html" ]]; then
    echo "ERROR: MultiQC finished, but no multiqc_report.html was created."
    exit 1
  fi

  echo "MultiQC report generated at: $MULTIQC_DIR/multiqc_report.html"
fi

echo "Done."
echo "Trimmed reads:  $FASTP_DIR"
echo "FastQC reports: $FASTQC_DIR"
if [[ "$RUN_MULTIQC" == "1" ]]; then
  echo "MultiQC report: $MULTIQC_DIR/multiqc_report.html"
fi