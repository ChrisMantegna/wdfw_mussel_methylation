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

# Optional: set to 1 to force re-run fastp even if outputs exist
FORCE_FASTP="${FORCE_FASTP:-0}"

# Optional: set to 1 to force re-run FastQC even if outputs exist
FORCE_FASTQC="${FORCE_FASTQC:-0}"

# -----------------------------
# Find tools
# -----------------------------
FASTP_BIN="${FASTP_BIN:-$(command -v fastp || true)}"
FASTQC_BIN="${FASTQC_BIN:-$(command -v fastqc || true)}"
MULTIQC_BIN="${MULTIQC_BIN:-$(command -v multiqc || true)}"

# Optional fallback search in shared space
if [[ -z "$FASTP_BIN" ]]; then
  FASTP_BIN="$(find /home/shared -type f -iname fastp -perm -111 2>/dev/null | head -n1 || true)"
fi

if [[ -z "$FASTQC_BIN" ]]; then
  FASTQC_BIN="$(find /home/shared -type f -iname fastqc -perm -111 2>/dev/null | head -n1 || true)"
fi

if [[ -z "$MULTIQC_BIN" ]]; then
  MULTIQC_BIN="$(find /home/shared -type f \( -iname 'multiqc' -o -iname 'multiqc*' \) -perm -111 2>/dev/null | head -n1 || true)"
fi

# -----------------------------
# Check required tools
# -----------------------------
if [[ -z "$FASTP_BIN" || ! -f "$FASTP_BIN" || ! -x "$FASTP_BIN" ]]; then
  echo "Error: fastp executable not found." >&2
  exit 1
fi

if [[ -z "$FASTQC_BIN" || ! -f "$FASTQC_BIN" || ! -x "$FASTQC_BIN" ]]; then
  echo "Error: fastqc executable not found." >&2
  exit 1
fi

if [[ "$RUN_MULTIQC" == "1" ]]; then
  if [[ -z "$MULTIQC_BIN" || ! -f "$MULTIQC_BIN" || ! -x "$MULTIQC_BIN" ]]; then
    echo "Warning: MultiQC requested, but executable not found. Skipping MultiQC." >&2
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
R1_FILES=("$RAW_DIR"/*_1.fastq.gz)

if [[ ${#R1_FILES[@]} -eq 0 ]]; then
  echo "Error: no *_1.fastq.gz files found in $RAW_DIR" >&2
  exit 1
fi

# fresh stderr log for this run
: > "$FASTP_DIR/fastp.stderr"

# -----------------------------
# Trim paired reads with skip logic
# -----------------------------
echo "=== Trimming paired reads with fastp ==="

trim_run_count=0
trim_skip_count=0

for R1 in "${R1_FILES[@]}"; do
  base="$(basename "$R1")"
  sample="${base%_1.fastq.gz}"
  R2="$RAW_DIR/${sample}_2.fastq.gz"

  if [[ ! -f "$R2" ]]; then
    echo "Skipping $sample: missing pair $R2" >&2
    continue
  fi

  out1="$FASTP_DIR/${sample}_R1.fastp-trim.fq.gz"
  out2="$FASTP_DIR/${sample}_R2.fastp-trim.fq.gz"
  md5_1="$FASTP_DIR/${sample}_R1.fastp-trim.fq.gz.md5"
  md5_2="$FASTP_DIR/${sample}_R2.fastp-trim.fq.gz.md5"
  html="$FASTP_DIR/${sample}.fastp-trim.report.html"
  json="$FASTP_DIR/${sample}.fastp-trim.report.json"

  if [[ "$FORCE_FASTP" -eq 0 && -f "$out1" && -f "$out2" && -f "$md5_1" && -f "$md5_2" && -f "$html" && -f "$json" ]]; then
    echo "Skipping fastp for $sample (trim outputs already exist)"
    ((trim_skip_count+=1))
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
    --html "$html" \
    --json "$json" \
    --out1 "$out1" \
    --out2 "$out2" \
    2>> "$FASTP_DIR/fastp.stderr"

  md5sum "$out1" > "$md5_1"
  md5sum "$out2" > "$md5_2"

  ((trim_run_count+=1))
done

echo "fastp complete."
echo "fastp run on:  $trim_run_count sample(s)"
echo "fastp skipped: $trim_skip_count sample(s)"

# -----------------------------
# FastQC on trimmed reads with skip logic
# -----------------------------
TRIMMED_FASTQ=("$FASTP_DIR"/*.fastp-trim.fq.gz)

if [[ ${#TRIMMED_FASTQ[@]} -eq 0 ]]; then
  echo "Error: no trimmed fastq files were produced or found." >&2
  exit 1
fi

echo "=== Running FastQC on trimmed reads ==="

qc_run_count=0
qc_skip_count=0

for fq in "${TRIMMED_FASTQ[@]}"; do
  base="$(basename "$fq")"
  sample="${base%.fq.gz}"

  html_out="$FASTQC_DIR/${sample}_fastqc.html"
  zip_out="$FASTQC_DIR/${sample}_fastqc.zip"

  if [[ "$FORCE_FASTQC" -eq 0 && -f "$html_out" && -f "$zip_out" ]]; then
    echo "Skipping FastQC for $base (outputs already exist)"
    ((qc_skip_count+=1))
    continue
  fi

  echo "Running FastQC on $base"
  "$FASTQC_BIN" -o "$FASTQC_DIR" "$fq"
  ((qc_run_count+=1))
done

echo "FastQC complete."
echo "FastQC run on:  $qc_run_count file(s)"
echo "FastQC skipped: $qc_skip_count file(s)"

# -----------------------------
# Optional MultiQC
# -----------------------------
if [[ "$RUN_MULTIQC" == "1" ]]; then
  mkdir -p "$MULTIQC_DIR"
  echo "=== Running MultiQC ==="
  "$MULTIQC_BIN" "$FASTQC_DIR" -m fastqc -o "$MULTIQC_DIR"
  echo "MultiQC report generated at: $MULTIQC_DIR"
fi

echo "Done."
echo "Trimmed reads:  $FASTP_DIR"
echo "FastQC reports: $FASTQC_DIR"