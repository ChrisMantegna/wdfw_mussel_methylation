#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# -----------------------------
# Configuration
# -----------------------------
RAW_URL="https://owl.fish.washington.edu/nightingales/M_trossulus/"
RAW_DIR="../data/raw"
OUT_DIR="../output/qc_reports"

# Local venv-installed MultiQC
MULTIQC_BIN="../.multiqc_env/bin/multiqc"

# -----------------------------
# Make output folders
# -----------------------------
mkdir -p "$RAW_DIR" "$OUT_DIR"

# -----------------------------
# Download FASTQ + checksum files
# -----------------------------
echo "Downloading sequence and checksum files..."
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
  "$RAW_URL"

echo "Download step complete."

# -----------------------------
# Combine checksum files
# -----------------------------
echo "Combining checksum files..."
md5_files=( "$RAW_DIR"/*.fastq.gz.md5sum )

if (( ${#md5_files[@]} > 0 )); then
  cat "${md5_files[@]}" > "$RAW_DIR/combined_checksums.txt"
  echo "Combined checksums written to: $RAW_DIR/combined_checksums.txt"
else
  echo "Warning: no checksum files found matching *.fastq.gz.md5sum"
fi

# -----------------------------
# Find FastQC
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

FASTQC_BIN="${FASTQC_BIN:-}"

if [[ -z "$FASTQC_BIN" ]]; then
  FASTQC_BIN="$(command -v fastqc 2>/dev/null || true)"
fi

FASTQC_BIN="$(find_executable \
  "$FASTQC_BIN" \
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
  echo "Make sure your .multiqc_env was built successfully."
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
# Run FastQC
# -----------------------------
echo "Running FastQC..."
"$FASTQC_BIN" -o "$OUT_DIR" "${fastq_files[@]}"
echo "FastQC complete."

# -----------------------------
# Run MultiQC
# -----------------------------
echo "Running MultiQC..."
"$MULTIQC_BIN" "$OUT_DIR" -o "$OUT_DIR" --force

if [[ ! -f "$OUT_DIR/multiqc_report.html" ]]; then
  echo "ERROR: MultiQC finished, but no multiqc_report.html was created."
  exit 1
fi

echo "MultiQC complete."
echo "Report created at: $OUT_DIR/multiqc_report.html"

# -----------------------------
# Done
# -----------------------------
echo "QC workflow complete."
echo "FastQC outputs: $OUT_DIR"
echo "MultiQC report: $OUT_DIR/multiqc_report.html"