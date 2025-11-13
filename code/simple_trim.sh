#!/usr/bin/env bash
set -euo pipefail

# paths
RAW_DIR="../data/raw"
TRIM_DIR="../data/trimmed"

mkdir -p "$TRIM_DIR"

# Tool lookup
FASTP_BIN="${FASTP_BIN:-$(command -v fastp || true)}"
FASTQC_BIN="${FASTQC_BIN:-$(command -v fastqc || true)}"
MULTIQC_BIN="${MULTIQC_BIN:-$(command -v multiqc || true)}"

# If not in PATH, search shared tool area
if [ -z "$FASTP_BIN" ]; then
  FASTP_BIN="$(find /home/shared -type f -iname fastp -perm -111 2>/dev/null | head -n1 || true)"
fi

if [ -z "$FASTQC_BIN" ]; then
  FASTQC_BIN="$(find /home/shared -type f -iname fastqc -perm -111 2>/dev/null | head -n1 || true)"
fi

if [ -z "$MULTIQC_BIN" ]; then
  for p in \
    /home/shared/MultiQC*/bin/multiqc \
    /home/shared/*/bin/multiqc \
    /home/shared/*/MultiQC*/bin/multiqc \
    /home/shared/*/MultiQC*/multiqc
  do
    if [ -f "$p" ] && [ -x "$p" ]; then   # <-- require a regular file
      MULTIQC_BIN="$p"
      break
    fi
  done
fi

# final broad search fallback (already OK: -type f)
if [ -z "$MULTIQC_BIN" ]; then
  MULTIQC_BIN="$(find /home/shared -type f \( -iname 'multiqc' -o -iname 'multiqc*' \) -perm -111 2>/dev/null | head -n1 || true)"
fi

echo "=== Trimming paired reads with fastp ==="
for R1 in "$RAW_DIR"/*_1.fastq.gz; do
  # Skip if no R1 files exist
  [ -e "$R1" ] || { echo "No *_1.fastq.gz files found in $RAW_DIR"; break; }

  base=$(basename "$R1")
  sample="${base%_1.fastq.gz}"
  R2="$RAW_DIR/${sample}_2.fastq.gz"

  if [ ! -f "$R2" ]; then
    echo "⚠️  Skipping $sample: missing pair $R2"
    continue
  fi

  echo "--- Sample: $sample ---"

  "$FASTP_BIN" \
    --in1 "$R1" \
    --in2 "$R2" \
    --detect_adapter_for_pe \
    --trim_front1 20 \
    --trim_front2 20 \
    --html "$TRIM_DIR/${sample}.fastp-trim.report.html" \
    --json "$TRIM_DIR/${sample}.fastp-trim.report.json" \
    --out1 "$TRIM_DIR/${sample}_R1.fastp-trim.fq.gz" \
    --out2 "$TRIM_DIR/${sample}_R2.fastp-trim.fq.gz" \
    2>> "$TRIM_DIR/fastp.stderr"

  # MD5s for trimmed files
  md5sum "$TRIM_DIR/${sample}_R1.fastp-trim.fq.gz" > "$TRIM_DIR/${sample}_R1.fastp-trim.fq.gz.md5"
  md5sum "$TRIM_DIR/${sample}_R2.fastp-trim.fq.gz" > "$TRIM_DIR/${sample}_R2.fastp-trim.fq.gz.md5"
done

echo "=== FastQC on trimmed reads ==="
"$FASTQC_BIN" -o "$TRIM_DIR" "$TRIM_DIR"/*.fastp-trim.fq.gz

echo "=== MultiQC summary ==="
"$MULTIQC_BIN" "$TRIM_DIR" -o "$TRIM_DIR"

echo "Done. Trimmed reads + reports are in: $TRIM_DIR"
