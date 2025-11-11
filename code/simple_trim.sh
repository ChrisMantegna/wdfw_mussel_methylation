#!/usr/bin/env bash
set -euo pipefail

# paths
RAW_DIR="../data/raw"
TRIM_DIR="../data/trimmed"

mkdir -p "$TRIM_DIR"

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

  fastp \
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
fastqc -o "$TRIM_DIR" "$TRIM_DIR"/*.fastp-trim.fq.gz

echo "=== MultiQC summary ==="
multiqc "$TRIM_DIR" -o "$TRIM_DIR"

# Optional: remove FastQC zip bundles (HTMLs are kept; zips are usually redundant)
# rm "$TRIM_DIR"/*.zip 2>/dev/null || true

echo "🎉 Done. Trimmed reads + reports are in: $TRIM_DIR"
