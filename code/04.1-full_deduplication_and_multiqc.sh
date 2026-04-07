#!/usr/bin/env bash
set -euo pipefail

# ----------------------------------
# Paths
# ----------------------------------
ALIGN_DIR="../output/aligned_full"
DEDUP_DIR="../output/deduplicated"

# Local venv-installed MultiQC
MULTIQC_BIN="../.multiqc_env/bin/multiqc"

mkdir -p "$DEDUP_DIR"

# ----------------------------------
# Tool discovery
# ----------------------------------
BISMARK_DEDUP_BIN="${BISMARK_DEDUP_BIN:-$(command -v deduplicate_bismark || true)}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-$(command -v samtools || true)}"
MULTIQC_BIN="${MULTIQC_BIN:-$(command -v multiqc || true)}"

# Search /home/shared if needed
if [ -z "$BISMARK_DEDUP_BIN" ]; then
  BISMARK_DEDUP_BIN="$(find /home/shared -type f -iname 'deduplicate_bismark' -perm -111 2>/dev/null | head -n1 || true)"
fi

if [ -z "$SAMTOOLS_BIN" ]; then
  SAMTOOLS_BIN="$(find /home/shared /usr/local /usr -type f -iname 'samtools' -perm -111 2>/dev/null | head -n1 || true)"
fi

if [ -z "$MULTIQC_BIN" ]; then
  MULTIQC_BIN="$(find /home/shared /usr/local /usr -type f -iname 'multiqc' -perm -111 2>/dev/null | head -n1 || true)"
fi

# ----------------------------------
# Guards
# ----------------------------------
if [ -z "$BISMARK_DEDUP_BIN" ] || [ ! -x "$BISMARK_DEDUP_BIN" ]; then
  echo "ERROR: deduplicate_bismark not found." >&2
  exit 1
fi

if [ -z "$SAMTOOLS_BIN" ] || [ ! -x "$SAMTOOLS_BIN" ]; then
  echo "ERROR: samtools not found." >&2
  exit 1
fi

echo "Using deduplicate_bismark: $BISMARK_DEDUP_BIN"
echo "Using samtools:           $SAMTOOLS_BIN"
echo "Using multiqc:            ${MULTIQC_BIN:-not found}"
echo "Alignment directory:      $ALIGN_DIR"
echo "Dedup output directory:   $DEDUP_DIR"

# ----------------------------------
# Check input BAMs
# ----------------------------------
shopt -s nullglob
BAMS=("$ALIGN_DIR"/*_bismark_bt2_pe.bam)

if [ ${#BAMS[@]} -eq 0 ]; then
  echo "ERROR: No BAM files found in $ALIGN_DIR" >&2
  exit 1
fi
shopt -u nullglob

# ----------------------------------
# [1/4] Deduplication
# ----------------------------------
echo "=== [1/4] Deduplication ==="

for BAM in "$ALIGN_DIR"/*_bismark_bt2_pe.bam; do
  sample="$(basename "$BAM" _bismark_bt2_pe.bam)"

  OUT_DEDUP="$DEDUP_DIR/${sample}_bismark_bt2_pe.deduplicated.bam"

  if [ -f "$OUT_DEDUP" ]; then
    echo "$sample already deduplicated. Skipping."
    continue
  fi

  echo "--- Deduplicating $sample ---"
  "$BISMARK_DEDUP_BIN" \
    --paired \
    --bam \
    --output_dir "$DEDUP_DIR" \
    "$BAM" \
    > "$DEDUP_DIR/${sample}_dedup_stdout.log" \
    2> "$DEDUP_DIR/${sample}_dedup_stderr.log"
done

# ----------------------------------
# [2/4] Sorting
# ----------------------------------
echo "=== [2/4] Sorting BAMs ==="

for BAM in "$DEDUP_DIR"/*.deduplicated.bam; do
  sample="$(basename "$BAM" .deduplicated.bam)"
  SORTED="$DEDUP_DIR/${sample}_dedup.sorted.bam"

  if [ -f "$SORTED" ]; then
    echo "$sample already sorted. Skipping."
    continue
  fi

  echo "--- Sorting $sample ---"
  "$SAMTOOLS_BIN" sort \
    -o "$SORTED" \
    "$BAM"
done

# ----------------------------------
# [3/4] Indexing
# ----------------------------------
echo "=== [3/4] Indexing BAMs ==="

for BAM in "$DEDUP_DIR"/*_dedup.sorted.bam; do
  if [ -f "${BAM}.bai" ]; then
    echo "$(basename "$BAM") already indexed. Skipping."
    continue
  fi

  echo "--- Indexing $(basename "$BAM") ---"
  "$SAMTOOLS_BIN" index "$BAM"
done

# ----------------------------------
# [4/4] MultiQC
# ----------------------------------
echo "=== [4/4] MultiQC ==="

if [ -n "${MULTIQC_BIN:-}" ] && [ -x "$MULTIQC_BIN" ]; then
  "$MULTIQC_BIN" "$DEDUP_DIR"
else
  echo "MultiQC not found — skipping report."
fi

# ----------------------------------
# Checksums
# ----------------------------------
echo "=== Generating checksums ==="

cd "$DEDUP_DIR"

: > checksums.md5
for file in *; do
  if [ "$file" != "checksums.md5" ]; then
    md5sum "$file" >> checksums.md5
  fi
done

echo "=== Done. Outputs in: $DEDUP_DIR ==="