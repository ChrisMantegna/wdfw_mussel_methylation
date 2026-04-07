#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob


# Paths
TRIM_ROOT="../data/trimmed"
TRIM_DIR="$TRIM_ROOT/fastp"
REF_DIR="../data/ref"
ALIGN_DIR="../output/aligned_full"

mkdir -p "$ALIGN_DIR"

# Tool discovery
BISMARK_BIN="${BISMARK_BIN:-$(command -v bismark || true)}"
BOWTIE2_BIN="${BOWTIE2_BIN:-$(command -v bowtie2 || true)}"

# Search /home/shared if not in PATH
if [ -z "$BISMARK_BIN" ]; then
  BISMARK_BIN="$(find /home/shared -type f \( -iname 'bismark' -o -iname 'bismark*' \) -perm -111 2>/dev/null | head -n1 || true)"
fi

if [ -z "$BOWTIE2_BIN" ]; then
  BOWTIE2_BIN="$(find /home/shared /opt /usr/local /usr -type f -iname 'bowtie2' -perm -111 2>/dev/null | head -n1 || true)"
fi

# Guards
if [ -z "$BISMARK_BIN" ] || [ ! -x "$BISMARK_BIN" ]; then
  echo "ERROR: bismark not found." >&2
  exit 1
fi

if [ -z "$BOWTIE2_BIN" ] || [ ! -x "$BOWTIE2_BIN" ]; then
  echo "ERROR: bowtie2 not found." >&2
  exit 1
fi

# Ensure bowtie2 is in PATH for bismark
export PATH="$(dirname "$BOWTIE2_BIN"):$PATH"

echo "Using bismark:              $BISMARK_BIN"
echo "Using bowtie2:              $BOWTIE2_BIN"
echo "Trimmed reads directory:    $TRIM_DIR"
echo "Reference directory:        $REF_DIR"
echo "Alignment output directory: $ALIGN_DIR"

# ----------------------------------
# Check reference is prepared
# ----------------------------------
if [ ! -d "$REF_DIR/Bisulfite_Genome" ]; then
  echo "ERROR: Bisulfite_Genome not found in $REF_DIR"
  echo "Run genome preparation step first."
  exit 1
fi

# Align with skip logic
echo "=== Aligning trimmed reads with Bismark ==="

for R1 in "${TRIMMED_R1[@]}"; do
  sample="$(basename "$R1" _R1.fastp-trim.fq.gz)"
  R2="$TRIM_DIR/${sample}_R2.fastp-trim.fq.gz"

  if [ ! -f "$R2" ]; then
    echo "Skipping $sample: missing mate file"
    continue
  fi

  # Detect completed alignment
  BAM_PATTERN="$ALIGN_DIR/${sample}"*_bismark_bt2_pe.bam

  echo "--- Aligning $sample ---"

  "$BISMARK_BIN" \
  --genome "$REF_DIR" \
  --non_directional \
  --score_min L,0,-0.6 \
  --parallel 4 \
  -1 "$R1" \
  -2 "$R2" \
  -o "$ALIGN_DIR" \
    > "$ALIGN_DIR/${sample}_stdout.log" \
    2> "$ALIGN_DIR/${sample}_stderr.log"

  ((align_run+=1))
done

echo "Alignment complete."
echo "Alignments run:  $align_run"
echo "Alignments skipped: $align_skip"

echo "Outputs in: $ALIGN_DIR"