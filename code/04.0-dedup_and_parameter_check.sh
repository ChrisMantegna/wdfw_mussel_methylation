#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# ----------------------------------
# Paths
# ----------------------------------
TRIM_ROOT="../data/trimmed"
TRIM_DIR="$TRIM_ROOT/fastp"
REF_DIR="../data/ref"

TEST_OUT_DIR="../output/all_bismark_min_score"
FULL_ALIGN_DIR="../output/aligned_full"

mkdir -p "$TRIM_ROOT" "$TRIM_DIR" "$REF_DIR" "$TEST_OUT_DIR" "$FULL_ALIGN_DIR"

# ----------------------------------
# Settings
# ----------------------------------
THREADS="${THREADS:-8}"

SCORES=(
  "L,0,-0.4"
  "L,0,-0.6"
  "L,0,-0.8"
  "L,0,-1.0"
  "L,-1,-0.6"
)

# ----------------------------------
# Helper: find executable
# ----------------------------------
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

score_to_label() {
  local score="$1"
  echo "$score" | sed 's/,/_/g; s/-/neg/g'
}

# ----------------------------------
# Tool discovery
# ----------------------------------
BISMARK_BIN="${BISMARK_BIN:-$(command -v bismark || true)}"
DEDUP_BIN="${DEDUP_BIN:-$(command -v deduplicate_bismark || true)}"
BOWTIE2_BIN="${BOWTIE2_BIN:-$(command -v bowtie2 || true)}"

BISMARK_BIN="$(find_executable \
  "$BISMARK_BIN" \
  /home/shared/Bismark*/bismark \
  /home/shared/*/bin/bismark \
  /home/shared/*/Bismark*/bismark \
  /usr/local/bin/bismark \
  /usr/bin/bismark \
)" || {
  echo "ERROR: bismark not found. Set BISMARK_BIN=/full/path/to/bismark" >&2
  exit 1
}

DEDUP_BIN="$(find_executable \
  "$DEDUP_BIN" \
  /home/shared/Bismark*/deduplicate_bismark \
  /home/shared/*/bin/deduplicate_bismark \
  /home/shared/*/Bismark*/deduplicate_bismark \
  /usr/local/bin/deduplicate_bismark \
  /usr/bin/deduplicate_bismark \
)" || {
  echo "ERROR: deduplicate_bismark not found. Set DEDUP_BIN=/full/path/to/deduplicate_bismark" >&2
  exit 1
}

BOWTIE2_BIN="$(find_executable \
  "$BOWTIE2_BIN" \
  /home/shared/*/bin/bowtie2 \
  /home/shared/bowtie2*/bowtie2 \
  /home/shared/*/bowtie2*/bowtie2 \
  /opt/*/bowtie2 \
  /opt/*/bowtie2*/bowtie2 \
  /usr/local/bin/bowtie2 \
  /usr/bin/bowtie2 \
)" || {
  echo "ERROR: bowtie2 not found. Set BOWTIE2_BIN=/full/path/to/bowtie2" >&2
  exit 1
}

if [[ ! -d "$REF_DIR" ]]; then
  echo "ERROR: reference directory not found at $REF_DIR" >&2
  exit 1
fi

if [[ ! -d "$REF_DIR/Bisulfite_Genome" ]]; then
  echo "ERROR: Bismark genome does not appear prepared in $REF_DIR" >&2
  echo "Expected: $REF_DIR/Bisulfite_Genome" >&2
  exit 1
fi

export BOWTIE2_BIN
export PATH="$(dirname "$BOWTIE2_BIN"):$PATH"

echo "Using bismark:              $BISMARK_BIN"
echo "Using deduplicate_bismark:  $DEDUP_BIN"
echo "Using bowtie2:              $BOWTIE2_BIN"
echo "Trimmed reads directory:    $TRIM_DIR"
echo "Reference directory:        $REF_DIR"
echo "Test output directory:      $TEST_OUT_DIR"
echo "Full BAM directory:         $FULL_ALIGN_DIR"

# ----------------------------------
# Check trimmed reads exist
# ----------------------------------
TRIMMED_R1=( "$TRIM_DIR"/*_R1.fastp-trim.fq.gz )

if (( ${#TRIMMED_R1[@]} == 0 )); then
  echo "ERROR: no trimmed R1 files found in $TRIM_DIR" >&2
  echo "Expected files like: SAMPLE_R1.fastp-trim.fq.gz" >&2
  exit 1
fi

# ----------------------------------
# [1/3] Align 10k-read subsets across score_min values
# ----------------------------------
echo "=== [1/3] Test alignment parameter sweep on 10k reads ==="

for R1 in "$TRIM_DIR"/*_R1.fastp-trim.fq.gz; do
  sample="$(basename "$R1" _R1.fastp-trim.fq.gz)"
  R2="$TRIM_DIR/${sample}_R2.fastp-trim.fq.gz"

  if [[ ! -f "$R2" ]]; then
    echo "Skipping $sample: missing mate file $R2"
    continue
  fi

  for score_min in "${SCORES[@]}"; do
    score_label="$(score_to_label "$score_min")"
    SAMPLE_TEST_DIR="$TEST_OUT_DIR/${sample}_${score_label}"

    mkdir -p "$SAMPLE_TEST_DIR"

    OUT_BAM_PATTERN="$SAMPLE_TEST_DIR/${sample}_${score_label}"*_bismark_bt2_pe.bam
    if compgen -G "$OUT_BAM_PATTERN" > /dev/null; then
      echo "$sample with $score_min already aligned (BAM exists). Skipping."
      continue
    fi

    echo "--- Aligning $sample with --score_min $score_min ---"
    "$BISMARK_BIN" \
      --genome "$REF_DIR" \
      --non_directional \
      -1 "$R1" \
      -2 "$R2" \
      -u 10000 \
      -p "$THREADS" \
      --score_min "$score_min" \
      -o "$SAMPLE_TEST_DIR" \
      --basename "${sample}_${score_label}" \
      > "$SAMPLE_TEST_DIR/${sample}_${score_label}_stdout.log" \
      2> "$SAMPLE_TEST_DIR/${sample}_${score_label}_stderr.log"

    echo "$sample with $score_min alignment complete."
  done
done

# ----------------------------------
# [2/3] Deduplicate test BAMs
# ----------------------------------
echo "=== [2/3] Deduplicate test BAMs ==="

TEST_BAMS=( "$TEST_OUT_DIR"/*/*_bismark_bt2_pe.bam )

if (( ${#TEST_BAMS[@]} == 0 )); then
  echo "No test BAMs found in $TEST_OUT_DIR"
else
  for bam in "${TEST_BAMS[@]}"; do
    [[ "$bam" == *.deduplicated.bam ]] && continue

    bam_dir="$(dirname "$bam")"
    bam_base="$(basename "$bam" .bam)"
    dedup_bam="$bam_dir/${bam_base}.deduplicated.bam"

    if [[ -f "$dedup_bam" ]]; then
      echo "$(basename "$bam") already deduplicated. Skipping."
      continue
    fi

    echo "--- Deduplicating test BAM: $(basename "$bam") ---"
    "$DEDUP_BIN" \
      --paired \
      --bam \
      "$bam" \
      > "$bam_dir/${bam_base}_dedup_stdout.log" \
      2> "$bam_dir/${bam_base}_dedup_stderr.log"

    echo "$(basename "$bam") deduplication complete."
  done
fi

# ----------------------------------
# [3/3] Deduplicate existing full BAMs
# ----------------------------------
echo "=== [3/3] Deduplicate existing full BAMs ==="

FULL_BAMS=( "$FULL_ALIGN_DIR"/*_bismark_bt2_pe.bam )

if (( ${#FULL_BAMS[@]} == 0 )); then
  echo "No existing full BAMs found in $FULL_ALIGN_DIR"
else
  for bam in "${FULL_BAMS[@]}"; do
    [[ "$bam" == *.deduplicated.bam ]] && continue

    bam_dir="$(dirname "$bam")"
    bam_base="$(basename "$bam" .bam)"
    dedup_bam="$bam_dir/${bam_base}.deduplicated.bam"

    if [[ -f "$dedup_bam" ]]; then
      echo "$(basename "$bam") already deduplicated. Skipping."
      continue
    fi

    echo "--- Deduplicating full BAM: $(basename "$bam") ---"
    "$DEDUP_BIN" \
      --paired \
      --bam \
      "$bam" \
      > "$bam_dir/${bam_base}_dedup_stdout.log" \
      2> "$bam_dir/${bam_base}_dedup_stderr.log"

    echo "$(basename "$bam") deduplication complete."
  done
fi

echo "=== Done. Test outputs are in: $TEST_OUT_DIR ; full deduplicated BAMs are in: $FULL_ALIGN_DIR ==="