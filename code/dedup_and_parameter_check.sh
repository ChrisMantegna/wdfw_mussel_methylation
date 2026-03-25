#!/usr/bin/env bash
set -euo pipefail


# Paths

READS_DIR="../data/trimmed/fastp"
GENOME_DIR="../data/ref"
TEST_OUT_DIR="../output/bismark_min_score"
FULL_OUT_DIR="../output/full_bismark"

# Tools
BISMARK_BIN="${BISMARK_BIN:-$(command -v bismark || true)}"
DEDUP_BIN="/home/shared/Bismark-0.24.0/deduplicate_bismark"
BOWTIE2_BIN="${BOWTIE2_BIN:-$(command -v bowtie2 || true)}"

# Threads
THREADS=8

# Full-run chosen score_min
FULL_SCORE_MIN="L,0,-1.0"

# Test-only score_min values
SCORES=(
  "L,0,-0.4"
  "L,0,-0.6"
  "L,0,-0.8"
  "L,0,-1.0"
  "L,-1,-0.6"
)

# Checks
mkdir -p "${TEST_OUT_DIR}" "${FULL_OUT_DIR}"

if [ -z "$BISMARK_BIN" ]; then
  for p in \
    /home/shared/Bismark*/bismark \
    /home/shared/*/bin/bismark \
    /home/shared/*/Bismark*/bismark
  do
    if [ -f "$p" ] && [ -x "$p" ]; then
      BISMARK_BIN="$p"
      break
    fi
  done
fi

if [ -z "$BISMARK_BIN" ]; then
  BISMARK_BIN="$(find /home/shared -type f \( -iname 'bismark' -o -iname 'bismark*' \) -perm -111 2>/dev/null | head -n1 || true)"
fi

if [[ ! -x "${BISMARK_BIN}" ]]; then
  echo "ERROR: bismark not found or not executable at ${BISMARK_BIN}"
  exit 1
fi

if [[ ! -x "${DEDUP_BIN}" ]]; then
  echo "ERROR: deduplicate_bismark not found or not executable at ${DEDUP_BIN}"
  exit 1
fi

if [ -z "$BOWTIE2_BIN" ]; then
  for p in \
    /home/shared/*/bin/bowtie2 \
    /home/shared/bowtie2*/bowtie2 \
    /home/shared/*/bowtie2*/bowtie2 \
    /opt/*/bowtie2 \
    /opt/*/bowtie2*/bowtie2 \
    /usr/local/bin/bowtie2 \
    /usr/bin/bowtie2
  do
    if [ -f "$p" ] && [ -x "$p" ]; then
      BOWTIE2_BIN="$p"
      break
    fi
  done
fi

if [ -z "$BOWTIE2_BIN" ]; then
  BOWTIE2_BIN="$(find /home/shared /opt /usr/local /usr -type f -iname 'bowtie2' -perm -111 2>/dev/null | head -n1 || true)"
fi

if [[ ! -d "${BOWTIE2_BIN}" ]]; then
  echo "ERROR: Bowtie2 directory not found at ${BOWTIE2_BIN}"
  exit 1
fi

if [[ ! -d "${GENOME_DIR}" ]]; then
  echo "ERROR: Genome directory not found at ${GENOME_DIR}"
  exit 1
fi

# Helper function

score_to_label() {
  local score="$1"
  echo "${score}" | sed 's/,/_/g; s/-/neg/g'
}

# Step 1: Test parameter sweep on 10k reads

echo "1/2: Test alignments on 10k reads for test samples"

for r1 in "${READS_DIR}"/*_R1.fastp-trim.fq.gz; do
  [[ -e "$r1" ]] || continue

  sample_name=$(basename "$r1" "_R1.fastp-trim.fq.gz")
  r2="${READS_DIR}/${sample_name}_R2.fastp-trim.fq.gz"

  if [[ ! -f "${r2}" ]]; then
    echo "Missing R2 for ${sample_name}; skipping."
    continue
  fi

  for score_min in "${SCORES[@]}"; do
    score_label="$(score_to_label "${score_min}")"
    param_dir="${TEST_OUT_DIR}/${sample_name}_${score_label}"
    bam="${param_dir}/${sample_name}_${score_label}_pe.bam"
    dedup_bam="${param_dir}/${sample_name}_${score_label}_pe.deduplicated.bam"

    mkdir -p "${param_dir}"

    if [[ -f "${bam}" ]]; then
      echo "[TEST ALIGN] Exists already: ${bam} -- skipping alignment"
    else
      echo "[TEST ALIGN] ${sample_name} with --score_min ${score_min}"

      "${BISMARK_BIN}" \
        --path_to_bowtie2 "${BOWTIE2_DIR}" \
        --genome "${GENOME_DIR}" \
        -1 "${r1}" \
        -2 "${r2}" \
        -u 10000 \
        -p "${THREADS}" \
        --score_min "${score_min}" \
        -o "${param_dir}" \
        --basename "${sample_name}_${score_label}" \
        2> "${param_dir}/${sample_name}_${score_label}_bismark_summary.txt"
    fi

    if [[ -f "${dedup_bam}" ]]; then
      echo "[TEST DEDUP] Exists already: ${dedup_bam} -- skipping deduplication"
    elif [[ -f "${bam}" ]]; then
      echo "[TEST DEDUP] Deduplicating ${bam}"

      "${DEDUP_BIN}" \
        --paired \
        --bam \
        "${bam}" \
        2> "${param_dir}/${sample_name}_${score_label}_dedup_summary.txt"
    else
      echo "[TEST DEDUP] No BAM found for ${sample_name}; skipping dedup"
    fi
  done
done

# Step 2: Deduplicate existing full BAMs

echo "2/2 Deduplicate existing full BAMs"


shopt -s nullglob
full_bams=( "${FULL_OUT_DIR}"/*_pe.bam )

if [[ ${#full_bams[@]} -eq 0 ]]; then
  echo "No existing full BAMs found in ${FULL_OUT_DIR}"
else
  for bam in "${full_bams[@]}"; do
    base="${bam%.bam}"
    dedup_bam="${base}.deduplicated.bam"

    if [[ -f "${dedup_bam}" ]]; then
      echo "[FULL DEDUP] Exists already: ${dedup_bam} -- skipping"
      continue
    fi

    echo "[FULL DEDUP] Deduplicating ${bam}"

    "${DEDUP_BIN}" \
      --paired \
      --bam \
      "${bam}" \
      2> "${base}_dedup_summary.txt"
  done
fi

echo
echo "Done."