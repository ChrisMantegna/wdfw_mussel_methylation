#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# Configuration
# -----------------------------
DEDUP_DIR="../output/deduplicated"
METHYL_DIR="../output/methylation_extraction"
REPORT_DIR="$METHYL_DIR/reports"
FILTER_DIR="$METHYL_DIR/filtered_cov"
LOG_DIR="$METHYL_DIR/logs"

# Genome folder used for --cytosine_report
GENOME_DIR="../data/ref"

# Processing controls
MIN_COV="${MIN_COV:-5}"
BUFFER_SIZE="${BUFFER_SIZE:-50%}"

mkdir -p "$METHYL_DIR" "$REPORT_DIR" "$FILTER_DIR" "$LOG_DIR"

# -----------------------------
# Tool discovery
# -----------------------------
BISMARK_METHYLATION_EXTRACTOR_BIN="${BISMARK_METHYLATION_EXTRACTOR_BIN:-$(command -v bismark_methylation_extractor || true)}"
BISMARK2REPORT_BIN="${BISMARK2REPORT_BIN:-$(command -v bismark2report || true)}"
BISMARK2SUMMARY_BIN="${BISMARK2SUMMARY_BIN:-$(command -v bismark2summary || true)}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-$(command -v samtools || true)}"
MULTIQC_BIN="${MULTIQC_BIN:-$(command -v multiqc || true)}"
GZIP_BIN="${GZIP_BIN:-$(command -v gzip || true)}"
ZCAT_BIN="${ZCAT_BIN:-$(command -v zcat || command -v gzcat || true)}"

if [ -z "$BISMARK_METHYLATION_EXTRACTOR_BIN" ]; then
  BISMARK_METHYLATION_EXTRACTOR_BIN="$(find /home/shared /usr/local /usr -type f -iname 'bismark_methylation_extractor' -perm -111 2>/dev/null | head -n1 || true)"
fi

if [ -z "$BISMARK2REPORT_BIN" ]; then
  BISMARK2REPORT_BIN="$(find /home/shared /usr/local /usr -type f -iname 'bismark2report' -perm -111 2>/dev/null | head -n1 || true)"
fi

if [ -z "$BISMARK2SUMMARY_BIN" ]; then
  BISMARK2SUMMARY_BIN="$(find /home/shared /usr/local /usr -type f -iname 'bismark2summary' -perm -111 2>/dev/null | head -n1 || true)"
fi

if [ -z "$SAMTOOLS_BIN" ]; then
  SAMTOOLS_BIN="$(find /home/shared /usr/local /usr -type f -iname 'samtools' -perm -111 2>/dev/null | head -n1 || true)"
fi

if [ -z "$MULTIQC_BIN" ]; then
  MULTIQC_BIN="$(find /home/shared /usr/local /usr -type f -iname 'multiqc' -perm -111 2>/dev/null | head -n1 || true)"
fi

# -----------------------------
# Guards
# -----------------------------
if [ -z "$BISMARK_METHYLATION_EXTRACTOR_BIN" ] || [ ! -x "$BISMARK_METHYLATION_EXTRACTOR_BIN" ]; then
  echo "ERROR: bismark_methylation_extractor not found." >&2
  exit 1
fi

if [ -z "$SAMTOOLS_BIN" ] || [ ! -x "$SAMTOOLS_BIN" ]; then
  echo "ERROR: samtools not found." >&2
  exit 1
fi

if [ -z "$GZIP_BIN" ] || [ ! -x "$GZIP_BIN" ]; then
  echo "ERROR: gzip not found." >&2
  exit 1
fi

if [ -z "$ZCAT_BIN" ] || ! command -v "$ZCAT_BIN" >/dev/null 2>&1; then
  echo "ERROR: zcat/gzcat not found." >&2
  exit 1
fi

if [ ! -d "$DEDUP_DIR" ]; then
  echo "ERROR: deduplicated BAM directory not found: $DEDUP_DIR" >&2
  exit 1
fi

if [ ! -d "$GENOME_DIR" ]; then
  echo "ERROR: genome directory not found: $GENOME_DIR" >&2
  echo "Update GENOME_DIR in this script before running." >&2
  exit 1
fi


echo "Using bismark_methylation_extractor: $BISMARK_METHYLATION_EXTRACTOR_BIN"
echo "Using bismark2report:              ${BISMARK2REPORT_BIN:-not found}"
echo "Using bismark2summary:             ${BISMARK2SUMMARY_BIN:-not found}"
echo "Using samtools:                    $SAMTOOLS_BIN"
echo "Using multiqc:                     ${MULTIQC_BIN:-not found}"
echo "Dedup directory:                   $DEDUP_DIR"
echo "Methylation output directory:      $METHYL_DIR"
echo "Genome directory:                  $GENOME_DIR"
echo "Minimum coverage filter:           $MIN_COV"

# -----------------------------
# Input BAMs
# -----------------------------
BAMS=(
  "$DEDUP_DIR"/*.deduplicated.bam
  "$DEDUP_DIR"/*_dedup.sorted.bam
  "$DEDUP_DIR"/*.sorted.bam
)

# de-duplicate any repeated matches
if [ ${#BAMS[@]} -eq 0 ]; then
  echo "ERROR: No deduplicated or sorted BAM files found in $DEDUP_DIR" >&2
  exit 1
fi

mapfile -t UNIQUE_BAMS < <(printf '%s\n' "${BAMS[@]}" | awk '!seen[$0]++')

# -----------------------------
# [1/6] Methylation extraction
# -----------------------------
echo "=== [1/6] Methylation extraction ==="

for BAM in "${UNIQUE_BAMS[@]}"; do
  [ -f "$BAM" ] || continue

  BAM_BASE="$(basename "$BAM")"
  SAMPLE_BASE="${BAM_BASE%.bam}"
  SPLIT_REPORT="$METHYL_DIR/${SAMPLE_BASE}_splitting_report.txt"

  if ls "$METHYL_DIR"/${SAMPLE_BASE}*.bismark.cov.gz >/dev/null 2>&1 || [ -f "$SPLIT_REPORT" ]; then
    echo "$SAMPLE_BASE already appears extracted. Skipping."
    continue
  fi

  echo "--- Extracting methylation for $SAMPLE_BASE ---"
  "$BISMARK_METHYLATION_EXTRACTOR_BIN" \
    --bedGraph \
    --counts \
    --scaffolds \
    --output "$METHYL_DIR" \
    --samtools_path "$(dirname "$SAMTOOLS_BIN")" \
    --gzip \
    --paired-end \
    --zero_based \
    --buffer_size "$BUFFER_SIZE" \
    --cytosine_report \
    --genome_folder "$GENOME_DIR" \
    "$BAM" \
    > "$LOG_DIR/${SAMPLE_BASE}_methylation_extractor.stdout.log" \
    2> "$LOG_DIR/${SAMPLE_BASE}_methylation_extractor.stderr.log"
done

# -----------------------------
# [2/6] Bismark report + summary
# -----------------------------
echo "=== [2/6] Bismark report and summary ==="

cd "$METHYL_DIR"

if [ -n "${BISMARK2SUMMARY_BIN:-}" ] && [ -x "$BISMARK2SUMMARY_BIN" ]; then
  if [ ! -f "$REPORT_DIR/bismark_summary_report.html" ]; then
    echo "--- Running bismark2summary ---"
    "$BISMARK2SUMMARY_BIN" ./*.bam ./*_report.txt ./*_splitting_report.txt  \
      > "$LOG_DIR/bismark2summary.stdout.log" \
      2> "$LOG_DIR/bismark2summary.stderr.log" || true

    shopt -s nullglob
    for f in bismark_summary_report.* bismark_summary_*; do
      mv "$f" "$REPORT_DIR/" 2>/dev/null || true
    done
    shopt -u nullglob
  else
    echo "bismark2summary output already present. Skipping."
  fi
else
  echo "bismark2summary not found — skipping."
fi

if [ -n "${BISMARK2REPORT_BIN:-}" ] && [ -x "$BISMARK2REPORT_BIN" ]; then
  shopt -s nullglob
  SPLIT_REPORTS=(./*_splitting_report.txt)
  ALIGN_REPORTS=(./*_report.txt)
  MBIAS_FILES=(./*.M-bias.txt)
  shopt -u nullglob

  if [ ${#SPLIT_REPORTS[@]} -gt 0 ]; then
    echo "--- Running per-sample bismark2report ---"

    for SPLIT in "${SPLIT_REPORTS[@]}"; do
      BASE="$(basename "$SPLIT" _splitting_report.txt)"
      SAMPLE_HTML="$REPORT_DIR/${BASE}.bismark_report.html"

      if [ -f "$SAMPLE_HTML" ]; then
        echo "Report for $BASE already exists. Skipping."
        continue
      fi

      ALIGN_MATCH=""
      MBIAS_MATCH=""

      for A in "${ALIGN_REPORTS[@]}"; do
        if [[ "$A" == *"$BASE"* ]]; then
          ALIGN_MATCH="$A"
          break
        fi
      done

      for M in "${MBIAS_FILES[@]}"; do
        if [[ "$M" == *"$BASE"* ]]; then
          MBIAS_MATCH="$M"
          break
        fi
      done

      echo "Generating report for $BASE"
      if [ -n "$ALIGN_MATCH" ] && [ -n "$MBIAS_MATCH" ]; then
        "$BISMARK2REPORT_BIN" \
          --alignment_report "$ALIGN_MATCH" \
          --splitting_report "$SPLIT" \
          --mbias_report "$MBIAS_MATCH" \
          > "$LOG_DIR/${BASE}_bismark2report.stdout.log" \
          2> "$LOG_DIR/${BASE}_bismark2report.stderr.log" || true
      elif [ -n "$ALIGN_MATCH" ]; then
        "$BISMARK2REPORT_BIN" \
          --alignment_report "$ALIGN_MATCH" \
          --splitting_report "$SPLIT" \
          > "$LOG_DIR/${BASE}_bismark2report.stdout.log" \
          2> "$LOG_DIR/${BASE}_bismark2report.stderr.log" || true
      else
        "$BISMARK2REPORT_BIN" \
          --splitting_report "$SPLIT" \
          > "$LOG_DIR/${BASE}_bismark2report.stdout.log" \
          2> "$LOG_DIR/${BASE}_bismark2report.stderr.log" || true
      fi

      shopt -s nullglob
      for HTML in ./*report.html; do
        if [ -f "$HTML" ]; then
          NEW_NAME="$REPORT_DIR/${BASE}.bismark_report.html"
          mv "$HTML" "$NEW_NAME" 2>/dev/null || true
        fi
      done
      shopt -u nullglob
    done
  else
    echo "No splitting reports found for bismark2report."
  fi
else
  echo "bismark2report not found — skipping."
fi

# -----------------------------
# [3/6] MultiQC
# -----------------------------
echo "=== [3/6] MultiQC ==="

if [ -n "${MULTIQC_BIN:-}" ] && [ -x "$MULTIQC_BIN" ]; then
  if [ ! -f "$REPORT_DIR/multiqc_report.html" ]; then
    "$MULTIQC_BIN" "$METHYL_DIR" --outdir "$REPORT_DIR" \
      > "$LOG_DIR/multiqc.stdout.log" \
      2> "$LOG_DIR/multiqc.stderr.log" || true
  else
    echo "MultiQC report already exists. Skipping."
  fi
else
  echo "MultiQC not found — skipping report."
fi

# -----------------------------
# [4/6] Coverage filtering
# -----------------------------
echo "=== [4/6] Coverage filtering of *.bismark.cov.gz ==="

COV_FILES=("$METHYL_DIR"/*.bismark.cov.gz)
if [ ${#COV_FILES[@]} -eq 0 ]; then
  echo "No *.bismark.cov.gz files found — skipping coverage filtering."
else
  for COV in "${COV_FILES[@]}"; do
    [ -f "$COV" ] || continue

    BASE="$(basename "$COV" .bismark.cov.gz)"
    FILTERED_OUT="$FILTER_DIR/${BASE}.min${MIN_COV}x.bismark.cov.gz"
    FILTER_SUMMARY="$FILTER_DIR/${BASE}.min${MIN_COV}x.filter_summary.txt"

    if [ -f "$FILTERED_OUT" ] && [ -f "$FILTER_SUMMARY" ]; then
      echo "$BASE already filtered at ${MIN_COV}x. Skipping."
      continue
    fi

    echo "--- Filtering $BASE at coverage >= $MIN_COV ---"

    TOTAL_LOCI="$($ZCAT_BIN "$COV" | wc -l | awk '{print $1}')"

    $ZCAT_BIN "$COV" \
      | awk -v min_cov="$MIN_COV" 'BEGIN{OFS="\t"} (($5 + $6) >= min_cov) {print $0}' \
      | "$GZIP_BIN" -c > "$FILTERED_OUT"

    KEPT_LOCI="$($ZCAT_BIN "$FILTERED_OUT" | wc -l | awk '{print $1}')"
    REMOVED_LOCI="$((TOTAL_LOCI - KEPT_LOCI))"

    {
      echo -e "sample\t$BASE"
      echo -e "input_file\t$(basename "$COV")"
      echo -e "min_coverage\t$MIN_COV"
      echo -e "total_loci\t$TOTAL_LOCI"
      echo -e "kept_loci\t$KEPT_LOCI"
      echo -e "removed_loci\t$REMOVED_LOCI"
      if [ "$TOTAL_LOCI" -gt 0 ]; then
        awk -v kept="$KEPT_LOCI" -v total="$TOTAL_LOCI" 'BEGIN{printf "percent_kept\t%.4f\n", (kept/total)*100}'
      else
        echo -e "percent_kept\t0"
      fi
    } > "$FILTER_SUMMARY"
  done
fi

# -----------------------------
# [5/6] Simple QC summaries from extractor outputs
# -----------------------------
echo "=== [5/6] Collating simple QC summaries ==="

SUMMARY_TSV="$REPORT_DIR/methylation_extraction_file_manifest.tsv"
MBIAS_MANIFEST="$REPORT_DIR/mbias_and_context_files.tsv"
FILTER_MANIFEST="$REPORT_DIR/filtering_manifest.tsv"

{
  echo -e "sample_base\tbedgraph_gz\tzero_cov\tcov_gz\tCpG_report_gz\tcytosine_context_summary\tmbias_txt\tsplitting_report"

  shopt -s nullglob
  for SPLIT in "$METHYL_DIR"/*_splitting_report.txt; do
    BASE="$(basename "$SPLIT" _splitting_report.txt)"
    BEDGRAPH="$(compgen -G "$METHYL_DIR/${BASE}*.bedGraph.gz" | head -n1 || true)"
    ZERO_COV="$(compgen -G "$METHYL_DIR/${BASE}*.bedGraph.gz.bismark.zero.cov" | head -n1 || true)"
    COV="$(compgen -G "$METHYL_DIR/${BASE}*.bismark.cov.gz" | head -n1 || true)"
    CPG="$(compgen -G "$METHYL_DIR/${BASE}*.CpG_report.txt.gz" | head -n1 || true)"
    CONTEXT="$(compgen -G "$METHYL_DIR/${BASE}*.cytosine_context_summary.txt" | head -n1 || true)"
    MBIAS="$(compgen -G "$METHYL_DIR/${BASE}*.M-bias.txt" | head -n1 || true)"

    echo -e "${BASE}\t${BEDGRAPH:-NA}\t${ZERO_COV:-NA}\t${COV:-NA}\t${CPG:-NA}\t${CONTEXT:-NA}\t${MBIAS:-NA}\t${SPLIT:-NA}"
  done
  shopt -u nullglob
} > "$SUMMARY_TSV"

{
  echo -e "sample_base\tmbias_txt\tcytosine_context_summary"
  shopt -s nullglob
  for MBIAS in "$METHYL_DIR"/*.M-bias.txt; do
    BASE="$(basename "$MBIAS" .M-bias.txt)"
    CONTEXT="$(compgen -G "$METHYL_DIR/${BASE}*.cytosine_context_summary.txt" | head -n1 || true)"
    echo -e "${BASE}\t${MBIAS}\t${CONTEXT:-NA}"
  done
  shopt -u nullglob
} > "$MBIAS_MANIFEST"

{
  echo -e "sample\tmin_coverage\ttotal_loci\tkept_loci\tremoved_loci\tpercent_kept"
  shopt -s nullglob
  for F in "$FILTER_DIR"/*.filter_summary.txt; do
    awk -F '\t' '
      BEGIN{sample=""; min=""; total=""; kept=""; removed=""; pct=""}
      $1=="sample"{sample=$2}
      $1=="min_coverage"{min=$2}
      $1=="total_loci"{total=$2}
      $1=="kept_loci"{kept=$2}
      $1=="removed_loci"{removed=$2}
      $1=="percent_kept"{pct=$2}
      END{print sample"\t"min"\t"total"\t"kept"\t"removed"\t"pct}
    ' "$F"
  done
  shopt -u nullglob
} > "$FILTER_MANIFEST"

# -----------------------------
# [6/6] Checksums
# -----------------------------
echo "=== [6/6] Generating checksums ==="

cd "$METHYL_DIR"
: > checksums.md5
for file in * reports/* filtered_cov/* logs/*; do
  [ -f "$file" ] || continue
  if [ "$(basename "$file")" != "checksums.md5" ]; then
    md5sum "$file" >> checksums.md5
  fi
done

echo "=== Done. Outputs in: $METHYL_DIR ==="
echo "Key report directory: $REPORT_DIR"
echo "Filtered coverage directory: $FILTER_DIR"
