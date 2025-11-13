#!/usr/bin/env bash
set -euo pipefail

# Configuration - loacation of sequences
RAW_URL="https://owl.fish.washington.edu/nightingales/M_trossulus/"

# Local folders (adjust if needed)
RAW_DIR="../data/raw"
OUT_DIR="../output/qc_reports"

# Find fastqc and multiqc
FASTQC_BIN="${FASTQC_BIN:-$(command -v fastqc || true)}"
MULTIQC_BIN="${MULTIQC_BIN:-$(command -v multiqc || true)}"

# If not in PATH, search shared tool area
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

# MultiQC
echo "Running MultiQC..."
"$MULTIQC_BIN" "$OUT_DIR" -o "$OUT_DIR"
echo "MultiQC report generated at $OUT_DIR"
