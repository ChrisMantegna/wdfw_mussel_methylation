#!/usr/bin/env bash
set -euo pipefail

# ----------------------------------
# Paths
# ----------------------------------
REF_URL="https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_036588685.1/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"

TRIM_ROOT="../data/trimmed"
TRIM_DIR="$TRIM_ROOT/fastp"      # matches the updated fastp script
REF_DIR="../data/ref"
ALIGN_DIR="../output/aligned"

mkdir -p "$TRIM_ROOT" "$TRIM_DIR" "$REF_DIR" "$ALIGN_DIR"

# ----------------------------------
# Tool discovery
# ----------------------------------
BISMARK_BIN="${BISMARK_BIN:-$(command -v bismark || true)}"
BISMARK_GENOME_PREP_BIN="${BISMARK_GENOME_PREP_BIN:-$(command -v bismark_genome_preparation || true)}"
BOWTIE2_BIN="${BOWTIE2_BIN:-$(command -v bowtie2 || true)}"
CURL_BIN="${CURL_BIN:-$(command -v curl || true)}"
UNZIP_BIN="${UNZIP_BIN:-$(command -v unzip || true)}"

# Search /home/shared if not in PATH
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

if [ -z "$BISMARK_GENOME_PREP_BIN" ]; then
  for p in \
    /home/shared/Bismark*/bismark_genome_preparation \
    /home/shared/*/bin/bismark_genome_preparation \
    /home/shared/*/Bismark*/bismark_genome_preparation
  do
    if [ -f "$p" ] && [ -x "$p" ]; then
      BISMARK_GENOME_PREP_BIN="$p"
      break
    fi
  done
fi

if [ -z "$BISMARK_GENOME_PREP_BIN" ]; then
  BISMARK_GENOME_PREP_BIN="$(find /home/shared -type f -iname 'bismark_genome_preparation*' -perm -111 2>/dev/null | head -n1 || true)"
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

# ----------------------------------
# Guards
# ----------------------------------
if [ -z "$CURL_BIN" ] || [ ! -x "$CURL_BIN" ]; then
  echo "ERROR: curl not found." >&2
  exit 1
fi

if [ -z "$UNZIP_BIN" ] || [ ! -x "$UNZIP_BIN" ]; then
  echo "ERROR: unzip not found." >&2
  exit 1
fi

if [ -z "$BISMARK_BIN" ] || [ ! -x "$BISMARK_BIN" ]; then
  echo "ERROR: bismark not found. Set BISMARK_BIN=/full/path/to/bismark" >&2
  exit 1
fi

if [ -z "$BISMARK_GENOME_PREP_BIN" ] || [ ! -x "$BISMARK_GENOME_PREP_BIN" ]; then
  echo "ERROR: bismark_genome_preparation not found. Set BISMARK_GENOME_PREP_BIN=/full/path/to/bismark_genome_preparation" >&2
  exit 1
fi

if [ -z "$BOWTIE2_BIN" ] || [ ! -x "$BOWTIE2_BIN" ]; then
  echo "ERROR: bowtie2 not found. Set BOWTIE2_BIN=/full/path/to/bowtie2" >&2
  exit 1
fi

export BOWTIE2_BIN
export PATH="$(dirname "$BOWTIE2_BIN"):$PATH"

echo "Using bismark:                  $BISMARK_BIN"
echo "Using bismark_genome_prep:      $BISMARK_GENOME_PREP_BIN"
echo "Using bowtie2:                  $BOWTIE2_BIN"
echo "Trimmed reads directory:        $TRIM_DIR"
echo "Reference directory:            $REF_DIR"
echo "Alignment output directory:     $ALIGN_DIR"

# ----------------------------------
# Check trimmed reads exist
# ----------------------------------
shopt -s nullglob
TRIMMED_R1=("$TRIM_DIR"/*_R1.fastp-trim.fq.gz)

if [ ${#TRIMMED_R1[@]} -eq 0 ]; then
  echo "ERROR: no trimmed R1 files found in $TRIM_DIR" >&2
  echo "Expected files like: SAMPLE_R1.fastp-trim.fq.gz" >&2
  exit 1
fi
shopt -u nullglob

# ----------------------------------
# [1/3] Download reference
# ----------------------------------
echo "=== [1/3] Download and prepare reference genome ==="

REF_ZIP="$REF_DIR/ref.zip"
UNPACK_DIR="$REF_DIR/unpacked"
GENOME_FA="$REF_DIR/genome.fa"

if [ ! -s "$REF_ZIP" ]; then
  echo "Downloading reference dataset ZIP..."
  "$CURL_BIN" -L "$REF_URL" -o "$REF_ZIP"
else
  echo "Reference ZIP already exists at $REF_ZIP"
fi

if [ ! -d "$UNPACK_DIR" ]; then
  echo "Unpacking reference..."
  mkdir -p "$UNPACK_DIR"
  "$UNZIP_BIN" -q "$REF_ZIP" -d "$UNPACK_DIR"
else
  echo "Reference already unpacked at $UNPACK_DIR (skipping unzip)."
fi

if [ ! -s "$GENOME_FA" ]; then
  echo "Combining FASTA files into $GENOME_FA ..."
  mapfile -d '' FASTA_FILES < <(find "$UNPACK_DIR" -type f \( -iname "*.fna" -o -iname "*.fa" -o -iname "*.fasta" \) -print0 | sort -z)

  if [ ${#FASTA_FILES[@]} -eq 0 ]; then
    echo "ERROR: no FASTA files found under $UNPACK_DIR" >&2
    exit 1
  fi

  : > "$GENOME_FA"
  for fa in "${FASTA_FILES[@]}"; do
    cat "$fa" >> "$GENOME_FA"
  done
else
  echo "FASTA already combined at $GENOME_FA (skipping combine)."
fi

# ----------------------------------
# [2/3] Prepare Bismark genome
# ----------------------------------
echo "=== [2/3] Prepare Bismark genome ==="

if [ ! -d "$REF_DIR/Bisulfite_Genome" ]; then
  echo "Running bismark_genome_preparation ..."
  "$BISMARK_GENOME_PREP_BIN" "$REF_DIR"
else
  echo "Bismark genome already prepared in $REF_DIR (skipping)."
fi

# ----------------------------------
# [3/3] Align trimmed reads
# ----------------------------------
echo "=== [3/3] Align trimmed reads with Bismark ==="

shopt -s nullglob
for R1 in "$TRIM_DIR"/*_R1.fastp-trim.fq.gz; do
  sample="$(basename "$R1" _R1.fastp-trim.fq.gz)"
  R2="$TRIM_DIR/${sample}_R2.fastp-trim.fq.gz"

  if [ ! -f "$R2" ]; then
    echo "Skipping $sample: missing mate file $R2"
    continue
  fi

  OUT_BAM_PATTERN="$ALIGN_DIR/${sample}"*_bismark_bt2_pe.bam
  if compgen -G "$OUT_BAM_PATTERN" > /dev/null; then
    echo "$sample already aligned (BAM exists). Skipping."
    continue
  fi

  echo "--- Aligning $sample ---"
  "$BISMARK_BIN" \
    --genome "$REF_DIR" \
    --non_directional \
    -1 "$R1" \
    -2 "$R2" \
    -o "$ALIGN_DIR" \
    > "$ALIGN_DIR/${sample}_stdout.log" \
    2> "$ALIGN_DIR/${sample}_stderr.log"

  echo "$sample alignment complete."
done
shopt -u nullglob

echo "=== Done. Aligned outputs are in: $ALIGN_DIR ==="