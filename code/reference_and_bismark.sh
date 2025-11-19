#!/usr/bin/env bash
set -euo pipefail

# paths
REF_URL="https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_036588685.1/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"

TRIM_DIR="../data/trimmed"        # from simple_trim.sh outputs
REF_DIR="../data/ref"             # where reference FASTA + Bismark index will live
ALIGN_DIR="../output/aligned"     # where Bismark results go

mkdir -p "$REF_DIR" "$ALIGN_DIR"

# Tool discovery
BISMARK_BIN="${BISMARK_BIN:-$(command -v bismark || true)}"
BISMARK_GENOME_PREP_BIN="${BISMARK_GENOME_PREP_BIN:-$(command -v bismark_genome_preparation || true)}"
BOWTIE2_BIN="${BOWTIE2_BIN:-$(command -v bowtie2 || true)}"

# If not in PATH, search common /home/shared locations 
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

# Final broad search fallback
if [ -z "$BISMARK_BIN" ]; then
  BISMARK_BIN="$(find /home/shared -type f \( -iname 'bismark' -o -iname 'bismark*' \) -perm -111 2>/dev/null | head -n1 || true)"
fi

# bismark genome preparation

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
  BISMARK_GENOME_PREP_BIN="$(
    find /home/shared -type f -iname 'bismark_genome_preparation*' -perm -111 2>/dev/null | head -n1 || true
  )"
fi

# Bowtie 2 
# Search common module locations
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

# Broad find() fallback (slow but thorough)
if [ -z "$BOWTIE2_BIN" ]; then
  BOWTIE2_BIN="$(
    find /home/shared /opt /usr/local /usr -type f -iname 'bowtie2' -perm -111 2>/dev/null | head -n1 || true
  )"
fi

# Final guard
if [ -z "$BOWTIE2_BIN" ] || [ ! -x "$BOWTIE2_BIN" ]; then
  echo "ERROR: bowtie2 not found. Set BOWTIE2_BIN=/full/path/to/bowtie2"
  exit 1
fi

# Export so Bismark knows where to find it
export BOWTIE2_BIN
export PATH="$(dirname "$BOWTIE2_BIN"):$PATH"

echo "Using Bowtie2 at: $BOWTIE2_BIN"

echo "=== [1/3] Download and prepare reference genome ==="

# Download the NCBI Datasets ZIP if missing
REF_ZIP="$REF_DIR/ref.zip"
if [ ! -s "$REF_ZIP" ]; then
  echo "Downloading reference dataset ZIP..."
  curl -L "$REF_URL" -o "$REF_ZIP"
else
  echo "Reference ZIP already exists at $REF_ZIP"
fi

# Unpack once
UNPACK_DIR="$REF_DIR/unpacked"
if [ ! -d "$UNPACK_DIR" ]; then
  echo "Unpacking reference..."
  mkdir -p "$UNPACK_DIR"
  unzip -q "$REF_ZIP" -d "$UNPACK_DIR"
else
  echo "Reference already unpacked at $UNPACK_DIR (skipping unzip)."
fi

# Collect FASTA(s) into a single genome.fa (Bismark expects fasta(s) in the genome folder)
GENOME_FA="$REF_DIR/genome.fa"
if [ ! -s "$GENOME_FA" ]; then
  echo "Combining FASTA files into $GENOME_FA ..."
  find "$UNPACK_DIR" -type f \( -iname "*.fna" -o -iname "*.fa" -o -iname "*.fasta" \) -print0 \
    | sort -z \
    | xargs -0 cat > "$GENOME_FA"
else
  echo "FASTA already combined at $GENOME_FA (skipping combine)."
fi

# Prepare Bismark genome (creates Bisulfite_Genome and Bowtie2 index); skip if already present
if [ ! -d "$REF_DIR/Bisulfite_Genome" ]; then
  echo "Running bismark_genome_preparation ..."
  "$BISMARK_GENOME_PREP_BIN" "$REF_DIR"
else
  echo "Bismark genome already prepared in $REF_DIR (skipping)."
fi

echo "=== [2/3] Align trimmed reads with Bismark ==="

# For each trimmed pair like SAMPLE_R1.fastp-trim.fq.gz / SAMPLE_R2.fastp-trim.fq.gz
shopt -s nullglob
for R1 in "$TRIM_DIR"/*_R1.fastp-trim.fq.gz; do
  sample="$(basename "$R1" _R1.fastp-trim.fq.gz)"
  R2="$TRIM_DIR/${sample}_R2.fastp-trim.fq.gz"

  if [ ! -f "$R2" ]; then
    echo "Skipping $sample: missing $R2"
    continue
  fi

  # Skip if already aligned (Bismark writes SAMPLE_R1.fastp-trim_bismark_bt2_pe.bam by default)
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
    > "$ALIGN_DIR/${sample}_stdout.log" 2> "$ALIGN_DIR/${sample}_stderr.log"

  echo "$sample alignment complete."
done
shopt -u nullglob

echo "=== [3/3] Done. Aligned outputs are in: $ALIGN_DIR ==="
