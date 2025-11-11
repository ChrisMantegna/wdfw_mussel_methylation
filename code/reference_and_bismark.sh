#!/usr/bin/env bash
set -euo pipefail

# paths
REF_URL="https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_036588685.1/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"

TRIM_DIR="../data/trimmed"        # from simple_trim.sh outputs
REF_DIR="../data/ref"             # where reference FASTA + Bismark index will live
ALIGN_DIR="../output/aligned"       # where Bismark results go

mkdir -p "$REF_DIR" "$ALIGN_DIR"

echo "=== [1/3] Download and prepare reference genome ==="

# Download the NCBI Datasets ZIP if missing
REF_ZIP="$REF_DIR/ref.zip"
if [ ! -s "$REF_ZIP" ]; then
  echo "Downloading reference dataset ZIP..."
  curl -L "$REF_URL" -o "$REF_ZIP"
else
  echo "Reference ZIP already exists at $REF_ZIP (skipping download)."
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
  # NCBI Datasets typically stores sequences under: unpacked/ncbi_dataset/data/*/*.fna (or .fa/.fasta)
  find "$UNPACK_DIR" -type f \( -iname "*.fna" -o -iname "*.fa" -o -iname "*.fasta" \) -print0 \
    | sort -z \
    | xargs -0 cat > "$GENOME_FA"
else
  echo "FASTA already combined at $GENOME_FA (skipping combine)."
fi

# Prepare Bismark genome (creates Bisulfite_Genome and Bowtie2 index); skip if already present
if [ ! -d "$REF_DIR/Bisulfite_Genome" ]; then
  echo "Running bismark_genome_preparation ..."
  bismark_genome_preparation "$REF_DIR"
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
    echo "⚠️  Skipping $sample: missing $R2"
    continue
  fi

  # Skip if already aligned (Bismark writes SAMPLE_R1.fastp-trim_bismark_bt2_pe.bam by default)
  OUT_BAM_PATTERN="$ALIGN_DIR/${sample}"*_bismark_bt2_pe.bam
  if compgen -G "$OUT_BAM_PATTERN" > /dev/null; then
    echo "✅ $sample already aligned (BAM exists). Skipping."
    continue
  fi

  echo "--- Aligning $sample ---"
  # Non-directional as in your example; remove flag if your library is directional
  bismark \
    --genome "$REF_DIR" \
    --non_directional \
    -1 "$R1" \
    -2 "$R2" \
    -o "$ALIGN_DIR" \
    > "$ALIGN_DIR/${sample}_stdout.log" 2> "$ALIGN_DIR/${sample}_stderr.log"

  echo "✅ $sample alignment complete."
done
shopt -u nullglob

echo "=== [3/3] Done. Aligned outputs are in: $ALIGN_DIR ==="
