#!/bin/bash
set -euo pipefail

GENOME_DIR="./reference_genome"
GENOME_NAME="Methylocystis_sp_ref"
FASTA="$GENOME_DIR/${GENOME_NAME}.fna"
READS_DIR="trimmed"
OUT_DIR="aligned"
THREADS=8

HISAT2="/home/nwezejus/miniconda3/envs/rnaseq_qc/bin/hisat2"
SAMTOOLS="/home/nwezejus/miniconda3/envs/rnaseq_qc/bin/samtools"

mkdir -p "$OUT_DIR"

# Build HISAT2 index if missing
if [ ! -f "$GENOME_DIR/${GENOME_NAME}.1.ht2" ]; then
    echo "Building HISAT2 index..."
    $HISAT2-build -p "$THREADS" "$FASTA" "$GENOME_DIR/$GENOME_NAME"
fi

# Align all R1 files
for R1 in "$READS_DIR"/*_1_trimmed.fq.gz; do
    [ -e "$R1" ] || continue
    SAMPLE=$(basename "$R1" "_1_trimmed.fq.gz")
    R2="$READS_DIR/${SAMPLE}_2_trimmed.fq.gz"

    if [[ ! -f "$R2" ]]; then
        echo "ERROR: Missing R2 for sample $SAMPLE"
        exit 1
    fi

    echo "Aligning sample $SAMPLE..."
    $HISAT2 -x "$GENOME_DIR/$GENOME_NAME" \
            -1 "$R1" -2 "$R2" \
            --rna-strandness RF \
            --threads "$THREADS" \
            --summary-file "$OUT_DIR/${SAMPLE}_alignment_summary.txt" \
            2> "$OUT_DIR/${SAMPLE}_hisat2.err" | \
    $SAMTOOLS view -@ "$THREADS" -b -o "$OUT_DIR/${SAMPLE}.bam"

done

# Sort and index all BAMs
for bam in "$OUT_DIR"/*.bam; do
    sorted="${bam%.bam}_sorted.bam"
    $SAMTOOLS sort -@ "$THREADS" -o "$sorted" "$bam"
    $SAMTOOLS index "$sorted"
    rm "$bam"  # optional: remove unsorted
done

echo "✅ All samples aligned and BAMs sorted & indexed!"

