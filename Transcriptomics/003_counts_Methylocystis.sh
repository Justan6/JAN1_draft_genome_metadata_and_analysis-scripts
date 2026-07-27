#!/bin/bash
set -euo pipefail

# ==========================
# Configuration
# ==========================
GFF="./reference_genome/Methylocystis_sp_ref.gff"
BAM_DIR="aligned"   # Output from hisat2_Methylocystis.sh
OUT_FILE="results/counts/counts_Methylocystis.txt"
THREADS=${SLURM_CPUS_PER_TASK:-10}

mkdir -p results/counts

# ==========================
# Sanity checks
# ==========================
if [ ! -f "$GFF" ]; then
    echo "ERROR: GFF file not found: $GFF"
    exit 1
fi

BAM_FILES=("$BAM_DIR"/*_sorted.bam)
if [ ${#BAM_FILES[@]} -eq 0 ]; then
    echo "ERROR: No BAM files found in $BAM_DIR"
    exit 1
fi

# ==========================
# Run featureCounts
# ==========================
echo "Running featureCounts on ${#BAM_FILES[@]} BAM files..."

featureCounts \
    -T "$THREADS" \
    -p -B -C \
    -a "$GFF" \
    -t CDS \
    -g Parent \
    -s 2 \
    -o "$OUT_FILE" \
    "${BAM_FILES[@]}"

echo "✅ Gene counts table saved to: $OUT_FILE"

# ==========================
# Alignment statistics
# ==========================
echo ""
echo "Generating flagstat summaries..."
for bam in "${BAM_FILES[@]}"; do
    sample=$(basename "$bam" .bam)
    echo "=== $sample ==="
    samtools flagstat "$bam" | grep "mapped"
    echo ""
done

