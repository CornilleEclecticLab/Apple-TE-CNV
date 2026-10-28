#!/bin/bash
#SBATCH --job-name=get_pileup
#SBATCH --output=../../logs/get_pileup.out
#SBATCH --error=../../logs/get_pileup.err
#SBATCH --time=7-00:00:00
#SBATCH -p compute
##SBATCH -c 64
##SBATCH --mem=32G

# Generate pileup files from BAM files that are already sorted (id.sorted.bam).

module load all gencore/2
module load samtools/1.9

# Exit on the first failing command
set -e

# Directories and files
BAM_DIR="../../data/BAM/filtered" # from 158,331 sequences to 3,780 seq
PILEUP_DIR="../../output/pileups"
REF_FASTA="../../data/TElibrary/Final_V2_DUAN_DACCORD_ATIP_DATABASE_non-redundant_150.reformat.reformat.reformat.fasta" #3,640 sequences

# Output directory for the pileups
mkdir -p "$PILEUP_DIR"

# Index the reference if needed
if [ ! -f "${REF_FASTA}.fai" ]; then
    echo "Indexing reference..."
    samtools faidx "$REF_FASTA"
fi

for SORTED in "$BAM_DIR"/*sorted.bam; do
    NAME=$(basename "$SORTED")
    OUT="${PILEUP_DIR}/${NAME}.pileup"

    # Skip if the pileup already exists
    if [ -s "$OUT" ]; then
        echo "Skipping $NAME: pileup already exists ($OUT)"
        continue
    fi

    echo "Processing $NAME"

    # Check and index the sorted BAM
    if [ ! -f "${SORTED}.bai" ]; then
        samtools index "$SORTED"
    fi

    # Build the pileup
    samtools mpileup -f "$REF_FASTA" "$SORTED" > "$OUT"

    echo " -> $OUT created"
done

echo "All pileups are in $PILEUP_DIR"
