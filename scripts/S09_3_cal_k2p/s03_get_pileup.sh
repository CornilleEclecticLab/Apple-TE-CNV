#!/bin/bash
#SBATCH --job-name=get_pileup
#SBATCH --output=../../logs/get_pileup.out
#SBATCH --error=../../logs/get_pileup.err
#SBATCH --time=7-00:00:00
#SBATCH -p compute
##SBATCH -c 64
##SBATCH --mem=32G

# Script pour générer des fichiers pileup à partir de BAM déjà triés (id.sorted.bam).

module load all gencore/2
module load samtools/1.9

# Exit si une commande échoue
set -e

# Répertoires et fichiers
BAM_DIR="../../data/BAM/filtered" # from 158,331 sequences to 3,780 seq
PILEUP_DIR="../../output/pileups"
REF_FASTA="../../data/TElibrary/Final_V2_DUAN_DACCORD_ATIP_DATABASE_non-redundant_150.reformat.reformat.reformat.fasta" #3,640 sequences

# Dossier de sortie pour les pileups
mkdir -p "$PILEUP_DIR"

# Index du génome si nécessaire
if [ ! -f "${REF_FASTA}.fai" ]; then
    echo "Indexing reference..."
    samtools faidx "$REF_FASTA"
fi

for SORTED in "$BAM_DIR"/*sorted.bam; do
    NAME=$(basename "$SORTED")
    OUT="${PILEUP_DIR}/${NAME}.pileup"

    # Vérifie si le pileup existe déjà
    if [ -s "$OUT" ]; then
        echo "Skipping $NAME: pileup already exists ($OUT)"
        continue
    fi

    echo "Processing $NAME"

    # Vérifie et indexe le BAM trié
    if [ ! -f "${SORTED}.bai" ]; then
        samtools index "$SORTED"
    fi

    # Crée le pileup
    samtools mpileup -f "$REF_FASTA" "$SORTED" > "$OUT"

    echo " -> $OUT created"
done

echo "All pileups are in $PILEUP_DIR"
