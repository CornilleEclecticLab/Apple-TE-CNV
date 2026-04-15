#!/bin/bash
#SBATCH --job-name=Duan_filter_BAM
#SBATCH --output=../../logs/filter_BAMs_Duan-out.txt
#SBATCH --error=../../logs/filter_BAMs_Duan-err.txt
#SBATCH --time=24:00:00
#SBATCH -p compute
##SBATCH -c 32
##SBATCH --mem=64G

# ---------------------------------------------------------------------------
# Minimal BAM filtering pipeline (variable input version)
# - process all BAMs from $INDIR matching "*sorting.bam"
# - keep all contigs from the reference
# - rebuild header from reference (-T)
# - outputs written to $OUTDIR
# ---------------------------------------------------------------------------

module load all gencore/2
module load samtools/1.9

# Directories
INDIR="../../data/BAM/originel/Duan"            # répertoire contenant les BAM
OUTDIR="../../data/BAM/filtered"  # répertoire de sortie
mkdir -p "$OUTDIR"

# Reference fasta (3,780 sequences)
REF="../../data/TElibrary/Final_V2_DUAN_DACCORD_ATIP_DATABASE_non-redundant_150.reformat.reformat.reformat.fasta"

# 1) make sure reference has an index (.fai)
samtools faidx "$REF"

# 2) build BED directly from the .fai (all contigs from REF)
BED="$OUTDIR/fasta_contigs.from_fai.bed"
awk '{print $1"\t0\t"$2}' "${REF}.fai" > "$BED"

# 3) loop over all BAMs ending with "sorting.bam"
for BAM in "$INDIR"/*sorting.bam; do
    base=$(basename "$BAM" .bam)

    FILTERED="$OUTDIR/${base}.filtered.bam"
    SORTED="$OUTDIR/${base}.filtered.sorted.bam"

    # if the file already exists, pass
    if [[ -f "$SORTED" ]]; then
        echo "[SKIP] $SORTED already exists, skipping."
        continue
    fi
    
    echo "[INFO] Processing $BAM"

    # filter alignments to contigs in BED, rebuild header from REF
    samtools view -b -L "$BED" -T "$REF" "$BAM" > "$FILTERED"

    # sort and index
    samtools sort -o "$SORTED" "$FILTERED"
    samtools index "$SORTED"

    # remove intermediate file
    rm -f "$FILTERED"

    echo "[OK] -> $SORTED"
done
