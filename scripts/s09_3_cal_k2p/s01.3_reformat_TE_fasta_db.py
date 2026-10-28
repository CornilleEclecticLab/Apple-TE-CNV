from Bio import SeqIO

# Aim: to check

# Input files
fasta_out = "../../data/TElibrary/Final_V2_DUAN_DACCORD_ATIP_DATABASE_non-redundant_150.reformat.reformat.reformat.fasta" # 3,780 seq
bam_names_file = "../../data/BAM/seqNames_from_BAM.txt" # 158,331 seq

# ---- Step 1: read contigs from FASTA ----
fasta_contigs = {record.id for record in SeqIO.parse(fasta_out, "fasta")}
print(f"Total contigs in FASTA: {len(fasta_contigs)}")

# ---- Step 2: read contigs from BAM names ----
with open(bam_names_file) as f:
    bam_contigs = {line.strip() for line in f}
print(f"Total contigs in BAM list: {len(bam_contigs)}")

# ---- Step 3: check which BAM contigs are missing in FASTA ----
missing = bam_contigs - fasta_contigs
found = bam_contigs & fasta_contigs

print(f"\n✅ Contigs from BAM found in FASTA: {len(found)}")
print(f"⚠️ Contigs from BAM missing in FASTA: {len(missing)}")

# if missing:
#     print("\nList of missing contigs:")
#     for c in sorted(missing):
#         print(" -", c)
# else:
#     print("\n❌ No BAM contigs are missing from the FASTA (all are present)")