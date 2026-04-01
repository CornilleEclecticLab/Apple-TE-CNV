from Bio import SeqIO

#Aim: BAM files contain all consensus but I just want the subset under studying in Antho's paper

# Input files
fasta_in = "../../data/TElibrary/Final_V2_DUAN_DACCORD_ATIP_DATABASE_non-redundant_150.fasta"
txt_files = ["../../data/TElists/250923_TE_structure.txt", "../../data/TElists/250923_TE_residuals.txt"]

# Output file
fasta_out = "../../data/TElibrary/Final_V2_DUAN_DACCORD_ATIP_DATABASE_non-redundant_150.reformat.fasta"

# ---- Step 1: read contig names to keep ----
contigs_to_keep = set()
for txt in txt_files:
    with open(txt) as f:
        for line in f:
            contigs_to_keep.add(line.strip())

    print(f"Reading {txt}... {len(contigs_to_keep)} contigs to search for in the FASTA file")

# ---- Step 2: filter the FASTA and track found contigs ----
found_contigs = set()
records_to_write = []
for record in SeqIO.parse(fasta_in, "fasta"):
    base_id = record.id.split("_PIF_Harbinger")[0] # Because I got ⚠️ Contigs not found in the FASTA file:
#  - TE_039238_LTR
#  - TE_051960_LTR
#  - TE_093261_LTR
#  - TE_118194_LTR
    if base_id in contigs_to_keep:
        records_to_write.append(record)
        found_contigs.add(record.id)

# ---- Step 3: write the new FASTA ----
count = SeqIO.write(records_to_write, fasta_out, "fasta")
print(f"{count} sequences written to {fasta_out}")

# ---- Step 4: report missing contigs ----
missing = contigs_to_keep - found_contigs
if missing:
    print("\n⚠️ Contigs not found in the FASTA file:")
    for c in sorted(missing):
        print(" -", c)
else:
    print("\n✅ All contigs were found in the FASTA file")
