#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import sys

# Input files
fasta_file = "../../data/TElibrary/Final_V2_DUAN_DACCORD_ATIP_DATABASE_non-redundant_150.reformat.fasta"
tsv_file = "../../data/TElibrary/V2_Pastec_TE_sorter_Annotation_final_V1.tsv" # TE classif

# Output file
out_fasta = "../../data/TElibrary/Final_V2_DUAN_DACCORD_ATIP_DATABASE_non-redundant_150.reformat.reformat.fasta"

# ---------------------------------------------------------------
# Step 1: Load mapping from TSV
#   Key   = contig_name_fasta_def
#   Value = ID from 2nd column
# ---------------------------------------------------------------
mapping = {}
errors = []

with open(tsv_file, newline='', encoding="utf-8") as f:
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)

    try:
        idx_contig = header.index("contig_name_fasta_def")
    except ValueError:
        sys.exit("ERROR: column 'contig_name_fasta_def' not found in TSV")

    for row in reader:
        if len(row) <= 1:
            continue

        contig = row[idx_contig].strip()
        final_id = row[1].strip()

        if contig.startswith(">"):
            contig = contig[1:]
        if final_id.startswith(">"):
            final_id = final_id[1:]

        if contig:
            if contig not in mapping:
                mapping[contig] = final_id

    #print(mapping)
    
# ---------------------------------------------------------------
# Step 2: Parse FASTA and replace headers
# ---------------------------------------------------------------
def fasta_records(path):
    """Yield (header, sequence) tuples from FASTA file"""
    header = None
    seq_chunks = []

    def flush():
        nonlocal header, seq_chunks
        if header is not None:
            yield (header, ''.join(seq_chunks))
        header = None
        seq_chunks = []

    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                for rec in flush():
                    yield rec
                header = line[1:].split()[0]
            else:
                seq_chunks.append(line)
        for rec in flush():
            yield rec

# ---------------------------------------------------------------
# Step 3: Write reformatted FASTA
# ---------------------------------------------------------------
total = 0
mapped = 0

with open(out_fasta, "w", encoding="utf-8") as w:
    for head, seq in fasta_records(fasta_file):
        total += 1
        new_id = mapping.get(head)

        if new_id:
            mapped += 1
            w.write(f">{new_id}\n")
        else:
            errors.append(head)
            w.write(f">{head}\n")

        for i in range(0, len(seq), 60):
            w.write(seq[i:i+60] + "\n")

# ---------------------------------------------------------------
# Step 4: Report
# ---------------------------------------------------------------
print(f"Total sequences: {total}", file=sys.stderr)
print(f"Retrieved sequences: {mapped}", file=sys.stderr)
print(f"Unretrieved sequences: {len(errors)}", file=sys.stderr)

if errors:
    print("Unretrieved contigs:", file=sys.stderr)
    for e in errors:
        print(" -", e, file=sys.stderr)

# BAM: _NA_NA Fasta: |NA|NA
# BAM: ) at the end to rm
# BAM: |LTR Fasta: _LTR
# sed 's/_NA_NA/|NA|NA/g' Final_V2_DUAN_DACCORD_ATIP_DATABASE_non-redundant_150.reformat.fasta > Final_V2_DUAN_DACCORD_ATIP_DATABASE_non-redundant_150.reformat.sed.fasta
