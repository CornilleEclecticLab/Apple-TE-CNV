#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Goal:
  Build a new FASTA whose headers are EXACTLY the names listed in contig_file,
  in the same order. For each BAM name, find the corresponding FASTA record
  using lenient matching rules, then write that record renamed to the BAM name.
  Sequences without a BAM counterpart are excluded from the output FASTA.

Matching rules (order-independent for matching only):
  - Strip any trailing ')'
  - Remove commas ','
  - Replace pipes '|' with underscores '_'
  - Drop a single trailing '_1'
Naming rule:
  - The written FASTA header is EXACTLY the BAM name (no normalization)
Outputs:
  - New FASTA with exact BAM names, same order as contig_file
  - TSV mapping: bam_name  <tab>  source_fasta_id
  - Unmatched BAM names
  - Unused FASTA records
"""

from Bio import SeqIO
from collections import defaultdict
from typing import List, Dict

# ---- Paths ----
fasta_file = "../../data/TElibrary/Final_V2_DUAN_DACCORD_ATIP_DATABASE_non-redundant_150.reformat.reformat.fasta"
contig_file = "../../data/BAM/seqNames_from_BAM.txt"
out_fasta   = "../../data/TElibrary/Final_V2_DUAN_DACCORD_ATIP_DATABASE_non-redundant_150.reformat.reformat.reformat.fasta"
map_tsv     = "../../data/BAM/rename_map_bam_to_fasta.tsv"
unmatched_bam_txt   = "../../data/BAM/unmatched_bam_names.txt"
unused_fasta_txt    = "../../data/BAM/unused_fasta_records.txt"
ambiguities_txt     = "../../data/BAM/ambiguous_matches.txt"

# ---- Normalization helpers for MATCHING (not for naming) ----
def strip_trailing_parens(s: str) -> str:
    while s.endswith(")"):
        s = s[:-1]
    return s

def normalize_for_match(s: str) -> str:
    """Order-independent normalization to compare BAM and FASTA ids."""
    s = s.strip()
    s = strip_trailing_parens(s)
    s = s.replace(",", "").replace("|", "_")
    if s.endswith("_1"):
        s = s[:-2]                    # drop only ONE trailing '_1'
    return s

# ---- Load FASTA and index by normalized id ----
fasta_by_norm: Dict[str, List[SeqIO.SeqRecord]] = defaultdict(list)
all_fasta_ids: List[str] = []

for rec in SeqIO.parse(fasta_file, "fasta"):
    all_fasta_ids.append(rec.id)
    key = normalize_for_match(rec.id)
    fasta_by_norm[key].append(rec)

# ---- Read BAM names in order ----
with open(contig_file) as f:
    bam_names = [line.strip() for line in f if line.strip()]

used_fasta_ids = set()
ambiguous_lines = []
mapping_rows = []
unmatched_bam = []

def pick_best_candidate(bam: str, candidates: List[SeqIO.SeqRecord]) -> SeqIO.SeqRecord:
    """
    Tie-breaker if multiple FASTA records share the same normalized key.
    Preference:
      1) Exact after simple pipe/comma removal (without dropping '_1')
      2) Exact after stripping trailing ')'
      3) Shortest original id (heuristic, deterministic)
    """
    simple = strip_trailing_parens(bam).replace(",", "").replace("|", "_")
    # 1) direct equality on transformed-but-without '_1' drop
    for r in candidates:
        if r.id == simple:
            return r
    # 2) equality after only stripping trailing ')'
    only_paren = strip_trailing_parens(bam)
    for r in candidates:
        if r.id == only_paren:
            return r
    # 3) fallback: shortest id
    return sorted(candidates, key=lambda r: (len(r.id), r.id))[0]

# ---- Build the output FASTA in BAM order ----
written = 0
with open(out_fasta, "w") as out:
    for bam in bam_names:
        key = normalize_for_match(bam)
        # Filter out already-used records to keep 1-to-1 mapping
        candidates = [r for r in fasta_by_norm.get(key, []) if r.id not in used_fasta_ids]

        if not candidates:
            unmatched_bam.append(bam)
            continue

        chosen = candidates[0] if len(candidates) == 1 else pick_best_candidate(bam, candidates)
        if len(candidates) > 1:
            ambiguous_lines.append(f"BAM '{bam}' matched multiple FASTA ids: "
                                   + ", ".join(sorted(r.id for r in candidates))
                                   + f" -> chosen: {chosen.id}")

        # Keep original source id for mapping before we rename the record
        source_id = chosen.id

        # Rename to EXACT BAM name
        chosen.id = bam
        chosen.name = bam
        chosen.description = bam

        SeqIO.write(chosen, out, "fasta")
        used_fasta_ids.add(source_id)
        mapping_rows.append((bam, source_id))
        written += 1

# ---- Write mapping and reports ----
with open(map_tsv, "w") as w:
    w.write("bam_name\tsource_fasta_id\n")
    for bam, src in mapping_rows:
        w.write(f"{bam}\t{src}\n")

with open(unmatched_bam_txt, "w") as w:
    for bam in unmatched_bam:
        w.write(bam + "\n")

with open(unused_fasta_txt, "w") as w:
    for fid in all_fasta_ids:
        if fid not in used_fasta_ids:
            w.write(fid + "\n")

if ambiguous_lines:
    with open(ambiguities_txt, "w") as w:
        for line in ambiguous_lines:
            w.write(line + "\n")

# ---- Console summary ----
print(f"Requested BAM names: {len(bam_names)}")
print(f"Written records    : {written}")
print(f"Unmatched BAM names: {len(unmatched_bam)} (see {unmatched_bam_txt})")
print(f"Unused FASTA ids   : {len([fid for fid in all_fasta_ids if fid not in used_fasta_ids])} (see {unused_fasta_txt})")
print(f"Mapping TSV        : {map_tsv}")
if ambiguous_lines:
    print(f"Ambiguities        : {len(ambiguous_lines)} (see {ambiguities_txt})")
print(f"Output FASTA       : {out_fasta}")
