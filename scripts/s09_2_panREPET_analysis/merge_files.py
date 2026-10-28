import pandas as pd
import re
from pathlib import Path
from typing import Optional, List

# ===============================================================
# === CONFIGURATION =============================================
# ===============================================================
BASE_DIR = Path("../../output")

INPUT_DIRS = {
    "Structure": BASE_DIR / "TE_structure",
    "Residuals": BASE_DIR / "TE_residuals"
}

FINAL_OUTPUT = BASE_DIR / "raw_data_Figure4.tsv"

# ===============================================================
# === HELPER FUNCTIONS ==========================================
# ===============================================================

def load_bed(path: str) -> pd.DataFrame:
    """Read BED (9 columns), type columns and drop only 100% identical rows."""
    df = pd.read_csv(path, sep="\t", header=None, dtype=str)
    df = df.iloc[:, :9]
    df.columns = ['chr','start','end','clique_id','accession','copy_name','consensus','clique_size','order']
    for c in ['start','end','clique_id','clique_size']:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    df = df.drop_duplicates(keep='first')
    return df

def load_inter(path: str) -> pd.DataFrame:
    """Read TE-vs-gene intersections (optionally with GO columns)."""
    df = pd.read_csv(path, sep="\t", header=0, dtype=str)
    df.columns = [c.strip() if isinstance(c, str) else c for c in df.columns]
    df = df.rename(columns={'#chr': 'chr'})
    for c in ['start','end','clique_id','clique_size','start_gene','end_gene']:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors='coerce')
    return df

def load_info(path: str) -> pd.DataFrame:
    """Read TE identities; normalize column names."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = [c.strip() if isinstance(c, str) else c for c in df.columns]
    df = df.rename(columns={
        'chromosome': 'chr',
        '#chr': 'chr',
        'TE insertions': 'clique_id'
    })
    for c in ['start','end','clique_size','identity','clique_id']:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors='coerce')
    return df

def extract_attr(attr: Optional[str], key: str = "ID") -> Optional[str]:
    """Extract 'key=value' from a semicolon-separated attribute string."""
    if not isinstance(attr, str):
        return None
    m = re.search(rf'{re.escape(key)}=([^;]+)', attr)
    return m.group(1) if m else None

def unique_join(values: List[str], sep: str = ";") -> Optional[str]:
    """
    Join unique, non-empty strings with 'sep', sorted for reproducibility.
    NOTE: This loses per-gene alignment with other columns by design.
    """
    vals = [v for v in values if isinstance(v, str) and v.strip() != ""]
    return sep.join(sorted(set(vals))) if vals else None

def gene_coord_str(row) -> Optional[str]:
    """Format 'chr:start-end(strand)' when gene coordinates are available."""
    if pd.notna(row.get('chr_gene')) and pd.notna(row.get('start_gene')) and pd.notna(row.get('end_gene')):
        try:
            s = f"{row['chr_gene']}:{int(row['start_gene'])}-{int(row['end_gene'])}"
        except Exception:
            return None
        if isinstance(row.get('strand_gene'), str):
            s += f"({row['strand_gene']})"
        return s
    return None

# ------------------- NEW: GO helpers (English comments) -------------------

def split_go_list(s: Optional[str]) -> List[str]:
    """
    Split a raw GO string like 'GO:0008150,GO:0009888' into a list.
    Ignore '-', 'NA', empty values. Returns [] when no valid codes.
    """
    if not isinstance(s, str):
        return []
    s = s.strip()
    if s in ("", "-", "NA", "None"):
        return []
    return [tok.strip() for tok in s.split(",") if tok.strip()]

def build_gene_GOs_per_TE(group: pd.DataFrame) -> pd.Series:
    """
    For a given TE group (KEYS fixed), build a per-gene GO string aligned to gene_ids.
    - For each gene_id: union of all GO codes across its rows (deduplicated, sorted).
    - Output 'gene_GOs' uses ';' between genes, and ',' inside each gene.
    - If a gene has no GO: output '-' for that gene.
    - Gene order matches unique_join(): sorted unique gene_ids.
    """
    # Collect GO codes per gene_id
    mapping = {}  # gene_id -> Set[GO]
    for _, r in group.iterrows():
        gid = r.get('gene_id')
        if not isinstance(gid, str) or gid.strip() == "":
            continue
        codes = split_go_list(r.get('GOs'))
        if gid not in mapping:
            mapping[gid] = set()
        mapping[gid].update(codes)

    if not mapping:
        return pd.Series({'gene_GOs': None})

    gene_list = sorted(mapping.keys())  # same ordering logic as unique_join
    per_gene = []
    for gid in gene_list:
        codes_sorted = sorted(mapping[gid])
        per_gene.append(",".join(codes_sorted) if codes_sorted else "-")

    return pd.Series({'gene_GOs': ";".join(per_gene)})

# ===============================================================
# === MAIN PIPELINE PER DIRECTORY ===============================
# ===============================================================

def process_te_directory(te_type: str, folder: Path) -> pd.DataFrame:
    print(f"\n➡ Processing folder: {te_type} ({folder})")

    if te_type.lower() == "structure":
        F_INFO = folder / "panREPET_identities_TE_structure_concatUniques.tsv"
    else:
        F_INFO = folder / "panREPET_identities_TE_residuals_concatUniques.tsv"

    F_BED      = folder / "panREPET.bed"
    F_INTER    = folder / "panREPET.bed.intersectVSgenes.sameAss"
    # Prefer GO-augmented intersections if present
    F_INTER_GO = folder / "panREPET.bed.intersectVSgenes.sameAss.GO"

    bed   = load_bed(F_BED)
    inter = load_inter(F_INTER_GO if F_INTER_GO.exists() else F_INTER)
    info  = load_info(F_INFO)

    # 🔑 Keys include clique_id to avoid collapsing distinct clusters
    KEYS = ['accession','chr','start','end','copy_name','consensus','clique_size','clique_id']

    # Keep all info rows; stable sort for reproducibility
    info_sorted = info.sort_values([k for k in KEYS if k in info.columns], kind="mergesort")

    # --- Gene annotation (extract IDs, names, coords) ---
    inter_annot = inter.copy()
    inter_annot['gene_id']    = inter_annot['features_gene'].apply(lambda s: extract_attr(s, "ID"))
    maybe_name                = inter_annot['features_gene'].apply(lambda s: extract_attr(s, "Name"))
    inter_annot['gene_name']  = maybe_name.where(maybe_name.notna(), inter_annot['gene_id'])
    inter_annot['gene_coord'] = inter_annot.apply(gene_coord_str, axis=1)

    # --- Aggregate per TE (same keys), keeping legacy columns ---
    agg_spec = {
        'n_genes_overlapping': ('gene_id',     lambda x: len(set([v for v in x if isinstance(v, str)]))),
        'gene_ids'           : ('gene_id',     lambda x: unique_join(x)),
        'gene_names'         : ('gene_name',   lambda x: unique_join(x)),
        'gene_types'         : ('type_gene',   lambda x: unique_join(x)),
        'gene_strands'       : ('strand_gene', lambda x: unique_join(x)),
        'gene_coords'        : ('gene_coord',  lambda x: unique_join(x)),
    }
    agg = (inter_annot
           .groupby([k for k in KEYS if k in inter_annot.columns], dropna=False)
           .agg(**agg_spec)
           .reset_index()
    )

    # --- NEW: build per-gene GO lists aligned with gene_ids ---
    if 'GOs' in inter_annot.columns:
        gene_GOs_df = (inter_annot
                       .groupby([k for k in KEYS if k in inter_annot.columns], dropna=False)
                       .apply(build_gene_GOs_per_TE)
                       .reset_index())
        # Merge gene_GOs into the aggregated table
        agg = pd.merge(agg, gene_GOs_df, how='left',
                       on=[k for k in KEYS if k in gene_GOs_df.columns])

    # --- Join TE BED + identities + aggregated gene annotations ---
    bed_info = pd.merge(bed, info_sorted, how='left',
                        on=[k for k in KEYS if k in bed.columns and k in info_sorted.columns])
    summary_df = pd.merge(bed_info, agg, how='left',
                          on=[k for k in KEYS if k in bed_info.columns and k in agg.columns])

    # Ensure the column exists even if GO file was absent
    if 'gene_GOs' not in summary_df.columns:
        summary_df['gene_GOs'] = pd.NA

    summary_df["TE_type"] = te_type
    print(f"  ✔ {summary_df.shape[0]} rows processed for {te_type}")
    return summary_df

# ===============================================================
# === RUN AND MERGE =============================================
# ===============================================================

all_summaries = [process_te_directory(label, folder) for label, folder in INPUT_DIRS.items()]
merged_all = pd.concat(all_summaries, ignore_index=True)

# Cast numeric columns cleanly (avoid 656.0)
for c in ['start','end','clique_id','clique_size']:
    if c in merged_all.columns:
        merged_all[c] = pd.to_numeric(merged_all[c], errors='coerce').astype('Int64')

# ===============================================================
# === FINAL FORMATTING & EXPORT =================================
# ===============================================================

final_order = [
    "TE_type", "accession", "clique_size", "chr", "start", "end",
    "consensus", "copy_name",
    "clique_id", "order", "identity",
    "n_genes_overlapping", "gene_ids",
    "gene_types", "gene_strands", "gene_coords",
    # NEW: per-gene GO lists (aligned to gene_ids; ';' between genes, ',' within a gene)
    "gene_GOs",
]

existing_cols = [c for c in final_order if c in merged_all.columns]
merged_all = merged_all[existing_cols]

# Stable final sort
merged_all = merged_all.sort_values(by=existing_cols, kind="mergesort")

merged_all.to_csv(FINAL_OUTPUT, sep="\t", index=False)

print("\n🎯 Clean merge completed successfully!")
print(f"Final file: {FINAL_OUTPUT}")
print(f"Total rows: {merged_all.shape[0]}")
print(f"Columns: {len(merged_all.columns)}")
print("Column order:")
print(", ".join(merged_all.columns))
