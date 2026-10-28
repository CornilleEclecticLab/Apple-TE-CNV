#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import glob
import math
import csv
from collections import defaultdict

import matplotlib.pyplot as plt

# Use SciPy if available for Welch's t-test; otherwise fallback to a permutation test
USE_SCIPY = True
try:
    from scipy import stats
except Exception:
    USE_SCIPY = False

# -------------- Helpers --------------

TE_ID_RE = re.compile(r"(TE_\d+)", re.IGNORECASE)

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path

def extract_te_id(s: str):
    """Extract TE_###### from strings like '>TE_000001_noCat' or 'TE_000001_LTR_Copia'."""
    if not s:
        return None
    s = s.strip().lstrip(">")
    m = TE_ID_RE.search(s)
    return m.group(1) if m else None

def load_te_list(path):
    """Load one TE per line; keep only 'TE_######' core IDs."""
    ids = set()
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            te = extract_te_id(line)
            if te:
                ids.add(te)
    return ids

def read_annotation_build_maps(mapping_tsv):
    """
    Read V2_Pastec_TE_sorter_Annotation_final_V1.tsv and build:
      - te_id -> old_contig_name_bis
      - old_contig_name_bis -> te_id
    Priority is on 'old_contig_name_bis' as the primary key.
    """
    te_to_oldbis = {}
    oldbis_to_te = {}

    with open(mapping_tsv, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        col_te = "contig_name_fasta_def"
        col_oldbis = "old_contig_name_bis"
        col_old = "old_contig_name"  # kept for optional fallback in mapping function (not used for grouping)

        if not reader.fieldnames or col_te not in reader.fieldnames or col_oldbis not in reader.fieldnames:
            raise ValueError(
                "Expected columns missing in mapping TSV. Need at least: "
                f"'{col_te}' and '{col_oldbis}'. Found: {reader.fieldnames}"
            )

        for row in reader:
            te_id = extract_te_id(row.get(col_te, ""))
            old_bis = row.get(col_oldbis, "").strip().lstrip(">")
            old = row.get(col_old, "").strip().lstrip(">")

            if te_id and old_bis:
                # one TE may appear multiple times with the same old_bis; first come first kept
                te_to_oldbis.setdefault(te_id, old_bis)
                # one old_bis should normally map to one TE; if multiple, keep first
                oldbis_to_te.setdefault(old_bis, te_id)

            # We do not build an old->te map for grouping, but we will keep 'old' handy for a soft fallback.

    return te_to_oldbis, oldbis_to_te

def read_all_k2p_rows(k2p_dir="../output/k2p"):
    """
    Concatenate all K2P TSV files with columns 'Consensus' and 'K2P'.
    Return a list of dicts: {file, consensus_raw, k2p}
    """
    rows = []
    for path in glob.glob(os.path.join(k2p_dir, "*.pileup.k2p.tsv")):
        with open(path, encoding="utf-8") as f:
            header = f.readline().strip().split("\t")
            try:
                idx_cons = header.index("Consensus")
                idx_k2p = header.index("K2P")
            except ValueError:
                if len(header) >= 2:
                    idx_cons, idx_k2p = 0, 1
                else:
                    continue
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) <= max(idx_cons, idx_k2p):
                    continue
                consensus_raw = parts[idx_cons]
                try:
                    k2p = float(parts[idx_k2p])
                except ValueError:
                    continue
                rows.append({
                    "file": os.path.basename(path),
                    "consensus_raw": consensus_raw,
                    "k2p": k2p
                })
    return rows

def consensus_to_oldbis_key(consensus_raw: str):
    """
    Normalize consensus names to match 'old_contig_name_bis' style.
    DO NOT TRIM. Only normalize the separator so it matches annotation style:
      'f2_NEW_CONTIG_MERGE_147|NA|NA'  -> 'f2_NEW_CONTIG_MERGE_147_NA_NA'
      'f2_NEW_CONTIG_MERGE_147'       -> 'f2_NEW_CONTIG_MERGE_147'
    Removes trailing ')' and replaces '|' by '_'.
    """
    if consensus_raw is None:
        return None
    # Keep the full string; replace '|' by '_' to harmonize with old_contig_name_bis.
    key = consensus_raw.strip().lstrip(">").replace("|", "_")
    # remove ')' at the end
    key = key.rstrip(")")
    return key

def build_oldbis_to_group(te_structure_ids, te_residual_ids, te_to_oldbis):
    """
    Build mapping: old_contig_name_bis -> group ('structure' or 'residual').
    This uses TE IDs from the lists, then maps to the corresponding old_contig_name_bis.
    If the same old_bis is assigned by two different TEs to different groups, last wins
    (but this should be rare if the mapping is clean).
    """
    oldbis_to_group = {}

    for te in te_structure_ids:
        old_bis = te_to_oldbis.get(te)
        if old_bis:
            oldbis_to_group[old_bis] = "structure"

    for te in te_residual_ids:
        old_bis = te_to_oldbis.get(te)
        if old_bis:
            # If an overlap occurs, you can choose a policy. Here, we keep 'residual' if it overwrites.
            # Change to a guarded policy if you want to detect/report conflicts instead.
            if old_bis in oldbis_to_group and oldbis_to_group[old_bis] != "residual":
                # report, but overwrite (or choose to skip)
                # print(f"Conflict on {old_bis}: was {oldbis_to_group[old_bis]}, now residual")
                pass
            oldbis_to_group[old_bis] = "residual"

    return oldbis_to_group

def summarize(values):
    import statistics
    clean = [v for v in values if v is not None and not math.isnan(v)]
    if not clean:
        return {"n": 0, "mean": float("nan"), "median": float("nan")}
    return {"n": len(clean), "mean": statistics.mean(clean), "median": statistics.median(clean)}

def welch_ttest(a, b):
    if not USE_SCIPY or len(a) < 2 or len(b) < 2:
        return (None, None)
    t, p = stats.ttest_ind(a, b, equal_var=False, nan_policy='omit')
    return (t, p)

def permutation_test_mean(a, b, n_iter=10000, random_state=42):
    import random, statistics as st
    rnd = random.Random(random_state)
    a = [x for x in a if x is not None and not math.isnan(x)]
    b = [x for x in b if x is not None and not math.isnan(x)]
    if not a or not b:
        return (float("nan"), float("nan"))
    diff_obs = abs(st.mean(a) - st.mean(b))
    pooled = a + b
    n_a = len(a)
    count = 0
    for _ in range(n_iter):
        rnd.shuffle(pooled)
        a_sim = pooled[:n_a]
        b_sim = pooled[n_a:]
        diff_sim = abs(st.mean(a_sim) - st.mean(b_sim))
        if diff_sim >= diff_obs:
            count += 1
    p_val = (count + 1) / (n_iter + 1)
    return (diff_obs, p_val)

def annotate_pvalue(ax, x1, x2, y, text, h=0.02):
    """
    Draws a p-value annotation above two groups.
    """
    # Compute the relative height according to the scale (linear or log)
    if ax.get_yscale() == "log":
        return
    else:
        y_line = y + h
        y_text = y + 1.5*h

    # Draw the connecting line
    ax.plot([x1, x1, x2, x2], [y, y_line, y_line, y], lw=1, c="k")

    # Place the label with a reasonable margin
    ax.text((x1 + x2) * 0.5, y_text, text, ha='center', va='bottom', fontsize=10)
    
def plot_violin_with_stats(groups_dict, out_png, test_label, p_value, means=None):
    labels = []
    data = []
    for label in ("structure", "residual"):
        if label in groups_dict and groups_dict[label]:
            labels.append(label)
            data.append(groups_dict[label])

    if len(data) < 2:
        print("Not enough data to plot both groups.")
        return

    plt.figure(figsize=(6, 7), dpi=160)
    ax = plt.gca()
    ax.violinplot(data, showmeans=False, showmedians=True, showextrema=True)
    #for i, d in enumerate(data, start=1):
    #    plt.scatter([i]*len(d), d, color='black', alpha=0.2, s=5)
    ax.set_yscale('log')
    ax.set_ylabel("K2P (log scale)")
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_title("K2P distribution by group (TE_structure vs TE_residuals)")

    if means is None:
        means = [sum(d)/len(d) for d in data]
    ax.plot([1, 2], means, marker="o", linestyle="-", linewidth=1)

    ymax = max(max(d) for d in data if d)
    ymin = min(min(d) for d in data if d)
    yspan = ymax - ymin if ymax > ymin else 1.0
    y_annot = ymax + 0.06 * yspan

    if p_value is not None and not math.isnan(p_value):
        text = f"{test_label}: p = {p_value:.8g}"
    else:
        text = f"{test_label}: p = NA"

    annotate_pvalue(ax, 1, 2, y_annot, text, h=0.03 * yspan)

    ax.set_ylim(bottom=max(0.0, ymin - 0.05*yspan), top=y_annot + 0.08*yspan)
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()
    print(f"Figure saved: {out_png}")

# -------------- Main --------------

def main():
    # Paths
    k2p_dir = "../../output/k2p"
    residuals_path = "../../data/TElists/250923_TE_residuals.txt"
    structure_path = "../../data/TElists/250923_TE_structure.txt"
    mapping_tsv = "../../data/TElibrary/V2_Pastec_TE_sorter_Annotation_final_V1.tsv"
    figures_dir = ensure_dir("../Figures")
    out_dir = ensure_dir(k2p_dir)

    # Load TE lists (as TE_###### IDs)
    te_residual_ids = load_te_list(residuals_path)
    te_structure_ids = load_te_list(structure_path)
    print(f"- TE residuals: {len(te_residual_ids)} IDs")
    print(f"- TE structure: {len(te_structure_ids)} IDs")

    # Read annotation and build maps
    te_to_oldbis, oldbis_to_te = read_annotation_build_maps(mapping_tsv)
    print(f"- TE mapped to old_contig_name_bis: {len(te_to_oldbis)}")

    # Build oldbis -> group using TE lists
    oldbis_to_group = build_oldbis_to_group(te_structure_ids, te_residual_ids, te_to_oldbis)
    print(f"- old_contig_name_bis assigned to a group: {len(oldbis_to_group)}")

    # Read K2P rows
    k2p_rows = read_all_k2p_rows(k2p_dir)
    print(f"- K2P rows loaded: {len(k2p_rows)}")

    # Merge and group WITHOUT trimming; only '|' -> '_' normalization
    merged = []
    groups = defaultdict(list)

    # For QA: count unmatched consensus
    unmatched = 0
    unmatched_list = []

    for row in k2p_rows:
        consensus_raw = row["consensus_raw"]
        k2p = row["k2p"]
        if math.isnan(k2p):
            continue
        k2p = max(0.000, k2p)

        # Build the exact key expected by 'old_contig_name_bis' style
        key_oldbis_style = consensus_to_oldbis_key(consensus_raw)

        # Assign group via old_contig_name_bis
        group = oldbis_to_group.get(key_oldbis_style)

        if not group:
            unmatched += 1
            unmatched_list.append({
                "file": row["file"],
                "consensus": consensus_raw,
                "key_oldbis_style": key_oldbis_style
            })
            continue

        merged.append({
            "file": row["file"],
            "consensus": consensus_raw,
            "old_contig_name_bis": key_oldbis_style,
            "group": group,
            "k2p": k2p
        })

        groups[group].append(k2p)

    print(f"- Merged rows used: {len(merged)}; unmatched (skipped): {unmatched}")


    # Write merged table (sorted)
    merged_tsv = os.path.join(out_dir, "k2p_sample_x_consensus.tsv")

    # Sort by: group, file, consensus, old_contig_name_bis, k2p
    merged_sorted = sorted(
        merged,
        key=lambda r: (r["group"], r["file"], r["consensus"], r["old_contig_name_bis"], r["k2p"])
    )
    
    with open(merged_tsv, "w", encoding="utf-8") as out:
        out.write("file\tconsensus\told_contig_name_bis\tgroup\tk2p\n")
        for r in merged_sorted:
            out.write(f"{r['file']}\t{r['consensus']}\t{r['old_contig_name_bis']}\t{r['group']}\t{r['k2p']}\n")
    
    print(f"Merged table: {merged_tsv} ({len(merged_sorted)} rows)")
    
    
    # Summary
    def write_summary(groups_dict, out_path):
        import statistics
    
        def stats_extended(values):
            clean = [v for v in values if v is not None and not math.isnan(v)]
            if not clean:
                return {"n": 0, "mean": float("nan"), "median": float("nan"),
                        "var": float("nan"), "min": float("nan"), "max": float("nan")}
            return {
                "n": len(clean),
                "mean": statistics.mean(clean),
                "median": statistics.median(clean),
                "var": statistics.variance(clean) if len(clean) > 1 else 0.0,
                "min": min(clean),
                "max": max(clean)
            }
    
        sum_struct = stats_extended(groups_dict.get("structure", []))
        sum_resid  = stats_extended(groups_dict.get("residual", []))
    
        with open(out_path, "w", encoding="utf-8") as out:
            out.write("group\tn\tmean\tmedian\tvariance\tmin\tmax\n")
            out.write(f"structure\t{sum_struct['n']}\t{sum_struct['mean']}\t{sum_struct['median']}\t"
                      f"{sum_struct['var']}\t{sum_struct['min']}\t{sum_struct['max']}\n")
            out.write(f"residual\t{sum_resid['n']}\t{sum_resid['mean']}\t{sum_resid['median']}\t"
                      f"{sum_resid['var']}\t{sum_resid['min']}\t{sum_resid['max']}\n")
    
        print(f"Summary saved: {out_path}")

    summary_tsv = os.path.join(out_dir, "summary_k2p.tsv")
    write_summary(groups, summary_tsv)

    # Statistical test
    a = groups.get("structure", [])
    b = groups.get("residual", [])
    p_value = float("nan")
    test_label = "Mean comparison"

    if a and b:
        means = [sum(a)/len(a), sum(b)/len(b)]
        if USE_SCIPY:
            t_stat, p_val = welch_ttest(a, b)
            test_label = "Welch's t-test"
            p_value = p_val if p_val is not None else float("nan")
            print(f"{test_label}: p = {p_value:.8e}")
        else:
            diff, p_val = permutation_test_mean(a, b, n_iter=10000, random_state=42)
            test_label = "Permutation test (means)"
            p_value = p_val
            print(f"{test_label}: p = {p_value} (10,000 shuffles; SciPy not found)")
    else:
        means = None
        print("Not enough data in both groups for a mean-comparison test.")

    # Plot
    fig_path = os.path.join(ensure_dir("../../Figures"), "violin_TE_groups.png")
    plot_violin_with_stats(groups, fig_path, test_label, p_value, means=means)

    print("Done.")


    #unmatched, why ?
    unmatched_tsv = os.path.join(out_dir, "unmatched_consensus.tsv")
    with open(unmatched_tsv, "w", encoding="utf-8") as out:
        out.write("file\tconsensus\tkey_oldbis_style\n")
        for r in unmatched_list:
            out.write(f"{r['file']}\t{r['consensus']}\t{r['key_oldbis_style']}\n")
    print(f"Unmatched details saved: {unmatched_tsv} ({len(unmatched_list)} rows)")
    

if __name__ == "__main__":
    main()
