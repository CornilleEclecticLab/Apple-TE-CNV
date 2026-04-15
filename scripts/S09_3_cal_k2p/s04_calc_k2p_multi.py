#!/usr/bin/env python3
import os, glob, math
import statistics

def is_transition(ref, alt):
    return (ref, alt) in [('A','G'), ('G','A'), ('C','T'), ('T','C')]

def is_transversion(ref, alt):
    return ref != alt and ref in "ACGT" and alt in "ACGT" and not is_transition(ref, alt)

def compute_k2p(pileup_file):
    results = {}
    counters = {}

    with open(pileup_file) as f:
        for line in f:
            # chrom, pos, ref, depth, bases, qual = line.strip().split("\t")
            cutting = line.strip().split("\t")
            if len(cutting) < 6:
                continue  # ligne incomplète
            chrom = cutting[0]
            pos = cutting[1]
            ref = cutting[2].upper()
            depth_field = cutting[3]
            bases = cutting[4]
            # qual = cutting[5]

            if ref not in "ACGT":
                continue

            try:
                depth = int(depth_field)
            except ValueError:
                continue

            if depth == 0:
                continue

            if chrom not in counters:
                counters[chrom] = {"ti": 0, "tv": 0, "sites": 0}

            counts = {b: bases.upper().count(b) for b in "ACGT"}
            for alt, n in counts.items():
                if alt == ref or n == 0:
                    continue
                if is_transition(ref, alt):
                    counters[chrom]["ti"] += n
                elif is_transversion(ref, alt):
                    counters[chrom]["tv"] += n
            counters[chrom]["sites"] += depth

    for chrom, vals in counters.items():
        sites = vals["sites"]
        if sites == 0:
            results[chrom] = float("nan")
            continue
        P = vals["ti"] / sites
        Q = vals["tv"] / sites
        try:
            K2P = -0.5 * math.log((1 - 2*P - Q) * math.sqrt(1 - 2*Q))
        except ValueError:
            K2P = float("nan")
        results[chrom] = K2P

    return results

# --- I/O ---
INPUT_PATTERN = "../../output/pileups/*.pileup"  # adapte si besoin ("*.pileup" pour le dossier courant)
OUTPUT_DIR = "../../output/k2p"
os.makedirs(OUTPUT_DIR, exist_ok=True)

pileups = glob.glob(INPUT_PATTERN)
print(pileups)

#print("Individu\tConsensus\tK2P")
all_results = {}

for pileup in pileups:
    indiv = os.path.basename(pileup).replace(".pileup", "")
    results = compute_k2p(pileup)

    # sortie par échantillon : *.pileup.k2p.tsv
    out_path = os.path.join(OUTPUT_DIR, f"{indiv}.pileup.k2p.tsv")
    with open(out_path, "w") as out:
        out.write("Consensus\tK2P\n")
        for consensus, k2p in results.items():
            out.write(f"{consensus}\t{k2p}\n")

    # impression console comme avant
    # print(indiv)
    for consensus, k2p in results.items():
        # print(f"{indiv}\t{consensus}\t{k2p}")
        all_results.setdefault(consensus, []).append(k2p)

# Résumé global par consensus (moyenne sur individus) – affichage console
# print("\n# Résumé par consensus (moyenne K2P sur tous les individus)")
# print("Consensus\tMoyenneK2P\tN_individus")
for consensus, values in all_results.items():
    clean_vals = [v for v in values if not math.isnan(v)]
    if clean_vals:
        mean_val = statistics.mean(clean_vals)
        # print(f"{consensus}\t{mean_val}\t{len(clean_vals)}")
