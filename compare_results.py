# This file compares S-PrediXcan results between:
# 1. Project's elastic net vs Kho et al. 2021 (Table 1)
# 2. Project's elastic net vs Project's mashr

import pandas as pd
import numpy as np
import os

# list of six GTEx v8 tissues used in the analysis
TISSUES = [
    "Adipose_Subcutaneous",
    "Adipose_Visceral_Omentum",
    "Whole_Blood",
    "Uterus",
    "Ovary",
    "Vagina",
]

# Genes reported as significant in Kho et al. (2021) Table 1
PAPER_GENES = ["GLDN", "CYP19A1", "HEY2", "SKAP1", "SNX11", "EVI2A", "EEFSEC"]

# Tissue abbreviations for table output
TISSUE_SHORT = {
    "Adipose_Subcutaneous": "Adi.SC",
    "Adipose_Visceral_Omentum": "Adi.Vis",
    "Whole_Blood": "WB",
    "Uterus": "Uterus",
    "Ovary": "Ovary",
    "Vagina": "Vagina",
}


# Load all results
def load_spredixcan(folder, prefix):
    results = {}
    for tissue in TISSUES:
        # Build the expected file path for this tissue's S-PrediXcan output
        path = f"{folder}/{prefix}{tissue}.csv"
        if os.path.exists(path):
            t = pd.read_csv(path)
            # Tag each row with its source tissue for later filtering
            t["tissue"] = tissue
            # Compute per-tissue Bonferroni threshold based on number of unique genes tested
            bonf = 0.05 / t["gene"].nunique()
            # Store as tuple: (dataframe, bonferroni threshold)
            results[tissue] = (t, bonf)
    return results


# Load elastic net S-PrediXcan results (from aim 2)
en = load_spredixcan("results/aim2", "spredixcan_allEC_")
# Load mashr S-PrediXcan results (Aim 3)
mashr = load_spredixcan("results/aim3/mashr", "spredixcan_mashr_")
# Load S-MultiXcan joint-tissue results
smx = pd.read_csv("results/aim2/smultixcan_allEC.csv", sep="\t")
# Compute Bonferroni threshold across all genes tested by S-MultiXcan
smx_bonf = 0.05 / len(smx)


# TABLE 1: S-MultiXcan significant genes
print("\n\nTABLE 1: S-MultiXcan significant genes (elastic net)")
print(f"Bonferroni threshold: {smx_bonf:.2e}")
print("-" * 40)

# Filter to genes passing Bonferroni threshold, sorted by p-value
sig_smx = smx[smx["pvalue"] < smx_bonf].sort_values("pvalue")
# Select and rename columns for display
t1 = sig_smx[["gene_name", "pvalue", "n", "z_min", "z_max", "t_i_best"]].copy()
t1.columns = ["Gene", "P-value", "N tissues", "Z min", "Z max", "Best tissue"]
# Format p-values and Z-scores for readability
t1["P-value"] = t1["P-value"].apply(lambda x: f"{x:.2e}")
t1["Z min"]   = t1["Z min"].apply(lambda x: f"{x:.2f}")
t1["Z max"]   = t1["Z max"].apply(lambda x: f"{x:.2f}")
# Flag whether each gene was reported in Kho et al. (2021) or is a novel finding
t1["In paper"] = t1["Gene"].apply(lambda g: "YES" if g in PAPER_GENES else "NEW")
print(t1.to_string(index=False))


# TABLE 2: Per-tissue Z-scores EN vs Mashr for paper genes
print("\n\nTABLE 2: Per-tissue Z-scores — Elastic Net vs Mashr (paper genes)")
print("Asterisk (*) = Bonferroni significant in that tissue/model")
print("-" * 40)

rows = []
for gene in PAPER_GENES:
    for tissue in TISSUES:
        row = {"Gene": gene, "Tissue": TISSUE_SHORT[tissue]}

        # Elastic net Z-score for this gene-tissue pair
        if tissue in en:
            t, bonf = en[tissue]
            # Look up this gene in the tissue dataframe
            g = t[t["gene_name"] == gene]
            if len(g) > 0:
                z = g["zscore"].values[0]
                p = g["pvalue"].values[0]
                # Append asterisk if gene passes Bonferroni threshold in this tissue
                sig = "*" if p < bonf else ""
                row["EN Z"] = f"{z:.2f}{sig}"
            else:
                # Gene has no trained prediction model in this tissue
                row["EN Z"] = "—"
        else:
            # Tissue file not found
            row["EN Z"] = "—"

        # Mashr Z-score for this gene-tissue pair
        if tissue in mashr:
            t, bonf = mashr[tissue]
            g = t[t["gene_name"] == gene]
            if len(g) > 0:
                z = g["zscore"].values[0]
                p = g["pvalue"].values[0]
                sig = "*" if p < bonf else ""
                row["Mashr Z"] = f"{z:.2f}{sig}"
            else:
                row["Mashr Z"] = "—"
        else:
            row["Mashr Z"] = "—"

        rows.append(row)

t2 = pd.DataFrame(rows)
print(t2.to_string(index=False))


# TABLE 3: Replication vs Kho et al. 2021
print("\n\nTABLE 3: Replication vs Kho et al. 2021")
print("-" * 40)

# Hardcoded reference values from Kho et al. (2021) Table 1
paper_smultixcan = {
    "EEFSEC":  {"p": 1.10e-6,  "best_tissue": "Multiple",   "coloc": "Yes"},
    "HEY2":    {"p": 9.94e-9,  "best_tissue": "Ovary",       "coloc": "Yes"},
    "GLDN":    {"p": 1.34e-12, "best_tissue": "Multiple",    "coloc": "No"},
    "CYP19A1": {"p": 9.52e-12, "best_tissue": "Adipose SC",  "coloc": "Yes"},
    "EVI2A":   {"p": 1.50e-6,  "best_tissue": "Multiple",    "coloc": "Yes"},
    "SKAP1":   {"p": 7.27e-9,  "best_tissue": "Whole Blood", "coloc": "Yes"},
    "SNX11":   {"p": 5.41e-7,  "best_tissue": "Multiple",    "coloc": "No"},
}

rep_rows = []
for gene in PAPER_GENES:
    pub = paper_smultixcan[gene]
    # Look up this gene in our S-MultiXcan results
    our = smx[smx["gene_name"] == gene]
    if len(our) > 0:
        our_p    = our["pvalue"].values[0]
        our_best = our["t_i_best"].values[0]
        # Mark as replicated if it passes our Bonferroni threshold
        our_sig  = "YES" if our_p < smx_bonf else "NO"
    else:
        our_p    = None
        our_best = "—"
        our_sig  = "NOT FOUND"

    rep_rows.append({
        "Gene":             gene,
        "Paper P":          f"{pub['p']:.2e}",
        "Paper best tissue": pub["best_tissue"],
        "Our P":            f"{our_p:.2e}" if our_p else "—",
        # Abbreviate tissue names for compact display
        "Our best tissue":  our_best.replace("Adipose_Subcutaneous", "Adi.SC")
                                    .replace("Adipose_Visceral_Omentum", "Adi.Vis")
                                    .replace("Whole_Blood", "WB") if our_best else "—",
        "Replicated":       our_sig,
        # Direction match only meaningful if gene was replicated
        "Direction match":  "YES" if our_p and our_p < smx_bonf else "—",
    })

t3 = pd.DataFrame(rep_rows)
print(t3.to_string(index=False))


# TABLE 4: EN vs Mashr significant gene overlap per tissue
print("\n\nTABLE 4: Significant genes — Elastic Net vs Mashr per tissue")
print("-" * 40)

for tissue in TISSUES:
    en_sig, mashr_sig = set(), set()

    # Get set of Bonferroni-significant genes under elastic net for this tissue
    if tissue in en:
        t, bonf = en[tissue]
        en_sig = set(t[t["pvalue"] < bonf]["gene_name"].tolist())

    # Get set of Bonferroni-significant genes under mashr for this tissue
    if tissue in mashr:
        t, bonf = mashr[tissue]
        mashr_sig = set(t[t["pvalue"] < bonf]["gene_name"].tolist())

    # Compute set intersections for overlap categorization
    both       = en_sig & mashr_sig       # significant in both models
    en_only    = en_sig - mashr_sig       # significant in elastic net only
    mashr_only = mashr_sig - en_sig       # significant in mashr only

    print(f"\n{TISSUE_SHORT[tissue]}:")
    print(f"  EN only    ({len(en_only):2d}): {', '.join(sorted(en_only)) or '—'}")
    print(f"  Both       ({len(both):2d}): {', '.join(sorted(both)) or '—'}")
    print(f"  Mashr only ({len(mashr_only):2d}): {', '.join(sorted(mashr_only)) or '—'}")