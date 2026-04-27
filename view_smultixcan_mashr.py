# S-MultiXcan Mashr Model

import pandas as pd
import numpy as np
print("S-MULTIXCAN MASHR MODEL RESULTS\n")


df = pd.read_csv("results/aim3/smultixcan_mashr.csv", sep="\t")

n_genes = len(df)
threshold = 0.05 / n_genes
print(f"Total genes tested: {n_genes:,}")
print(f"Bonferroni threshold: {threshold:.2e}")

sig = df[df["pvalue"] < threshold].sort_values("pvalue")
print(f"Significant genes: {len(sig)}\n")
print(sig[["gene", "gene_name", "pvalue", "n", "n_indep", "z_min", "z_max", "t_i_best"]].to_string(index=False))

paper_genes = ["GLDN", "CYP19A1", "HEY2", "SKAP1", "SNX11", "EVI2A", "EEFSEC"]
print("\nPaper genes in mashr S-MultiXcan:")
paper = df[df["gene_name"].isin(paper_genes)].sort_values("pvalue")
print(paper[["gene_name", "pvalue", "t_i_best"]].to_string(index=False))