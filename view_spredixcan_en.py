import pandas as pd
import os

print("S-PREDIXCAN ELASTIC NET RESULTS\n")

TISSUES = [
    "Adipose_Subcutaneous",
    "Adipose_Visceral_Omentum",
    "Whole_Blood",
    "Uterus",
    "Ovary",
    "Vagina",
] # list of tissues analyzed

paper_genes = ["GLDN", "CYP19A1", "HEY2", "SKAP1", "SNX11", "EVI2A", "EEFSEC"] # genes highlighted in Kho et al. 2021
all_tissue_results = [] # list to store results from all tissues

for tissue in TISSUES: # loop through each tissue and read the corresponding S-PrediXcan results
    path = f"results/aim2/spredixcan_allEC_{tissue}.csv" # path to the S-PrediXcan results for this tissue
    if not os.path.exists(path): # check if the file exists
        print(f"\n  {tissue}: file not found") # if the file is missing, print a message and skip to the next tissue
        continue

    t = pd.read_csv(path) # read the S-PrediXcan results into a DataFrame
    t.columns = t.columns.str.strip() # remove any leading/trailing whitespace from column names
    t["tissue"] = tissue # add a column for tissue name

    bonf_t = 0.05 / t["gene"].nunique() # calculate the Bonferroni threshold for this tissue based on the number of unique genes tested
    sig_t  = t[t["pvalue"] < bonf_t].sort_values("pvalue") # filter for significant genes and sort by p-value

    print(f"\n{tissue}")
    print(f" Genes tested : {t['gene'].nunique():,}") # print the number of unique genes tested with comma as thousands separator
    print(f" Bonferroni threshold: {bonf_t:.2e}") # print the Bonferroni threshold in scientific notation
    print(f" Significant genes : {len(sig_t)}") # print the number of significant genes
    if len(sig_t) > 0: # if there are significant genes, print the top hits
        print(f" Top hits:")
        print(sig_t[["gene_name", "zscore", "pvalue", "effect_size"]].head(10).to_string(index=False)) # print the top 10 significant genes with their z-score, p-value, and effect size, without the index

    all_tissue_results.append(t) # add the results for this tissue to the list

# Per-tissue Z-scores for paper genes
print("\nPer-tissue Z-scores for Kho et al. 2021 genes\n")
all_df = pd.concat(all_tissue_results, ignore_index=True) # combine results from all tissues into a single DataFrame
for gene in paper_genes: # loop through each gene highlighted in the paper and print its Z-scores and p-values across tissues
    g = all_df[all_df["gene_name"] == gene][["tissue", "zscore", "pvalue"]] # filter the combined DataFrame for rows corresponding to this gene and select the tissue, z-score, and p-value columns
    if len(g) > 0: # if the gene is found in any tissue, print its Z-scores and p-values across tissues
        print(f"\n {gene}:") # print the gene name
        print(g.to_string(index=False)) # print the tissue, z-score, and p-value for this gene across all tissues, without the index
    else:
        print(f"\n {gene}: not found in any tissue") # if the gene is not found in any tissue, print a message indicating that

        