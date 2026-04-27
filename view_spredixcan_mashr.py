import pandas as pd
import os

TISSUES = [
    "Adipose_Subcutaneous",
    "Adipose_Visceral_Omentum",
    "Whole_Blood",
    "Uterus",
    "Ovary",
    "Vagina",
] # list of tissues analyzed

print("VIEW S-PREDIXCAN MASHR RESULTS\n")

paper_genes = ["GLDN", "CYP19A1", "HEY2", "SKAP1", "SNX11", "EVI2A", "EEFSEC"] # genes identified in Kho et al. 2021 as associated with endometrial cancer risk
all_results = [] # empty list to store results from all tissues

for tissue in TISSUES: # for loop through each tissue and read the corresponding S-PrediXcan results
    path = f"results/aim3/mashr/spredixcan_mashr_{tissue}.csv" # path to the S-PrediXcan results for this tissue
    t = pd.read_csv(path) # read the S-PrediXcan results into a DataFrame
    t["tissue"] = tissue # add a column for tissue name
    bonf = 0.05 / t["gene"].nunique() # calculate the Bonferroni threshold for this tissue based on the number of unique genes tested
    sig = t[t["pvalue"] < bonf] # filter for significant genes based on the Bonferroni threshold
    print(f"\n{tissue}: {t['gene'].nunique():,} genes, {len(sig)} significant") # print the number of genes tested and the number of significant genes for this tissue
    if len(sig) > 0: # if there are significant genes, print the top hits
        print(sig[["gene_name","zscore","pvalue"]].to_string(index=False)) # print the gene name, z-score, and p-value for the significant genes, without the index
    all_results.append(t) # add the results for this tissue to the list

print("\nPaper genes across mashr tissues:")
all_df = pd.concat(all_results) # combine results from all tissues into a single DataFrame
for gene in paper_genes: # for loop through each gene highlighted in the paper and print its Z-scores and p-values across tissues
    g = all_df[all_df["gene_name"] == gene][["tissue","zscore","pvalue"]] # filter the combined DataFrame for rows corresponding to this gene and select the tissue, z-score, and p-value columns
    if len(g) > 0: # if the gene is found in any tissue, print its Z-scores and p-values across tissues
        print(f"\n{gene}:") # print the gene name
        print(g.to_string(index=False)) # print the tissue, z-score, and p-value for this gene across all tissues, without the index
