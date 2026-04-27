#S-MultiXcan Elastic Net Model
import pandas as pd

print("S-MULTIXCAN ELASTIC NET MODEL RESULTS\n")

df = pd.read_csv("results/aim2/smultixcan_allEC.csv", sep="\t") # read the S-MULTIXCAN results into a DataFrame
df.columns = df.columns.str.strip() # remove any leading/trailing whitespace from column names

bonferroni = 0.05 / len(df) # calculate the Bonferroni threshold for significance based on the total number of genes tested
sig = df[df["pvalue"] < bonferroni].sort_values("pvalue") # filter for significant genes and sort by p-value

print(f"Total genes tested : {len(df):,}") # print the total number of genes tested with comma as thousands separator
print(f"Bonferroni threshold : {bonferroni:.2e}") # print the Bonferroni threshold in scientific notation
print(f"Significant genes : {len(sig)}") # print the number of significant genes
print(f"\nSignificant genes (sorted by p-value):") # print a header for the significant genes table
print(sig[["gene", "gene_name", "pvalue", "n", "n_indep",
           "z_min", "z_max", "z_mean", "t_i_best"]].to_string(index=False)) # print the significant genes with the corresponding columns, without the index

# Replication of Kho et al. 2021 Table 1
print("\n\nReplication vs Kho et al. 2021 (Table 1)")
paper_genes = ["GLDN", "CYP19A1", "HEY2", "SKAP1", "SNX11", "EVI2A", "EEFSEC"] # genes highlighted in Kho et al. 2021
for gene in paper_genes: # loop through each gene highlighted in the paper and check if it is significant in our S-MULTIXCAN results
    row = df[df["gene_name"] == gene] # filter the S-MULTIXCAN results for the row corresponding to this gene
    if len(row) > 0: # if the gene is found in our results, check if it is significant and print the p-value and replication status
        p = row["pvalue"].values[0] # get the p-value for this gene
        status = "REPLICATED" if p < bonferroni else "not significant" # determine if the gene is replicated based on the Bonferroni threshold
        print(f"  {gene:<12}: p={p:.2e}  {status}") # print the gene name, p-value in scientific notation, and replication status
    else:
        print(f"  {gene:<12}: not found") # if the gene is not found in our results, print a message indicating that
        