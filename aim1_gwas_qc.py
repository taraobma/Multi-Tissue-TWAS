# Aim 1: GWAS Processing and QC for Endometrial Cancer

# This script must be run in two stages:
# STAGE 1 — Run on local first:
#   1. Load harmonized GWAS summary statistics from the GWAS Catalog
#   2. Apply QC filters (autosomes, biallelic SNPs, palindromic removal, MHC exclusion)
#   3. Construct GTEx v8 variant IDs for mashr model matching
#   6. Prepare input file for LDSC heritability estimation (data/gwas/allEC_ldsc.gz)
#
# Then copy outputs to BU SCC:
#   scp data/gwas/allEC_ldsc.gz nobma@scc1.bu.edu:.../final_proj/
#   scp run_ldsc.sh nobma@scc1.bu.edu:.../final_proj/
#
# STAGE 2 — Run on BU SCC:
#   - LD clumping using PLINK v1.90b6.21 with 1000 Genomes EUR reference panel
#     (run via run_plink_clump.sh; output: results/aim1/clumped_loci.clumped)
#   - SNP heritability using LDSC with Pan-UKBB EUR LD scores
#     (run via run_ldsc.sh; output: results/aim1/allEC_h2.log)
#
# Then copy results back to local:
#   scp "nobma@scc1.bu.edu:.../results/aim1/*" results/aim1/
#
# STAGE 3 — Run on local (macOS) again:
#   4. Load PLINK LD clumping results
#   5. Check replication of published EC GWAS loci (Kho et al. 2021)
#   6. Generate Manhattan and QQ plots for QC visualization


# Outputs after a succesful run of Aim 1:
# data/gwas/allEC_aim1.gz - QC-passed SNPs (for S-PrediXcan)
# data/gwas/allEC_ldsc.gz - LDSC-formatted sumstats (for heritability)
# results/aim1/top_loci.csv - independent genome-wide significant loci
# results/aim1/loci_replication.csv - Replication of Kho et al. 2021 genes
# results/aim1/manhattan_qq.png - Manhattan and QQ plots
# results/aim1/qc_summary.txt - QC filter counts and summary statistics


import os
import sqlite3
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2 as scipy_chi2

os.makedirs("results/aim1", exist_ok=True)


# Step 1: Load harmonized GWAS summary statistics
print("Step 1: Loading GWAS summary statistics")


df = pd.read_csv("data/gwas/GCST006464.h.tsv.gz", sep="\t", low_memory=False) # read the harmonized GWAS summary stats into a df
print(f" Loaded {len(df):,} rows") # print out the number of rows loaded with comma as thousands separator
print(f" Columns: {list(df.columns)}\n") # print out the column names to check what fields are available



# Step 2: QC filters
print("Step 2: Applying QC filters")


qc_log = {"raw": len(df)} # print out the initial number of rows before QC filtering


# Convert chromosome and position to numeric so we can filter on them
df["hm_chrom"] = pd.to_numeric(df["hm_chrom"], errors="coerce")
df["hm_pos"] = pd.to_numeric(df["hm_pos"], errors="coerce")


# Drop rows missing any of the essential fields
required_cols = ["hm_chrom", "hm_pos", "hm_effect_allele", "hm_other_allele",
                 "hm_beta", "standard_error", "p_value", "hm_rsid"] # fields needed for downstream analyses
df = df.dropna(subset=required_cols).copy() # drop rows with missing values in the defined column above
qc_log["missing_dropped"] = len(df) # log the number of rows remaining after dropping missing values
print(f" After dropping rows with missing fields: {len(df):,}") # print out the number of rows remaining after dropping missing values with comma as thousands separator


# Convert to integers now that NaN rows are gone
df["hm_chrom"] = df["hm_chrom"].astype(int) # convert chromosome column to integer type
df["hm_pos"] = df["hm_pos"].astype(int) # convert position column to integer type


# Keep autosomes only (chromosomes 1-22)
df = df[df["hm_chrom"].between(1, 22)].copy() # filter to keep only rows where chromosome is between 1 and 22 (inclusive)
qc_log["autosomes"] = len(df) # log the number of rows remaining after keeping autosomes only
print(f"  After keeping autosomes only (chr 1-22): {len(df):,}") # print out the number of rows remaining after keeping autosomes only


# Keep biallelic SNPs - both alleles must be a single A/T/C/G nucleotide
valid_bases = {"A", "T", "C", "G"}
df = df[df["hm_effect_allele"].isin(valid_bases) &
        df["hm_other_allele"].isin(valid_bases)].copy() # filter to keep rows with biallelic SNPs where both effect and other alleles are in the set of valid bases
qc_log["biallelic"] = len(df) # log the number of rows remaining after keeping biallelic SNPs only
print(f"  After keeping biallelic SNPs only: {len(df):,}") # print out the number of rows after filtering


# Remove palindromic SNPs (A/T and C/G pairs) - these are strand-ambiguous
# and cannot be reliably harmonized without allele frequency information
palindromic = (
    ((df["hm_effect_allele"] == "A") & (df["hm_other_allele"] == "T")) |
    ((df["hm_effect_allele"] == "T") & (df["hm_other_allele"] == "A")) |
    ((df["hm_effect_allele"] == "C") & (df["hm_other_allele"] == "G")) |
    ((df["hm_effect_allele"] == "G") & (df["hm_other_allele"] == "C"))
) # boolean mask to identify palindromic SNPs where effect and other alleles form A/T or C/G pairs
df = df[~palindromic].copy() # filter to keep only non-palindromic SNPs by negating the palindromic mask
qc_log["no_palindromic"] = len(df) # log the number of rows remaining after removing palindromic SNPs
print(f"  After removing palindromic SNPs: {len(df):,}") # print out the number of rows after filtering


# Exclude the MHC region (chr6:26-34 Mb) as this region has complex LD
# was also excluded by Kho et al. 2021 in their TWAS analysis


# define the chromosome and position boundaries for the MHC region
mhc_mask = (
    (df["hm_chrom"] == 6) &
    (df["hm_pos"] >= 26_000_000) &
    (df["hm_pos"] <= 34_000_000)
)
df = df[~mhc_mask].copy() # filter to keep only rows that are not in the MHC region by negating the mhc_mask
qc_log["no_mhc"] = len(df) # log the number of rows remaining after fitlering out the MHC region
print(f"  After removing MHC region (chr6:26-34 Mb): {len(df):,}") # print out the number of rows after filtering


# Remove duplicate positions
df = df.drop_duplicates(subset=["hm_chrom", "hm_pos",
                                 "hm_effect_allele", "hm_other_allele"]).copy() # drop dup rowsbased on chrom, pos, effect allele, and other allele to ensure unique variants
qc_log["no_duplicates"] = len(df) # log the rows remaining after removing duplicates
print(f"  After removing duplicate chr:pos:allele combinations: {len(df):,}") # print out the number of rows after filtering


# Sanity checks on SE and p-value
df = df[(df["standard_error"] > 0) &
        (df["p_value"] > 0) &
        (df["p_value"] <= 1)].copy() # filter to keep only rows with positive standard error and valid p-values between 0 and 1
qc_log["valid_stats"] = len(df) # log the number of rows remaining after sanity checks
print(f"  After SE and p-value sanity checks: {len(df):,}\n") # print out the number of rows after filtering



# Step 3: Rename columns and build variant IDs
print("Step 3: Renaming columns and constructing variant IDs")


df = df.rename(columns={
    "hm_rsid": "SNPID",
    "hm_effect_allele": "EA",
    "hm_other_allele": "OA",
    "hm_beta": "BETA",
    "standard_error": "SE",
    "p_value": "P",
    "hm_chrom": "CHR",
    "hm_pos": "POS",
    "hm_effect_allele_frequency": "EAF",
}) # rename columns for downstream analyses


# Build GTEx v8 PredictDB variant ID format: chr{N}_{pos}_{ref}_{alt}_b38
# The non-effect allele (OA) is the reference, effect allele (EA) is the alternate
# This format is required by mashr models; elastic net models use rsIDs
df["varID"] = ("chr" + df["CHR"].astype(str) + "_" +
               df["POS"].astype(str) + "_" +
               df["OA"] + "_" + df["EA"] + "_b38") # chr{N}_{pos}_{ref}_{alt}_b38


# Build the flipped version in case the model stores alleles in the other order
df["varID_flip"] = ("chr" + df["CHR"].astype(str) + "_" +
                    df["POS"].astype(str) + "_" +
                    df["EA"] + "_" + df["OA"] + "_b38")


# Select output columns
out_cols = ["varID", "varID_flip", "SNPID", "CHR", "POS", "EA", "OA", "BETA", "SE", "P"] # columns to include in the output file
if "EAF" in df.columns: #if effect allele freq is available
    out_cols.append("EAF") # add effect allele frequency to output if available
df_out = df[out_cols].copy() # create a new DataFrame with only the selected output columns


print(f" Final QC-passed SNPs: {len(df_out):,}") # print out a message that includes the number of QC-passed SNPs
print(f" Sample output:") # show the sample output from the code above
print(df_out.head(3).to_string(index=False)) # print out the top 3 rows of the final QC-passed DataFrame with all columns, no index


# Compute lambda GC from the full QC-passed dataset
# Formula: lambda_GC = median(chi2) / 0.4549
# where chi2 is derived from p-values via the chi-squared inverse survival function (df=1)
# 0.4549 is the theoretical median of chi2(1) under the null
all_pvals = df_out["P"].values.clip(1e-300) # clip to avoid log(0)
chi2_stats = scipy_chi2.isf(all_pvals, df=1) # convert p-values to chi-squared statistics
lambda_gc  = np.median(chi2_stats) / 0.4549  # divide by expected median under null

# Step 4: Check SNP match rate against GTEx v8 elastic net model on the Uterus tissue 
print("\nStep 4: SNP match rate check against GTEx v8 elastic net (Uterus)")


db_path = "data/model/eqtl/elastic_net/en_Uterus.db" # path to the tissue from GTEx v8 elastic net models
# use sql to read the model SNPs from the weights table and check how many of our QC-passed SNPs match the model SNPs based on both rsID and varID formats
if os.path.exists(db_path): 
    conn = sqlite3.connect(db_path)
    model_snps = set(pd.read_sql("SELECT rsid FROM weights", conn)["rsid"])
    conn.close() 


    n = len(df_out) # total number of QC-passed SNPs
    # rsID for 
    print(f" rsID match: {df_out['SNPID'].isin(model_snps).sum():,} / {n:,}") # count the matched SNPs based on rsID and print the count
    print(f"({100 * df_out['SNPID'].isin(model_snps).mean():.1f}%)") # print out the percentage of SNPs that match the model based on rsID
    print(f" varID match: {df_out['varID'].isin(model_snps).sum():,} / {n:,}") # count the matched SNPs based on varID and print the count
else: # if the model file is not found
    print(f" Model file not found at {db_path}") # print out a message



# Step 5: Save QC-passed file for S-PrediXcan and LDSC
print("\nStep 5: Saving QC-passed GWAS file")


out_path = "data/gwas/allEC_aim1.gz" # save the QC-passed to a file called allEC_aim1, with gzip compressed
df_out.to_csv(out_path, sep="\t", index=False, compression="gzip") 
print(f" Saved: {out_path}") # once saved, print out a message that the file is save followed by the path of the file
print(f" Lambda GC: {lambda_gc:.3f}") # print the lambda GC computed from the full QC-passed dataset


# Step 6: Load PLINK LD clumping results
# LD clumping was run on BU SCC using PLINK v1.90b6.21 - run it prior to this step using ld_clumping.py
# Reference panel: 1000 Genomes EUR (N=503, GRCh38, duplicates removed with PLINK2)
# Parameters: P1=5e-8, P2=1e-5, r2<0.1, window=500kb


print("\nStep 6: Loading PLINK LD clumping results")
clump_path = "results/aim1/clumped_loci.clumped" # path to PLINK clumping output file
# if the clumping results file exists, read it in and extract the independent loci
# otherwise print a message that the file is not found
if os.path.exists(clump_path):
    loci_df = pd.read_csv(clump_path, sep=r"\s+")
    loci_df = loci_df[["CHR", "BP", "SNP", "P"]].copy()
    loci_df.columns = ["CHR", "POS", "SNPID", "P"]
    loci_df = loci_df.sort_values(["CHR", "POS"]).reset_index(drop=True)
    loci_df["locus_id"] = range(1, len(loci_df) + 1)


# print out the following params used in the PLINK command for LD clumping, along with the number of indep loci identified
    print(f" Independent loci identified: {len(loci_df)}")
    print(f" (P < 5e-8, r2 < 0.1, window = 500 kb, reference = 1000G EUR)\n") 
    
    # header for LD clumping results
    print(loci_df[["locus_id", "CHR", "POS", "SNPID", "P"]].to_string(index=False))


    loci_df.to_csv("results/aim1/top_loci.csv", index=False) # save the result to a CSV file
    print("\n Saved: results/aim1/top_loci.csv") # once saved, print out a message


else: # is clumping results not found print out the following message and create an empty df with the identified columns
    print(f" {clump_path} not found.") 
    print(" Run PLINK clumping on SCC first — use ld_clumping.py.")
    loci_df = pd.DataFrame(columns=["CHR", "POS", "SNPID", "P", "locus_id"])



# Step 7: Replication check against Kho et al. 2021 Table 1 genes
# These are the 7 genes identified by S-MultiXcan in the original paper.
# We check whether our GWAS data shows a signal at each underlying locus.


print("\nStep 7: Replication of Kho et al. 2021 S-MultiXcan genes (Table 1)")


# UCSC Genome Browser on Human (GRCh38/hg38) band boundaries
# Dict of band coordinate
band_coords = {
    "3q21.3": (3, 126_100_000, 129_500_000),
    "6q22.31": (6, 118_100_000, 125_800_000),
    "15q21.2": (15, 49_200_000, 52_600_000),
    "17q11.2": (17, 27_400_000, 33_500_000),
    "17q21.32": (17, 46_800_000, 49_300_000),
}


# DF of the article's Table 1
paper_loci = pd.DataFrame([
    {"locus": "3q21.3", "gene": "EEFSEC", "paper_p": 1.10e-6},
    {"locus": "6q22.31", "gene": "HEY2", "paper_p": 9.94e-9},
    {"locus": "15q21.2", "gene": "GLDN", "paper_p": 1.34e-12},
    {"locus": "15q21.2", "gene": "CYP19A1", "paper_p": 9.52e-12},
    {"locus": "17q11.2", "gene": "EVI2A", "paper_p": 1.50e-6},
    {"locus": "17q21.32", "gene": "SKAP1", "paper_p": 7.27e-9},
    {"locus": "17q21.32", "gene": "SNX11", "paper_p": 5.41e-7},
])


# parse chromosome number from the locus string (e.g., "15q21.2" to 15)
# Add band coordinates from the lookup
paper_loci["chr"] = paper_loci["locus"].map(lambda x: band_coords[x][0])
paper_loci["band_start"] = paper_loci["locus"].map(lambda x: band_coords[x][1])
paper_loci["band_end"] = paper_loci["locus"].map(lambda x: band_coords[x][2])


# For each paper locus, find the best clumped locus on the same chromosome
rep_results = [] # empty list to store the rep results


for _, row in paper_loci.iterrows(): # loop over each paper gene
    in_band = df_out[(df_out["CHR"] == row["chr"]) &
                     (df_out["POS"] >= row["band_start"]) &
                     (df_out["POS"] <= row["band_end"])] # subset SNPs to only those within the cytoband boundaries
    if len(in_band) == 0:
        rep_results.append({**row, "proj_snp": None, "proj_p": None,
                            "gw_sig": False, "suggestive": False}) # append a no signal record so the gene still appears in the output table
        continue
    best = in_band.loc[in_band["P"].idxmin()] # get the smallest p-val within the cytoband


    # Record the paper gene's metadata alongside our best matching locus and its significance flags
    rep_results.append({
        **row, # unpack all paper_loci columns for this gene
        "proj_snp": best["SNPID"], # rsID of our most significant SNP within that cytoband
        "proj_p": best["P"], # p-value
        "gw_sig": best["P"] < 5e-8, # True if genome-wide significant (std threshold)
        "suggestive": best["P"] < 1e-5, # True if at least suggestive (relaxed threshold)
    })


rep_df = pd.DataFrame(rep_results) # convert list of dicts to a tidy DataFrame
rep_df.to_csv("results/aim1/loci_replication.csv", index=False) # save to csv


n_gw = rep_df["gw_sig"].sum() # count genes replicated at genome-wide significance
n_sug = rep_df["suggestive"].sum() # count genes with at least a suggestive signal


# Print summary counts and the full replication table
print(f" Genes checked (Kho et al. 2021 Table 1): {len(rep_df)}")
print(f" Replicated at genome-wide significance (P < 5e-8): {n_gw}/{len(rep_df)}")
print(f" Replicated at suggestive significance (P < 1e-5): {n_sug}/{len(rep_df)}\n")
print(rep_df[["locus", "gene", "paper_p", "proj_snp",
              "proj_p", "gw_sig", "suggestive"]].to_string(index=False))
# Select only the display columns and print without the pandas row index



# Step 8: Manhattan and QQ plots
print("\nStep 8: Generating Manhattan and QQ plots")


sig_snps = df_out[df_out["P"] < 1e-5] # keep all SNPs below P < 1e-5


# Randomly sample up to 200,000 non-significant SNPs to reduce plotting time and memory with fixed seed for reproducibility
other_snps = df_out[df_out["P"] >= 1e-5].sample(n=min(200_000, len(df_out[df_out["P"] >= 1e-5])), random_state=42)


# Combine the two subsets and sort by genomic position for correct left-to-right ordering
plot_df = pd.concat([sig_snps, other_snps]).sort_values(["CHR", "POS"]).reset_index(drop=True)


# Compute –log10(p) for plotting; clip at 1e-300 to avoid log(0) = -inf
plot_df["neglog10p"] = -np.log10(plot_df["P"].clip(lower=1e-300))


# Build cumulative chromosome x-axis offsets
chrom_offset = {} # dict mapping chromosome number to cumulative base-pair offset
offset = 0 # running total of base pairs placed so far
for chrom in range(1, 23): # iterate chromosomes 1 through 22
    chrom_offset[chrom] = offset # store current offset before adding this chromosome's width
    chrom_snps = plot_df[plot_df["CHR"] == chrom] # all SNPs on this chromosome
    if len(chrom_snps) > 0:
        offset += chrom_snps["POS"].max() # advance offset by the max position on this chromosome


# Absolute genomic position = within-chromosome position + chromosome's cumulative offset
plot_df["abs_pos"] = plot_df.apply(lambda r: r["POS"] + chrom_offset[r["CHR"]], axis=1)


# Create a 2-panel figure: Manhattan on top, QQ on bottom
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 10))
fig.suptitle("Endometrial Cancer GWAS QC — GCST006464", fontsize=14, fontweight="bold")


# Alternating blue shades to visually separate adjacent chromosomes
colors = ["#3B6CB5", "#8EB4E3"]
for i, chrom in enumerate(range(1, 23)): # loop through each chromosome
    c_df = plot_df[plot_df["CHR"] == chrom] # SNPs belonging to this chromosome
    if len(c_df) > 0:
        ax1.scatter(c_df["abs_pos"], c_df["neglog10p"],
                    s=2, # tiny dot size avoids overplotting with millions of SNPs
                    c=colors[i % 2], # alternate between two blues for even/odd chromosomes
                    alpha=0.5, # semi-transparency helps show overlapping density
                    rasterized=True) # render as bitmap to keep SVG/PDF file size manageable


# Draw horizontal reference lines at the two standard thresholds
ax1.axhline(-np.log10(5e-8), color="red", linestyle="--", linewidth=1, label="P = 5×10⁻⁸") # GW significance
ax1.axhline(-np.log10(1e-5), color="orange", linestyle=":", linewidth=1, label="P = 1×10⁻⁵") # suggestive threshold


# Annotate each independent clumped locus with its rsID above its point
for _, row in loci_df.iterrows():
    xpos = row["POS"] + chrom_offset[row["CHR"]] # absolute x-position on the plot
    ypos = -np.log10(row["P"]) # y-position on the –log10 scale
    ax1.annotate(row["SNPID"], xy=(xpos, ypos), xytext=(0, 5),
                 textcoords="offset points", # offset 5 pts upward so label doesn't overlap the dot
                 fontsize=6, ha="center", color="darkred")


# Set x-axis tick marks at the midpoint of each chromosome band
xticks, xlabels = [], []
for chrom in range(1, 23):
    c_df = plot_df[plot_df["CHR"] == chrom]
    if len(c_df) > 0:
        xticks.append(c_df["abs_pos"].median()) # median abs_pos = visual center of chromosome band
        xlabels.append(str(chrom)) # label is just the chromosome number


ax1.set_xticks(xticks)
ax1.set_xticklabels(xlabels, fontsize=7) # small font to prevent overlapping chromosome labels
ax1.set_xlabel("Chromosome", fontsize=10)
ax1.set_ylabel("−log₁₀(P)", fontsize=10)
ax1.set_title("Manhattan plot", fontsize=11)
ax1.legend(fontsize=8)

# QQ plot 
# QQ plot — use a uniform random sample across ALL p-values
# to avoid the subsampling artifact where keeping all significant SNPs
# distorts the expected vs observed ratio
qq_sample = df_out.sample(n=min(500_000, len(df_out)), random_state=42)

p_sorted = np.sort(qq_sample["P"].values)   # ascending: index 0 = smallest p-value
n_qq     = len(p_sorted)

# Rank-based expected quantiles: rank i/(n+1) under uniform null
# Result is descending (-log10 of small ranks = large values)
# observed is also descending because p_sorted[0] is smallest p → largest -log10
# So both arrays are paired correctly at every index with no flipping needed
expected = -np.log10(np.arange(1, n_qq + 1) / (n_qq + 1))  # descending
observed = -np.log10(p_sorted.clip(1e-300))                  # descending (-log10 of ascending p)

step = max(1, n_qq // 50_000)  # thin to ~50k points so scatter renders quickly
ax2.scatter(expected[::step], observed[::step],
            s=2, alpha=0.5, color="#3B6CB5", rasterized=True)
ax2.plot([0, expected[0]], [0, expected[0]],
         "r--", linewidth=1, label="Expected (null)")  # diagonal null line; departure = inflation
ax2.set_ylim(0, min(observed.max(), 7))   # cap at wherever null signal ends
ax2.set_xlim(0, expected[0] * 1.05)
ax2.set_xlabel("Expected −log₁₀(P)", fontsize=10)
ax2.set_ylabel("Observed −log₁₀(P)", fontsize=10)
ax2.set_title(f"QQ plot (λ_GC = {lambda_gc:.3f})", fontsize=11)
ax2.legend(fontsize=8)

plt.tight_layout()
plt.savefig("results/aim1/manhattan_qq.png", dpi=150, bbox_inches="tight")
plt.close()
print(" Saved: results/aim1/manhattan_qq.png")

# Step 9: Prepare LDSC input file
# LDSC munge_sumstats.py expects columns: SNP, A1, A2, BETA, SE, P, N
# File is used on BU SCC to estimate SNP heritability


print("\nStep 9: Preparing LDSC-ready sumstats file")


df_ldsc = df_out[["SNPID", "EA", "OA", "BETA", "SE", "P"]].copy() # columns required by LDSC's munge_sumstats.py
df_ldsc.columns = ["SNP", "A1", "A2", "BETA", "SE", "P"] # rename to the exact column names expected by munge_sumstats.py
df_ldsc["N"] = 121_885 # sample size column for LDSC


# Write as a gzip-compressed TSV; gzip compression to save space for large GWAS files
df_ldsc.to_csv("data/gwas/allEC_ldsc.gz", sep="\t", index=False, compression="gzip")


# Once done print out a message and reminder to run LDSC on the SCC
print(f" Saved: data/gwas/allEC_ldsc.gz ({len(df_ldsc):,} SNPs)")
print(" Copy to SCC and run run_ldsc.sh to estimate SNP heritability")


# Step 10: Write QC summary report
print("\nStep 10: Writing QC summary")


with open("results/aim1/qc_summary.txt", "w") as f: # open file for writing
    f.write("GWAS QC Summary\n")
    f.write("Study: GCST006464 (O'Mara et al. 2018)\n") # document the source study
    f.write("N = 121,885 (12,906 cases and 108,979 controls) - all European ancestry\n\n") # record cohort size


    f.write("QC Filters\n")
    # List of (qc_log key, easy readable label) pairs in pipeline order
    steps = [
        ("raw", "Raw SNPs"),
        ("missing_dropped", "After dropping missing fields"),
        ("autosomes", "After keeping autosomes (chr 1-22)"),
        ("biallelic", "After keeping biallelic SNPs"),
        ("no_palindromic", "After removing palindromic SNPs"),
        ("no_mhc", "After removing MHC region (chr6:26-34 Mb)"),
        ("no_duplicates", "After removing duplicates"),
        ("valid_stats", "After SE/p-value sanity checks"),
    ]

    prev = None # track the previous step's SNP count to compute how many were removed
    for key, label in steps:
        n = qc_log.get(key, "N/A") # look up SNP count for this step; default "N/A" if step wasn't logged

        # Compute how many SNPs were dropped relative to the previous step
        removed = f" (-{prev - n:,})" if prev and isinstance(n, int) else ""
        f.write(f" {label:<45}: {n:>10,}{removed}\n") # left-align label, right-align count, append removed count
        prev = n if isinstance(n, int) else prev # update prev only if we got a real integer count


    f.write(f"\nFinal SNPs for TWAS: {len(df_out):,}\n") # SNP count passed on to downstream analysis
    f.write("-" * 40 + "\n")
    f.write(f"\nLD Clumping (PLINK v1.90b6.21 on BU SCC)\n")
    f.write(f" Reference panel : 1000 Genomes EUR (N=503, GRCh38)\n") # LD reference used for clumping
    f.write(f" Significance (P1) : P < 5e-8\n") # only SNPs below this are index SNPs
    f.write(f" Secondary (P2) : P < 1e-5\n") # SNPs above P2 are excluded from clump windows
    f.write(f" LD threshold (r2) : < 0.1\n") # SNPs in LD with index SNP below r² 0.1 are pruned
    f.write(f" Window : 500 kb\n") # physical window around each index SNP
    f.write(f" Independent loci : {len(loci_df)}\n") # number of loci surviving clumping

    f.write(f"\nReplication vs Kho et al. 2021 Table 1\n")
    f.write("-" * 40 + "\n")
    f.write(f" Genes checked : {len(rep_df)}\n") # total paper genes tested
    f.write(f" Replicated (P<5e-8) : {n_gw}/{len(rep_df)}\n") # genome-wide replications
    f.write(f" Replicated (P<1e-5) : {n_sug}/{len(rep_df)}\n") # suggestive replications

    f.write(f"\nGenomic Inflation\n")
    f.write("-" * 40 + "\n")
    f.write(f" Lambda GC : {lambda_gc:.3f}\n") # λ_GC

print(" Saved: results/aim1/qc_summary.txt") # once saved, print out a message