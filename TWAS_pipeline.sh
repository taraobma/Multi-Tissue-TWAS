# Shell Script used for this project

## Download the Data
wget "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006464/harmonised/30093612-GCST006464-EFO_1001512.h.tsv.gz" \
  -O data/gwas/GCST006464.h.tsv.gz

## md5sum file
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006464/harmonised/md5sum.txt

## Show the study metadata
wget "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006464/harmonised/30093612-GCST006464-EFO_1001512.h.tsv.gz-meta.yaml" -O - | cat

## Download GTEx v8 MASHR S-PrediXcan models
# Download MASHR eQTL models
wget https://zenodo.org/record/3518299/files/mashr_eqtl.tar?download=1 -O mashr_eqtl.tar
tar -xvf mashr_eqtl.tar

# S-MultiXcan LD covariance file (MASHR)
wget https://zenodo.org/record/3518299/files/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz?download=1 \
  -O gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz

## Download MetaXcan software
git clone https://github.com/hakyimlab/MetaXcan

# Create and activate the pre-built conda environment
conda env create -f MetaXcan/software/conda_env.yaml
conda activate imlabtools

# Install additional dependencies not in imlabtools by default
pip install matplotlib # for Manhattan and QQ plots


# Test both tools work
/opt/anaconda3/envs/imlabtools/bin/python MetaXcan/software/SPrediXcan.py --help
/opt/anaconda3/envs/imlabtools/bin/python MetaXcan/software/SMulTiXcan.py --help

## Download elastic net models
wget "https://zenodo.org/record/3519321/files/elastic_net_eqtl.tar?download=1" \
  -O data/model/elastic_net_eqtl.tar

# S-MultiXcan SNP covariance (elastic net)
wget "https://zenodo.org/record/3519321/files/gtex_v8_expression_elastic_net_snp_smultixcan_covariance.txt.gz?download=1" \
  -O data/gtex_v8_expression_elastic_net_snp_smultixcan_covariance.txt.gz

# Extract elastic net models
mkdir -p data/model/eqtl/elastic_net
tar -xf data/model/elastic_net_eqtl.tar -C data/model/eqtl/elastic_net

# Verify 6 tissues are present
ls data/model/eqtl/elastic_net/ | grep -E "Uterus|Ovary|Vagina|Adipose|Whole_Blood"

# S-MultiXcan SNP covariance (mashr)
wget "https://zenodo.org/records/3518299/files/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz" \
  -O data/model/eqtl/mashr/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz

# Get the UCSC genome band coordinate for each gene and chromosome location from Table 1 - for step 7 
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz \
  -O data/cytoband/cytoBand.txt.gzgunzip cytoBand.txt.gz

gzcat data/cytoband/cytoBand.txt.gz | awk '
  ($1=="chr3"  && $4=="q21.3")  ||
  ($1=="chr6"  && $4=="q22.31") ||
  ($1=="chr15" && $4=="q21.2")  ||
  ($1=="chr17" && $4=="q11.2")  ||
  ($1=="chr17" && $4=="q21.32") ' > data/cytoband/chr_band.txt

# == AIM 1 – GWAS QC ==========
# Filters: autosomes only, biallelic SNPs, palindromic SNP removal, 
# MHC exclusion (chr6:26-34Mb), deduplication, SE/p sanity checks
# Constructs varID (chr_pos_OA_EA_b38) and varID_flip columns for mashr

# Outputs:  data/gwas/allEC_aim1.gz — QC-passed SNPs for Aim 2 (SNPID=rsID)
#           data/gwas/allEC_ldsc.gz — LDSC-ready munged sumstats
#           results/aim1/top_loci.csv — independent GW-significant loci
#           results/aim1/loci_replication.csv — replication vs Kho et al. 2021
#           results/aim1/manhattan_qq.png — Manhattan and QQ plots
#           results/aim1/qc_summary.txt — QC filter counts and lambda GC

python aim1_gwas_qc.py 2>&1 | tee results/aim1/aim1_run.log

# View the QC-passed GWAS file (first 6 rows)
gzcat data/gwas/allEC_aim1.gz | head -6

# View QC summary
cat results/aim1/qc_summary.txt

#LDSC SNP heritability (run on BU SCC)
# Copy files to SCC
scp data/gwas/allEC_ldsc.gz nobma@scc1.bu.edu:/projectnb/bs859/students/nobma/final_proj/
scp run_ldsc.sh nobma@scc1.bu.edu:/projectnb/bs859/students/nobma/final_proj/

# On BU SCC run:
./run_ldsc.sh 2>&1 | tee results/aim1/ldsc_run.log

# Copy results back to local
scp "nobma@scc1.bu.edu:/projectnb/bs859/students/nobma/final_proj/results/aim1/*" results/aim1/


# == AIM 2 – S-PrediXcan (elastic net, 6 tissues) ==========
# Outputs:  results/aim2/spredixcan_allEC_{TISSUE}.csv for each tissue

# make a directory for aim 2's results 
mkdir -p results/aim2

# run S-PrediXcan Elastic Net Model
./run_spredixcan_en.sh 2>&1 | tee results/aim2/spredixcan_en_run.log

# view the s-predixcan elastic net results
python view_spredixcan.py

# AIM 2 – S-MultiXcan (joint 6-tissue, elastic net)
# Output:   results/aim2/smultixcan_allEC.csv, results/aim2/smultixcan_allEC.log
./run_smultixcan.sh

# check the output files
head -5 results/aim2/smultixcan_allEC.csv
wc -l results/aim2/smultixcan_allEC.csv

# view the the s-multixcan results of the elastic net 
python view_smultixcan.py


# AIM 3 – MASHR S-PrediXcan (6 tissues)
# Note: mashr models use varID as key (not rsID like elastic net)
# Note: use --keep_non_rsid to retain chr_pos format SNPs
# Output:   results/aim3/mashr/spredixcan_mashr_{TISSUE}.csv for each tissue

# make a directory for aim 3's results
mkdir -p results/aim3/mashr

# run S-PrediXcan Mashr Model
./run_spredixcan_mashr.sh 2>&1 | tee results/aim3/mashr/spredixcan_mashr_run.log


# AIM 3 – Compare elastic net vs mashr results
# Produces 4 comparison tables:
#   Table 1: S-MultiXcan significant genes
#   Table 2: Per-tissue Z-scores EN vs Mashr for paper genes
#   Table 3: Replication vs Kho et al. 2021
#   Table 4: Significant gene overlap EN vs Mashr per tissue

python compare_results.py
