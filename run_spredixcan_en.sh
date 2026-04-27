mkdir -p results/aim2

TISSUES=(
  "Adipose_Subcutaneous"
  "Adipose_Visceral_Omentum"
  "Whole_Blood"
  "Uterus"
  "Ovary"
  "Vagina"
)
# list of tissue to analyze

#for each tissue, run S-PrediXcan with the corresponding elastic net model and covariance file, using SNPID (rsID) for SNP matching
for TISSUE in "${TISSUES[@]}"; do
  echo "Running S-PrediXcan for ${TISSUE}..."
  python MetaXcan/software/SPrediXcan.py \
    --model_db_path data/model/eqtl/elastic_net/en_${TISSUE}.db \
    --covariance data/model/eqtl/elastic_net/en_${TISSUE}.txt.gz \
    --gwas_folder  data/gwas/ \
    --gwas_file_pattern "allEC_aim1.gz" \
    --snp_column SNPID \
    --effect_allele_column EA \
    --non_effect_allele_column OA \
    --beta_column BETA \
    --se_column SE \
    --pvalue_column P \
    --gwas_N 121885 \
    --remove_ens_version \
    --throw \
    --output_file results/aim2/spredixcan_allEC_${TISSUE}.csv
  echo "Done: ${TISSUE}"
done


# model_db_path specifies the path to the S-PrediXcan model database for this tissue
# covariance specifies the path to the SNP covariance file for this tissue
# gwas_folder specifies the folder containing the GWAS summary statistics files
# gwas_file_pattern specifies the pattern to match the GWAS summary statistics files in the data/gwas/ folder
# gwas_N is the sample size
# remove_ens_version removes the version number from Ensembl gene IDs (e.g. ENSG000001234.5 -> ENSG000001234) to match the gene IDs in the model db
# throw will raise an error if any issues are encountered (e.g. missing columns, mismatched SNP IDs) instead of silently skipping them - importantt for QC and troubleshooting

