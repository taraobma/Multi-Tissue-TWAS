#!/bin/bash
set -e

mkdir -p results/aim3/mashr

TISSUES=(
  "Adipose_Subcutaneous"
  "Adipose_Visceral_Omentum"
  "Whole_Blood"
  "Uterus"
  "Ovary"
  "Vagina"
)

for TISSUE in "${TISSUES[@]}"; do
  echo "Running S-PrediXcan (mashr) for ${TISSUE}..."

  /opt/anaconda3/envs/imlabtools/bin/python MetaXcan/software/SPrediXcan.py \
    --model_db_path data/model/eqtl/mashr/mashr_${TISSUE}.db \
    --covariance    data/model/eqtl/mashr/mashr_${TISSUE}.txt.gz \
    --gwas_folder   data/gwas/ \
    --gwas_file_pattern "allEC_aim1.gz" \
    --model_db_snp_key         varID \
    --snp_column               varID \
    --effect_allele_column     EA \
    --non_effect_allele_column OA \
    --beta_column              BETA \
    --se_column                SE \
    --pvalue_column            P \
    --gwas_N 121885 \
    --separator $'\t' \
    --keep_non_rsid \
    --remove_ens_version \
    --throw \
    --output_file results/aim3/mashr/spredixcan_mashr_${TISSUE}.csv

  echo "Done: ${TISSUE}"
done

echo "All mashr S-PrediXcan runs completed."