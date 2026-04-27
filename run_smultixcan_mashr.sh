# s-multixcan for mashr model

python MetaXcan/software/SMulTiXcan.py \
  --models_folder data/model/eqtl/mashr/ \
  --models_name_filter "mashr_(Adipose_Subcutaneous|Adipose_Visceral_Omentum|Ovary|Uterus|Vagina|Whole_Blood).db" \
  --models_name_pattern "mashr_(.*).db" \
  --snp_covariance data/model/eqtl/mashr/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz \
  --metaxcan_folder results/aim3/mashr/ \
  --metaxcan_filter "spredixcan_mashr_(.*).csv" \
  --metaxcan_file_name_parse_pattern "spredixcan_mashr_(.*).csv" \
  --gwas_file data/gwas/allEC_aim1.gz \
  --snp_column varID \
  --non_effect_allele_column OA \
  --effect_allele_column EA \
  --beta_column BETA \
  --pvalue_column P \
  --se_column SE \
  --keep_non_rsid \
  --cutoff_condition_number 30 \
  --trimmed_ensemble_id \
  --verbosity 2 \
  --throw \
  --output results/aim3/smultixcan_mashr.csv \
  2>&1 | tee results/aim3/smultixcan_mashr.log