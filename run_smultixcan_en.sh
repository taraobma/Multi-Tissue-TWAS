# s-multixcan for elastic net model

python MetaXcan/software/SMulTiXcan.py \
  --models_folder data/model/eqtl/elastic_net/ \
  --models_name_filter "en_(Adipose_Subcutaneous|Adipose_Visceral_Omentum|Ovary|Uterus|Vagina|Whole_Blood).db" \
  --models_name_pattern "en_(.*).db" \
  --snp_covariance data/model/gtex_v8_expression_elastic_net_snp_smultixcan_covariance.txt.gz \
  --metaxcan_folder results/aim2/ \
  --metaxcan_filter "spredixcan_allEC_(.*).csv" \
  --metaxcan_file_name_parse_pattern "spredixcan_allEC_(.*).csv" \
  --gwas_file data/gwas/allEC_aim1.gz \
  --snp_column SNPID \
  --non_effect_allele_column OA \
  --effect_allele_column EA \
  --beta_column BETA \
  --pvalue_column P \
  --se_column SE \
  --cutoff_condition_number 30 \
  --trimmed_ensemble_id \
  --verbosity 2 \
  --throw \
  --output results/aim2/smultixcan_allEC.csv \
  2>&1 | tee results/aim2/smultixcan_allEC.log

  # trimmed_ensemble_id is to strip everything after the dot in the ensemble id for matching gene IDs across file\
  # GTEx v8 elastic net - uses unversioned IDs
  # ENSG00000139618.15 to ENSG00000139618