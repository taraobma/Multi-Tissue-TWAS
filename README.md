# Multi-Tissue TWAS of Endometrial Cancer

A variant-to-gene pipeline for endometrial cancer that integrates GWAS summary statistics with GTEx v8 multi-tissue eQTL expression prediction models. The analysis replicates and extends Kho et al. (2021) using elastic net and mashr prediction models across six biologically relevant tissues.

---

## Overview

| Aim | Description | Tools |
|-----|-------------|-------|
| Aim 1 | GWAS QC, LD clumping, SNP heritability | Python, PLINK, LDSC |
| Aim 2 | Single-tissue and joint-tissue TWAS (elastic net) | S-PrediXcan, S-MultiXcan |
| Aim 3 | Model comparison: elastic net vs mashr | S-PrediXcan, Python |

**GWAS:** GCST006464 (O'Mara et al. 2018) — 12,906 cases, 108,979 controls, European ancestry  
**Tissues:** Subcutaneous adipose, visceral omentum adipose, whole blood, uterus, ovary, vagina  
**Models:** GTEx v8 elastic net and mashr (PredictDB)

---

## Repository Structure

Only source files are tracked. Data, results, logs, and model files are excluded via `.gitignore`.

```
859_final_proj/
├── aim1_gwas_qc.py            # Aim 1: GWAS QC pipeline (runs locally)
├── ld_clumping.sh             # Aim 1: PLINK LD clumping (run on BU SCC)
├── run_ldsc.sh                # Aim 1: LDSC heritability (run on BU SCC)
├── run_spredixcan_en.sh       # Aim 2: S-PrediXcan elastic net (runs locally)
├── run_smultixcan_en.sh       # Aim 2: S-MultiXcan elastic net (runs locally)
├── run_spredixcan_mashr.sh    # Aim 3: S-PrediXcan mashr (runs locally)
├── run_smultixcan_mashr.sh    # Aim 3: S-MultiXcan mashr (runs locally)
├── view_spredixcan_en.py      # View elastic net S-PrediXcan results
├── view_spredixcan_mashr.py   # View mashr S-PrediXcan results
├── view_smultixcan_en.py      # View elastic net S-MultiXcan results
├── view_smultixcan_mashr.py   # View mashr S-MultiXcan results
├── compare_results.py         # Compare elastic net vs mashr results
├── TWAS_pipeline.sh           # Master pipeline script (all steps in order)
├── final_report.ipynb         # Project write-up notebook
├── .gitignore
└── README.md
```

---

## How to Run

This pipeline runs in three stages. Stages 1 and 3 run locally; Stage 2 runs on BU SCC.

See `TWAS_pipeline.sh` for the full pipeline in order.

### Stage 1 — Local: GWAS QC

```bash
# Run QC pipeline (produces allEC_aim1.gz and allEC_ldsc.gz)
/opt/anaconda3/envs/imlabtools/bin/python aim1_gwas_qc.py

# Copy outputs to SCC
scp data/gwas/allEC_ldsc.gz nobma@scc1.bu.edu:/projectnb/bs859/students/nobma/final_proj/
scp run_ldsc.sh nobma@scc1.bu.edu:/projectnb/bs859/students/nobma/final_proj/
scp ld_clumping.sh nobma@scc1.bu.edu:/projectnb/bs859/students/nobma/final_proj/
```

### Stage 2 — BU SCC: LD Clumping and Heritability

```bash
module load plink2/2.00a2.3
module load plink/1.90b6.21
module load python2/2.7.16
module load ldsc/2020-05-05_github3d0c446

bash ld_clumping.sh
bash run_ldsc.sh

# Copy results back to local
scp "nobma@scc1.bu.edu:/projectnb/bs859/students/nobma/final_proj/results/aim1/*" results/aim1/
```

### Stage 3 — Local: TWAS and Model Comparison

```bash
# Aim 2: S-PrediXcan and S-MultiXcan (elastic net)
bash run_spredixcan_en.sh
bash run_smultixcan_en.sh

# View Aim 2 results
/opt/anaconda3/envs/imlabtools/bin/python view_spredixcan_en.py
/opt/anaconda3/envs/imlabtools/bin/python view_smultixcan_en.py

# Aim 3: S-PrediXcan and S-MultiXcan (mashr)
bash run_spredixcan_mashr.sh
bash run_smultixcan_mashr.sh

# View Aim 3 results and compare
/opt/anaconda3/envs/imlabtools/bin/python view_spredixcan_mashr.py
/opt/anaconda3/envs/imlabtools/bin/python view_smultixcan_mashr.py
/opt/anaconda3/envs/imlabtools/bin/python compare_results.py
```

---

## Key Results

### Aim 1
- **7,112,203** SNPs retained after QC
- **16** independent genome-wide significant loci (PLINK clumping, r² < 0.1, 500 kb)
- **h²_SNP = 0.025** (SE = 0.005), LDSC intercept = 1.066, λ_GC = 1.097
- **6/7** Kho et al. (2021) loci replicated at genome-wide significance

### Aim 2
- **9** Bonferroni-significant genes in S-MultiXcan (P < 3.79×10⁻⁶)
- All **7** Kho et al. (2021) genes replicated with concordant tissue assignments
- **2** additional genes: SRP14 and TSEN2 (present in paper at FDR < 0.01)
- Key tissue-specific findings: CYP19A1 (adipose), HEY2 (ovary), SKAP1 (whole blood)

### Aim 3
- Elastic net and mashr showed substantial differences in tissue assignments
- Whole blood was the most model-consistent tissue (CYP19A1, EEFSEC, EVI2A significant in both)
- HEY2 ovary signal present in elastic net but absent in mashr (no mashr model for HEY2 in ovary)
- Elastic net better preserved tissue-specific signals consistent with published colocalization results