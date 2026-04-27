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

See `TWAS_pipeline.sh` for the full pipeline with all commands in order, including data download, environment setup, and all three aims.

> **Note:**  Stage 2 (LD clumping and LDSC heritability) requires BU SCC access. All other steps run locally on macOS.
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