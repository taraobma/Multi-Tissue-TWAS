#!/bin/bash

module load python2/2.7.16
module load ldsc/2020-05-05_github3d0c446

mkdir -p results/aim1

# Munge sumstats
munge_sumstats.py \
    --sumstats allEC_ldsc.gz \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --signed-sumstats BETA,0 \
    --p P \
    --N-col N \
    --merge-alleles /projectnb/bs859/data/ldscore_files/w_hm3.snplist \
    --out results/aim1/allEC_munged

# SNP heritability using Pan-UKBB EUR LD scores
ldsc.py \
    --h2 results/aim1/allEC_munged.sumstats.gz \
    --ref-ld /projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscore/UKBB.EUR.rsid \
    --w-ld   /projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscore/UKBB.EUR.rsid \
    --out results/aim1/allEC_h2

cat results/aim1/allEC_h2.log