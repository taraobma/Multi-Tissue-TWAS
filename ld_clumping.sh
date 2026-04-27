# This was run on the BU SCC

mkdir -p results/aim1
mkdir -p data/ref

# use plink2 to remove duplicates from the 1000G EUR reference panel
# then use plink1.9 to perform LD clumping on the allEC_ldsc.gz file

module load plink2/2.00a2.3

plink2 \
    --bfile /projectnb/bs859/data/1000G/plinkformat/1000G_EUR \
    --rm-dup force-first \
    --make-bed \
    --out data/ref/1000G_EUR_nodup \
    2>&1 | tee results/aim1/plink_nodup.log


module load plink/1.90b6.21

plink \
    --bfile data/ref/1000G_EUR_nodup \
    --clump allEC_ldsc.gz \
    --clump-snp-field SNP \
    --clump-p1 5e-8 \
    --clump-p2 1e-5 \
    --clump-r2 0.1 \
    --clump-kb 500 \
    --out results/aim1/clumped_loci \
    2>&1 | tee results/aim1/plink_clump.log