#!/bin/bash
# This script filters genotype data for GWAS.
# Run with: bash /scripts/bash/prep-gwas-files.sh

DATA_PREFIX="spVL.shcs.rs"
FAM_W_SEX="spVL.shcs.rs_with_sex_SN.fam"

# Filter data to clean dataset
./plink2 \
--bfile data/$DATA_PREFIX \
--fam data/$FAM_W_SEX \
--keep output/gwas_individuals_to_keep.txt \
--king-cutoff 0.09375 \
--maf 0.01 \
--geno 0.05 \
--hwe 0.00005 \
--make-bed \
--memory 10000 \
--out output/${DATA_PREFIX}.filtered

# Summarize the QC-passed genome dataset
./plink2 \
--bfile output/${DATA_PREFIX}.filtered \
--freq \
--missing \
--hardy \
--memory 10000 \
--out output/${DATA_PREFIX}.filtered