# This script runs GWAS with both phenotypes, sex and top 5 PCs of genetic variation as co-variates.
# Run as: bash scripts/bash/run-gwas.sh

DATA_PREFIX="spVL.shcs.rs.filtered"
MEM=10000

# Get variants in LD to exclude from PCA calculation
./plink2 \
--bfile output/$DATA_PREFIX \
--indep-pairwise 50 5 0.2 \
--memory $MEM \
--out output/$DATA_PREFIX

# Get top 5 PCs for covariates based on SNPs in linkage equilibrium
./plink2 \
--bfile output/$DATA_PREFIX \
--exclude output/${DATA_PREFIX}.prune.out \
--pca \
--out output/$DATA_PREFIX \
--memory $MEM
awk '{print $1,$2,$3,$4,$5,$6,$7}' output/$DATA_PREFIX.eigenvec >  output/${DATA_PREFIX}.pc5.txt

# Run GWAS
./plink2 \
--bfile output/$DATA_PREFIX \
--pheno output/gwas_phenotypes.txt \
--glm sex \
--covar output/spVL.shcs.rs.filtered.pc5.txt \
--memory $MEM \
--out output/$DATA_PREFIX

# Extract SNP-only results, without covariates
head -1 output/${DATA_PREFIX}.trait.glm.linear > output/${DATA_PREFIX}.trait.glm.linear.nocovariates
grep -h 'ADD' output/${DATA_PREFIX}.trait.glm.linear >> output/${DATA_PREFIX}.trait.glm.linear.nocovariates

head -1 output/${DATA_PREFIX}.h.MWA.glm.linear > output/${DATA_PREFIX}.h.MWA.glm.linear.nocovariates
grep -h 'ADD' output/${DATA_PREFIX}.h.MWA.glm.linear >> output/${DATA_PREFIX}.h.MWA.glm.linear.nocovariates

# Paste results from both trait values together
paste "output/${DATA_PREFIX}.trait.glm.linear.nocovariates" "output/${DATA_PREFIX}.h.MWA.glm.linear.nocovariates" | awk '{print $1, $2, $3, $9, $12, $21, $24}' > output/gwas_results.txt
HEADER="CHROM POS ID BETA_standard P_standard BETA_corrected P_corrected"
sed -i.bak "1 s/^.*$/$HEADER/" output/gwas_results.txt

# Add p-value column to results
awk '
function abs(v) { return v < 0 ? -v : v }
BEGIN { OFS = " " }
NR == 1 {
$8 = "neg_log10_P_standard";
$9 = "neg_log10_P_corrected";
$10 = "neg_log10_P_standard_minus_corrected";
$11 = "abs_neg_log10_P_standard_minus_corrected" }
NR >= 2 {
$8 = -log($5)/log(10);
$9 = -log($7)/log(10);
$10 = $8 - $9
$11 = abs($10)} 1' < output/gwas_results.txt > output/gwas_results.pvals.txt

# Add beta (effect size) column to results
awk '
function abs(v) { return v < 0 ? -v : v }
BEGIN { OFS = " " }
NR == 1 {
$12 = "abs_BETA_standard_minus_corrected" }
NR >= 2 {
$12 = abs($4 - $6) } 1' < output/gwas_results.pvals.txt > output/gwas_results.pvals.effectsize.txt


