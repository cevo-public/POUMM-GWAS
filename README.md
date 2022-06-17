# Infectious disease GWAS project

This project applies the Phylogenetic Ornstein-Uhlenback Mixed Model (POUMM) to estimate and remove pathogen effects from infectious disease trait data prior to Genome-Wide Association Study (GWAS).

## Run simuations
Simulate data under the POUMM using a randomly generated, HIV-like phylogeny and specified parameter values/ranges. Fit the POUMM to estimate back the parameters and check accuracy against the true (simulated) values.

### Requirements
* Docker
* Specified parameter values/ranges in [config-simulation.yaml](config-simulation-full.yaml)

### Instructions

* Build the containers specified in the Dockerfiles:
  ```
  docker build -t simulate-poumm-accuracy -f Dockerfile-simulate-poumm-accuracy .
  docker build -t simulate-poumm-gwas -f Dockerfile-simulate-poumm-gwas .
  ```
* To run locally:
  ```
  mkdir -p output
  docker run \
  -v /Users/nadeaus/Repos/poumm-gwas/output:/output \
  -v /Users/nadeaus/Repos/poumm-gwas/config-simulation.yaml:/config-simulation.yaml \
  simulate-poumm-accuracy
  docker run \
  -v /Users/nadeaus/Repos/poumm-gwas/output:/output \
  -v /Users/nadeaus/Repos/poumm-gwas/config-simulation.yaml:/config-simulation.yaml \
  simulate-poumm-gwas
  ```
* To run on Euler:
  * Push the images to the ETHZ repository:
  ```
  docker build -t registry.ethz.ch/nadeaus/poumm-gwas/simulate-poumm-accuracy -f Dockerfile-simulate-poumm-accuracy .
  docker build -t registry.ethz.ch/nadeaus/poumm-gwas/simulate-poumm-gwas -f Dockerfile-simulate-poumm-gwas .
  docker push registry.ethz.ch/nadeaus/poumm-gwas/simulate-poumm-accuracy
  docker push registry.ethz.ch/nadeaus/poumm-gwas/simulate-poumm-gwas
  ```
  * Log in to Euler, clone this repository, navigate to top-level directory
  * Pull the images from gitlab and convert them to singularity images in interactive jobs:
  ```
  module load eth_proxy
  bsub -I -W 0:05 -o singularity_build_%J.log "singularity build --docker-login simulate-poumm-accuracy.sif docker://registry.ethz.ch/nadeaus/poumm-gwas/simulate-poumm-accuracy:latest"
  # Enter username and then password for ETH gitlab (no prompt will come, just enter them and wait until job ends)
  bsub -I -W 0:05 -o singularity_build_%J.log "singularity build --docker-login simulate-poumm-gwas.sif docker://registry.ethz.ch/nadeaus/poumm-gwas/simulate-poumm-gwas:latest"
  ```
  * Run the singularity images in Euler jobs:
  ```
  mkdir -p output
  bsub -N -o simulate-poumm-accuracy_%J.log \
  "singularity run --bind config-simulation.yaml:/config-simulation.yaml --bind output:/output simulate-poumm-accuracy.sif"
  bsub -N -o simulate-poumm-gwas_%J.log \
  "singularity run --bind config-simulation.yaml:/config-simulation.yaml --bind output:/output simulate-poumm-gwas.sif"
  ```
  
## Apply method to GWAS for host genetic determinants of HIV spVL 

### Data
Data cannot be published due to privacy protections. See the [Swiss HIV Cohort Study (SHCS)](http://www.shcs.ch/) for more information.

* Host genetic data provided by the SHCS.

  * genotypes (spVL.shcs.rs.bed, pVL.shcs.rs.bim, and pVL.shcs.rs.fam)

  * HLA genotypes imputed with SNP2HLA (spVL.shcs.rs.hla.bed, spVL.shcs.rs.hla.bim, spVL.shcs.rs.hla.fam)

* SHCS cohort scores along top principal components from genotype matrix of merged SHCS and HapMap data. Provided by Christian Thorball.

* Viral load measurements and other clinical data provided by the SHCS.

* Viral genetic data provided by the SHCS.

  * pol gene sequences (newfasta2019-10-25.fas)

  * Viral subtypes

### Calculate spVL and prepare sequence data, metadata
* ``` Rscript scripts/R/calculate_spvl.R``` produces spVL values calculated a few different ways. They are compared in `output/spvl_calculation_comparison.png`. I use the mean of viral load measurements taken before treatment start for further analysis.
* ```Rscript scripts/R/filter_sequences.R``` attaches spVL and subtype data to the pol sequence data. I filter to only subtype B sequences with spVL values (focal) and A sequences (outgroup) of at least 750 non-gap, non-N characters. The sequence header format is `<patient id>_<collection date %Y-%m-%d>_<'outgroup' or 'focal'>_<spVL value>`.

### Build pathogen phylogeny (IQ-TREE)
* Build container, mount volume with prepared sequence data, run container.
* This is done locally because I don't have sudo permissions to run Docker on the server where the data lives.
* The script generates an alignment, trims characters after position 1505, and constructs and approximate maximum-likelihood tree.
``` 
docker build -t build-tree -f Dockerfile-build-tree .
# Connect to smb://d.ethz.ch/groups/bsse/stadler/
docker run \
--volume=/Volumes/stadler/SHCSData/data/newfasta2019-10-25.fas:/sequences.fasta:ro \
--volume=`pwd`/output:/output build-tree
```

<!-- ### Build pathogen phylogeny (BEAST2)
* Using BEAST2 version 2.6.3
* Used same alignment as produced for IQ-TREE tree building (output/pathogen_trimmed.fasta) except without 5 type-A outgroup sequences (output/pathogen_trimmed_no_outgroup.fasta)
* Used tip dates
* GTR substitution model, 4 gamma-distributed rate categories, empirical base frequencies
* Strict clock with rate fixed to 0.00079 subs / site / year as in Stadler et al, PNAS, 2013
* Birth-death tree prior with serial sampling implemented in BDMM package
* Birth, death, sampling rate priors taken from Stadler et al, PNAS, 2013, table shows non-default priors
* Ran 3 independent chains for 6 million samples each
* Using un-bounded sampling proportion meant sampling  proportion highly correlated with become-uninfectius rate
* Sampling proportion in all of SHCS >= 45% (Swiss HIV Cohort Study et al, International Journal of Epidemiology, 2009)
* Sampling starts 1994 in analyzed dataset

| Parameter | Prior distribution | Notes |
| --- | --- | --- |
| reproductive number | LogN(0.5,1) | 95% interquartile range: 0.2 - 11.7 |
| become-uninfectious rate | LogN(-1, 1) | 95% interquartile range: 0.05 - 2.61 (20 years - 140 days) | 
| sampling proportion | Uniform(0.35, 0.75) | Note: fixed to 0 until time 23.4 (one month prior to first sample) |
| time of outbreak origin | LogN(3.3, 0.1)| 95% interquartile range: 22.3 - 33 years before first sample |

* I had a really hard time getting this to converge, after 6 million steps the posterior was still increasing almost linearly -->


### Root the phylogeny

* ```Rscript scripts/R/root_tree.R``` roots the phylogeny with type "A" sequences as the outgroup, then removes the outgroup.

### Fit the POUMM, apply phylogenetic correction to spVL trait
* ```fit_poumm.R``` fits the POUMM to the phylogeny and calculated spVL values, generating maximum-likelihood parameter estimates.
* ```correct_trait.R``` generates estimates for individual-specific viral, environmental parts of trait using POUMM parameters and trait values.

Note: results in `output_revisions` are from fitting the POUMM and correcting trait values based on a different approximate ML tree output by IQ-TREE using the `-wt` parameter.

### Prepare human genotype data
* `scripts/R/filter_gwas_individuals.R` generates a list of SHCS individuals of European descent carrying subtype B HIV.
* Filter the host genotype files based on individuals to keep, variant thresholds.
* Summarize allele frequencies, missingness in filtered human genotype data.
```
docker build -t prep-gwas-files -f Dockerfile-prep-gwas-files .
# Connect to smb://d.ethz.ch/groups/bsse/stadler/
docker run \
--volume=/Volumes/stadler/SHCSData/data:/data:ro \
--volume=`pwd`/output:/output prep-gwas-files
```
* ```make_gwas_phenotypes_file.R``` generates a phenotype file for GWAS with raw and estaimated environmental-only trait values.

### Run comparative GWAS (PLINK)
* Get top 5 principal components (PCs) of host genetic variation.
* Run GWAS using PLINK with sex, top 5 PCs as covariates.
* Filter PLINK results to SNPs only, not covariates and add p-value, effect size columns.
```
docker build -t run-gwas -f Dockerfile-run-gwas .
docker run --volume=`pwd`/output:/output run-gwas
```

<!-- ### Run comparative GWAS (ATOMM)
* Format host genotype data into haplotype genotype matrices.
```
docker build -t plink -f Dockerfile-plink .
# Run image interactively:
docker run -it --rm --mount type=bind,src=$PWD/scripts,dst=/scripts --mount type=bind,src=$PWD/output,dst=/output plink
# Run the commands in scripts/bash/get_atomm_host_genotypes.sh
```
* Format pathogen genotype data into haplotype genotype matrices.
```
Rscript get_atomm_pathogen_genotypes.R
```
* Format phenotype data into required tabular format
* Run program as in demo.m, testing only marginal effects on the host side 

I didn't continue this analysis since I got access to the data used in the ATOMM study to run our method on.
Anyways I was having trouble getting the awk-based replacement to work to transform the host genotypes. -->

## Apply method to GWAS for A. thaliana genetic determinants of QDR against X. arboricola

### Data

* A. thaliana host genetic data from 1001 genomes project (https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz)

* QDR (Quantitative Disease Resistance) measurements provided by ATOMM study authors (full_phenotype_data.txt).

* X. arboricola pathogen genetic data from ATOMM study (pathogen_alleles.txt)

### Build pathogen phylogeny (Neighbor-Joining) and calculate average QDR trait

* `scripts/R/make_nj_tree_xanthamonas.R`

### Fit the POUMM, apply phylogenetic correction to QDR trait

* `scripts/R/fit_poumm_xanthamonas.R` fits the POUMM to the tree and the mean QDR score across all hosts & leaves
* `scripts/R/correct_trait_xanthamonas.R` uses the maximum posterior POUMM parameters to estiamtes the pathogen part of the mean QDR score for each pathogen strain, then randomly selects one pathogen/host pairing per host type and subtracts the pathogen part for the respective pathogen from the mean QDR score for that host

On Euler cluster:
```
cd $SCRATCH/arabidopsis_gwas
wget https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz

# Transfer phenotype file and list of samples to filter host VCF to Euler (gwas_host_ids.txt and gwas_phenotypes.txt)

# Filter host VCF file to only samples with phenotypes
env2lmod
module load bcftools/1.12 htslib/1.12
bsub "tabix -p vcf 1001genomes_snp-short-indel_only_ACGTN.vcf.gz"
bsub "bcftools view --samples-file gwas_host_ids.txt 1001genomes_snp-short-indel_only_ACGTN.vcf.gz > host_genotypes.vcf"

# Get host genotypes in PLINK bed format
bsub -N "$HOME/programs/plink2 --vcf host_genotypes.vcf --max-alleles 2 --make-bed --out arabidopsis"

# Filter data to clean dataset
bsub -N "$HOME/programs/plink2 \
--bfile arabidopsis \
--maf 0.1 \
--make-bed \
--out arabidopsis.filtered"

# Get top 5 PCs for covariates based on all SNPs
bsub -N "$HOME/programs/plink2 \
--bfile arabidopsis.filtered \
--pca \
--out arabidopsis.filtered"

awk '{print $1,$2,$3,$4,$5,$6,$7}' arabidopsis.filtered.eigenvec >  arabidopsis.filtered.pc5.txt

# Run GWAS
bsub -N "$HOME/programs/plink2 \
--bfile arabidopsis \
--pheno gwas_phenotypes.txt \
--glm \
--covar arabidopsis.filtered.pc5.txt \
--out arabidopsis.filtered"

# Paste results from both trait values together
paste "arabidopsis.filtered.trait.glm.linear" "arabidopsis.filtered.h.MWA.glm.linear" | awk '{print $1, $2, $3, $9, $12, $21, $24}' > gwas_results.txt
HEADER="CHROM POS ID BETA_standard P_standard BETA_corrected P_corrected"
sed -i.bak "1 s/^.*$/$HEADER/" gwas_results.txt

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
$11 = abs($10)} 1' < gwas_results.txt > gwas_results.pvals.txt

# Add beta (effect size) column to results
awk '
function abs(v) { return v < 0 ? -v : v }
BEGIN { OFS = " " }
NR == 1 {
$12 = "abs_BETA_standard_minus_corrected" }
NR >= 2 {
$12 = abs($4 - $6) } 1' < gwas_results.pvals.txt > gwas_results.pvals.effectsize.txt

```
