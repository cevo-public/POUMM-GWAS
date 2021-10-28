# Infectious disease GWAS project

Written by Sarah Nadeau \
Last updated 30.08.2021

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

### Build pathogen phylogeny
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

* ```Rscript scripts/R/root_tree.R``` roots the phylogeny with type "A" sequences as the outgroup, then removes the outgroup.

### Fit the POUMM, apply phylogenetic correction to spVL trait
* ```fit_poumm.R``` fits the POUMM to the phylogeny and calculated spVL values, generating maximum-likelihood parameter estimates.
* ```correct_trait.R``` generates estimates for individual-specific viral, environmental parts of trait using POUMM parameters and trait values.
* ```make_gwas_phenotypes_file.R``` generates a phenotype file for GWAS with raw and estaimated environmental-only trait values.

### Prepare human genotype data
* `scripts/R/filter_gwas_individuals.R` generates a list of SHCS individuals of European descent carrying subtype B HIV.
* Filter the host genotype files based on individuals to keep, variant thresholds.
```
docker build -t prep-gwas-files -f Dockerfile-prep-gwas-files .
# Connect to smb://d.ethz.ch/groups/bsse/stadler/
docker run \
--volume=/Volumes/stadler/SHCSData/data:/data:ro \
--volume=`pwd`/output:/output prep-gwas-files
```
* Adjust the .fam file: add POUMM-calcuated GWAS trait values.
* Get top 5 principal components of host genetic variation.

### Run comparative GWAS
* Run GWAS using PLINK.
* Filter PLINK results to SNPs only, not covariates and add p-value, beta columns.

[comment]: <> (* add_sex_to_fam_SN.R: Add sex information to .fam file.)

[comment]: <> (* filter_genotype_data_and_add_sex: Filter and QC host genetic data using PLINK.)

[comment]: <> (* summarize_data.sh: Get SNP frequencies and missingness reports for filtered data.)

[comment]: <> (* generate_pc_covariates.sh: Get top 5 principal components of host genetic variation.)

[comment]: <> (* run_gwas.sh: Run GWAS using PLINK.)

[comment]: <> (* manipulate_plink_results_for_plotting.sh: Filter PLINK results to SNPs only, not covariates and add p-value, beta columns.)



