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
  
[comment]: <> (## Application)

[comment]: <> (### Data)

[comment]: <> (* Host genetic data provided by the Swiss HIV Cohort Study &#40;SHCS&#41;.)

[comment]: <> (  * genotypes &#40;spVL.shcs.rs.bed, pVL.shcs.rs.bim, and pVL.shcs.rs.fam&#41;)

[comment]: <> (  * HLA genotypes imputed with SNP2HLA &#40;spVL.shcs.rs.hla.bed, spVL.shcs.rs.hla.bim, spVL.shcs.rs.hla.fam&#41;)

[comment]: <> (* SHCS cohort scores along top principle components from genotype matrix of merged SHCS and HapMap data. Provided by Christian Thorball.)

[comment]: <> (* Viral load measurements provided by the SHCS.)

[comment]: <> (* Viral genetic data provided by the SHCS.)

[comment]: <> (  * pol gene sequences &#40;newfasta2019-10-25.fas&#41;)

[comment]: <> (  * Viral subtypes)

[comment]: <> (### Analysis)

[comment]: <> (* generate_alignment.sh: Align pol sequences with MAFFT.)

[comment]: <> (* build_tree.sh: Construct pol phylogeny with IQTREE.)

[comment]: <> (* root_tree.R: Root the phylogeny and remove the outgroup.)

[comment]: <> (* calculate_spvl.R: Calculate spVL based on viral load measurements.)

[comment]: <> (* fit_poumm.R: Fit the POUMM to the phylogeny and spVL trait values.)

[comment]: <> (* correct_spvl.R: Generate alternate GWAS endpoint using POUMM estimates and spVL trait values.)

[comment]: <> (* get_europeans_to_keep.R: Get list of individuals of SHCS individuals of European ancestry based on principle component scores compared to HapMap individuals.)

[comment]: <> (* generate_no_b_hosts_list.R: Get list of individuals carrying subtype B HIV.)

[comment]: <> (* add_sex_to_fam_SN.R: Add sex information to .fam file.)

[comment]: <> (* filter_genotype_data_and_add_sex: Filter and QC host genetic data using PLINK.)

[comment]: <> (* summarize_data.sh: Get SNP frequencies and missingness reports for filtered data.)

[comment]: <> (* make_pheno_file.R: Add alternate GWAS endpoint to .fam file.)

[comment]: <> (* generate_pc_covariates.sh: Get top 5 principal components of host genetic variation.)

[comment]: <> (* run_gwas.sh: Run GWAS using PLINK.)

[comment]: <> (* manipulate_plink_results_for_plotting.sh: Filter PLINK results to SNPs only, not covariates and add p-value, beta columns.)



