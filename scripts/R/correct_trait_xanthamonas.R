# This script is to apply phylogenetic correction to QDR using POUMM parameters.
# Run as: Run as: Rscript scripts/R/correct_trait.R

FIGDIR <- "output_arabidopsis_xanthomonas"
OUTDIR <- "output_arabidopsis_xanthomonas"

require(ape)
require(POUMM)
require(dplyr)
require(ggtree)
require(ggplot2)

source("scripts/R/functions/POUMM_utility_functions.R")

# Load the tree
print("Loading POUMM fit and tree.")
tree <- ape::read.tree(file = paste(OUTDIR, "pathogen_nj.newick", sep = "/"))

# Parse trait data from tip labels (this ensures tips and trait values are in same order)
tip_data <- data.frame(
  label = tree$tip.label) %>%
  tidyr::separate(label, into = c("id", "sampledate", "is_outgroup", "trait"), sep = "\\|", remove = F) %>%
  mutate(trait = as.numeric(trait))

# Load the estimated POUMM parameters
POUMM_summary <- read.csv(file = paste(OUTDIR, "poumm_parameter_estimates.csv", sep = "/"), row.names = 1)

# Make data structure for inferred POUMM params
POUMM_params <- data.frame(
  alpha  = POUMM_summary['alpha', 'PostMean'],
  theta  = POUMM_summary['theta', 'PostMean'],
  sigma  = POUMM_summary['sigma', 'PostMean'],
  sigmae = POUMM_summary['sigmae', 'PostMean'],
  g0     = POUMM_summary['g0', 'PostMean'])  # using posterior means

# Apply phylo-correction
paste("Inferring maximum likelihood pathogen and non-pathogen parts of trait.")
inference_results <- do.call(
  what = inferVHfromPOUMMfit,
  args = c(POUMM_params, list(tree = tree, z = tip_data$trait)))

inference_results <- cbind(tip_data, inference_results)

# Write out phylo-corrected trait values
write.csv(
  x = inference_results,
  file = paste(OUTDIR, "poumm_corrected_traits.csv", sep = "/"),
  row.names = F)

# Plot tree with phylo-corrected trait
p <- ggtree(tr = tree) %<+% inference_results +
  geom_tippoint(aes(color = as.numeric(v.MWA))) +
  scale_color_gradient(low = "green", high = "red", name = "Phylo-estimated\npathogen trait value")

ggsave(plot = p, paste(FIGDIR, "pathogen_tree_with_phylo-estimated_pathogen_trait.png", sep = "/"), height = 10, width = 7, units = "in")

p2 <- ggtree(tr = tree) %<+% inference_results +
  geom_tippoint(aes(color = as.numeric(h.MWA))) +
  scale_color_gradient(low = "green", high = "red", name = "Phylo-estimated\nenvironmental trait value")

ggsave(plot = p2, paste(FIGDIR, "pathogen_tree_with_phylo-estimated_env_trait.png", sep = "/"), height = 10, width = 7, units = "in")

# Load full phenotype data
phenotypes <- read.table("data/arabidopsis_xanthomonas/full_phenotype_data.txt", header = T, sep = "\t") %>%
    mutate(
        Strain = recode(
            Strain,
            "ME_DV_A37" = "MEDV_A37",
            "ME_DV_P25" = "MEDV_P25",
            "ME_DV_A40" = "MEDV_A40",
            "ME_DV_P39" = "MEDV_P39",
        ),
        Leaf_1 = as.integer(Leaf_1),
        Leaf_2 = as.integer(Leaf_2),
        Leaf_3 = as.integer(Leaf_3),
        Leaf_4 = as.integer(Leaf_4)
    ) %>%
    tidyr::pivot_longer(
      cols = c("Leaf_1", "Leaf_2", "Leaf_3", "Leaf_4"), 
      values_to = "QDR_score", 
      names_to = "Leaf", 
      names_prefix = "Leaf_"
    )

# Load host genotype accession IDs from 1001 genomes project
host_genotype_accession_ids <- read.csv("data/arabidopsis_xanthomonas/host_accession_metadata.csv", header = F) %>%
    select(c(V1, V3)) %>%
    rename("id" = "V1", "Name_acc" = "V3")

# Filter to accesssions & strains with genotype data, scored leaves
phenotypes_filtered <- phenotypes %>%
  filter(!is.na(QDR_score)) %>%
  filter(Strain %in% inference_results$id) %>%
  filter(Name_acc %in% host_genotype_accession_ids$Name_acc)

# Calculate mean QDR score across 3 replicates for each plant/pathogen pairing
mean_phenotypes <- phenotypes_filtered %>%
  group_by(Strain, Name_acc) %>%
  summarize(n_leaves_scored = n(), mean_QDR_score = mean(QDR_score), .groups = "drop")

# Pick a random host/pathogen pair for each host plant (prioritizing more leaves scored)
set.seed(54321)
phenotypes_for_gwas <- mean_phenotypes %>%
    group_by(Name_acc) %>%
    slice_sample(n = 1, weight_by = n_leaves_scored)
phenotypes_for_gwas_summary <- phenotypes_for_gwas %>%
  group_by(Strain) %>%
  summarize(n_plant_accessions_paired = n())
write.table(
  x = phenotypes_for_gwas_summary,
  file = "output_arabidopsis_xanthomonas/host_pathogen_pairings_for_gwas.txt"
)

# Generate phenotypes file for PLINK (with accession ids matching vcf file from 1001 genomes)
phenotypes_for_gwas_with_id <- phenotypes_for_gwas %>% left_join(host_genotype_accession_ids)

phenotypes_for_gwas_with_correction <- phenotypes_for_gwas_with_id %>%
    full_join(inference_results, by = c("Strain" = "id")) %>%
    mutate(non_pathogen_trait = mean_QDR_score - v.MWA) %>% 
    ungroup()

phenotypes <- phenotypes_for_gwas_with_correction %>% 
  select(id, non_pathogen_trait, mean_QDR_score) %>%
  rename("IID" = "id", "h.MWA" = "non_pathogen_trait", "trait" = "mean_QDR_score")

write.table(
  x = phenotypes,
  file = "output_arabidopsis_xanthomonas/gwas_phenotypes.txt",
  row.names = F,
  col.names = T,
  quote = F)

# Write out list of samples to filter host VCF to
write.table(
  x = phenotypes$IID,
  file = "output_arabidopsis_xanthomonas/gwas_host_ids.txt",
  quote = F,
  row.names = F,
  col.names = F
)
