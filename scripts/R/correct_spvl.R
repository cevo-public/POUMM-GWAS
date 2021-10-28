# This script is to apply phylogenetic correction to spVL using POUMM parameters.
# Run as: Run as: /home/nadeaus/programs/R-3.6.2/bin/Rscript /links/groups/stadler/SHCSData/analysis/200207_redo_whole_pipeline/correct_spvl.R

OUTDIR <- "output"
DATADIR <- "/Volumes/stadler/SHCSData/"

require(ape)
require(POUMM)

source("scripts/R/functions/POUMM_utility_functions.R")

# Load the tree
print("Loading POUMM fit and tree.")
tree <- ape::read.tree(file = paste(OUTDIR, "pathogen_no_outgroup.newick", sep = "/"))

# Parse trait data from tip labels (this ensures tips and trait values are in same order)
tip_data <- data.frame(
  label = tree$tip.label) %>%
  tidyr::separate(label, into = c("id", "sampledate", "is_outgroup", "trait"), sep = "_", remove = F) %>%
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
paste("Inferring maximum likelihood viral and non-viral parts of trait.")
inference_results <- do.call(
  what = inferVHfromPOUMMfit,
  args = c(POUMM_params, list(tree = tree, z = tip_data$trait)))

inference_results <- cbind(tip_data, inference_results)

# Write out phylo-corrected trait values
write.csv(
  x = inference_results,
  file = paste(OUTDIR, "poumm_corrected_traits.csv", sep = "/"),
  row.names = F)

