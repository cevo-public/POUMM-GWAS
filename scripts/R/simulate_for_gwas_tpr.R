# ==============================================================================
# Test GWAS TPR improvement using MWA estimation method.
# ==============================================================================

# To run: Rscript simulate_for_gwas_tpr.R"

# Source necessary functions
source("scripts/R/functions/simulate_epidemic.R")
source("scripts/R/functions/simulation_control_functions.R")
source("scripts/R/functions/utility_functions.R")
source("scripts/R/functions/plot_functions.R")
source("scripts/R/functions/PLINK_utility_functions.R")
source("scripts/R/functions/POUMM_utility_functions.R")
source("scripts/R/functions/GWAS_utility_functions.R")

# Set seed for replicability of simulation
set.seed(1)

# Read parameters for data simulation from config file
# config_values <- yaml::read_yaml(file = "../../config-simulation.yaml")
config_values <- yaml::read_yaml(file = "config-simulation.yaml")
N <- config_values$N

# Designate output file
outfile <- generateOutfileName(outdir = "output", description = "GWAS_TPR_MLE")

# Simulate HIV phylogeny
tree <- simulateHIVTreeExpBL(N)

# Set constant parameters for data simulation
param.list <- list(
  N = N,
  g0 = config_values$g0,
  alpha = config_values$alpha,
  theta = config_values$theta,
  M = config_values$M,
  K = config_values$K,
  d = config_values$d,
  delta = config_values$delta,
  H2 = config_values$H2,
  H2.h = config_values$H2.h,
  var.z = config_values$var.z,
  N.reps = config_values$N.reps)

param.df <- permuteParameters(
  N.reps = param.list$N.reps,
  parameterlist = param.list)

# Simulate data and calculate estimator error
data.list <- mapply(
  simulateGWASPValues,
  N = param.df$N,
  g0 = param.df$g0,
  alpha = param.df$alpha,
  theta = param.df$theta,
  M = param.df$M,
  K = param.df$K,
  d = param.df$d,
  delta = param.df$delta,
  H2 = param.df$H2,
  H2.h = param.df$H2.h,
  var.z = param.df$var.z,
  MoreArgs = list(tree = tree,
                  is.MLE = config_values$is.MLE))

# Clean up and save simulation results
data <- makeDataFrameFromMapplyResults(mapply.results = data.list)

write.table(
  data,
  file = outfile,
  quote = F,
  sep = "\t",
  row.names = F)

# Run comparative PGLS simulations ---------------------------------------------

# Designate output file
outfile <- generateOutfileName(outdir = "output", description = "GWAS_TPR_PGLS_MLE")

# Simulate data and calculate estimator error
data.list <- mapply(
  simulatePGLSPValues,
  N = param.df$N,
  g0 = param.df$g0,
  alpha = param.df$alpha,
  theta = param.df$theta,
  M = param.df$M,
  K = param.df$K,
  d = param.df$d,
  delta = param.df$delta,
  H2 = param.df$H2,
  H2.h = param.df$H2.h,
  var.z = param.df$var.z,
  MoreArgs = list(tree = tree,
                  is.MLE = config_values$is.MLE,
                  is.estimatePOUMMparams = T))

# Clean up and save simulation results
data <- makeDataFrameFromMapplyResults(mapply.results = data.list)

write.table(
  data,
  file = outfile,
  quote = F,
  sep = "\t",
  row.names = F)