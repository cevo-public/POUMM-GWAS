# This script is to fit the POUMM.
# Run as: Rscript scripts/R/fit_poumm.R

FIGDIR <- "output_arabidopsis_xanthomonas"
OUTDIR <- "output_arabidopsis_xanthomonas"

N_SAMPLES_MCMC <- 4E6
N_ADAPT_MCMC <- 2E5
THIN_MCMC <- 1000
ACCEPTANCE_RATE_MCMC <- 0.01
TREE_FN <- paste(OUTDIR, "pathogen_nj.newick", sep = "/")
SEED <- 1500

require(parallel)
require(ape)
require(tidyr)
require(POUMM)
require(bayesmeta)
require(dplyr)
require(ggplot2)
# require(la)

set.seed(SEED)  # Seed for POUMM parameter inference

# Write out log file
con <- file(paste(OUTDIR, "fit_poumm.log", sep = "/"))
open(con, "w")
writeLines(
  con = con,
  text = paste(
    "n_samples_mcmc", N_SAMPLES_MCMC, "\n",
    "n_adapt_mcmc", N_ADAPT_MCMC, "\n",
    "thin_mcmc", THIN_MCMC, "\n",
    "acceptance_rate_mcmc:", ACCEPTANCE_RATE_MCMC, "\n"))
close(con)

# Load the tree (Important: tree must be without outgroup)
print("loading the tree")
tree <- ape::read.tree(file = TREE_FN)

# Parse trait data from tip labels
tip_data <- data.frame(
  label = tree$tip.label) %>%
  tidyr::separate(label, into = c("id", "sampledate", "is_outgroup", "trait"), sep = "\\|", remove = F) %>%
  mutate(trait = as.numeric(trait))

# Fit the POUMM
print("Inferring the POUMM parameters")
POUMM_spec <- POUMM::specifyPOUMM_ATH2tMeanSeG0(
  # par vector is c(alpha, theta, H2tMean, sigmae, g0)
  z = tip_data$trait,
  tree = tree,
  nSamplesMCMC = N_SAMPLES_MCMC,
  parPriorMCMC = function(par) {
      # alpha
      dexp(par[1], rate = 0.02, log = T) +
      # g0 or theta
      dnorm(par[2], mean = 0.4, sd = 0.2, log = T) + 
      # H2
      dunif(par[3],  min = 0, max = 1, log = T) + 
      # SigmaE
      bayesmeta::drayleigh(par[4], scale = 1/sqrt(2*0.02), log = T) + 
      # g0 or theta
      dnorm(par[5], mean = 0.4, sd = 0.2, log = T)
  },
  nAdaptMCMC = N_ADAPT_MCMC,
  accRateMCMC = ACCEPTANCE_RATE_MCMC,
  thinMCMC = THIN_MCMC)

POUMM_fit <- POUMM::POUMM(
  z = tip_data$trait,
  tree = tree,
  spec = POUMM_spec)

filename <- paste(OUTDIR, "POUMM_fit.RData", sep = "/")
save(POUMM_fit, file = filename)

# Generate figures, write out estimates
plots <- plot(POUMM_fit, doZoomIn = TRUE, doPlot = F)
p <- plots$densplot +
  scale_fill_manual(
    values = c("grey", "red", "blue"),
    aesthetics = c("fill", "color"),
    labels = c("Prior", "Chain 1", "Chain 2"),
    name = element_blank()) +
  theme_bw() +
  labs(x = "Value", y = "Density") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = "bottom")
ggsave(plot = p, filename = paste(FIGDIR, "poumm_parameter_estimates.png", sep = "/"), width = 6, height = 4, units = "in")

plots2 <- plot(POUMM_fit, doZoomIn = F, doPlot = F)
p2 <- plots2$traceplot + 
    scale_fill_manual(
        values = c("grey", "red", "blue"),
        aesthetics = c("fill", "color"),
        labels = c("Prior", "Chain 1", "Chain 2"),
        name = element_blank()) +
  theme_bw() +
  labs(x = "Iteration", y = "Density") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = "bottom")
ggsave(plot = p2, filename = paste(FIGDIR, "poumm_traceplot.png", sep = "/"), width = 6, height = 4, units = "in")  

estimates <- summary(POUMM_fit)[stat %in% c("alpha", "theta", "sigma", "sigmae", "g0", "H2tMean")]
estimates_pretty <- estimates %>%
  mutate(
    MLE = as.character(round(MLE, 2)),
    HPD_round = lapply(FUN = round, X = HPD, 2),
    HPD = unlist(lapply(FUN = paste, X = HPD_round, collapse=', '))) %>%
  select(-HPD_round)

write.csv(x = estimates_pretty, row.names = F, file = paste(OUTDIR, "poumm_parameter_estimates.csv", sep = "/"))
