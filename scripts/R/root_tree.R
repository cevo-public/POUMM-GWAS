# This script is to root the HIV tree, then write out the tree w/o outgroup.
# Run as: Rscript scripts/R/root_tree.R
library(dplyr)
library(ggplot2)
library(ggtree)
library(tidyr)

OUTDIR <- "output"
DATADIR <- "/Volumes/stadler/SHCSData/data"
TREE_FN <- "pathogen.treefile"

paste("loading tree")
tree <- ape::read.tree(file = paste(OUTDIR, TREE_FN, sep = "/"))

# Parse outgroup info from tip labels
tip_data <- data.frame(
  label = tree$tip.label) %>%
  tidyr::separate(label, into = c("id", "sampledate", "is_outgroup", "trait"), sep = "_", remove = F)

# Plot tree with outgroup
p <- ggtree(tr = tree) %<+% tip_data +
  geom_tippoint(aes(color = is_outgroup))

ggsave(plot = p, paste(OUTDIR, "pathogen_tree_with_outgroup.png", sep = "/"), height = 10, width = 7, units = "in")

# Plot tree with trait
p2 <- ggtree(tr = tree) %<+% tip_data +
  geom_tippoint(aes(color = as.numeric(trait))) +
  scale_color_gradient(low = "green", high = "red", name = "Trait value")

ggsave(plot = p2, paste(OUTDIR, "pathogen_tree_with_trait.png", sep = "/"), height = 10, width = 7, units = "in")

# Root the tree
paste("Rooting tree according to outgroup and writing to file")
outgroup_ids <- as.character(tip_data[tip_data$is_outgroup == "outgroup", "label"])
rooted_tree <- ape::root(
  phy = tree,
  outgroup = outgroup_ids)

# Write out the rooted tree with only B seqeunces for POUMM fitting
rooted_tree_no_outgroup <- ape::drop.tip(
  phy = rooted_tree,
  tip = outgroup_ids)

treename <- "pathogen_no_outgroup.newick"
ape::write.tree(
  phy = rooted_tree_no_outgroup,
  file = paste(OUTDIR, treename, sep = "/"))