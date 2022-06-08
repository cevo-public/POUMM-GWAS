# This script is to build a Neighbor-Joining tree from allele data by way of a similarity matrix

require(dplyr)

# Load data
infile <- "data/arabidopsis_xanthomonas/pathogen_alleles.txt"
metadatafile <- "data/arabidopsis_xanthomonas/full_phenotype_data.txt"

alleles <- read.table(infile, header = T)
metadata <- read.table(metadatafile, header = T, sep = "\t")

# Generate genetic relatedness matrix
alleles_clean <- alleles %>% 
    select(-c(count, freq1, freq2, stat_pro, df, p)) %>%
    rename_with(~ gsub("_allele$", "", .x))
G <- t(as.matrix(alleles_clean))
dist <- ape::dist.gene(G, method = "pairwise")

# Build neighbor-joining tree
tree_nj <- ape::nj(dist)

# Plot tree
png("output_arabidopsis_xanthomonas/tree.png")
plot(tree_nj)
dev.off()

# Calculate phenotype values
metadata_clean <- metadata %>%
    mutate(
        Strain = recode(
            Strain,
            "ME_DV_A37" = "MEDV_A37",
            "ME_DV_P25" = "MEDV_P25",
            "ME_DV_A40" = "MEDV_A40",
            "ME_DV_P39" = "MEDV_P39",
        )
    )

metadata_long <- metadata_clean %>% pivot_longer(
    cols = c("Leaf_1", "Leaf_2", "Leaf_3", "Leaf_4"),
    names_to = "Leaf",
    names_prefix = "Leaf_",
    values_to = "QDR_score"
)

average_phenotype <- metadata_long %>%
    mutate(QDR_score = as.integer(QDR_score)) %>%
    group_by(Strain) %>%
    summarize(
        n_leaves = sum(!is.na(QDR_score)),
        mean_QDR_score = mean(QDR_score, na.rm = T),
        sd_QDR_score = sd(QDR_score, na.rm = T)
    )

# Plot phenotype value distribution
trait_dist <- ggplot(data = average_phenotype, aes(x = mean_QDR_score)) + 
    geom_histogram(bins = 10) + 
    labs(x = "Mean QDR score")
ggsave(
    plot = trait_dist,
    file = "output_arabidopsis_xanthomonas/trait_dist.png",
    width = 5, height = 5, units = "in"
)

# Annotate tips with phenotype values
tip_label_data <- data.frame(
    id = tree_nj$tip.label,
    sampledate = NA,
    is_outgroup = F
) %>% 
    left_join(average_phenotype, by = c("id" = "Strain")) %>%
    mutate(label = paste(id, sampledate, is_outgroup, mean_QDR_score, sep = "|"))

tree_nj$tip.label <- tip_label_data$label
ape::write.tree(tree_nj, "output_arabidopsis_xanthomonas/pathogen_nj.newick")

# Write out phenotype data
write.table(x = average_phenotype, file = "output_arabidopsis_xanthomonas/trait_data.txt")
