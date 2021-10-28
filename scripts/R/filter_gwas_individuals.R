# This script is to filter individuals for GWAS to those of European ancestry, carrying subtype B HIV.
# Run as: Rscript scripts/R/filter_gwas_individuals.R

library(dplyr)
library(gridExtra)
library(ggplot2)

OUTDIR <- "output"
DATADIR <- "/Volumes/stadler/SHCSData/data"

# Load principle components data
pc_scores <- read.delim(
  file = paste(DATADIR, "spVL.shcs.rs.hm3.merge.SNPqc.mds", sep = "/"),
  header = T,
  sep = "") %>%
  mutate(Population = case_when(
    Pop == "Test" ~ "Swiss HIV Cohort Study",
    Pop == "YRI" ~ "Yoruba in Ibadan, Nigeria (YRI)",
    Pop == "JPT+CHB" ~ "Japanese in Tokyo, Japan (JPT) +\nHan Chinese in Beijing, China (CHB)",
    Pop == "CEU" ~ "CEPH/Utah Collection (CEU) ",
    Pop == "MKK" ~ "Maasai in Kinyawa, Kenya (MKK)",
    Pop == "LWK" ~ "Luhya in Webuye, Kenya (LWK)",
    Pop == "CHD" ~ "Chinese in Metropolitan Denver, CO, USA (CHD)",
    Pop == "GIH" ~ "Gujarati Indians in Houston, TX, USA (GIH)",
    Pop == "TSI" ~ "Toscani in Italia (TSI)",
    Pop == "MEX" ~ "Mexican Ancestry in LA, CA, USA (MXL)",
    Pop == "ASW" ~ "African Ancestry in SW USA (ASW)",
    T ~ Pop),
  is_european = PC1 > 0.01 & PC1 < 0.025
    & PC2 > 0.035 & PC2 < 0.05
    & PC3 > -0.03 & PC3 < 0.03)

# Load subtype data from sequence metadata
seq_metadata <- read.csv(file = paste(OUTDIR, "seq_metadata.csv", sep = "/"))

gwas_metadata <- merge(
  x = seq_metadata, y = pc_scores,
  by.x = "id", by.y = "IID",
  all = T) %>%
  filter(subtype == "B" | Pop != "Test")

# Plot combined SHCS and HapMap samples along the first few PCs
pop_colors <- RColorBrewer::brewer.pal(n = 11, name = "Set3")
names(pop_colors) <- unique(gwas_metadata$Population)
swap_col <- pop_colors[9]
pop_colors[9] <- pop_colors[1]
pop_colors[1] <- swap_col

p1 <- ggplot(
  data = gwas_metadata %>% pivot_longer(cols = c(PC2, PC3), names_to = "PC_number", values_to = "PC_value"),
  aes(x = PC1, y = PC_value, fill = Population, color = is_european)) +
  geom_point(pch=21, size=3) +
  theme_bw() +
  scale_fill_manual(values = pop_colors) +
  scale_color_manual(values = c("transparent", "black")) +
  labs(y = element_blank()) +
  guides(color = "none") +
  facet_grid(PC_number ~ .)

ggsave(plot = p1, filename = paste(OUTDIR, "host_genotype_pca.png", sep = "/"))

# Get list of samples filtered based on European ancestry
gwas_individual_filtering_summary <- gwas_metadata %>%
  filter(Pop == "Test") %>%
  group_by(is_european) %>%
  summarise(n_individuals = n(),
            n_subtype_b_individuals = sum(subtype == "B"))

write.csv(x = gwas_individual_filtering_summary,
          file = paste(OUTDIR, "gwas_individual_filtering_summary.csv", sep = "/"),
          row.names = F)

gwas_individuals_to_keep <- gwas_metadata %>% filter(is_european, Pop == "Test")
print(paste("Keeping", length(unique(gwas_individuals_to_keep$id)), "individuals for GWAS."))

write.table(
  x = gwas_individuals_to_keep %>% select(id),
  file = paste(OUTDIR, "gwas_individuals_to_keep.txt", sep = "/"),
  row.names = F,
  col.names = F,
  sep = "\t",
  quote = F)