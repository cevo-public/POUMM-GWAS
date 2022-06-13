# env2lmod
# module load r/4.1.3
# bsub -I -R "rusage[mem=40960]" "Rscript plot_manhattan_euler.R"

print("Loading libraries")
library(ggplot2)
library(dplyr)
library(tidyr)

OUTDIR <- "output_arabidopsis_xanthomonas"

print("Sourcing functions")
source("scripts/R/functions/plot_functions.R")
source("scripts/R/functions/utility_functions.R")

standard_gwas_name <- "GWAS with standard trait value"
corrected_gwas_name <- "GWAS with estimated non-pathogen part of trait"
shared_theme <- theme_bw()
linewidth <- 16.5  # cm
height <- 11  # cm

# Load data
print("Loading data")
gwas_results_xanthomonas <- read.delim(
    file = paste(OUTDIR, "gwas_results.pvals.effectsize.txt", sep = "/"), 
    sep = "", 
    header = T
)

pt_size <- 1

# Add genome position to data
chr_lengths <- gwas_results_xanthomonas %>%
        group_by(CHROM) %>%
        summarise(chr_len = max(POS))
chr_starts <- chr_lengths %>%
        mutate(CHR_start = cumsum(as.numeric(chr_len)) - chr_len) %>%
        select(-chr_len)
gwas_results_annotated <- left_join(gwas_results_xanthomonas, chr_starts, by = "CHROM")
gwas_results_annotated <- gwas_results_annotated %>%
        arrange(CHROM, POS) %>%
        mutate(BP_position = CHR_start + POS)

# This is maybe not the best way, but makes sure tick marks and labels are at the center of each chromosome
x_axis_ticks <- chr_starts
max_BP_position <- max(gwas_results_annotated[gwas_results_annotated$CHROM == 5, "BP_position"])
x_axis_ticks <- rbind(
        x_axis_ticks,
        data.frame(
                CHROM = c(5, 6),
                CHR_start = c(max_BP_position, max_BP_position + 1000)))
x_axis_ticks$midpoint <- x_axis_ticks$CHR_start + (lead(x_axis_ticks$CHR_start) - x_axis_ticks$CHR_start)/2

# Format data for plotting
results_long <- tidyr::pivot_longer(
        data = gwas_results_annotated,
        cols = c("P_standard", "P_corrected"),
        names_to = "gwas_type",
        names_prefix = "P_",
        values_to = "P")

# Fix facet labels for plotting
results_long$gwas_type <- factor(
        results_long$gwas_type,
        levels = c("standard", "corrected"),
        labels = c(standard_gwas_name, corrected_gwas_name))

# Make comparative GWAS overview plot
print("Making plot")
bonferroni_signifiance_threshold <- 0.05 / nrow(gwas_results_annotated)
overview_comparison <- ggplot(
        data = results_long %>% sample_frac(0.1),  # TODO: plot all of results, not sub-sample of them
        aes(x = BP_position, y = -log10(P))) +
        geom_point(aes(color = CHROM %% 2 == 0), size = pt_size) +
        scale_color_manual(values = c("grey", "dark grey")) +
        scale_x_continuous(
                breaks = x_axis_ticks$midpoint,
                labels = x_axis_ticks$CHROM,
                expand = c(0.01, 0.01)) +
        scale_y_continuous(limits = c(0, 16), expand = c(0, 0)) +
        labs(x = "Chromosome", y = expression('-Log'[10]*'(p)')) +
        geom_hline(
                yintercept = -log10(bonferroni_signifiance_threshold),
                linetype = "dashed") +
        shared_theme +
        theme(
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                legend.position = "none"
        ) +
        facet_grid(. ~ gwas_type)

print("Saving plot")
ggsave(
        plot = overview_comparison,
        filename = paste(OUTDIR, "gwas_results_xanthomonas.png", sep = "/"),
        width = linewidth, height = height / 2, unit = "cm"
)