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
    file = paste("gwas_results.maxmaf.pvals.txt", sep = "/"), 
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
                CHROM = c(6),
                CHR_start = c(max_BP_position)))
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
        data = results_long,
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
        filename = paste("gwas_results_xanthomonas.png", sep = "/"),
        width = linewidth, height = height / 2, unit = "cm"
)

# Make qq-plot
# Fiter to variants with p-value calculated
gwas_results_xanthomonas_filtered <- gwas_results_xanthomonas %>% filter(!is.na(P_standard))

qq_data_standard <- data.frame(
        observed_log10p = -log10(sort(gwas_results_xanthomonas_filtered$P_standard, decreasing = F)),
        expected_log10p = -log10(ppoints(nrow(gwas_results_xanthomonas_filtered))),
        gwas_type = standard_gwas_name)
qq_data_phylo <- data.frame(
        observed_log10p = -log10(sort(gwas_results_xanthomonas_filtered$P_corrected, decreasing = F)),
        expected_log10p = -log10(ppoints(nrow(gwas_results_xanthomonas_filtered))),
        gwas_type = corrected_gwas_name)
qq_results <- rbind(qq_data_standard, qq_data_phylo)

p <- ggplot(data = qq_results, aes(x = expected_log10p, y = observed_log10p)) +
        geom_point() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                    color = "blue") +
        facet_wrap(. ~ gwas_type) +
        labs(x = expression("expected log"[10]*"(p)"), y = expression("observed log"[10]*"(p)")) +
        shared_theme
ggsave(plot = p, filename = paste("qq_plots.png", sep = "/"),
       width = linewidth, height = height * 2/3, units = "cm")