# This script is to calcuate spVL for the SHCS data in several ways.
#   * liberal filter = exclude meaurements after treatment
#   * strict filter = exclude measurements possibly < 6 mo. after infection and after treatment or AIDS
#   * liberal 1: liberal filter --> take first of remaining measurements 
#   * liberal 2: liberal filter --> take mean of remaining measurements
#   * strict 1: strict filter  --> take mean of remaining measurements
# Note that spVL measured to be 0 is replaced with log10(25)

# Run as: Rscript scripts/R/calculate_spvl.R

require(dplyr)
require(ggplot2)
require(lubridate)
require(tidyr)

OUTDIR <- "output"
DATADIR <- "/Volumes/stadler/SHCSData/data"
LAB_FN <- "LAB_rna_measurement.csv"
TAILOR_1_FN <- "TAILOR_rna_measurement.csv"
TAILOR_2_FN <- "RNA_measurement_additional_info.csv"

# Load data --------------------------------------------------------------------
print("loading data")
lab_data_raw <- read.table(paste(DATADIR, LAB_FN, sep = "/"), sep = ",", header = T)
tailor_data1 <- read.table(paste(DATADIR, TAILOR_1_FN, sep = "/"), sep = ",", header = T)
tailor_data2 <- read.table(paste(DATADIR, TAILOR_2_FN, sep = "/"), sep = ",", header = T)
tailor_data_raw <- merge(x = tailor_data1, y = tailor_data2)

# Parse dates -----------------------------------------------------------------
print("Parsing dates.")

tailor_data <- tailor_data_raw %>%
  mutate(seroc_date = as.Date(seroc_date, format = "%d.%m.%y"),
         haart_start_date = as.Date(haart_start_date, format = "%d.%m.%y"),
         tri_start_date = as.Date(tri_start_date, format = "%d.%m.%y"),
         first_c_date = as.Date(first_c_date, format = "%d.%m.%y"),
         art_start_date = as.Date(art_start_date, format = "%d.%m.%y"),
         cd4_200_fd = as.Date(cd4_200_fd, format = "%d.%m.%y"),
         regdate = as.Date(regdate, format = "%d.%m.%y"),
         hiv_posdate = as.Date(hiv_posdate, format = "%d.%m.%y"),
         hiv_posdocdate = as.Date(hiv_posdocdate, format = "%d.%m.%y"))

lab_data <- lab_data_raw %>%
  mutate(labdate = as.Date(labdate, format = "%d.%m.%y"),
         rna = as.numeric(rna),
         cd4 = as.numeric(cd4))

patient_data_all <- merge(
  x = tailor_data, y = lab_data,
  by = "id")

# Apply liberal filter ---------------------------------------------------------
print("Applying liberal filter to RNA measurements")

patient_data <- patient_data_all %>% mutate(
  measuresment_status = case_when(
    is.na(rna) ~ "no_rna",
    labdate > art_start_date ~ "after_art",
    labdate > haart_start_date ~ "after_haart",
    labdate > tri_start_date ~ "after_tri",
    T ~ "valid_measurement"),
   log10rna = ifelse(test = rna == 0, yes = log10(25), no = log10(rna)))

print("Measurement failure reasons for liberal filtering:")
table(patient_data$measuresment_status)

# Filter lab data based on these exclusion criteria
patient_data_filtered <- patient_data %>% filter(measuresment_status == "valid_measurement")

# Get liberal spVL 1 -----------------------------------------------------------
spvl_liberal_1 <- patient_data_filtered %>%
  group_by(id) %>%
  arrange(labdate) %>%
  summarize(spVL = log10rna[1])

# Get liberal spVL 2 -----------------------------------------------------------
spvl_liberal_2 <- patient_data_filtered %>%
  group_by(id) %>%
  arrange(labdate) %>%
  summarize(spVL = mean(log10rna))

# Filter lab measurements for strict spVL definition ---------------------------
print("Applying strict filter to RNA measurements")
time_cutoff_days <- 6 * 30  # roughly 6 months
patient_data <- patient_data_all %>% mutate(
  measuresment_status = case_when(
    is.na(rna) ~ "no_rna",
    labdate < regdate + time_cutoff_days ~ "close_to_regdate",
    labdate < hiv_posdocdate + time_cutoff_days ~ "close_to_hiv_posdocdate",
    labdate < hiv_posdate + time_cutoff_days ~ "close_to_hiv_posdate",
    labdate < seroc_date + time_cutoff_days ~ "close_to_seroc",
    cd4 < 200 ~ "after_cd4_200",
    labdate > cd4_200_fd ~ "after_cd4_200",
    labdate > first_c_date ~ "after_aids_defining_event",
    labdate > art_start_date ~ "after_art",
    labdate > haart_start_date ~ "after_haart",
    labdate > tri_start_date ~ "after_tri",
    T ~ "valid_measurement"),
  log10rna = ifelse(test = rna == 0, yes = log10(25), no = log10(rna)))

print("Measurement failure reasons for strict filtering:")
table(patient_data$measuresment_status)

# Filter lab data based on these exclusion criteria
patient_data_filtered <- patient_data %>% filter(measuresment_status == "valid_measurement")

# Get strict spVL version 1: take average --------------------------------------
spvl_strict_1 <- patient_data_filtered %>%
  group_by(id) %>%
  summarise(spVL = mean(log10rna))

# Combine and write out all calculated spVL values
spvl_values <- rbind(
  spvl_liberal_1 %>% mutate(spvl_method = "spvl_lenient_filter_first_measurement"),
  spvl_liberal_2 %>% mutate(spvl_method = "spvl_lenient_filter_mean_measurement"),
  spvl_strict_1 %>% mutate(spvl_method = "spvl_strict_filter_mean_measurement")) %>% pivot_wider(
names_from = spvl_method,
values_from = spVL)

write.csv(x = spvl_values, file = paste(OUTDIR, "spvl_values.csv", sep = "/"))

p <- ggplot(data = spvl_values,
       aes(x = spvl_lenient_filter_first_measurement,
           y = spvl_lenient_filter_mean_measurement,
           color = spvl_strict_filter_mean_measurement)) +
  scale_color_gradient(low = "green", high = "red") +
  coord_equal() +
  geom_point() +
  theme_bw()

ggsave(plot = p, filename = paste(OUTDIR, "spvl_calculation_comparison.png", sep = "/"))
