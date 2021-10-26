# This script is to filter the HIV pol sequence alignment to focal and outgroup subtypes only, annotate sequences with date, trait value and focal/outgroup status.
# Data file is newfasta2019-10-25.fas which contains 2012 pol sequences corresponding to 2012/2125 genotyped samples (the remaining host genotypes are not associated with a pol sequence)
# The sequences are not aligned (some are longer than others) even though gap characters have been introduced --> what does this mean??
# Dates in the fasta file correspond to "sampledate" in the subtyping data.
# Run as: Rscript scripts/R/filter_alignment.R

library(dplyr)
library(ape)
library(lubridate)
library(tidyr)

OUTDIR <- "../../output"
DATADIR <- "/Volumes/stadler/SHCSData/data"
FASTA_FN <- "newfasta2019-10-25.fas"
SPVL_FN <- "spvl_values.csv"
OUTGROUP_TYPE <- "A"

paste("Loading subtype data.")
subtypes_raw <- read.delim(
  file = paste(DATADIR, "subtype.csv", sep = "/"),
  header = T,
  sep = ",",
  colClasses = rep("character", 6))
subtypes <- subtypes_raw %>%
  mutate(sampledate = as.Date(sampledate, format = "%d.%m.%y"),
         seq_dt = as.Date(seq_dt, format = "%d.%m.%y")) %>%
  select(id, sampledate, subtype)

paste("Loading spVL data.")
spvl_values <- read.csv(file = paste(OUTDIR, SPVL_FN, sep = "/"))

paste("Loading sequence file.")
sequences <- ape::read.FASTA(file = paste(DATADIR, FASTA_FN, sep = "/"))

# Parse sequence names to get patient IDs, collection dates
temp <- strsplit(x = names(sequences), split = "\\|")
temp2 <- matrix(unlist(temp), ncol = 3, byrow = T)
seq_data <- data.frame(orig_label = names(sequences), id = temp2[, 2], sampledate = temp2[, 3])
seq_data$sampledate <- as.Date(seq_data$sampledate)

# Match sequence data to subtype data
tip_data <- merge(
  x = seq_data, y = subtypes,
  by = c("id", "sampledate"),
  all.x = T, all.y = F) %>%
  filter(subtype != "No sequence") %>%  # don't consider failed sequences
  group_by(orig_label, id, sampledate) %>%  # consider all subtypes when samples have more than one entry in subtype
  summarise(subtypes = paste(unique(subtype), collapse = "; "),
            n_different_subtypes_recorded = length(unique(subtype)),
            .groups = "drop")

# Match sequence data to spVL data
tip_data_spvl <- merge(
  x = tip_data, y = spvl_values,
  by = "id",
  all.x = T, all.y = F)

# Write out summary of sequence data
seq_subtype_summary <- tip_data_spvl %>%
  group_by(subtypes) %>%
  summarize(n_seqs = n(), n_seqs_with_spvl = sum(!is.na(spvl_lenient_filter_mean_measurement))) %>%
  arrange(desc(n_seqs))

seq_date_summary <- tip_data_spvl %>%
  mutate(year = lubridate::year(sampledate)) %>%
  group_by(year) %>%
  summarise(n_seqs = n(),
            n_subtype_B_seqs = sum(subtypes == "B"),
            n_subtype_B_seqs_with_spvl = sum(subtypes == "B" & !is.na(spvl_lenient_filter_mean_measurement))) %>%
  arrange(year)

write.csv(x = seq_subtype_summary, file = paste(OUTDIR, "seq_subtype_summary.csv", sep = "/"), row.names = F)
write.csv(x = seq_date_summary, file = paste(OUTDIR, "seq_date_summary.csv", sep = "/"), row.names = F)

# Filter and write out sequence data, metadata
seq_metadata <- tip_data_spvl %>%
  filter((!is.na(spvl_lenient_filter_mean_measurement) & subtypes == "B") |  # focal sequences
         subtypes == OUTGROUP_TYPE) %>%  # outgroup sequences
  select(-c("spvl_lenient_filter_first_measurement", "spvl_strict_filter_mean_measurement", "n_different_subtypes_recorded")) %>%
  rename(subtype = subtypes) %>%
  mutate(is_outgroup = case_when(subtype == OUTGROUP_TYPE ~ "outgroup", T ~ "focal")) %>%
  tidyr::unite(col = "label", id, sampledate, is_outgroup, spvl_lenient_filter_mean_measurement, sep = "_", remove = F))

write.csv(x = seq_metadata, file = paste(OUTDIR, "seq_metadata.csv", sep = "/"), row.names = F)

sequences_filtered <- subset(
  x = sequences,
  subset = names(sequences) %in% seq_metadata$orig_label)

new_labels <- seq_metadata$label[match(names(sequences_filtered), seq_metadata$orig_label)]
names(sequences_filtered) <- new_labels

ape::write.FASTA(x = sequences_filtered, file = paste(OUTDIR, "sequences.fasta", sep = "/"))

