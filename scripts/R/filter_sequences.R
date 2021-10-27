# This script is to filter the HIV pol sequence alignment to focal and outgroup subtypes only, annotate sequences with date, trait value and focal/outgroup status.
# Data file is newfasta2019-10-25.fas which contains 2012 pol sequences corresponding to 2012/2125 genotyped samples (the remaining host genotypes are not associated with a pol sequence)
# The sequences are not aligned (some are longer than others) even though gap characters have been introduced --> what does this mean??
# Dates in the fasta file correspond to "sampledate" in the subtyping data.
# Run as: Rscript scripts/R/filter_sequences.R

library(dplyr)
library(ape)
library(lubridate)
library(tidyr)

OUTDIR <- "output"
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

# Get length of sequences (non-gap characters)
get_seq_length <- function(seq) {
  non_gap_seq <- seq[!(seq %in% c("n", "N", "-"))]
  return(length(non_gap_seq))
}

n_chars_in_seqs <- unlist(lapply(FUN = get_seq_length, X = as.character(sequences)))
seq_length_data <- data.frame(n_chars_in_seq = n_chars_in_seqs, orig_label = names(n_chars_in_seqs))
tip_data_spvl_with_length <- merge(x = tip_data_spvl, y = seq_length_data, by = "orig_label", all = T)

# Write out summary of sequence data
seq_subtype_summary <- tip_data_spvl_with_length %>%
  group_by(subtypes) %>%
  summarize(n_seqs = n(),
            n_seqs_with_spvl = sum(!is.na(spvl_lenient_filter_mean_measurement)),
            n_long_enough_seqs_with_spvl = sum(!is.na(spvl_lenient_filter_mean_measurement) & n_chars_in_seq > 750)) %>%
  arrange(desc(n_seqs))

seq_date_summary <- tip_data_spvl_with_length %>%
  mutate(year = lubridate::year(sampledate)) %>%
  group_by(year) %>%
  summarise(n_seqs = n(),
            n_subtype_B_seqs = sum(subtypes == "B"),
            n_subtype_B_seqs_with_spvl = sum(subtypes == "B" & !is.na(spvl_lenient_filter_mean_measurement)),
            n_subtype_B_seqs_with_spvl_long_enough = sum(subtypes == "B" &
                                                           !is.na(spvl_lenient_filter_mean_measurement) &
                                                           n_chars_in_seq > 750)) %>%
  arrange(year)

write.csv(x = seq_subtype_summary, file = paste(OUTDIR, "seq_subtype_summary.csv", sep = "/"), row.names = F)
write.csv(x = seq_date_summary, file = paste(OUTDIR, "seq_date_summary.csv", sep = "/"), row.names = F)

# Filter and write out sequence data, metadata
seq_metadata <- tip_data_spvl_with_length %>%
  mutate(is_outgroup = case_when(subtypes == OUTGROUP_TYPE ~ "outgroup", subtypes == "B" ~ "focal", T ~ "other")) %>%
  filter((!is.na(spvl_lenient_filter_mean_measurement) & is_outgroup == "focal") |  # focal sequences
           is_outgroup == "outgroup",  # outgroup sequences
           n_chars_in_seq > 750) %>%  # seq length
  select(-c("spvl_lenient_filter_first_measurement", "spvl_strict_filter_mean_measurement", "n_different_subtypes_recorded")) %>%
  rename(subtype = subtypes) %>%
  tidyr::unite(col = "label", id, sampledate, is_outgroup, spvl_lenient_filter_mean_measurement, sep = "_", remove = F) %>%
  group_by(is_outgroup) %>%
  mutate(idx = sample(1:n(), size = n(), replace = F)) %>%
  filter(is_outgroup == "focal" | idx <= 5) %>%  # Take up to 5 randomly chosen sequences from outgroup taxa for outgroup
  select(-idx)

write.csv(x = seq_metadata, file = paste(OUTDIR, "seq_metadata.csv", sep = "/"), row.names = F)

sequences_filtered <- subset(
  x = sequences,
  subset = names(sequences) %in% seq_metadata$orig_label)

new_labels <- seq_metadata$label[match(names(sequences_filtered), seq_metadata$orig_label)]
names(sequences_filtered) <- new_labels

ape::write.FASTA(x = sequences_filtered, file = paste(OUTDIR, "sequences.fasta", sep = "/"))

