# ==============================================================================
# Turn pathogen genome alignment into ATOMM-formatted genotype matrix
# Get phenotype file for ATOMM
# ==============================================================================

require(ape)
require(tidyr)
require(dplyr)
require(readr)

aln_file <- "output/pathogen_trimmed.fasta"  # prepared alignment

# Read in alignment, convert to matrix
aln <- ape::read.dna(file = aln_file, format = "fasta", as.matrix = T, as.character = T)

# Get pathogen info, including pathogen-host-phenotype correspondance
path_info <- data.frame(
    seq_name = rownames(aln)
) %>% tidyr::separate(
    col = seq_name, 
    into = c('host_id', 'date', 'is_focal', 'spvl'), 
    sep = "_",
    remove = F
)

# Get host order in host genotype file
hosts_str <- gsub("\n", "", read_file("output/atomm_sequence_host_header.txt"))
hosts <- strsplit(hosts_str, split = "\t")[[1]]
hosts <- hosts[7:length(hosts)]
host_info <- data.frame(
    atomm_sequence_host_col_index = seq_len(length(hosts)),
    host_id = hosts
) %>% tidyr::separate(
    col = host_id,
    into = c("host_id", "tmp"),
    sep = "_"
) %>% select(-tmp)

# Remove sequences not included in GWAS (hosts filtered out, outgroup)
path_host_info <- left_join(host_info, path_info, by = "host_id")

rows_to_keep <- which(rownames(aln_raw) %in% path_host_info$seq_name)
gwas_aln <- aln_raw[rows_to_keep, ]

# Get consensus base at each position (excluding missing '-')
get_mode <- function(c) {
   uniqc <- unique(c[c != '-'])
   mode <- uniqc[which.max(tabulate(match(c, uniqc)))]
   return(mode)
}
consensus <- apply(X = gwas_aln, MARGIN = 2, FUN = get_mode)

# Code SNPs as consensus (0) or non-consensus (1) ('-' assumed to be reference)
get_snp_code <- function(consensus, r) {
    non_ref <- r != consensus
    non_del <- r != '-'
    is_snp <- non_ref & non_del
    return(as.numeric(is_snp))
}
snps <- apply(
    X = gwas_aln,
    MARGIN = 1,
    FUN = get_snp_code,
    consensus = consensus
)

# Get ATOMM pathogen ID (order) info
path_order_atomm <- data.frame(
    seq_name = colnames(snps),
    path_id = seq_len(ncol(snps))
)

# Add ATOMM first 2 columns (chromosome ID, SNP ID)
atomm_genotypes <- cbind(
    rep(1, nrow(snps)),
    seq_len(nrow(snps)),
    snps
)

# Write out genotype matrix
write.table(
    x = atomm_genotypes,
    file = "output/atomm_sequence_pathogen.txt",
    row.names = F,
    col.names = F,
    sep = "\t",
    quote = F
)

# Get host sex information for covariate
fam <- read.delim(file =  "output/spVL.shcs.rs.filtered.fam", sep = "", header = F)
colnames(fam) <- c("FID", "IID", "IID_father", "IID_mother", "sex", "phenotype")  # https://www.cog-genomics.org/plink/1.9/formats#fam
fam$IID <- as.character(fam$IID)

# Write out phenotype file
path_host_info_w_orders <- merge(path_host_info, path_order_atomm)
path_host_info_w_sex <- left_join(
    path_host_info_w_orders,
    fam %>% select("IID", "sex"),
    by = c("host_id" = "IID")
)

write.table(
    x = path_host_info_w_sex,
    file = "output/atomm_metadata.txt",
    sep = "\t"
)

atomm_phenotypes <- path_host_info_w_sex %>%
    mutate(intercept = 1) %>%
    select(atomm_sequence_host_col_index, path_id, intercept, sex, spvl)

write.table(
    x = atomm_phenotypes,
    file = "output/atomm_phenotype.txt",
    quote = F,
    sep = "\t",
    col.names = F,
    row.names = F
)
