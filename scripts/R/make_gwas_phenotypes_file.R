# This script is to make a .fam (phenotype & family relationships) file with corrected and total trait values for PLINK.
# Run as: Run as: Rscript scripts/R/make_gwas_phenotypes_file.R

require(tidyr)
require(dplyr)

OUTDIR <- "output"
FAM_FN <- "spVL.shcs.rs.filtered.fam"

# Load files
print("Loading fam and corrected trait data.")
fam <- read.delim(file =  paste(OUTDIR, FAM_FN, sep = "/"), sep = "", header = F)
colnames(fam) <- c("FID", "IID", "IID_father", "IID_mother", "sex", "phenotype")  # https://www.cog-genomics.org/plink/1.9/formats#fam

corrected_trait <- read.csv(file = paste(OUTDIR, "poumm_corrected_traits.csv", sep = "/"))

print("Merging corrected trait into fam file.")
phenotypes <- merge(
  fam, corrected_trait %>% select(id, h.MWA, trait),
  by.x = "IID", by.y = "id",
  all.x = T, all.y = F)  # may have used more individuals for POUMM fitting than have genotype data
if(any(is.na(phenotypes$h.MWA) | is.na(phenotypes$trait))) {
  stop("Some genomes after filtering don't have a phenotype!")
}
phenotypes <- phenotypes %>% select(FID, IID, h.MWA, trait)

write.table(
  x = phenotypes,
  file = paste(OUTDIR, "gwas_phenotypes.txt", sep = "/"),
  row.names = F,
  col.names = T,
  quote = F)
