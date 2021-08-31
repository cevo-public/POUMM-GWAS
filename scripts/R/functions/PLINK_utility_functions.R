generateMapFile <- function(variant.names, K, out, filename) {
  # generate dummy map dataframe for plink, write out to file
  map <- data.frame(
    CHR = rep(1, length.out = K),
    variantID = variant.names,
    POS = rep(0, length.out = K))
  write.table(
    map, file = paste(out, paste(filename, "map", sep = "."), sep = "/"),
    quote = F, sep = " ", row.names = F, col.names = F)
}

generatePedFile <- function(G, N, out, filename) {
  # generate ped dataframe with dummy phenotype
  ped <- list(
    FID = rep(1, length.out = N),
    IID = 1:N,
    IID_father = rep(0, length.out = N),
    IID_mother = rep(0, length.out = N),
    sex = rep(0, length.out = N),
    phenotype = rep(-9, length.out = N))
  
  # add causal variant calls to ped dataframe
  for (j in 1:dim(G)[2]) {
    allele.1 <- rep("N", N) 
    allele.1[G[, j] > 0] <- "Y"
    allele.2 <- rep("N", N)
    allele.2[G[, j] == 2] <- "Y"
    ped[[6 + 2 * j - 1]] <- allele.1
    ped[[6 + 2 * j]] <- allele.2
  }
  attributes(ped) <- list(row.names=c(NA_integer_, N),
                        class="data.frame",
                        names=make.names(names(ped),
                                         unique=TRUE))  # for memory efficient coversion to df
  write.table(
    ped, file = paste(out, paste(filename, "ped", sep = "."), sep = "/"),
    quote = F, sep = " ", row.names = F, col.names = F)
}

generatePhenoFile <- function(epidemic, N, out, filename) {
  # generate pheno file with all phenotypes (h estimators)
  pheno <- data.frame(
    FID = rep(1, length.out = N),
    IID = 1:N,
    z = epidemic$z - mean(epidemic$z),
    h.OU = epidemic$inferred.h.values$h.OU,
    h.MWA = epidemic$inferred.h.values$h.MWA,
    h.plus.e = epidemic$simulated.h.values$h + epidemic$simulated.e.values$e)
  write.table(
    pheno, file = paste(out, paste(filename, "pheno", sep = "."), sep = "/"),
    quote = F, sep = " ", row.names = F, col.names = T)
}

convertPLINKTextToBinary <- function(plink.path, out, filename) {
  # convert plink text to plink binary to save space
  system2(command = plink.path,
          args = c(paste("--file", paste(out, filename, sep = "/")),
                   paste("--out", paste(out, filename, sep = "/"))))
}

runPLINKAssociationTests <- function(plink.path, out, filename) {
  # gather all filenames 
  # TODO: run with --ped and --map rather than --bfile filename ### left off here
  # run PLINK association test with generated bfiles
  system2(command = plink.path, 
          args = c(paste("--bfile", paste(out, filename, sep = "/")),
                   paste("--pheno", paste(out, paste(filename, "pheno", sep = "."), sep = "/")),
                   paste("--out", paste(out, filename, sep = "/")),
                   "--all-pheno", "--assoc", "--allow-no-sex"))
}

summarizeQAssocResults <- function(filepath, filename) {
  # TODO: make this function genearalizable to load arbitrary qassoc files
  # load qassoc files
  h.MWA.p.values <- read.delim(
    file = paste(paste(filepath, filename, sep = "/"), "h.MWA", "qassoc", sep = "."),
    sep = "")
  h.OU.p.values <- read.delim(
    file = paste(paste(filepath, filename, sep = "/"), "h.OU", "qassoc", sep = "."),
    sep = "")
  h.plus.e.p.values <- read.delim(
    file = paste(paste(filepath, filename, sep = "/"), "h.plus.e", "qassoc", sep = "."),
    sep = "")
  norm.z.p.values <- read.delim(
    file = paste(paste(filepath, filename, sep = "/"), "z", "qassoc", sep = "."),
    sep = "")
  
  # make p-value dataframe
  is.causal <- getCausalVariants(variant.names = as.character(h.MWA.p.values$SNP))
  p.values <- data.frame(
    variant = h.MWA.p.values$SNP,
    is.causal.variant = is.causal,
    h.MWA.p.value = h.MWA.p.values$P,
    h.OU.p.value = h.OU.p.values$P,
    h.plus.e.p.value = h.plus.e.p.values$P,
    z.p.value = z.p.values$P)
  return(p.values)
}