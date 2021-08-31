# ==============================================================================
# This script has functions to simulate an epidemic where patients contribute
# h and environment contributes e to a trait value z, pathogens contribute v.
# v values are modeled using the POUMM, h values are modeled as contributions
# from causal SNPs, and and e values are drawn from a normal dist. with mean 0.
#
# Parameters:
# tree = viral phylogenetic tree
# N = # of samples = # pathogen strains = # patients
# g0, aplha, theta, sigma = POUMM parameters
# M = # causal host SNPs
# K = total # host SNPs
# d = ploidy of the host
# p = allele frequency of SNPs
# delta = effect size of a copy of a causal SNP
# var.z = phenotypic variance in the population
#
# author: Sarah Nadeau
# date created: 19.06.2019
#
# Remember: I fix g0 to mean(z) for v estimation by MLE g0 still estimated --> alpha can be way off.
# Remember: I implicitly assume perfect linkage equilibrium with the actual causal SNP.
# Remember: the h- and v-detrmined components of z are assumed to be independent and we know this is not true
# TODO: p-values dataframe has duplicated columns from merge
# ==============================================================================

simulateEpidemic <- function(tree, N, g0, alpha, theta, M, K, d, p = NA, delta = NA, H2, H2.h, var.z, is.MLE = T) {
  # Simulate data for GWAS
  sigma <- calculateSigmaFromParameterization(
    tree = tree,
    alpha = alpha,
    H2 = H2,
    var.z = var.z)
  simulated.v.values <- simulateVirusValues(
    tree = tree,
    N = N,
    g0 = g0,
    alpha = alpha,
    theta = theta,
    sigma = sigma)
  p.delta <- calculatePorDelta(
    p = p, 
    delta = delta, 
    M = M, 
    H2.h = H2.h, 
    var.z = var.z,
    d = d)
  p <- p.delta$p
  delta <- p.delta$delta
  simulated.h.values <- simulateHostValues(
    N = N,
    M = M,
    K = K,
    d = d,
    p = p,
    delta = delta)
  simulated.e.values <- simulateEnvironmentalValuesFromVarZ(
    N = N,
    H2 = H2,
    var.gh = delta^2 * M * d * p * (1 - p),
    var.z = var.z)
  z <- getTraitValues(
    v = simulated.v.values$v,
    h = simulated.h.values$h,
    e = simulated.e.values$e)
  if (is.MLE) {
    POUMM.fit.summary <- fitPOUMMtoData(
      z = z,
      tree = tree)
    inferred.v.values <- InferVValues(
      z = z,
      tree = tree,
      POUMM.fit.summary = POUMM.fit.summary)
  } else {
    POUMM.fit.summary <- NULL
    inferred.v.values <- InferVValuesUsingTruePOUMMParams(
      z = z,
      tree = tree,
      g0 = g0,
      alpha = alpha,
      theta = theta,
      sigma = sigma,
      sigmae = sqrt((1 - H2) * var.z))
  }
  inferred.h.values <- InferHValues(
    z = z,
    v.OU = inferred.v.values$v.OU,
    v.MWA = inferred.v.values$v.MWA)
  
  # Summarize simulated data 
  epidemic <- list(
    POUMM.fit.summary = POUMM.fit.summary,
    z = z,
    sigma = sigma,
    simulated.v.values = simulated.v.values,
    simulated.h.values = simulated.h.values,
    simulated.e.values = simulated.e.values,
    inferred.v.values = inferred.v.values,
    inferred.h.values = inferred.h.values,
    empirical.H2 = var(simulated.v.values$v)/var(z),
    empirical.H2.h = var(simulated.h.values$h)/var(z),
    p = p,
    delta = delta)
  return(epidemic)
}

simulateVirusValues <- function(tree, N, g0, alpha, theta, sigma) {
  v <- POUMM::rVNodesGivenTreePOUMM(
    tree = tree,
    z0 = g0,
    alpha = alpha,
    theta = theta,
    sigma = sigma,
    sigmae = 0)[1:N]
  return(list(tree = tree, v = v))
}

calculateP <- function(deltas, M, H2.h, var.z, d) {
  # TODO: implement randomly dist. effect sizes
}

calculatePorDelta <- function(p, delta, M, H2.h, var.z, d) {
  if (!(is.na(p)) & !(is.na(delta))) {
    warning("Both p and delta specified. Delta will be calculated based on p.")
  }
  if (is.na(p)) {
    if (is.na(delta)) {
      stop("Either p or delta must be fixed. The other will be calculated.")
    }
    p <- calculatePFromMandDelta(M = M, delta = delta, d = d, H2.h = H2.h, var.z = var.z)
  } else {
    delta <- calculateDeltaFromMandP(M = M, p = p, d = d, H2.h = H2.h, var.z = var.z)
  }
  return(list(p = p, delta = delta))
}

calculateDeltaFromMandP <- function(M, p, d, H2.h, var.z) {
  delta <- sqrt(H2.h * var.z / (M * d * p * (1 - p)))
  return(delta)
}

calculatePFromMandDelta <- function(M, delta, d, H2.h, var.z) {
  c <- H2.h * var.z / (delta^2 * M * d)
  p.min.root <- (1 - sqrt(1 - 4*c))/2  # return smaller possible p value
  return(p.min.root)
}

simulateHostValues <- function(N, M, K, d, p, delta) {
  genotypes <- rbinom(
    n = N*K,
    size = d,
    p = p)
  genotypes <- matrix(
    data = genotypes,
    nrow = N,
    ncol = K)
  if (K < M ) {
    stop("Total # SNPs K must be >= # causal SNPs M.")
  }
  variant.names <- paste("causal_variant", 1:M, sep="_")
  if (K > M) {
    variant.names <- c(variant.names, paste("variant", (M + 1):K, sep="_"))
  }
  colnames(genotypes) <- variant.names
  if ( M %% 2 == 1 ) {
    stop("# causal SNPs M must be even to ensure expectation of h is 0.")
  }
  deltas <- c(rep(-delta, M / 2), rep(delta, M / 2))
  h = apply(X = genotypes[, 1:M], MARGIN = 1,
            FUN = function(x, y) {return(sum(x * y))}, y = deltas)
  return(list(h = h, genotypes = genotypes, variant.names = variant.names))
}

simulateEnvironmentalValuesFromVarZ <- function(N, H2, var.gh, var.z) {
  var.gv <- H2 * var.z
  var.e <- var.z - var.gh - var.gv
  if (abs(var.e) < 1E-12) {
    var.e <- 0  # deal with numerical imprecision which makes var.e slightly negative when it should be 0
  }
  if (var.e < 0) {
    stop(paste("variance in h:", format(var.gh, digits = 2),
               "+ variance in v:", format(var.gv, digits = 2),
               "is greater than specified variance in z:", var.z))
  }
  e <- rnorm(
    n = N,
    mean = 0,
    sd = sqrt(var.e))
  return(list(e = e, sigmae = sqrt(var.e)))
}

getTraitValues <- function(v, h, e) {
  z <- v + h + e
  return(z)
}

# TODO: figure out when adding a prior for g0 to the true value helps
fitPOUMMtoData <- function(z, tree) {
  POUMM.fit <- POUMM::POUMM(
    z = z,
    tree = tree,
    spec = POUMM::specifyPOUMM(
      # g0Prior = list("mean" = mean(z), "var" = 0),
      nSamplesMCMC = 0))
  POUMM.fit.summary <- summary(POUMM.fit)
  return(POUMM.fit.summary)
}

InferVValues <- function(z, tree, POUMM.fit.summary) {
  inference <- POUMM:::gPOUMM(
    z = z,
    tree = tree,
    g0 = mean(z),
    alpha = POUMM.fit.summary[[1, 3]],
    theta = POUMM.fit.summary[[2, 3]],
    sigma = POUMM.fit.summary[[3, 3]],
    sigmae = POUMM.fit.summary[[4, 3]])
  v.OU <- inference$mu.g
  v.MWA <- inference$mu.g.poumm
  return(list(v.OU = v.OU, v.MWA = v.MWA))
}

InferVValuesUsingTruePOUMMParams <- function(z, tree, g0, alpha, theta, sigma, sigmae) {
  inference <- POUMM:::gPOUMM(
    z = z,
    tree = tree,
    g0 = g0,
    alpha = alpha,
    theta = theta,
    sigma = sigma,
    sigmae = sigmae)
  v.OU <- inference$mu.g
  v.MWA <- inference$mu.g.poumm
  return(list(v.OU = v.OU, v.MWA = v.MWA))
}

InferHValues <- function(z, v.OU, v.MWA) {
  h.OU <- z - v.OU
  h.MWA <- z - v.MWA
  return(list(h.OU = h.OU, h.MWA = h.MWA))
}
