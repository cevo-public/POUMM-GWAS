# ==============================================================================
# Control functions for testing GWAS methods on simulated data
# ==============================================================================

addSimulationParamsToResults <- function(df, epidemic, g0, alpha, theta, K, M, H2, H2.h, var.z, is.MLE) {
  # If POUMM parameters were inferred using MLE, add to results dataframe
  if (is.MLE) {
    df <- addParamsToDataFrame(
      df = df,
      param.list = list(
        g0.AVG = mean(epidemic$z),
        alpha.MLE = epidemic$POUMM.fit.summary[[1, 3]],
        sigma.MLE = epidemic$POUMM.fit.summary[[3, 3]],
        theta.MLE = epidemic$POUMM.fit.summary[[2, 3]],
        sigmae.MLE = epidemic$POUMM.fit.summary[[4, 3]]))
  }
  df <- addParamsToDataFrame(
    df = df,
    param.list = list(
      # Actual POUMM parameters
      g0 = g0,
      alpha = alpha,
      sigma = epidemic$sigma,
      theta = theta,
      sigmae = epidemic$simulated.e.values$sigmae,
      # Actual other fixed parameters
      K = K,
      M = M,
      p = epidemic$p,
      delta = epidemic$delta,
      H2.h = H2.h,
      H2 = H2,
      var.z = var.z,
      # Empirical variables from simulation
      empirical.H2 = epidemic$empirical.H2,
      empirical.H2.h = epidemic$empirical.H2.h))
  return(df)
}

simulateFullDataset <- function(tree, N, g0, alpha, theta, M, K, d, p = NA, delta = NA, H2, H2.h, var.z, is.MLE = T) {
  # Simulate data
  epidemic <- simulateEpidemic(
    tree = tree,
    N = N,
    g0 = g0,
    alpha = alpha,
    theta = theta,
    M = M,
    K = K, 
    d = 2, 
    p = p,
    delta = delta,
    H2 = H2,
    H2.h = 0.25,
    var.z = var.z,
    is.MLE = T)
  
  # Make results table
  v.h.e.z.data <- data.frame(
    sample = names(epidemic$z),
    # Actual and estimated v, h, e, z values
    v = epidemic$simulated.v.values$v,
    h = epidemic$simulated.h.values$h,
    e = epidemic$simulated.e.values$e,
    z = epidemic$z,
    v.OU = epidemic$inferred.v.values$v.OU,
    v.MWA = epidemic$inferred.v.values$v.MWA,
    h.OU = epidemic$inferred.h.values$h.OU,
    h.MWA = epidemic$inferred.h.values$h.MWA)
  
  # Add parameter values to results table
  v.h.e.z.data <- addSimulationParamsToResults(
    df = v.h.e.z.data,
    epidemic = epidemic, 
    g0 = g0, alpha = alpha, theta = theta, K = K, M = M, H2 = H2, H2.h = H2.h, 
    var.z = var.z, is.MLE = is.MLE)
  genotype.data <- as.data.frame(epidemic$simulated.h.values$genotypes)
  genotype.data$sample <- names(epidemic$z)
  data <- merge(v.h.e.z.data, genotype.data, by = "sample")
  return(data)
}

simulateVHEZValues <- function(tree, N, g0, alpha, theta, M, K, d, p = NA, delta = NA, H2, H2.h, var.z, is.MLE = T) {
  # Simulate data
  epidemic <- simulateEpidemic(tree, N, g0, alpha, theta, M, K, d, p, delta, H2, H2.h, var.z, is.MLE)
  # Make results table
  v.h.e.z.values <- data.frame(
    sample = names(epidemic$z),
    # Actual and estimated v, h, e, z values
    v = epidemic$simulated.v.values$v,
    h = epidemic$simulated.h.values$h,
    e = epidemic$simulated.e.values$e,
    z = epidemic$z,
    v.OU = epidemic$inferred.v.values$v.OU,
    v.MWA = epidemic$inferred.v.values$v.MWA,
    h.OU = epidemic$inferred.h.values$h.OU,
    h.MWA = epidemic$inferred.h.values$h.MWA)
  # Add parameter values to results table
  v.h.e.z.values <- addSimulationParamsToResults(
    df = v.h.e.z.values,
    epidemic = epidemic, 
    g0 = g0, alpha = alpha, theta = theta, K = K, M = M, H2 = H2, H2.h = H2.h, 
    var.z = var.z, is.MLE = is.MLE)
  return(v.h.e.z.values)
}

simulateVEHZValuesError <- function(tree, N, g0, alpha, theta, M, K, d, p = NA, delta = NA, H2, H2.h, var.z, is.MLE = T) {
  print(paste0("Simulating for accuracy with these parameters: alpha=", alpha, ", H2=", H2, ", H2.h=", H2.h, ", var.z=", var.z))
  # Simulate data
  epidemic <- simulateEpidemic(
    tree = tree,
    N = N,
    g0 = g0,
    alpha = alpha,
    theta = theta,
    M = M,
    K = K, 
    d = 2, 
    p = p,
    delta = delta,
    H2 = H2,
    H2.h = 0.25,
    var.z = var.z,
    is.MLE = T)
  
  # Extract actual and estimated v and h values, actual z values
  h <- epidemic$simulated.h.values$h
  h.OU <- epidemic$inferred.h.values$h.OU
  h.MWA <- epidemic$inferred.h.values$h.MWA
  v <- epidemic$simulated.v.values$v
  v.OU <- epidemic$inferred.v.values$v.OU
  v.MWA <- epidemic$inferred.v.values$v.MWA
  z <- epidemic$z
  
  # Calculate error between actual and estimated values
  error.values <- data.frame(
    # v estimators
    RMSE.z.v = calculateRMSE(predicted = z, actual = v),  # error for non-normalized spVL 
    RMSE.v.MWA.v = calculateRMSE(predicted = v.MWA, actual = v),
    RMSE.v.OU.v = calculateRMSE(predicted = v.OU, actual = v),
    MAPE.z.v = calculateMAPE(predicted = z, actual = v),
    MAPE.v.MWA.v = calculateMAPE(predicted = v.MWA, actual = v),
    MAPE.v.OU.v = calculateMAPE(predicted = v.OU, actual = v),
    UMBRAE.v.MWA.z = calculateUMBRAE(predicted = v.MWA, actual = v, benchmark = z),
    UMBRAE.v.OU.z = calculateUMBRAE(predicted = v.OU, actual = v, benchmark = z),
    # h estimators
    RMSE.z.h = calculateRMSE(predicted = z - mean(z), actual = h),  # error for normalized spVL
    RMSE.h.MWA.h = calculateRMSE(predicted = h.MWA, actual = h),
    RMSE.h.OU.h = calculateRMSE(predicted = h.OU, actual = h),
    MAPE.z.h = calculateMAPE(predicted = h, actual = h),
    MAPE.h.MWA.h = calculateMAPE(predicted = h.MWA, actual = h),
    MAPE.h.OU.h = calculateMAPE(predicted = h.OU, actual = h),
    UMBRAE.h.MWA.z = calculateUMBRAE(predicted = h.MWA, actual = h, benchmark = z),
    UMBRAE.h.OU.z = calculateUMBRAE(predicted = h.OU, actual = h, benchmark = z))
  
  # Store calculated values
  error.values$p <- epidemic$p
  error.values$sigma <- epidemic$sigma
  
  # Add parameter values to results table
  error.values <- addSimulationParamsToResults(
    df = error.values,
    epidemic = epidemic, 
    g0 = g0, alpha = alpha, theta = theta, K = K, M = M, H2 = H2, H2.h = H2.h, 
    var.z = var.z, is.MLE = is.MLE)
  return(error.values)
}

simulateGWASPValuesWithPLINK <- function(tree, N, g0, alpha, theta, M, K, d, p = NA, delta = NA, H2, H2.h, var.z, is.MLE = T, out = ".", plink.path, out.filename = NULL, write.out = F) {
  # The write.out option writes results to a file during mapply rather than keeping
  # all mapply results in memory until all simulations complete, then writing out
  if (write.out) {
    runID <<- runID + 1  # increment counter outside function scope
  } else {
    runID <- "NA"
  }

  # Simulate data
  epidemic <- simulateEpidemic(
    tree = tree,
    N = N,
    g0 = g0,
    alpha = alpha,
    theta = theta,
    M = M,
    K = K, 
    d = 2, 
    p = p,
    delta = delta,
    H2 = H2,
    H2.h = 0.25,
    var.z = var.z,
    is.MLE = T)
  
  # Generate temporary PLINK input files
  filename <- "temp_simulation_file"
  generateMapFile(
    variant.names = epidemic$simulated.h.values$variant.names, 
    K = K, 
    out = out, 
    filename = filename)
  generatePedFile(
    G = epidemic$simulated.h.values$genotypes, 
    N = N, 
    out = out, 
    filename = filename)
  generatePhenoFile(
    epidemic = epidemic, 
    N = N, 
    out = out, 
    filename = filename)
  
  # Feed input files to PLINK for association tests
  convertPLINKTextToBinary(
    plink.path = plink.path, 
    out = out, 
    filename = filename)
  rmFiles(
    filepath = out, 
    filename = filename, 
    extensions = c("map", "ped"))
  runPLINKAssociationTests(
    plink.path = plink.path, 
    out = out, 
    filename = filename)
  rmFiles(
    filepath = out, 
    filename = filename, 
    extensions = c("bed", "bim", "fam", "nosex", "pheno", "log"))
  
  # Report p-values from PLINK
  p.values <- summarizeQAssocResults(filepath = out, filename = filename)
  
  # Add parameter values to results table
  p.values <- addSimulationParamsToResults(
    df = p.values,
    epidemic = epidemic, 
    g0 = g0, alpha = alpha, theta = theta, K = K, M = M, H2 = H2, H2.h = H2.h, 
    var.z = var.z, is.MLE = is.MLE)
  if (write.out) {
    p.values <- addParamsToDataFrame(df = p.values, param.list = list("runID" = runID))
  }
  
  # Write results for simulation to file
  if (write.out) {
    write.table(
      p.values,
      file = paste(out, paste(out.filename, "txt", sep = "."), sep = "/"),
      quote = F,
      sep = " ",
      row.names = F,
      col.names = F,
      append = T)
  }
  return(p.values)
}

simulateGWASPValues <- function(tree, N, g0, alpha, theta, M, K, d, p = NA, delta = NA, H2, H2.h, var.z, is.MLE = T) {
  print(paste0("Simulating for GWAS + POUMM TPR with these parameters: alpha=", alpha, ", H2=", H2, ", H2.h=", H2.h, ", var.z=", var.z))
  # Simulate data
  epidemic <- simulateEpidemic(
    tree = tree,
    N = N,
    g0 = g0,
    alpha = alpha,
    theta = theta,
    M = M,
    K = K, 
    d = 2, 
    p = p,
    delta = delta,
    H2 = H2,
    H2.h = 0.25,
    var.z = var.z,
    is.MLE = T)

  # Report p-values from PLINK
  p.values <- summarizeAssociationResults(epidemic = epidemic)
  
  # Add parameter values to results table
  p.values <- addSimulationParamsToResults(
    df = p.values,
    epidemic = epidemic, 
    g0 = g0, alpha = alpha, theta = theta, K = K, M = M, H2 = H2, H2.h = H2.h, 
    var.z = var.z, is.MLE = is.MLE)
  return(p.values)
}

simulatePGLSPValues <- function(tree, N, g0, alpha, theta, M, K, d, p = NA, delta = NA, H2, H2.h, var.z, is.MLE = T, is.estimatePOUMMparams = T) {
  print(paste0("Simulating for GWAS + PGLS TPR with these parameters: alpha=", alpha, ", H2=", H2, ", H2.h=", H2.h, ", var.z=", var.z))
  # Simulate data
  epidemic <- simulateEpidemic(
    tree = tree,
    N = N,
    g0 = g0,
    alpha = alpha,
    theta = theta,
    M = M,
    K = K, 
    d = 2, 
    p = p,
    delta = delta,
    H2 = H2,
    H2.h = 0.25,
    var.z = var.z,
    is.MLE = T)
  
  # Define correlation structure
  if (is.estimatePOUMMparams) {
    cov <- POUMM::covVTipsGivenTreePOUMM(
      tree = tree,
      alpha = epidemic$POUMM.fit.summary[[1, 3]],
      sigma = epidemic$POUMM.fit.summary[[3, 3]],
      sigmae = epidemic$POUMM.fit.summary[[4, 3]])
  } else {
  cov <- POUMM::covVTipsGivenTreePOUMM(
    tree = tree,
    alpha = alpha,
    sigma = epidemic$sigma,
    sigmae = epidemic$simulated.e.values$sigmae)
  }
  corr <- cov2cor(cov)
  
  # Calculate weights for heteroskedasticity due to non-ultrametricity
  weights <- diag(cov)
  
  # Translate corr matrix to corrPhyl class for nmle
  value <- corr[lower.tri(corr)]  # lower triangular entries of corr matrix stacked column-wise into a vector 
  correlation <- nlme::corSymm(
    value = value,
    form = ~ 1,
    fixed = T) 
  correlation <- nlme::Initialize(correlation, data = data.frame(epidemic$z))
  
  getPGLSSlopePValue <- function(variant.index, z, epidemic, correlation, weights) {
    x <- epidemic$simulated.h.values$genotypes[, variant.index]
    data <- data.frame(z = z, x = x, weights = weights)
    PGLS.fit <- nlme::gls(
      z ~ x,
      data = data, 
      correlation = correlation,
      weights = nlme::varFixed(~weights))
    PGLS.summary <- summary(PGLS.fit)
    p.value <- PGLS.summary$tTable["x", "p-value"]
    return(p.value)
  }
  
  PGLS.p.values <- unlist(lapply(
    X = 1:M,
    FUN = getPGLSSlopePValue,
    z = epidemic$z, epidemic = epidemic, correlation = correlation, weights = weights))

  # Build results dataframe
  is.causal <- getCausalVariants(variant.names = epidemic$simulated.h.values$variant.names)
  p.values <- data.frame(
    variant = epidemic$simulated.h.values$variant.names,
    is.causal.variant = is.causal,
    empirical.H2 = epidemic$empirical.H2,
    pgls.p.value = PGLS.p.values)
  p.values <- addSimulationParamsToResults(
    df = p.values,
    epidemic = epidemic,
    g0 = g0, alpha = alpha, theta = theta, K = K, M = M,
    H2 = H2, H2.h = H2.h, var.z = var.z, is.MLE = is.MLE)
  return(p.values)
}

# TODO: update with new parameterization

# simulateGWASInput <- function(tree, N, g0, alpha, theta, M, K, d, p, delta, H2, var.z) {
#   epidemic <- simulateEpidemic(tree, N, g0, alpha, theta, M, K, d, p, delta, H2, var.z)
#   genotypes <- epidemic$simulated.h.values$genotypes
#   GWAS.values <- data.frame(
#     h = epidemic$simulated.h.values$h,
#     z = epidemic$z,
#     h.OU = epidemic$inferred.h.values$h.OU,
#     h.MWA = epidemic$inferred.h.values$h.MWA)
#   GWAS.values <- cbind(
#     GWAS.values, genotypes)
#   return(GWAS.values)
# }
#
# 
# simulateFalsePositives <- function(tree, N, g0, alpha, theta, M, K, d, p, delta, H2,
#                                    var.z, is.MLE, out, plink.path, chunk.size, filename, temp.filename) {
#   runID <<- runID + 1  # increment counter outside function scope
#   source("functions/PLINK_utility_functions.R")
#   
#   # check that K - M is a multiple of chunk size
#   if (!((K - M) %% chunk.size == 0)) {
#     stop("# noncausal variants (K - M) must be a multiple of chunk.size.")
#   }
#   
#   # generate pheno PLINK input files from causal variants
#   epidemic <- simulateEpidemic(
#     tree, N, g0, alpha, theta, M, K = M, d, p, delta, H2, var.z, is.MLE)
#   generatePhenoFile(epidemic, N, out, temp.filename)
#   
#   # simulate non-causal genotypes and build ped file
#   N.noncausal.variants <- K - M
#   N.chunks <- floor(N.noncausal.variants / chunk.size)
#   
#   for (i in 1:N.chunks) {  # for each non-causal variant chunk
#     start <- chunk.size * i - (chunk.size - 1)
#     end <- chunk.size * i
#     print(paste("beginning chunk of SNPS from", start, ":", end))
#     generateMapFile(variant.names = start:end, chunk.size, out, temp.filename)
#     
#     G <- rbinom(n = N * chunk.size, size = 2, prob = p)
#     dim(G) <- c(N, chunk.size)
#     generatePedFile(G, N, out, temp.filename)
#     
#     # feed input files to PLINK for association tests
#     system2(
#       command = plink.path,
#       args = c(paste("--file", paste(out, temp.filename, sep = "/")),
#                paste("--pheno", paste(out, paste(temp.filename, "pheno", sep = "."), sep = "/")),
#                paste("--out", paste(out, paste(temp.filename, i, sep = "_"), sep = "/")),
#                "--all-pheno", "--assoc", "--allow-no-sex"))
#     
#     # report p-values from PLINK
#     p.values <- summarizeQAssocResults(
#       filepath = out,
#       filename = paste(temp.filename, i, sep = "_"))
#     p.values <- addParamsToDataFrame(
#       df = p.values,
#       param.list = list(
#         "g0" = g0,
#         "alpha" = alpha,
#         "theta" = theta,
#         "sigmae" = epidemic$sigma,
#         "sigma" = epidemic$sigma,
#         "H2" = H2,
#         "p" = p,
#         "delta" = delta,
#         "N" = N,
#         "M" = M,
#         "K" = K,
#         "is.MLE" = is.MLE,
#         "var.z" = var.z,
#         "runID" = runID))
#     
#     # write out p-values to file
#     write.table(
#       p.values,
#       file = paste(out, paste(filename, "txt", sep = "."), sep = "/"),
#       quote = F,
#       sep = " ",
#       row.names = F,
#       col.names = F,
#       append = T)
#     
#     # remove temporary files
#     rmFiles(filepath = out,
#             filename = paste(temp.filename, i, sep = "_"),
#             extensions = c("nosex", "log", "h.MWA.qassoc", "h.OU.qassoc", "h.plus.e.qassoc", "z.qassoc"))
#   }
#   rmFiles(filepath = out,
#           filename = temp.filename,
#           extensions = c("map", "ped", "pheno"))
# }
# 
# checkMLE <- function(tree, N, g0, alpha, theta, M, K, d, p, delta, H2, var.z) {
#   epidemic <- simulateEpidemic(tree, N, g0, alpha, theta, M, K, d, p, delta, H2, var.z)
#   mle.param.values <- data.frame(
#     g0.estimated = mean(epidemic$z),
#     g0.mle = epidemic$POUMM.fit.summary[[9, 3]],
#     alpha.mle = epidemic$POUMM.fit.summary[[1, 3]],
#     theta.mle = epidemic$POUMM.fit.summary[[2, 3]],
#     sigma.mle = epidemic$POUMM.fit.summary[[3, 3]],
#     sigmae.mle = epidemic$POUMM.fit.summary[[4, 3]])
#   mle.param.values <- addParamsToDataFrame(
#     df = mle.param.values,
#     param.list = list(
#       "g0" = g0,
#       "alpha" = alpha,
#       "theta" = theta,
#       "sigma" = epidemic$sigma,
#       "H2" = H2,
#       "p" = p,
#       "delta" = delta,
#       "N" = N,
#       "M" = M,
#       "K" = K,
#       "is.MLE" = is.MLE,
#       "var.z" = var.z))
#   return(mle.param.values)
# }