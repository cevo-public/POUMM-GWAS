# This function is incorrect! It was used for the initial manuscript submission
# calculateSigmaFromParameterization_old <- function(tree, alpha, H2, var.z) {
#   var.gv <- H2 * var.z
#   t.avg <- mean(POUMM::nodeTimes(tree = tree, tipsOnly = TRUE))
#   sigma <- sqrt((2 * alpha * var.gv) / (1 - exp(-2 * alpha * t.avg)))
#   return(sigma)
# }

# This function is correct! It was used for the manuscript revision
calculateSigmaFromParameterization <- function(tree, alpha, H2, var.z) {
  t.avg <- mean(POUMM::nodeTimes(tree = tree, tipsOnly = TRUE))
  sigma <- sqrt((-2 * alpha * H2 * (var.z * (0.75 - H2))) / ((H2 - 1) * (1 - exp(-2 * alpha * t.avg))))
  return(sigma)
}

permutePOUMMparameters <- function(N.simulations = 1, alphas, sigmas, thetas, 
                                   sigmaes = NA, H2s = NA, tree = NA, ...) {
  # Creates parameter dataframe with one row of parameter values per simulation. 
  # If H2 is specified, sigmae is determined from sigmae. 
  # If H2 is not specified and sigmae is, H2 is determined from sigmae.  
  #   Returns dataframe. 
  checkValidPOUMMParameters(alphas, sigmas, thetas, sigmaes, H2s, tree)
  if (isDoubleList(H2s)) {
    param.grid <- makeParamGridWithH2(alphas = alphas, 
                                      sigmas = sigmas, 
                                      thetas = thetas, 
                                      H2s = H2s, 
                                      tree = tree)
  } else if (isDoubleList(sigmaes)) {
    param.grid <- makeParamGridWithSigmae(alphas = alphas, 
                                          sigmas = sigmas, 
                                          thetas = thetas, 
                                          sigmaes = sigmaes, 
                                          tree = tree)
  } else {
    param.grid <- expand.grid(alphas = alphas, 
                              sigmas = sigmas, 
                              thetas = thetas,
                              sigmaes = sigmaes,
                              H2s = H2s)
  }
  param.grid <- do.call(what = "rbind", 
                        args = replicate(n = N.simulations, 
                                         expr = param.grid, 
                                         simplify = FALSE))
  return(param.grid)
}

checkValidPOUMMParameters <- function(alphas, sigmas, thetas, H2s, sigmaes, tree, ...) {
  if (!isDoubleList(alphas)) {
    stop("alphas must be a vector of doubles with at least one entry.")
  } else if (!isDoubleList(sigmas)) {
    stop("sigmas must be a vector of doubles with at least one entry.")
  } else if (!isDoubleList(thetas)) {
    stop("thetas must be a vector of doubles with at least one entry.")
  } else if (!(isDoubleList(sigmaes) || is.na(sigmaes))) {
    stop("POUMM sigmaes must be either NA or a vector of doubles with >= 1 entry.")
  } else if (!(isDoubleList(H2s) || is.na(H2s))) {
    stop("H2s must be either NA or a vector of doubles with >= 1 entry.")
  } else if (isDoubleList(H2s) || isDoubleList(sigmaes)) {
    if (!is.list(tree)) {
      stop("Tree must exist to calculate H2 or POUMM sigmae.")
    }
  }
}

makeParamGridWithH2 <- function(alphas, sigmas, thetas, H2s, tree) {
  param.grid <- expand.grid(alphas = alphas, 
                            sigmas = sigmas, 
                            thetas = thetas, 
                            H2s = H2s)
  param.grid$sigmaes <- mapply(calculateSigmaeFromH2, 
                               H2 = param.grid$H2s,
                               alpha = param.grid$alphas, 
                               sigma = param.grid$sigmas, 
                               MoreArgs = list(tree = tree))
  return(param.grid)
}

makeParamGridWithSigmae <- function(alphas, sigmas, thetas, sigmaes, tree) {
  param.grid <- expand.grid(alphas = alphas, 
                            sigmas = sigmas, 
                            thetas = thetas, 
                            sigmaes = sigmaes)
  param.grid$H2s <- mapply(calculateH2FromSigmae, 
                           sigmae = param.grid$sigmaes, 
                           alpha = param.grid$alphas, 
                           sigma = param.grid$sigmas, 
                           MoreArgs = list(tree = tree))
  return(param.grid)
}

calculateSigmaeFromH2 <- function(H2, alpha, sigma, tree) {
  # Computes sigmae given H2 at mean root-tip distance and other POUMM params.
  #   Returns sigmae (a scalar).
  if (H2 == 0) {
    stop("H2 must be > 0 because sigmae is undefined when H2 = 0.")
  }
  time <- mean(POUMM::nodeTimes(tree = tree, tipsOnly = T))
  var.v <- sigma^2 * (1 - exp(-2 * alpha * time)) / (2 * alpha)
  sigmae <- sqrt(var.v * (1 - H2) / H2)
  return(sigmae)
}

calculateH2FromSigmae <- function(sigmae, alpha, sigma, tree) {
  # Computes heritability at mean root-tip distance on tree given POUMM params.
  #   Returns H2 (a scalar).
  if (sigmae == 0) {
    stop("Sigma must be > 0 because sigma = 0 -> H2 = 0 -> sigmae undefined.")
  }
  time <- mean(POUMM::nodeTimes(tree = tree, tipsOnly = T))
  var.v <- sigma^2 * (1 - exp(-2 * alpha * time)) / (2 * alpha)
  H2 <- var.v/(var.v + sigmae^2)
  return(H2)
}

fitPOUMMMLE <- function(z, tree, is.fitg0 = T, ...) {
  # Fits the POUMM to a given tree and set of trait values.
  # If !is.fitg0, sets g0 to the mean trait value. 
  POUMM.fit <- POUMM::POUMM(z, tree, doMCMC = F)
  POUMM.summary <- summary(POUMM.fit)
  if (is.fitg0) {
    g0 = POUMM.summary[[9, 3]]
  } else {
    g0 = mean(z)
  }
  MLE.POUMM.params <- data.frame(
    alpha  = POUMM.summary[[1, 3]],
    theta  = POUMM.summary[[2, 3]],
    sigma  = POUMM.summary[[3, 3]],
    sigmae = POUMM.summary[[4, 3]],
    g0     = g0
  )
  return(MLE.POUMM.params)
}

inferVHfromPOUMMfit <- function(z, tree, g0, alpha, theta, sigma, sigmae, 
                                make_table = F, ...) {
  # Performs v and h inference based on fitted POUMM, tree and z values.
  #   Returns a dataframe: N samples x [v.OU = genotype value inferred from 
  #   POUMM parameters, v.MWA = genotype value inferred as matrix weighted 
  #   average of v.OU and z, h.OU = z - v.OU, h.MWA = z - v.MWA]
  #   Samples are reported in the same order as tree$tip.label
  gPOUMM.calculations <- POUMM:::gPOUMM(z = z, 
                                        tree = tree, 
                                        g0 = g0, 
                                        alpha = alpha, 
                                        theta = theta, 
                                        sigma = sigma, 
                                        sigmae = sigmae)
  v.OU <- gPOUMM.calculations$mu.g
  v.MWA <- c(gPOUMM.calculations$mu.g.poumm)
  h.OU <- c(z - v.OU)
  h.MWA <- c(z - v.MWA)
  v.h.estimates <- data.frame(v.OU = v.OU, 
                              v.MWA = v.MWA,
                              h.OU = h.OU,
                              h.MWA = h.MWA)
  return(v.h.estimates)
}
