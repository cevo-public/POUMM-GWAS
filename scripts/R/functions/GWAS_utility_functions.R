summarizeAssociationResults <- function(epidemic) {
  # Calculate p-values for association tests
  h.MWA.p.values = runAssociationTests(
    pred.vars = epidemic$simulated.h.values$genotypes,
    response.var = epidemic$inferred.h.values$h.MWA)
  h.OU.p.values = runAssociationTests(
    pred.vars = epidemic$simulated.h.values$genotypes,
    response.var = epidemic$inferred.h.values$h.OU)
  h.plus.e.p.values = runAssociationTests(
    pred.vars = epidemic$simulated.h.values$genotypes,
    response.var = epidemic$simulated.h.values$h + epidemic$simulated.e.values$e)
  norm.z.p.values = runAssociationTests(
    pred.vars = epidemic$simulated.h.values$genotypes,
    response.var = epidemic$z - mean(epidemic$z))

  # make p-value dataframe
  is.causal <- getCausalVariants(variant.names = epidemic$simulated.h.values$variant.names)
  p.values <- data.frame(
    variant = epidemic$simulated.h.values$variant.names,
    is.causal.variant = is.causal,
    h.MWA.p.value = h.MWA.p.values,
    h.OU.p.value = h.OU.p.values,
    h.plus.e.p.value = h.plus.e.p.values,
    norm.z.p.value = norm.z.p.values)
  return(p.values)
}

runAssociationTests <- function(pred.vars, response.var) {
  # Returns a vector of p-values, 1 for each column in the pred.vars dataframe.
  p.values <- apply(pred.vars, 2, fitLinearModel, response.var = response.var)
  return(as.vector(p.values))
}

fitLinearModel <- function(pred.var, response.var) {
  # Estimates the linear slope coefficient, returns p-value for it.
  linear.model <- lm(response.var ~ pred.var)
  linear.model.coefficients <- summary(linear.model)$coefficients
  if (dim(linear.model.coefficients)[1] == 1) {
    warning("No linear model p value estimated if response variable constant.")
    p.value <- NA
  } else {
    p.value <- linear.model.coefficients[2, "Pr(>|t|)"]
  }
  return(p.value)
}

calculateGWASpower <- function(p.values, causal.variants = NA, alpha) {
  if (!(is.na(causal.variants))) {
    p.values <- p.values[causal.variants]
  }
  false.negatives <- p.values > alpha
  type.II.error <- sum(false.negatives)/length(p.values)
  return(1 - type.II.error)
}

runPGLSTests <- function(pred.vars, response.var, corr.structure, weights) {
  # Returns a vector of p-values from PGLS, 1 for each col in the pred.vars df.
  p.values <- apply(pred.vars, 2, fitPGLS, 
                    response.var = response.var, 
                    corr.structure = corr.structure,
                    weights = weights)
  return(as.vector(p.values))
}

fitPGLS <- function(pred.var, response.var, corr.structure, weights) {
  geno.pheno.data <- data.frame(
    pred.var = pred.var, 
    response.var = response.var,
    weights = weights)
  
  # get GLS coefficient estimates
  # for non-ultrametric tree: http://blog.phytools.org/2012/04/using-nlmegls-for-phylogenetic.html
  pgls.fit <- nlme::gls(
    response.var ~ pred.var, 
    data = geno.pheno.data,
    correlation = corr.structure, 
    weights = nlme::varFixed(~ weights))
  pgls.summary <- summary(pgls.fit)
  p.value <- pgls.summary$tTable["pred.var", "p-value"]
  if (p.value > 0.8) {
    print(summary(pgls.fit)) ###
  }
  return(p.value)
}
