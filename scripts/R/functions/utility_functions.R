loadOrGenerateAndSaveTree <- function(N.samples = NA, tree.filename = NA, ...) {
  if (is.na(tree.filename)) {
    tree <- ape::rtree(N.samples)
  } else {
    warning(paste("Attempting to load", tree.filename, "from", getwd()))
    tree <- readRDS(paste(getwd(), tree.filename, sep = "/"))
  }
  warning(paste("Saving tree to", getwd()))
  saveRDS(tree, file = "./tree.RData")
  return(tree)
}

correctSignificanceLevelForMultTesting <- function(significance.level, 
                                                   significance.correction = NA, 
                                                   N.tests = NULL) {
  if (is.na(significance.correction)) {
    return(significance.level)
  } else if (is.null(N.tests)) {
    stop("Must specify the number of tests performed.")
  }
  if (significance.correction == "GenomeWideSignificance") {
    return(5E-8) 
  } else if (significance.correction == "Bonferroni") {
    return(significance.level/N.tests)
  } else {
    error("significance.correction must be one of [NA, 'GenomeWideSignificance', 'Bonferroni']")
  }
}

writeParametersToStdOut <- function(params) {
  names <- names(params)
  for (i in 1:length(params)) {
    param.name <- names[[i]]
    param.vals <- params[[i]]
    concat.param.vals <- paste(param.vals, collapse = " ")
    cat(paste(param.name, ": ", concat.param.vals, "\n", sep=""))
  }
}

makeDataFrameFromMapplyResults <- function(mapply.results) {
  dataframe <- as.data.frame(t(as.data.frame(mapply.results)))
  cols_to_unnest <- colnames(dataframe)
  dataframe$runID <- 1:dim(dataframe)[1]
  dataframe <- tidyr::unnest(dataframe, cols = all_of(cols_to_unnest))
  return(dataframe)
}

isDoubleList <- function(x) {
  # Returns T if x is a list of type double, F otherwise.
  if(is.vector(x) & all(is.numeric(x))) {
    return(T)
  } else {
    return(F)
  }
}

echoParameters <- function(parameterlist) {
  # TODO: have this handle non-constant parameters (value = a list)
  # Returns a list of parameter names and values sep
  outstring <- ""
  l <- length(parameterlist)
  i = 1
  while (i < l) {
    name <- names(parameterlist)[i]
    value <- parameterlist[[i]]
    if (length(value) > 1) {
      value <- paste(value, collapse=", ")
    }
    outstring <- paste(outstring, name, ": ", value, ", ", sep="")
    i <- i + 1
  }
  outstring <- paste(outstring, names(parameterlist)[l], ": ", parameterlist[[l]], sep="")
  return(outstring)
}

permuteParameters <- function(N.reps, parameterlist) {
  # Returns parameter dataframe with all parameter permulations.
  validateParameters(parameterlist)
  param.grid <- expand.grid(parameterlist)
  param.grid <- do.call(
    what = "rbind", 
    args = replicate(
      n = N.reps, 
      expr = param.grid, 
      simplify = FALSE))
  param.grid$runID <- 1:dim(param.grid)[1]
  return(param.grid)
}

validateParameters <- function(parameterlist) {
  # Check input is valid
  # TODO
}

wrapString <- function(string, width = 80) {
  return(paste(strwrap(x = string, width = width), collapse = "\n"))
}

# Credit here: https://gist.github.com/bshishov/5dc237f59f019b26145648e2124ca1c9
calculateMAPE <- function(predicted, actual, benchmark = NA, pad.by = 0) {
  # Returns mean absolute percent error between actual and predicted values 
  #   (scalar).
  return(100 * sum(abs((actual - predicted) / (actual + pad.by))) / length(actual))
}

calculateRMSE <- function(predicted, actual, benchmark = NA) {
  # Returns root mean square error between actual and predicted values (scalar).
  return (sqrt(sum((actual - predicted)^2) / length(actual)))
}

calculateBRAE <- function(predicted, actual, benchmark) {
  # Returns bounded relative absolute error (1 x N where N = number of values).
  abs.error <- abs(actual - predicted)
  abs.benchmark.error <- abs(actual - benchmark)
  denominator <- abs.error + abs.benchmark.error
  if (any(denominator == 0)) {
    warning("Predicted and benchmark estimates were both perfect for at least
            one sample. BRAE is defined to be 0.5 in this case.")
    return(0.5)
  }
  return(abs.error/denominator)
}

calculateMBRAE <- function(predicted, actual, benchmark) {
  # Returns mean bounded relative absolute error (scalar).
  BRAE <- calculateBRAE(predicted, actual, benchmark)
  return(mean(BRAE))
}

calculateUMBRAE <- function(predicted, actual, benchmark) {
  # Returns unscaled mean bounded relative absolute error (scalar).
  MBRAE <- calculateMBRAE(predicted, actual, benchmark)
  if (MBRAE == 1) {
    warning("Benchmark estimate was perfect for all samples. 
            Returning Inf for UMBRAE.")
  }
  return(MBRAE/(1 - MBRAE))
}

addParamsToDataFrame <- function(df, param.list) {
  param.df <- as.data.frame(param.list)
  df <- merge(df, param.df)
  return(df)
}

scaleBranchLengths <- function(tree, mean.root.tip.time) {
  require(POUMM)
  mean.t <- mean(POUMM::nodeTimes(tree = tree, tipsOnly = TRUE))
  tree.scale.factor <- mean.root.tip.time/mean.t  
  scaled.tree <- tree
  scaled.tree$edge.length <- tree$edge.length*tree.scale.factor
  return(scaled.tree)
}

rmFiles <- function(filepath, filename, extensions) {
  # Some example PLINK extensions are "map", "ped", "bim", "fam", "log", "nosex", "txt"
  for (extension in extensions) {
    full.filename <- paste(paste(filepath, filename, sep = "/"), extension, sep = ".")
    system2(command = "rm", args = full.filename)
  }
}

getCausalVariants <- function(variant.names) {
  # determine if a variant is causal or not based on variant names of the form
  # "causal_variant_X" or "variant_X"
  split.list <- strsplit(
    x = variant.names, 
    split = "_")
  is.causal <- lapply(split.list, function(l) l[[1]])
  return(is.causal == "causal")
}

simulateHIVTreeExpBL <- function(N) {
  require(ape)
  HIV_N <- 8483 # number of tips in UK tree
  base_tree <- ape::rtree(n = HIV_N, br = rexp, rate = 1)
  HIV_tree <- scaleBranchLengths(tree = base_tree, mean.root.tip.time = 0.14)  # based on UK HIV tree (see VM email 02.07.19)
  subsampled_HIV_tree <- ape::drop.tip(
    phy = HIV_tree,
    tip = sample(x = HIV_tree$tip.label, size = HIV_N - N))
  return(subsampled_HIV_tree)
}

simulateHIVTreeGeometricBL <- function(N) {
  require(ape)
  # Simulates an HIV tree but where branch lengths can take only limited values
  generate_branch_lengths <- function(nbr) {
    geom_draws <- rgeom(n = nbr, prob = 1 / 10)
    branch_lengths <- unlist(lapply(FUN = max, X = geom_draws, 1))
    return(branch_lengths / 10)  # minimum branch length is 0.1
  }
  HIV_N <- 8483 # number of tips in UK tree
  if (N > HIV_N) {
    stop("Cannot simulate tree of > 8,483 tips.")
  }
  base_tree <- ape::rtree(n = HIV_N, br = generate_branch_lengths)  
  HIV_tree <- scaleBranchLengths(tree = base_tree, mean.root.tip.time = 0.14)  # based on UK HIV tree (see VM email 02.07.19)
  subsampled_HIV_tree <- ape::drop.tip(
    phy = HIV_tree,
    tip = sample(x = HIV_tree$tip.label, size = HIV_N - N))
  return(subsampled_HIV_tree)
}

printExpectedH2FromParameterization <- function(tree = NULL, sigma = NULL, alpha = NULL, d = 2, delta, M, p, var.z) {
  # Use to determine M, p, delta values such that H2.h is reasonable value
  if (!(is.null(tree))) {
    mean.t <- mean(POUMM::nodeTimes(tree = tree, tipsOnly = T))
    expected.var.v <- sigma^2 * (1 - exp(-2 * alpha * mean.t))/(2 * alpha)
    expected.H2.v <- expected.var.v / var.z
    print(paste("Expected H2.v:", format(expected.H2.v, digits = 2)))
  }
  expected.var.h <- delta^2 * M * d * p * (1 - p)
  print(paste("Expected var.h:", format(expected.var.h, digits = 3)))
  expected.H2.h <- expected.var.h / var.z
  print(paste("Expected H2.h:", format(expected.H2.h, digits = 3)))
}

generateOutfileName <- function(outdir, description, suffix = "txt") {
  # Generate a results filename with full path of the form /path/date_description_i.suffix"
  date <- Sys.Date()
  filename <- paste(outdir, "/", date, "_", description, ".", suffix, sep = "")
  i <- 2
  while (file.exists(filename)) {
    filename <- paste(outdir, "/", date, "_", description, "_", i, ".", suffix, sep = "")
    i <- i + 1
  }
  return(filename)
}

write_table_safe <- function(
  x, file, row.names = F, col.names = T, overwrite = F, quote = F, sep = "\t") {
  # Writes out a txt file. If the file already exists, appends timestamp to file name
  if (file.exists(file) & !(overwrite)) {
    warning(paste("Appending date & time to", filename, "because it already exists in", workdir))
    timestamp <- gsub(x = Sys.time(), " |:|-", "_")
    datetime_file <- paste(file, timestamp, sep = "_")
    write.table(
      x = x,
      file = datetime_file,
      col.names = col.names, row.names = row.names, quote = quote, sep = sep)
  } else {
    write.table(
      x = x,
      file = file,
      col.names = col.names, row.names = row.names, quote = quote, sep = sep)
  }
}
