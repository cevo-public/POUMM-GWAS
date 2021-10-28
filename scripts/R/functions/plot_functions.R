chooseErrorFunction <- function(ERROR_METRIC) {
  # Returns function for error metric calculation.
  if (ERROR_METRIC == "MAPE__") {
    return(calculateMAPE)
  } else if (ERROR_METRIC == "RMSE") {
    return(calculateRMSE)
  } else if (ERROR_METRIC == "MBRAE") {
    return(calculateMBRAE)
  } else if (ERROR_METRIC == "UMBRAE") {
    return(calculateUMBRAE)
  } else {
    stop("Unknown error metric.")
  }
}  # TODO: implement a scale-invariant error metric (sMAPE?)

addNextHighestValueToDataFrame <- function(df, colname) {
  values <- unique(df[, colname])
  if (length(values) < 2) {
    stop("Attempting to add list of values with all the same value.")
  }
  nextcol <- paste(colname, "next", sep=".")
  df[, nextcol] <- unlist(lapply(df[, colname], getNextHighestValueInList, values))
  return(df)
}

getNextHighestValueInList <- function(value, list) {
  list <- sort(list, decreasing = F)
  next.highest.value <- list[which(list > value)[1]]
  if (is.na(next.highest.value)) {
    next.highest.value <- max(list) + max(list) - list[length(list) - 1]
  }
  return(next.highest.value)
}

plotHeatmap <- function(df, plotname = NULL, N.sim, is.MLE,
                        scale.limits, legend.name, fill.lo, fill.hi,
                        fill.mid = NA, midpoint = NA, x.var, y.var, fill.var,
                        facet.x = NA, facet.y = NA, is.oob.squish = T,
                        is.real.scale = T, labeller = plot_labeller) {
  if (!(is.na(facet.x)) | !(is.na(facet.y))) {
    facet.spec <- constructFacettingSpec(facet.x, facet.y, labeller)
    no_facet <- F
  } else {
    no_facet <- T
  }
  oob.choice <- getOOBFunction(is.oob.squish)
  color.scale <- constructColorScale(fill.lo, fill.hi, legend.name, scale.limits, oob.choice, fill.mid, midpoint)

  # get tile bounds
  if(is.real.scale) {
    df <- addNextHighestValueToDataFrame(df = df, colname = x.var)
    df <- addNextHighestValueToDataFrame(df = df, colname = y.var)
    tiles <- geom_rect(aes_string(xmin = x.var,
                                  xmax = paste(x.var, "next", sep="."),
                                  ymin = y.var,
                                  ymax = paste(y.var, "next", sep=".")))
    x.axis <- scale_x_continuous(expand = c(0, 0), limits = c(0, NA))
    y.axis <- scale_y_continuous(expand = c(0, 0))
  } else {
    x.levels.df <- generateVarLevelsDF(df = df, var = x.var, axis = "x")
    y.levels.df <- generateVarLevelsDF(df = df, var = y.var, axis = "y")
    df <- merge(df, x.levels.df)
    df <- merge(df, y.levels.df)
    tiles <- geom_rect(aes_string(xmin = "x.start",
                                  xmax = "x.end",
                                  ymin = "y.start",
                                  ymax = "y.end"))
    x.axis <- createAxis(var.levels.df = x.levels.df, var = x.var, is.x = T)
    y.axis <- createAxis(var.levels.df = y.levels.df, var = y.var, is.x = F)
  }

  # warn user to summarize replicates for each tile
  if (N.sim > 1) {
    warning("Have you already summarized your replicates?")
  }

  # construct heatmap
  if (no_facet) {
    heatmap <- ggplot(
      data = df,
      aes_string(x = x.var, y = y.var, fill = fill.var)) +
      tiles +
      x.axis +
      y.axis +
      color.scale +
      labs(title = plotname)
    return(heatmap)
  }
  heatmap <- ggplot(
    data = df,
    aes_string(x = x.var, y = y.var, fill = fill.var)) +
    tiles +
    x.axis +
    y.axis +
    color.scale +
    facet.spec +
    labs(title = plotname)
  return(heatmap)
}

constructFacettingSpec <- function(facet.x, facet.y, labeller) {
  if (!is.na(facet.x) & !is.na(facet.y)) {
    facet.spec <-facet_grid(reformulate(facet.x, facet.y),
                            labeller = labeller)
  } else if (!is.na(facet.y)) {
    facet.spec <- facet_grid(reformulate(".", facet.y),
                             labeller = labeller)
  } else {
    facet.spec <- facet_grid(reformulate(facet.x, "."),
                             labeller = labeller)
  }
  return(facet.spec)
}

getOOBFunction <- function(is.oob.squish) {
  # set whether to squish or censor out-of-bounds fill values
  if (is.oob.squish) {
    oob.choice <- scales::squish
  } else {
    oob.choice <- scales::censor
  }
  return(oob.choice)
}

constructColorScale <- function(fill.lo, fill.hi, legend.name, scale.limits, oob.choice, fill.mid = NA, midpoint = NA) {
  if (!is.na(fill.mid)) {
    if (is.na(midpoint)) {
      stop("Must specify midpoint value if fill.mid given.")
    }
    return(scale_fill_gradient2(
      low = fill.lo,
      high = fill.hi,
      name = legend.name,
      limits = scale.limits,
      oob = oob.choice,
      mid = fill.mid,
      midpoint = midpoint))
  } else {
    return(scale_fill_gradient(
      low = fill.lo,
      high = fill.hi,
      name = legend.name,
      limits = scale.limits,
      oob = oob.choice))
  }
  return(color.scale)
}

generateVarLevelsDF <- function(df, var, axis) {
  var.levels.df <- unique(df[var])
  n.breaks <- dim(as.data.frame(var.levels.df))[1]
  var.levels.df[paste(axis, "start", sep = ".")] <- 0:(n.breaks - 1)
  var.levels.df[paste(axis, "end", sep = ".")] <- 1:(n.breaks)
  return(var.levels.df)
}

createAxis <- function(var.levels.df, var, is.x) {
  n.breaks = dim(var.levels.df)[1]
  if (is.x) {
    axis <- scale_x_continuous(
      breaks = seq(0.5, n.breaks - 0.5, 1),
      labels = as.character(var.levels.df[[var]]),
      expand = c(0, 0),
      limits = c(0, n.breaks))
  } else {
    axis <- scale_y_continuous(
      breaks = seq(0.5, n.breaks - 0.5, 1),
      labels = as.character(var.levels.df[[var]]),
      expand = c(0, 0),
      limits = c(0, n.breaks))
  }
  return(axis)
}

plot_labeller <- function(variable, value) {
  # Returns facet label with parameter name prepended.
  if (is.character(value)) {
    return(value)
  } else {
    return(paste(variable, "=", value))
  }
}

plotROCCurve <- function(df, plotname = NULL, N.sim, facet.x, facet.y = NULL) {
  if (!("FPR" %in% colnames(df) & "TPR" %in% colnames(df))) {
    stop("ROC plot data must have true and false positive rate columns named 'TPR' and 'FPR'.")
  }
  ROC.dotplot <- ggplot(data = df, aes(x = FPR, y = TPR)) +
    geom_point() +
    facet_grid(reformulate(facet.x, facet.y),
               labeller = plot_labeller) +
    labs(title = paste(plotname,
                       "simulations:", N.sim))
  return(ROC.dotplot)
}

summarizeGWASResultsForPlotting <- function(df, sig.level) {
  df.long <- tidyr::gather(
    data = df,
    key = "inference.method",
    value = "p.value",
    h.plus.e.p.value, z.p.value, h.OU.p.value, h.MWA.p.value)

  df.long$TP <- ifelse(
    df.long$p.value < sig.level & df.long$is.causal.variant, T, F)
  df.long$FP <- ifelse(
    df.long$p.value < sig.level & !(df.long$is.causal.variant), T, F)

  # calculate TPR
  df.TPR <- as.data.frame(
    dplyr::summarise(
      dplyr::group_by(
        df.long,
        runID, inference.method, alpha, sigma, H2),
      TPR = sum(TP, na.rm = T)/sum(!is.na(TP) & is.causal.variant),
      FP = sum(FP),
      N.noncausal.variants = sum(!(is.causal.variant))))

  # summarize replicates
  df.summary <- as.data.frame(
    dplyr::summarise(
      dplyr::group_by(
        df.TPR,
        inference.method, alpha, sigma, H2),
      mean.TPR = mean(TPR),
      sd.TPR = sd(TPR),
      mean.FP = mean(FP),
      sd.FP = sd(FP),
      N.noncausal.variants = mean(N.noncausal.variants),
      N.sim = dplyr::n()))

  # set plot factor order
  df.summary$inference.method <- factor(
    df.summary$inference.method,
    levels = c("h.plus.e.p.value", "h.MWA.p.value", "h.OU.p.value", "z.p.value"),
    labels = c("h + e", "h_MWA", "h_OU", "z"))

  return(df.summary)
}

summarizePGLSResultsForPlotting <- function(df, sig.level) {
  df$TP <- ifelse(
    df$pgls.p.value < sig.level & df$is.causal.variant, T, F)
  df$FP <- ifelse(
    df$pgls.p.value < sig.level & !(df$is.causal.variant), T, F)

  # calculate TPR
  df.TPR <- as.data.frame(
    dplyr::summarise(
      dplyr::group_by(
        df,
        runID, alpha.x, sigma, H2.x, p),
      TPR = sum(TP, na.rm = T)/sum(!is.na(TP) & is.causal.variant),
      FP = sum(FP),
      N.noncausal.variants = sum(!(is.causal.variant))))

  # summarize replicates
  df.summary <- as.data.frame(
    dplyr::summarise(
      dplyr::group_by(
        df.TPR,
        alpha.x, sigma, H2.x),
      mean.TPR = mean(TPR),
      sd.TPR = sd(TPR),
      mean.FP = mean(FP),
      sd.FP = sd(FP),
      N.noncausal.variants = mean(N.noncausal.variants),
      N.sim = dplyr::n()))

  return(df.summary)
}

summarizeErrorResultsForPlotting <- function(df, error.cols) {
  df.long <- do.call(
    what = tidyr::gather,
    args = c(error.cols,
             list(data = df, key = "inference.method", value = "error")))
  df.summary <- as.data.frame(
    dplyr::summarise(
      dplyr::group_by(
        df.long,
        inference.method, alpha, H2),
      mean.error = mean(error),
      sd.error = sd(error),
      N.reps = dplyr::n()))
  return(df.summary)
}

makeTPRPlot <- function(df, is.MLE, scale.limits = c(0, 1),
                        fill.lo = "red", midpoint = NA, fill.mid = NA, fill.hi = "green",
                        x.var = "alpha", y.var = "H2", facet.y = "inference.method") {
  p <- plotHeatmap(
    df = df,
    plotname = element_blank(),
    N.sim = df$N.sim[1],
    is.MLE = is.MLE,
    scale.limits = scale.limits,
    legend.name = "TPR",
    fill.lo = fill.lo,
    fill.hi = fill.hi,
    fill.mid = fill.mid,
    midpoint = midpoint,
    x.var = x.var,
    y.var = y.var,
    fill.var = "mean.TPR",
    facet.y = facet.y,
    is.real.scale = F)
  return(p)
}

makeFPPlot <- function(df, is.MLE, scale.limits = c(0, 10)) {
  p <- plotHeatmap(
    df = df,
    plotname = element_blank(),
    N.sim = df$N.sim[1],
    is.MLE = is.MLE,
    scale.limits = scale.limits,
    legend.name = "mean # FP",
    fill.lo = "green",
    fill.hi = "red",
    x.var = "alpha",
    y.var = "H2",
    fill.var = "mean.FP",
    facet.y = "inference.method",
    is.real.scale = F)
  if (is.MLE) {
    p <- p + ggtitle("MLE")
  } else {
    p <- p + ggtitle("no MLE")
  }
  return(p)
}

summarizeGWASResultsForPValuePlotting <- function(df) {
  df.long <- tidyr::gather(
    data = df,
    key = "inference.method",
    value = "p.value",
    h.plus.e, z, h.OU, h.MWA)
  df.long$inference.method <- factor(df.long$inference.method, levels = c("h.plus.e", "h.MWA", "h.OU", "z"))
  return(df.long)
}

makePValuesPlot <- function(df, alpha) {
  p <- ggplot(data = df,
              aes(x = inference.method, y = -log10(p.value))) +
    geom_violin() +
    facet_grid(alpha ~ H2,
               labeller = plot_labeller) +
    geom_hline(yintercept = -log10(alpha))
  return(p)
}

extractParamSettingsFromDF <- function(df) {
  param.list <- list(
    g0 = unique(df$g0),
    alpha = unique(df$alpha),
    theta = unique(df$theta),
    M = unique(df$M),
    K = unique(df$K),
    d = unique(df$d),
    p = unique(df$p),
    delta = unique(df$delta),
    H2 = unique(df$H2),
    var.z = unique(df$var.z),
    N.reps = length(unique(df$runID)))
  return(echoParameters(param.list))
}
