library("lattice")
library("svd")
library("forecast")
library(ggplot2)
source("eossa_new.R", local = TRUE)

iossa_post_grouping_trend <- function(ssa_obj,
                                signal_rank,
                                auto_trend_freq = 1/ 24,
                                auto_threshold = 0.7,
                                iossa_maxiter = 500,
                                iossa_tolerance = 1e-3,
                                kappa = 2){
  ioss <- iossa(ssa_obj, nested.groups = as.list(c(1:signal_rank)), tol = iossa_tolerance, 
                maxiter = iossa_maxiter, kappa = kappa)
  g <- grouping.auto(ioss, base = "series", groups = list(c(1:signal_rank)),
                     freq.bins = list(Trend = auto_trend_freq), 
                     threshold = auto_threshold)
  
  trend_comp <- g$Trend
  signal_indices <- 1:signal_rank
  res_comp <- signal_indices[!signal_indices %in% trend_comp]
  rec <- reconstruct(ioss, groups = list(Trend = trend_comp, Periodics = res_comp))
  
  return(list(model_object = ioss, 
              grouping = list(Trend = trend_comp, Periodics = res_comp),
              iterations = ioss$iossa.result$iter,
              trend = rec$Trend,
              periodics = rec$Periodics,
              residuals = attr(rec, "residuals")))
}

iossa_manual_trend <- function(ssa_obj,
                         signal_rank,
                         auto_trend_freq = 1 / 24,
                         auto_threshold = 0.7,
                         iossa_maxiter = 500,
                         iossa_tolerance = 1e-3,
                         manual_grouping,
                         kappa = 2){
  ioss <- iossa(ssa_obj, nested.groups = manual_grouping, tol = iossa_tolerance,
                maxiter = iossa_maxiter, kappa = kappa)
  g <- grouping.auto(ioss, base = "series", groups = manual_grouping,
                     freq.bins = list(Trend = auto_trend_freq), 
                     threshold = auto_threshold)
  
  trend_comp <- g$Trend
  signal_indices <- 1:signal_rank
  res_comp <- signal_indices[!signal_indices %in% trend_comp]
  
  rec <- reconstruct(ioss, groups = list(Trend = trend_comp, Periodics = res_comp))
  return(list(model_object = ioss,
              grouping = list(Trend = trend_comp, Periodics = res_comp),
              iterations = ioss$iossa.result$iter,
              trend = rec$Trend,
              periodics = rec$Periodics,
              residuals = attr(rec, "residuals")))
}

iossa_auto_grouping_ssa_trend <- function(ssa_obj,
                                    signal_rank,
                                    auto_trend_freq = 1 / 24,
                                    auto_threshold = 0.7,
                                    auto_initial_threshold = 0.6,
                                    iossa_tolerance = 1e-3,
                                    iossa_maxiter = 1,
                                    kappa = 2){
  auto_grouping <- grouping.auto(ssa_obj, base = "series", groups = list(c(1:signal_rank)),
                                 freq.bins = list(Trend = auto_trend_freq), 
                                 threshold = auto_initial_threshold)
  
  trend_comp <- auto_grouping$Trend
  signal_indices <- 1:signal_rank
  res_comp <- signal_indices[!signal_indices %in% trend_comp]
  
  result <- list()
  
  if ((length(trend_comp) > 0) & (length(res_comp) > 0)) {
    ioss <- iossa(ssa_obj, nested.groups = list(trend_comp, res_comp), tol = iossa_tolerance,
                  maxiter = iossa_maxiter, kappa = kappa)
    g <- grouping.auto(ioss, base = "series", groups = list(c(1:signal_rank)),
                       freq.bins = list(Trend = auto_trend_freq), 
                       threshold = auto_threshold)
    
    trend_comp <- g$Trend
    signal_indices <- 1:signal_rank
    res_comp <- signal_indices[!signal_indices %in% trend_comp]
    
    result[["model_object"]] <- ioss
    result[["iterations"]] <- ioss$iossa.result$iter
  } else {
    result[["model_object"]] <- ssa_obj
    result[["iterations"]] <- 1
  }
  
  rec <- reconstruct(result[["model_object"]], groups = list(Trend = trend_comp, Periodics = res_comp))
  
  result[["grouping"]] <- list(Trend = trend_comp, Periodics = res_comp)
  result[["trend"]] <- rec$Trend
  result[["periodics"]] <- rec$Periodics
  result[["residuals"]] <- attr(rec, "residuals")
  
  return(result)
}

iossa_auto_grouping_fossa_trend <- function(ssa_obj,
                                      signal_rank,
                                      auto_trend_freq = 1 / 24,
                                      auto_threshold = 0.7,
                                      auto_initial_threshold = 0.6,
                                      iossa_tolerance = 1e-3,
                                      iossa_maxiter = 1,
                                      kappa = 2,
                                      ...){
  foss <- fossa(ssa_obj, nested.groups = 1:signal_rank)
  auto_grouping <- grouping.auto(foss, base = "series", groups = list(c(1:signal_rank)),
                                 freq.bins = list(Trend = auto_trend_freq), 
                                 threshold = auto_initial_threshold)
  
  trend_comp <- auto_grouping$Trend
  signal_indices <- 1:signal_rank
  res_comp <- signal_indices[!signal_indices %in% trend_comp]
  
  result <- list()
  
  if ((length(trend_comp) > 0) & (length(res_comp) > 0)) {
    ioss <- iossa(foss, nested.groups = list(trend_comp, res_comp), tol = iossa_tolerance,
                  maxiter = iossa_maxiter, kappa = kappa)
    g <- grouping.auto(ioss, base = "series", groups = list(c(1:signal_rank)),
                       freq.bins = list(Trend = auto_trend_freq), 
                       threshold = auto_threshold)
    
    trend_comp <- g$Trend
    signal_indices <- 1:signal_rank
    res_comp <- signal_indices[!signal_indices %in% trend_comp]
    
    result[["model_object"]] <- ioss
    result[["iterations"]] <- ioss$iossa.result$iter
  } else {
    result[["model_object"]] <- foss
    result[["iterations"]] <- 1
  }
  
  rec <- reconstruct(result[["model_object"]], groups = list(Trend = trend_comp, Periodics = res_comp))
  
  result[["grouping"]] <- list(Trend = trend_comp, Periodics = res_comp)
  result[["trend"]] <- rec$Trend
  result[["periodics"]] <- rec$Periodics
  result[["residuals"]] <- attr(rec, "residuals")
  
  return(result)
}

grouping.auto.ssa.custom <- function(x, groups,
                                       base = c("series", "eigen", "factor"),
                                       freq.bins = 2,
                                       threshold = 0,
                                       method = c("constant", "linear"),
                                       ...,
                                       drop = TRUE) {
  
  base <- match.arg(base)
  method <- match.arg(method)
  
  if (missing(groups))
    groups <- as.list(1:min(nsigma(x), nu(x)))
  
  # Continue decomposition, if necessary
  Rssa:::.maybe.continue(x, groups = groups, ...)
  
  groups <- sort(unique(unlist(groups)))
  
  if (identical(base, "eigen")) {
    Fs <- x$U[,groups,drop=FALSE] # Rssa:::.U(x)[, groups, drop = FALSE]
  } else if (identical(base, "factor")) {
    Fs <- calc.v(x, groups, ...)
  } else if (identical(base, "series")) {
    N <- x$length
    Fs <- matrix(unlist(reconstruct(x, groups = as.list(groups), ...)), nrow = N)
  }
  
  pgs <- Rssa:::pgram(Fs)
  
  if (!is.list(freq.bins)) {
    if (length(freq.bins) == 1 && freq.bins >= 2) {
      freq.bins <- seq(0, 0.5, length.out = freq.bins + 1)[-1]
    }
    
    freq.lower.bound <- c(-Inf, head(freq.bins, -1))
    freq.upper.bound <- freq.bins
  } else {
    freq.bins <- lapply(freq.bins, function(x) if (length(x) == 1) c(-Inf, x) else x)
    freq.lower.bound <- lapply(freq.bins, function(x) x[1])
    freq.upper.bound <- lapply(freq.bins, function(x) x[2])
  }
  
  nresult <- max(length(threshold),
                 length(freq.bins))
  
  if (length(freq.lower.bound) < nresult) {
    freq.lower.bound <- rep_len(freq.lower.bound, nresult)
  }
  
  if (length(freq.upper.bound) < nresult) {
    freq.upper.bound <- rep_len(freq.upper.bound, nresult)
  }
  
  if (length(threshold) < nresult) {
    threshold <- rep_len(threshold, nresult)
  }
  
  norms <- colSums(pgs$spec)
  contributions <- matrix(NA, length(groups), nresult)
  for (i in seq_len(nresult)) {
    contributions[, i] <- if (identical(method, "constant"))
      colSums(pgs$spec[pgs$freq < freq.upper.bound[i] & pgs$freq >= freq.lower.bound[i],, drop = FALSE]) / norms
    else if (identical(method, "linear"))
      sapply(pgs$cumspecfuns, function(f) diff(f(c(freq.lower.bound[i], freq.upper.bound[i])))) / norms
  }
  
  type <- if (all(threshold <= 0)) "splitting" else "independent"
  
  if (identical(type, "splitting")) {
    gi <- max.col(contributions, ties.method = "first")
    result <- lapply(seq_len(nresult), function(i) groups[gi == i])
  } else if (identical(type, "independent")) {
    result <- lapply(seq_len(nresult), function(i) groups[contributions[, i] >= threshold[i]])
  } else {
    stop("Unknown type for pgrouping.auto.pgram")
  }
  
  names(result) <- if (!is.null(names(freq.bins))) Rssa:::.group.names(freq.bins) else Rssa:::.group.names(threshold)
  colnames(contributions) <- names(result)
  rownames(contributions) <- as.character(groups)
  
  if (drop) {
    result <- result[sapply(result, length) > 0]
  }
  
  attr(result, "contributions") <- contributions
  attr(result, "type") <- type
  attr(result, "threshold") <- threshold
  
  class(result) <- "grouping.auto.pgram"
  
  result
}

eossa_auto_trend <- function(ssa_obj,
                       signal_rank,
                       clust_type = c("hierarchial", "distance", "frequency"),
                       k = 2,
                       delta = 1e-3,
                       auto_trend_freq = 1 / 24,
                       auto_threshold = 0.7,
                       base='series'){
    eoss <- eossa_new(ssa_obj, 
                      nested.groups = list(1:signal_rank), 
                      clust_type = clust_type, 
                      k = k, 
                      delta = delta,
                      auto_trend_freq = auto_trend_freq)

    g <- grouping.auto.ssa.custom(eoss, base = base, groups = eoss$iossa.groups,
                       freq.bins = list(Trend = auto_trend_freq), 
                       threshold = auto_threshold)
  
  trend_comp <- g$Trend
  signal_indices <- 1:signal_rank
  res_comp <- signal_indices[!signal_indices %in% trend_comp]
  
  rec <- reconstruct(eoss, groups = list(Trend = trend_comp, Periodics = res_comp))
  return(list(model_object = eoss,
              grouping = list(Trend = trend_comp, Periodics = res_comp),
              trend = rec$Trend,
              periodics = rec$Periodics,
              residuals = attr(rec, "residuals")))
}

basic_auto_trend <- function(ssa_obj,
                       signal_rank,
                       auto_trend_freq = 1 / 24,
                       auto_threshold = 0.7){
  g <- grouping.auto(ssa_obj, base = "series", groups = list(1:signal_rank),
                     freq.bins = list(Trend = auto_trend_freq), 
                     threshold = auto_threshold)
  
  trend_comp <- g$Trend
  signal_indices <- 1:signal_rank
  res_comp <- signal_indices[!signal_indices %in% trend_comp]
  
  rec <- reconstruct(ssa_obj, groups = list(Trend = trend_comp, Periodics = res_comp))
  
  return(list(model_object = ssa_obj,
              grouping = list(Trend = trend_comp, Periodics = res_comp),
              trend = rec$Trend,
              periodics = rec$Periodics,
              residuals = attr(rec, "residuals")))
}

fossa_auto_trend <- function(ssa_obj,
                       signal_rank,
                       auto_trend_freq = 1 / 24,
                       auto_threshold = 0.7){
  foss <- fossa(ssa_obj, nested.groups = 1:signal_rank)
  g <- grouping.auto(foss, base = "series", groups = list(1:signal_rank),
                     freq.bins = list(Trend = auto_trend_freq), 
                     threshold = auto_threshold)
  
  
  trend_comp <- g$Trend
  signal_indices <- 1:signal_rank
  res_comp <- signal_indices[!signal_indices %in% trend_comp]
  
  rec <- reconstruct(foss, groups = list(Trend = trend_comp, Periodics = res_comp))
  return(list(model_object = foss,
              grouping = list(Trend = trend_comp, Periodics = res_comp),
              trend = rec$Trend,
              periodics = rec$Periodics,
              residuals = attr(rec, "residuals")))
}

auto_trend_model <- function(ssa_obj,
                             method = c("iossa.post_grouping", "iossa.manual", "iossa.auto_grouping.ssa",
                                        "iossa.auto_grouping.fossa", "eossa.auto", "basic.auto", "fossa.auto"),
                             signal_rank,
                             base='series',
                             ...){
  # Input: methods --- vector of names of ssa methods:
  #           "iossa.post_grouping"       --- IOSSA with separate initial grouping, using 
  #                                           automated trend identification after IOSSA. 
  #                                           Tolerance parameter is  provided in 
  #                                           parameters["iossa_tolerance"]. Maxiter  parameter
  #                                           is in parameters["iossa_maxiter"].
  #           "iossa.manual"              --- IOSSA with manual initial grouping, grouping is 
  #                                           provided in parameters["manual_grouping"]. 
  #                                           Tolerance parameter is provided in 
  #                                           parameters["iossa_tolerance"]. Maxiter parameter 
  #                                           is in parameters["iossa_maxiter"].
  #           "iossa.auto_grouping.ssa"   --- IOSSA with auto initial grouping, Tolerance 
  #                                           parameter is provided in 
  #                                           parameters["iossa_tolerance"]. Maxiter
  #                                           parameter is in parameters["iossa_maxiter"].
  #           "iossa.auto_grouping.fossa" --- same as iossa.auto_grouping.ssa but using
  #                                           FOSSA for initial grouping.
  #           "eossa.auto"                --- EOSSA with auto grouping.
  #           "basic.auto"                --- SSA with auto grouping.
  #           "fossa.auto"                --- FOSSA with auto grouping.
  #       parameters --- list with method parameters:
  #           "series_length"             --- time series length.
  #           "signal_rank"               --- rank of signal for nested decomposition.
  #           "window_length"             --- window length in SSA.
  #           "auto_trend_freq"           --- trend max frequency for auto grouping.
  #           "auto_threshold"            --- threshold for auto grouping.
  #           "iossa_maxiter"             --- max num iterations in IOSSA.
  #           "iossa_tolerance"           --- tolerance in IOSSA.
  #           "manual_grouping"           --- nested grouping, list of vectors, in IOSSA.
  # Output: function creates and returns an object of class auto_trend,
  #         object contains result of decomposition and ssa information.
  
  if (!prod(method %in% c("iossa.post_grouping", "iossa.manual", "iossa.auto_grouping.ssa",
                          "iossa.auto_grouping.fossa", "eossa.auto", "basic.auto", "fossa.auto"))){
    stop("Wrong method parameter.")
  }
  method <- c(method)
  
  kwargs <- list(...)
  
  this <- list(method = method)
  attr(this, ".env") <- new.env();
  
  if (length(method) > 1){
    for (m in method){
      this[[m]] <- list()
    }
  }
  
  results <- list()
  
  class(this) <- c("auto_trend")
  if ("iossa.post_grouping" %in% method){
    res <- iossa_post_grouping_trend(ssa_obj = ssa_obj, 
                               signal_rank = signal_rank, 
                               auto_trend_freq = kwargs$auto_trend_freq, 
                               auto_threshold = kwargs$auto_threshold,
                               iossa_maxiter = kwargs$iossa_maxiter, 
                               iossa_tolerance = kwargs$iossa_tolerance,
                               kappa = kwargs$kappa)
    results[["iossa.post_grouping"]] <- res
    if (length(method) == 1)
      this$iterations <- res$iterations
    else {
      this[["iossa.post_grouping"]]$iterations = res$iterations
    }
  }
  if ("iossa.manual" %in% method){
    res <- iossa_manual_trend(ssa_obj = ssa_obj, 
                        signal_rank = signal_rank, 
                        auto_trend_freq = auto_trend_freq, 
                        auto_threshold = auto_threshold,
                        iossa_maxiter = iossa_maxiter, 
                        iossa_tolerance = iossa_tolerance, 
                        manual_grouping = manual_grouping,
                        kappa = kappa)
    results[["iossa.manual"]] <- res
    if (length(method) == 1)
      this$iterations <- res$iterations
    else{
      this[["iossa.manual"]]$iterations = res$iterations
    }
  }
  if ("iossa.auto_grouping.ssa" %in% method){
    res <- iossa_auto_grouping_ssa_trend(ssa_obj = ssa_obj, 
                                   signal_rank = signal_rank, 
                                   auto_trend_freq = auto_trend_freq, 
                                   auto_threshold = auto_threshold,
                                   auto_initial_threshold = auto_initial_threshold, 
                                   iossa_tolerance = iossa_tolerance,
                                   kappa = kappa)
    results[["iossa.auto_grouping.ssa"]] <- res
    if (length(method) == 1)
      this$iterations <- res$iterations
    else{
      this[["iossa.auto_grouping.ssa"]]$iterations = res$iterations
    }
  }
  if ("iossa.auto_grouping.fossa" %in% method){
    res <- iossa_auto_grouping_fossa_trend(ssa_obj = ssa_obj, 
                                     signal_rank = signal_rank, 
                                     auto_trend_freq = kwargs$auto_trend_freq, 
                                     auto_threshold = kwargs$auto_threshold,
                                     iossa_tolerance = kwargs$iossa_tolerance,
                                     iossa_maxiter = kwargs$iossa_maxiter,
                                     kappa = kwargs$iossa_kappa)
    results[["iossa.auto_grouping.fossa"]] <- res
    if (length(method) == 1)
      this$iterations <- res$iterations
    else{
      this[["iossa.auto_grouping.fossa"]]$iterations = res$iterations
    }
  }
  if ("eossa.auto" %in% method){
    res <- eossa_auto_trend(ssa_obj, 
                      signal_rank = signal_rank, 
                      clust_type = kwargs$clust_type,
                      auto_trend_freq = kwargs$auto_trend_freq, 
                      auto_threshold = kwargs$auto_threshold,
                      delta = kwargs$delta,
                      base=base)
    results[["eossa.auto"]] <- res
  }
  if ("basic.auto" %in% method){
    res <- basic_auto_trend(ssa_obj = ssa_obj, 
                      signal_rank = signal_rank, 
                      auto_trend_freq = kwargs$auto_trend_freq, 
                      auto_threshold = kwargs$auto_threshold)
    results[["basic.auto"]] <- res
  }
  if ("fossa.auto" %in% method){
    res <- fossa_auto_trend(ssa_obj, 
                      signal_rank = signal_rank, 
                      auto_trend_freq = auto_trend_freq, 
                      auto_threshold = auto_threshold)
    results[["fossa.auto"]] <- res
  }
  
  this$original <- ssa_obj$F
  this$signal_rank <- signal_rank
  
  if (length(method) == 1){
    this$model_object <- res$model_object
    this$grouping <- res$grouping
    this$trend <- res$trend
    this$periodics <- res$periodics
    this$residuals <- res$residuals
  }
  else{
    for (m in method){
      this[[m]]$model_object <- results[[m]]$model_object
      this[[m]]$grouping <- results[[m]]$grouping
      this[[m]]$trend <- results[[m]]$trend
      this[[m]]$periodics <- results[[m]]$periodics
      this[[m]]$residuals <- results[[m]]$residuals
    }
  }
  
  this
}

pgram_1d <- function(x) {
  if (!is.matrix(x)) x <- as.matrix(x)
  stopifnot(all(is.finite(x)))
  
  X <- mvfft(x)
  n <- nrow(x)
  
  N <- n %/% 2 + 1
  spec <- abs(X[seq_len(N),, drop = FALSE])^2
  
  if (n %% 2 == 0) {
    if (N > 2) spec[2:(N-1), ] <- 2 * spec[2:(N-1), ]
  } else {
    if (N >= 2) spec[2:N, ] <- 2 * spec[2:N, ]
  }
  
  freq <- seq(0, 1, length.out = n + 1)[seq_len(N)]
  
  cumspecfuns <- lapply(seq_len(ncol(x)),
                        function(j)
                          approxfun(c(0, freq[-N] + 1/(2*n), 0.5),
                                    c(0, cumsum(spec[, j])),
                                    rule = 2))
  
  list(spec = spec, freq = freq, cumspecfuns = cumspecfuns)
}

periodic_grouping_freq <- function(ssa_obj, 
                                   groups,  
                                   s_0 = 1, 
                                   rho_0 = 0.9,
                                   eps = 1e-9,
                                   ...){
  L <- ssa_obj$window
  max_k <- length((0:(L %/% 2)) / L)
  
  if (missing(groups))
    groups <- as.list(1:min(nsigma(ssa_obj),n))
  
  groups <- sort(unique(unlist(groups)))
  groups2 <- groups[-length(groups)]
  
  n <- length(groups)
  
  Fs <- ssa_obj$U[, groups, drop = FALSE]
  # periodogram
  pgs <- pgram_1d(Fs)
  pgs$spec <- pgs$spec / L
  
  ### part one
  max_pgram <- apply(pgs$spec, 2, function(x) pgs$freq[which.max(x)])
  
  I_1 <- groups2[L * abs(max_pgram[1:(n-1)] - max_pgram[2:n]) - s_0 < eps]
  I_2 <- groups[L * abs(max_pgram - 0.5) - s_0 < eps]
  
  I_1_freqs <- numeric(0)
  
  if (length(groups2) == 0)
    I_1 <- numeric(0)
  
  ### part two
  if (length(I_1) != 0){
    rho_I_1 <- (pgs$spec[,I_1]  + pgs$spec[,(I_1 + 1)])/2
    if (length(I_1) == 1)
      rho_I_1 <- matrix(rho_I_1, ncol = 1)
    
    prefinal_mask <- apply(rho_I_1, 2, function(x) max(x[1:(max_k-1)] + x[2:max_k]) - rho_0 > eps)
    
    rho_values <- apply(rho_I_1, 2, function(x) max(x[1:(max_k-1)] + x[2:max_k]))
    rho_values <- rho_values[rho_values - rho_0 > eps]
    I_1_prefinal <- I_1[prefinal_mask]
    I_1_freqs <- I_1_freqs[prefinal_mask]
    
    filter_mask <- rep(FALSE, length(I_1_prefinal))
    sorted_indices <- order(rho_values, decreasing = TRUE)
    if (length(rho_values) > 0){
      for (j in seq_along(sorted_indices)) {
        left_taken <- FALSE
        if (sorted_indices[j] > 1){
          if ((I_1_prefinal[sorted_indices[j]] == I_1_prefinal[sorted_indices[j] - 1]) & 
              (filter_mask[sorted_indices[j] - 1] == TRUE))
            left_taken <- TRUE 
        }
        
        right_taken <- FALSE
        if (sorted_indices[j] < length(sorted_indices)){
          if ((I_1_prefinal[sorted_indices[j]] == I_1_prefinal[sorted_indices[j] + 1]) & 
              (filter_mask[sorted_indices[j] + 1] == TRUE))
            right_taken <- TRUE
        }
          
        if (!left_taken & !right_taken)
          filter_mask[sorted_indices[j]] <- TRUE
      }
    
      I_1_prefinal <- I_1_prefinal[c(filter_mask, FALSE)]
      rho_I_1 <- (pgs$spec[,I_1_prefinal]  + pgs$spec[,(I_1_prefinal + 1)])/2
      if (length(I_1_prefinal) == 1)
        rho_I_1 <- matrix(rho_I_1, ncol = 1)
      
      argmax_rho <- apply(rho_I_1, 2, function(x) which.max(x[1:(max_k-1)] + x[2:max_k]))
    
    } else{
      argmax_rho <- numeric(0)
    }
    
    I_1_freqs <- numeric(length(argmax_rho))
    if (length(argmax_rho) > 0)
      for (i in 1:length(argmax_rho)){
        left_freq_rho <- rho_I_1[argmax_rho[i]]
        right_freq_rho <- rho_I_1[argmax_rho[i] + 1]
        if (left_freq_rho / (left_freq_rho + right_freq_rho) < 0.7){
          if (right_freq_rho / (left_freq_rho + right_freq_rho) < 0.7){
            I_1_freqs[i] <- (2 * argmax_rho[i] - 1) / 2 / L
          } else {
            I_1_freqs[i] <- argmax_rho[i] / L
          }
        } else {
          I_1_freqs[i] <- (argmax_rho[i] - 1) / L
        }
      }
    
    if (length(I_1_prefinal) != 0) {
      I_1_final <- sort(c(I_1_prefinal, I_1_prefinal+1))
    } else {
      I_1_final <- numeric(0)
    }
  } else {
    I_1_final <- numeric(0)
  }
  
  if (length(I_2) != 0) {
    rho_I_2 <- pgs$spec[,I_2] 
    if (length(I_2) == 1)
      rho_I_2 <- matrix(rho_I_2, ncol = 1)
    r1 <- rho_I_2[which(pgs$freq == (L %/% 2 - 1) / L),]
    r2 <- rho_I_2[which(pgs$freq == L %/% 2 / L),]
    
    I_2_final <- I_2[r1 + r2 - rho_0 > eps]
  } else {
    I_2_final <- numeric(0)
  }
  
  list(I_1 = I_1_final, 
       I_1_freqs = I_1_freqs,
       I_2 = I_2_final)
}

angle.fun <- function(P, Q) {
  angle <- function(P1, P2, Q1, Q2) {
    t <- acos((P1 * P2 + Q1 * Q2) / sqrt(
      P1 ^ 2 + Q1 ^ 2) / sqrt(P2 ^ 2 + Q2 ^ 2))
    t[!is.nan(t)]
  }
  angles <- angle(P1 = P[-length(P)],
                  P2 = P[-1],
                  Q1 = Q[-length(Q)],
                  Q2 = Q[-1])
  
  vv <- var(angles)
  mm <- mean(angles)
  tau_value <- vv / min(1, mm ^ 2)
  return(list(tau_value = tau_value,
              mean_angle = mm))
}

periodic_grouping_angle_reg <- function(ssa_obj, 
                                        groups,
                                        threshold = 0.01, 
                                        eps = 1e-9,
                                        s_0 = 1,
                                        rho_0 = 0.9,
                                        p_0 = 0.05,
                                        ...) {
  groups <- sort(unique(unlist(groups)))
  
  if (missing(groups)) {
    groups <- 1:nu(ssa_obj)
  }

  # Identification of periodics with frequency 0.5
  L <- ssa_obj$window
  Fs <- ssa_obj$U[, groups, drop = FALSE]
  if (length(groups) == 1){
    Fs <- matrix(Fs, ncol = 1)
    sums <- sum(abs(sign(Fs[1:(L - 1),]) + sign(Fs[2:L,])) / 2) 
  } else {
    sums <- colSums(abs(sign(Fs[1:(L - 1),]) + sign(Fs[2:L,])) / 2) 
  }
  one_el_gamma <- sums / (L - 1)
  one_el_gamma <- one_el_gamma[one_el_gamma < p_0]
  I_2_final <- groups[sums / (L - 1) < p_0]
  
  groups <- setdiff(groups, I_2_final)
  
  if (length(groups) >= 2){
    filter_mask <- rep(FALSE, length(groups) - 1)
    tau <- numeric(length(groups) - 1)
    I_1_freqs <- numeric(length(groups) - 1)
    for (j in seq_along(groups)[-length(groups)]) {
      res <- angle.fun(P = ssa_obj$U[, groups[j]], 
                       Q = ssa_obj$U[, groups[j + 1]])
      tau[j] <- res$tau_value
      I_1_freqs[j] <- res$mean_angle / (2 * pi)
    }
    sorted_indices <- order(tau)
    for (j in seq_along(sorted_indices)) {
      left_taken <- FALSE
      if (sorted_indices[j] > 1){
        if ((groups[sorted_indices[j]] == groups[sorted_indices[j] - 1]) &
            (filter_mask[sorted_indices[j] - 1] == TRUE))
          left_taken <- TRUE
      }
      
      right_taken <- FALSE
      if (sorted_indices[j] < length(sorted_indices)){
        if ((groups[sorted_indices[j]] == groups[sorted_indices[j] + 1]) &
            (filter_mask[sorted_indices[j] + 1] == TRUE))
          right_taken <- TRUE
      }
      
      if (!left_taken & !right_taken & (tau[sorted_indices[j]] - threshold < eps))
        filter_mask[sorted_indices[j]] <- TRUE
    }
    
    tau_raw <- tau
    two_el_tau <- tau[c(filter_mask, FALSE)]
    I_1_final <- groups[c(filter_mask, FALSE)]
    I_1_final <- sort(c(I_1_final, I_1_final + 1))
    I_1_freqs <- I_1_freqs[filter_mask]
  }else{
    I_1_final <- numeric(0)
    I_1_freqs <- numeric(0)
    two_el_tau <- numeric(0)
    tau_raw <- numeric(0)
  }
  
  list(I_1 = I_1_final,
       two_el_tau=two_el_tau,
       I_1_freqs = I_1_freqs,
       I_2 = I_2_final,
       tau_raw=tau_raw,
       one_el_gamma=one_el_gamma)
}

auto_periodic_model <- function(ssa_obj,
                                groups,
                                method = "frequency",
                                ...) {
  # Input:
  #       ssa_obj --- ssa object.
  #       groups  --- list with groups for periodic identification.
  #       method --- string with name of periodic identification method. Must be one of:
  #           "frequency"                 --- frequency method (default).
  #           "angle_reg"                 --- angle regularity method.
  #     *extra parameters*:
  #       s0                              --- parameter for frequency method 
  #                                           (frequencies proximity), default = 1
  #       rho0                            --- threshold for rho measure, default = 0.9
  
  allowed_methods <- c("frequency", "angle_reg")
  
  if (!(method %in% allowed_methods))
    stop(paste0("Wrong method parameter, must be one of: ", allowed_methods))
  
  this <- list(method = method)
  attr(this, ".env") <- new.env();
  class(this) <- c("auto_periodic_model")
  
  if (method == "frequency"){
    res <- periodic_grouping_freq(ssa_obj = ssa_obj,
                                  groups = groups,
                                  ...)
  }
  if (method == "angle_reg"){
    res <- periodic_grouping_angle_reg(ssa_obj = ssa_obj,
                                       groups = groups,
                                       ...)
    this$two_el_tau <- res$two_el_tau
    this$one_el_gamma <- res$one_el_gamma
  }
  
  this$series <- list()
  I_1_harmonics_num <- length(res$I_1) %/% 2
  
  index <- 1
  while (index <= I_1_harmonics_num){
    series_name <- paste0("P", index)
    this$series[[series_name]] <- list(series = reconstruct(ssa_obj, 
                                                            groups = list(c(res$I_1[index], 
                                                                            res$I_1[index + 1])))$F1,
                                       freq = res$I_1_freqs)
    index <- index + 1
  }
  
  index <- 1
  while (index <= length(res$I_2)){
    series_name <- paste0("P", index + I_1_harmonics_num)
    this$series[[series_name]] <- list(series = reconstruct(ssa_obj, 
                                                            groups = list(c(res$I_2[index])))$F1,
                                       freq = 0.5)
    index <- index + 1
  }
  
  
  this$two_el_indices <- res$I_1
  this$two_el_freqs <- res$I_1_freqs
  this$one_el_indices <- res$I_2
  this$undetermined_indices <- setdiff(unique(groups), 
                                       c(res$I_1, res$I_2))
  if (length(this$undetermined_indices) > 0){
    this[["undetermined_series"]] <- reconstruct(ssa_obj, 
                                              groups = list(this$undetermined_indices))$F1
  } else {
    this[["undetermined_series"]] <- numeric(0)
  }
  
  this
}

auto_ssa <- function(time_series,
                     L = length(time_series) %/% 2,
                     ssa.method = "eossa.auto",
                     rank_est.method = "svd",
                     periodic_id.method = "frequency",
                     ...
                     ){
  # Input:
  #       time_series --- real-value vector
  #       L           --- window length, default is length(time_series) %/% 2
  #       ssa.method --- string with name of ssa method:
  #           "iossa.auto_grouping.fossa" --- same as iossa.auto_grouping.ssa but using
  #                                           FOSSA for initial grouping.
  #           "eossa.auto"                --- EOSSA with auto grouping (default).
  #           "fossa.auto"                --- FOSSA with auto grouping.
  #       rank_est.method --- string with name of rank estimation method:
  #           "svd"                       --- svd method (default).
  #           "mgn"                       --- mgn method.
  #       periodic_id.method --- string with name of periodic identification method:
  #           "frequency"                 --- frequency method (default).
  #           "angle_reg"                 --- angle regularity method.
  #     *extra parameters*:
  #           "signal_rank"               --- rank of signal for nested decomposition.
  #                                           If not provided, rank estimation method will be used.
  #           "auto_trend_freq"           --- trend max frequency for auto grouping.
  #           "auto_trend_threshold"      --- threshold for trend auto grouping.
  #           "iossa_maxiter"             --- max num iterations in IOSSA.
  #           "iossa_tolerance"           --- tolerance in IOSSA.
  #           "eossa_clust_type"          --- clustering type for EOSSA, string. Must be one of:
  #                 "hierarchial"               --- hierarchial clustering.
  #                 "distance"                  --- distance clustering.
  #                 "frequency"                 --- frequency clustering.
  #           "max_rank"                  --- maximum rank for estimation.
  # Output: function creates and returns an object of class auto_ssa,
  #         object contains result of decomposition and ssa information.
  
  if (missing(ssa.method))
    ssa.method <- "eossa.auto"
  
  allowed_methods <- c("fossa.auto", 
                       "eossa.auto", "iossa.auto_grouping.fossa")
  
  if (!(ssa.method %in% allowed_methods))
    stop(paste0("Wrong ssa.method parameter. Must be one of ", 
                allowed_methods))
  
  this <- list(method = method,
               trend_model <- NA,
               periodic_model <- NA)
  attr(this, ".env") <- new.env();
  class(this) <- c("auto_ssa")
  
  #Rank estimation
  estimated_rank <- estimated_rank(time_series = time_series,
                                   L = L,
                                   method = rank_est.method,
                                   ...)
  #SSA object
  ssa_obj <- ssa(x = time_series, 
                 L = L,
                 ...)
  
  #Separability improving and trend extraction
  this$trend_model <- auto_trend_model(ssa_obj = ssa_obj,
                                       method = ssa.method,
                                       signal_rank = estimated_rank,
                                       ...)
  
  #Periodics extraction
  this$periodics_model <- auto_periodic_model(
    ssa_obj = this$trend_model$model_object,
    groups = this$trend_model$grouping$Periodics,
    method = periodic_id.method,
    ...)
  
  #Result components
  this$components <- list()
  this$components$Trend <- this$trend_model$trend
  this$components$Residuals <- this$trend_model$residuals
  this$components$Undetermined <- this$periodic_model$series$undetermined
  
  #for (i in 1:length)
}

get_weight <- function(index, N, window_length){
  # Input: index          --- index of series element
  #        N              --- series length
  #        window_length  --- length of window in SSA
  # Result: returns weight in weighted correlation
  K <- N - window_length + 1
  L_star <- min(window_length, K)
  K_star <- max(window_length, K)
  return(ifelse((index <= L_star), index, ifelse(index <= K_star, L_star, N - index + 1)))
}

wcor_series <- function(comp1, comp2, window_length) {
  # Input: comp1         --- series
  #        comp2         --- series
  #        window_length --- length of window in SSA
  if (length(comp1) != length(comp2))
    stop("Lengths of series are not equal.")
  N <- length(comp1)
  wcor_sum <- 0
  fr_norm1 <- 0
  fr_norm2 <- 0
  for (i in 1:N){
    wcor_sum <- wcor_sum + get_weight(i, N, window_length) * comp1[i] * comp2[i]
    fr_norm1 <- fr_norm1 + get_weight(i, N, window_length) * comp1[i] * comp1[i]
    fr_norm2 <- fr_norm2 + get_weight(i, N, window_length) * comp2[i] * comp2[i]
  }
  fr_norm1 <- sqrt(fr_norm1)
  fr_norm2 <- sqrt(fr_norm2)
  return(wcor_sum / fr_norm1 / fr_norm2)
}

plot.auto_trend <- function(x, type = c("series", "vectors", "values", "elementary"),
                            filename = NA,
                            file_width = 1000,
                            file_height = 600,
                            size = 0.7,
                            compare = FALSE,
                            actual_trend = numeric(0),
                            compare_only = FALSE,
                            plot.contrib = FALSE,
                            original_series = numeric(0),
                            color_list,
                            methods_map,
                            indices = NA
){
  if (is.na(indices))
    indices <- 1:x$signal_rank
  
  type = match.arg(type)
  if (type == "series"){
    if (compare){
      if (length(actual_trend) == 0)
        stop("Actual trend must be passed to plot.")
      if (length(x$method) == 1){
        if (!is.na(filename))
          pdf(file = filename, width = file_width, height = file_height)
        
        plot_df <- data.frame("Reconstructed trend" = x$trend, 
                              "Index" = 0:(length(x$original)-1), 
                              "Actual trend" = actual_trend)
        
        if (length(original_series) > 0)
          plot_df[,"Original series"] <- original_series
        
        p <- ggplot(plot_df, aes(Index)) + 
          geom_line(aes(y = Actual.trend, colour = "Actual trend"), size = size) + 
          geom_line(aes(y = Reconstructed.trend, colour = "Reconstructed trend"), size = size)
        
        if (length(original_series) > 0)
          p <- p  + geom_line(aes(y = Actual.trend, colour = "Original series"), size = size)
        
        p <- p + 
          labs(x = "Index", y = "Value", color = "Trends") + 
          scale_color_manual(values=c("#07973D", "#EE3D3D"))
        
        print(p + theme_classic())
        
        if (!is.na(filename))
          dev.off()
      }
      else{
        index <- 1
        if (!compare_only){
          for (m in x$method){
            if (!is.na(filename)){
              if (length(c(filename)) == 1)
                png(filename = paste0(filename, '.', m, '.png'), width = file_width, height = file_height)
              else
                png(filename = filename[index], width = file_width, height = file_height)
            }
            
            plot_df <- data.frame("Reconstructed trend" = x[[m]]$trend, 
                                  "Index" = 0:(length(x$original)-1), 
                                  "Actual trend" = actual_trend)
            
            p <- ggplot(plot_df, aes(Index)) + 
              geom_line(aes(y = Actual.trend, colour = "Actual trend"), size = size) + 
              geom_line(aes(y = Reconstructed.trend, colour = "Reconstructed trend"), size = size) + 
              labs(x = "Index", y = "Value", color = "Series") + 
              scale_color_manual(values=c("#EE3D3D", "#0805B5"))
            
            print(p)
            
            if (!is.na(filename))
              dev.off()
            
            index <- index + 1
          }
        }
        
        if (!is.na(filename)){
          if (length(c(filename)) == 1)
            pdf(file = paste0(filename, '.', 'compare', '.pdf'), width = file_width, height = file_height)
          else
            pdf(file = filename[index], width = file_width, height = file_height)
        }
        
        plot_df <- data.frame("Index" = 0:(length(x$original)-1), 
                              "Actual trend" = actual_trend)
        
        if (length(original_series) > 0)
          plot_df$Original.series <- original_series
        for (m in x$method){
          plot_df[m] <- x[[m]]$trend
        }
        
        index <- 1
        p <- ggplot(plot_df, aes(Index)) + 
          geom_line(aes_string(y = "Actual.trend", colour = as.factor(index)), size = size)
        index <- 2
        colors <- c(color_list$actual.trend)
        
        if (length(original_series) > 0) {
          p <- p + geom_line(aes_string(y = "Original.series", colour = as.factor(index)), size = size)
          index <- 3
        }
        labels <- c("Actual trend")
        for (m in x$method){
          p <- p + geom_line(aes_string(y = m, colour = as.factor(index)), size = size)
          colors <- c(colors, color_list[[m]])
          labels <- c(labels, methods_map[[m]])
          index <- index + 1
        }
        
        p <- p + scale_colour_manual(name = "Trends", labels = labels, values = colors) +
          labs(title = "", x = "Time", y = "Sales")
        
        print(p + theme_classic())
        
        if (!is.na(filename))
          dev.off()
      }
    }
    else{
      if (length(x$method) == 1){
        if (!is.na(filename))
          pdf(file = filename, width = file_width, height = file_height)
        
        plot_df <- data.frame("Original series" = x$original, "Reconstructed trend" = x$trend, 
                              "Reconstructed periodics" = x$periodics, "Index" = 0:(length(x$original)-1), 
                              "Reconstructed residuals" = x$residuals)
        
        p <- ggplot(plot_df, aes(Index)) + 
          geom_line(aes(y = Original.series, colour = "Original series"), size = size) + 
          geom_line(aes(y = Reconstructed.trend, colour = "Reconstructed trend"), size = size) + 
          geom_line(aes(y = Reconstructed.periodics, colour = "Reconstructed periodics"), size = size) + 
          geom_line(aes(y = Reconstructed.residuals, colour = "Reconstructed residuals"), size = size) + 
          labs(title = "", x = "Index", y = "Value", color = "Series") + 
          scale_color_manual(values=c(color_list[["original.series"]], 
                                      color_list[["reconstructed.periodics"]], 
                                      color_list[["reconstructed.residuals"]], 
                                      color_list[[x$method]]))
        
        print(p + theme_classic())
        
        if (!is.na(filename))
          dev.off()
      }
      else{
        index <- 1
        if (!compare_only){
          for (m in x$method){
            if (!is.na(filename)){
              if (length(filename) == 1)
                png(filename = paste0(filename, '.', m, '.png'), width = file_width, height = file_height)
              else
                png(filename = filename[index], width = file_width, height = file_height)
            }
            
            plot_df <- data.frame("Original series" = x$original, "Trend" = x[[m]]$trend, 
                                  "Periodics" = x[[m]]$periodics, "Index" = 0:(length(x$original)-1), 
                                  "Residuals" = x[[m]]$residuals)
            
            p <- ggplot(plot_df, aes(Index)) + 
              geom_line(aes(y = Original.series, colour = "Original series"), size = size) + 
              geom_line(aes(y = Trend, colour = "Trend"), size = size) + 
              geom_line(aes(y = Periodics, colour = "Periodics"), size = size) + 
              geom_line(aes(y = Residuals, colour = "Residuals"), size = size) + 
              labs(title = paste0(m, ' components'), x = "Index", y = "Value", color = "Series") + 
              scale_color_manual(values=c("#52FA5E", "#EE3D3D", "#645502", "#0805B5"))
            
            print(p)
            
            if (!is.na(filename))
              dev.off()
            
            index <- index + 1
          } 
        }
        
        if (!is.na(filename)){
          if (length(c(filename)) == 1)
            pdf(file = paste0(filename, '.', 'compare', '.pdf'), width = file_width, height = file_height)
          else
            pdf(file = filename[index], width = file_width, height = file_height)
        }
        
        plot_df <- data.frame("Index" = seq(from = as.Date('1980/1/1'), to = as.Date('1994/6/1'), by= "month"))
        
        if (length(original_series) > 0)
          plot_df$Original.series <- ts(as.vector(original_series), start = 1980, frequency = 12)
        
        for (m in x$method){
          plot_df[m] <- ts(as.vector(x[[m]]$trend), start = 1980, frequency = 12)
        }
        
        p <- ggplot(plot_df, aes(Index))
        
        index <- 0
        colors <- numeric(0)
        labels <- numeric(0)
        if (length(original_series) > 0) {
          p <- p + geom_line(aes_string(y = "Original.series", colour = as.factor(index)), size = size * 0.5)
          labels <- c("Original series")
          colors <- c(color_list$original.series)
        }
        
        index <- 1
        for (m in x$method){
          p <- p + geom_line(aes_string(y = m, colour = as.factor(index)), size = size)
          colors <- c(colors, color_list[[m]])
          labels <- c(labels, methods_map[[m]])
          index <- index + 1
        }
        
        
        p <- p + scale_colour_manual(labels = labels, values = colors) +
          labs(x = "Time", y = "Sales", color = "")
        
        print(p + theme_classic())
        
        if (!is.na(filename))
          dev.off()
      }
    }
  }
  if (type == "vectors"){
    if (length(x$method) == 1){
      if (!is.na(filename))
        png(filename = filename, width = file_width, height = file_height)
      
      print(plot(x$model_object, type = "vectors", idx = 1:x$signal_rank, main = paste0("Eigenvectors, ", x$method)))
      
      if (!is.na(filename))
        dev.off() 
    }
    else {
      index <- 1
      for (m in x$method){
        if (!is.na(filename)){
          if (length(filename) == 1)
            png(filename = paste0(filename,'.', m, '.png'), width = file_width, height = file_height)
          else
            png(filename = filename[index], width = file_width, height = file_height)
        }
        
        print(plot(x[[m]]$model_object, type = "vectors", idx = 1:x$signal_rank, main = paste0("Eigenvectors, ", m)))
        
        if (!is.na(filename))
          dev.off() 
        index <- index + 1
      }
    }
  }
  if (type == "values"){
    if (length(x$method) == 1){
      if (!is.na(filename))
        png(filename = filename, width = file_width, height = file_height)
      
      print(plot(x$model_object, type = "values", main = paste0("Component norms, ", x$method)))
      
      if (!is.na(filename))
        dev.off() 
    }
    else {
      index <- 1
      for (m in x$method){
        if (!is.na(filename)){
          if (length(filename) == 1)
            png(filename = paste0(filename,'.', m, '.png'), width = file_width, height = file_height)
          else
            png(filename = filename[index], width = file_width, height = file_height)
        }
        
        print(plot(x[[m]]$model_object, type = "values", main = paste0("Component norms, ", m)))
        
        if (!is.na(filename))
          dev.off() 
        index <- index + 1
      }
    }
  }
  if (type == "elementary"){
    if (length(x$method) == 1){
      if (!is.na(filename))
        pdf(file = filename, width = file_width, height = file_height)
      
      print(plot(x$model_object, type = "series", groups = indices, 
                 main = "",
                 plot.contrib = plot.contrib))
      
      if (!is.na(filename))
        dev.off() 
    }
    else {
      index <- 1
      for (m in x$method){
        if (!is.na(filename)){
          if (length(filename) == 1)
            png(filename = paste0(filename,'.', m, '.png'), width = file_width, height = file_height)
          else
            png(filename = filename[index], width = file_width, height = file_height)
        }
        
        print(plot(x[[m]]$model_object, type = "series", groups = 1:x$signal_rank, 
                   main = paste0("Reconstructed series, ", m),
                   plot.contrib = plot.contrib))
        
        if (!is.na(filename))
          dev.off() 
        index <- index + 1
      }
    }
    
  }
}