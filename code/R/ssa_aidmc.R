library("here")

source(here("R", "mc-mssa.R"), local = TRUE)
source(here("R", "scripts.R"), local = TRUE)

grouping.auto.freqs <- function(x,
                                groups,
                                freqs,
                                C_max = 5,
                                threshold = 0.5) {
  N <- x$length
  n_triples <- sapply(freqs, function(freq)
    ifelse(abs(freq - 0.5) < .Machine$double.eps^0.5, 1, 2))
  deltas <- sapply(freqs, function(freq)
    calculate_delta(N, freq, C_max))
  freq.bins <- lapply(seq_along(freqs), function(i)
    c(freqs[i] - deltas[i], freqs[i] + deltas[i] + 1e-9))
  g <- grouping.auto.pgram(x, groups, freq.bins = freq.bins, threshold = threshold)
  
  components <- integer()
  for (j in seq_along(freq.bins)) {
    gg <- setdiff(g[[paste0("F", j)]], components)
    len <- min(length(gg), n_triples[j])
    components <- c(components, gg[seq_len(len)])
  }
  components
}

auto_mcssa <- function(x,
                       ssa_obj,
                       L,
                       delta = NULL,
                       C_max = 5,
                       groups = 1:30,
                       conf.level = 0.95,
                       auto_trend_freq = 1 / 24,
                       threshold = 0.5,
                       method = c("adaptive", "semi-adaptive", "fixed"),
                       maxit = 10,
                       detrend = TRUE,
                       ...) {
  method <- match.arg(method)
  
  N <- length(x)
  
  args <- list(...)
  if (!is.null(args$fixed))
    fixed <- args$fixed
  else
    fixed <- list(phi = NA, d = NA)
  
  trend_components <- integer()
  trend_est <- rep(0, N)
  noise_est <- x
  if (detrend) {
    trend_model <- auto_trend_model(
      ssa_obj,
      method = "eossa.auto",
      signal_rank = length(groups),
      clust_type = "distance",
      auto_trend_freq = auto_trend_freq,
      auto_threshold = 0.5,
      delta = 0.001
    )
    trend_est <- trend_model$trend
    trend_components <- trend_model$grouping$Trend
    noise_est <- x - trend_est
    
    groups <- seq_len(length(groups) - length(trend_components))
    ssa_obj <- ssa(
      x = noise_est,
      L = ssa_obj$window,
      svd.method = "svd",
      neig = length(groups)
    )
  }
  
  m <- mcssa(
    x = noise_est,
    L = L,
    basis = "cos",
    fixed = fixed,
    conf.level = conf.level,
    freq.range = c(auto_trend_freq, 0.5)
  )
  
  freq.exclude <- list(c(0, auto_trend_freq))
  
  freq <- numeric()
  freq.bins <- list()
  n_triples <- integer()
  proj_vectors <- list()
  i <- 1
  
  groups_current <- groups
  periodics_components <- integer()
  periodics_est <- rep(0, N)
  while (m$reject && i <= maxit) {
    freq[i] <- get_signif_freq(m, n_periodics = 1)
    delta <- calculate_delta(N, freq[i], C_max)
    freq.bins[[i]] <- c(freq[i] - delta - 1e-9, freq[i] + delta + 1e-9)
    
    n_triples[i] <- ifelse(freq[i] == 0.5, 1, 2)
    g <- grouping.auto.pgram(ssa_obj, groups_current, "series", freq.bins[i], threshold)
    if (length(g)) {
      g <- g[[1]]
      len <- min(length(g), n_triples[i])
      components <- g[seq_len(len)]
      groups_current <- setdiff(groups_current, components)
    } else {
      components <- integer()
      freq.exclude[[length(freq.exclude) + 1]] <- freq.bins[[i]]
    }
    
    periodic <- reconstruct(ssa_obj, list(components))[[1]]
    
    noise_est <- noise_est - periodic
    periodics_est <- periodics_est + periodic
    periodics_components <- c(periodics_components, components)
    
    model <- as.list(arfima_whittle(noise_est, c(fixed$phi, fixed$d), freq.exclude))
    m <- mcssa(
      x = noise_est,
      L = L,
      basis = "cos",
      model = model,
      conf.level = conf.level,
      freq.range = NULL,
      freq.exclude = freq.exclude
    )
    
    i <- i + 1
  }
  
  list(
    trend = trend_est,
    trend_components = trend_components,
    periodics = periodics_est,
    periodics_components = periodics_components,
    freq.bins = freq.bins,
    freqs = freq
  )
}


ssa_aidmc <- function(x,
                      L1,
                      L2 = (N + 1) %/% 2,
                      delta = NULL,
                      C_max = 5,
                      rank = 30,
                      conf.level = 0.95,
                      auto_trend_freq = 1 / 24,
                      tau_threshold = 0.05,
                      threshold = 0.5,
                      method = c("adaptive", "semi-adaptive", "fixed"),
                      maxit = 10,
                      trace = FALSE,
                      ...) {
  args <- list(...)
  if (!is.null(args$fixed))
    fixed <- args$fixed
  else
    fixed <- list(phi = NA, d = NA)
  
  N <- length(x)
  
  dec <- ssa_aid(
    x,
    L = L2,
    signal_rank = rank,
    conf.level = conf.level,
    auto_trend_freq = auto_trend_freq,
    tau_threshold = tau_threshold,
    mcssa.arguments = list(L = L1, fixed = fixed),
    trace = trace
  )
  
  x_resid <- dec$residuals
  
  groups <- setdiff(1:rank, seq_len(dec$signal_rank))
  
  model_obj <- dec$model_obj
  if (dec$signal_rank == 0) {
    model_obj <- ssa(x, L2, svd.method = "svd", neig = rank)
  }
  
  m <- auto_mcssa(
    x_resid,
    model_obj,
    L1,
    delta,
    C_max,
    groups,
    conf.level,
    auto_trend_freq,
    threshold,
    method,
    maxit,
    detrend = FALSE,
    ...
  )
  
  j_max <- dec$signal_rank
  if (length(m$periodics_components)) {
    wcors <- suppressWarnings(wcor(model_obj, groups = seq_len(rank))[m$periodics_components, , drop = FALSE])
    for (i in seq_along(m$periodics_components)) {
      w <- wcors[i, ]
      j_max <- max(j_max, which(w >= 0.9) |> tail())
    }
  } else {
    return(dec)
  }
  
  dec2 <- ssa_aid(
    x,
    L = L2,
    signal_rank = j_max,
    conf.level = conf.level,
    auto_trend_freq = auto_trend_freq,
    tau_threshold = tau_threshold,
    mcssa.arguments = list(L = L1, fixed = fixed),
    nstages = 1,
    trace = trace
  )
  dec2
}