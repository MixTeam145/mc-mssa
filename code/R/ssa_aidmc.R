library("here")

source(here("R", "mc-mssa.R"), local = TRUE)
source(here("R", "scripts.R"), local = TRUE)

get_signif_freq <- function(x, n_periodics) {
  significance <- x$statistic$contribution - x$predint$upper
  idx <- which(significance > 0)
  
  if (missing(n_periodics))
    n_periodics <- length(idx)
  
  n_periodics <- min(n_periodics, length(idx))
  
  if (n_periodics == 1) {
    max_idx <- idx[which.max(significance[idx])]
    
    neighbors <- c()
    if ((max_idx - 1) %in% idx)
      neighbors <- c(neighbors, max_idx - 1)
    if ((max_idx + 1) %in% idx)
      neighbors <- c(neighbors, max_idx + 1)
    
    idx <- c(max_idx, neighbors)
    freq <- weighted.mean(x$statistic$freq[idx], significance[idx])
  } else {
    o <- order(-significance[idx])[1:n_periodics]
    idx <- idx[o]
    freq <- x$statistic$freq[idx]
  }
  
  freq
}

calculate_delta <- function(N,
                            omega0,
                            C,
                            eps = 5e-2,
                            threshold = 0.9) {
  omega <- seq(0, 0.5, 1 / N)
  
  if (abs(omega0 * N - round(omega0 * N)) < .Machine$double.eps^0.5)
    delta0 <- 0
  else {
    o <- order(abs(omega[omega <= 0.5] - omega0))
    delta0 <- abs(omega0 - omega[o[2]])
  }
  
  deltas <- delta0 + omega
  
  x <- em_harmonic(N, omega0, C)
  
  spec <- abs(fft(x)[seq_len(N %/% 2 + 1)])^2
  spec <- spec / sum(spec)
  S <- numeric(length(deltas))
  
  for (i in seq_along(deltas)) {
    S[i] <- sum(spec[abs(omega - omega0) <= deltas[i] + .Machine$double.eps^0.5])
    if (abs(S[i] - 1) <= .Machine$double.eps^0.5)
      break
  }
  idx <- which.max(diff(S) < eps & S[-1] > threshold) + 1
  
  deltas[idx]
}

auto_mcssa <- function(x,
                       ssa_obj,
                       L,
                       groups = 1:30,
                       conf.level = 0.95,
                       auto_trend_freq = 1 / 24,
                       auto_trend_threshold = 0.5,
                       auto_periodic_threshold = 0.5,
                       C_max = 5,
                       method = c("adaptive", "semi-adaptive", "fixed"),
                       maxit = 10,
                       detrend = TRUE,
                       fixed = list(phi = 0, d = 0),
                       ...) {
  method <- match.arg(method)
  
  N <- length(x)
  
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
      ...
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
    g <- grouping.auto.pgram(ssa_obj, groups_current, "series", freq.bins[i], auto_periodic_threshold)
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
                      L = (N + 1) %/% 2,
                      rank = 30,
                      conf.level = 0.95,
                      auto_trend_freq = 1 / 24,
                      auto_trend_threshold = 0.5,
                      tau_threshold = 0.05,
                      p_0 = 0.03,
                      auto_periodic_threshold = 0.5,
                      C_max = 5,
                      method = c("adaptive", "semi-adaptive", "fixed"),
                      maxit = 10,
                      trace = FALSE,
                      mcssa.arguments = list(L = L, fixed = list(phi = 0, d = 0)),
                      ...) {
  N <- length(x)
  
  dec <- ssa_aid(
    x,
    L = L,
    signal_rank = rank,
    conf.level = conf.level,
    auto_trend_freq = auto_trend_freq,
    auto_trend_threshold = auto_trend_threshold,
    tau_threshold = tau_threshold,
    p_0 = p_0,
    auto_periodic_threshold = auto_periodic_threshold,
    C_max = C_max,
    mcssa.arguments = mcssa.arguments,
    trace = trace,
    ...
  )
  
  x_resid <- dec$residuals
  
  groups <- setdiff(1:rank, seq_len(dec$signal_rank))
  
  model_obj <- dec$model_obj
  if (dec$signal_rank == 0) {
    model_obj <- ssa(x, L, svd.method = "svd", neig = rank)
  }
  
  m <- auto_mcssa(
    x_resid,
    ssa_obj = model_obj,
    L = mcssa.arguments$L,
    groups = groups,
    conf.level = conf.level,
    auto_trend_freq = auto_trend_freq,
    t_threshold = t_threshold,
    C_max = C_max,
    method = method,
    maxit = maxit,
    detrend = FALSE,
    fixed = mcssa.arguments$fixed,
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
    L = L,
    signal_rank = j_max,
    conf.level = conf.level,
    auto_trend_freq = auto_trend_freq,
    auto_trend_threshold = auto_trend_threshold,
    tau_threshold = tau_threshold,
    p_0 = p_0,
    auto_periodic_threshold = auto_periodic_threshold,
    C_max = C_max,
    mcssa.arguments = mcssa.arguments,
    nstages = 1,
    trace = trace
  )
  dec2
}