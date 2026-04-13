library("here")

source(here("R", "mc-mssa.R"), local = TRUE)
source(here("R", "scripts.R"), local = TRUE)


get_signif_freq <- function(x, n_periodics, weighted = TRUE) {
  significance <- x$statistic$contribution - x$predint$upper
  idx <- which(significance > 0)
  
  if (missing(n_periodics))
    n_periodics <- length(idx)
  
  n_periodics <- min(n_periodics, length(idx))
  
  if (n_periodics == 1 && weighted) {
    max_idx <- idx[which.max(significance[idx])]
    
    neighbors <- c()
    if ((max_idx - 1) %in% idx)
      neighbors <- c(neighbors, max_idx - 1)
    if ((max_idx + 1) %in% idx)
      neighbors <- c(neighbors, max_idx + 1)
    
    idx <- c(max_idx, neighbors)
    
    if (x$statistic$freq[max_idx] < 0.5)
      freq <- weighted.mean(x$statistic$freq[idx], significance[idx])
    else
      freq <- 0.5
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
  
  # if (abs(omega0 * N - round(omega0 * N)) < .Machine$double.eps^0.5)
  #   delta0 <- 0
  # else {
  #   o <- order(abs(omega[omega <= 0.5] - omega0))
  #   delta0 <- abs(omega0 - omega[o[2]])
  # }
  # 
  # deltas <- delta0 + omega
  
  x <- em_harmonic(N, omega0, C)
  
  spec <- abs(fft(x)[seq_len(N %/% 2 + 1)])^2
  spec <- spec / sum(spec)
  
  # S <- numeric(length(deltas))
  # for (i in seq_along(deltas)) {
  #   S[i] <- sum(spec[abs(omega - omega0) <= deltas[i] + .Machine$double.eps^0.5])
  #   # if (S[i] > threshold) {
  #   #   idx <- i
  #   #   break
  #   # }
  #   if (abs(S[i] - 1) <= .Machine$double.eps^0.5)
  #     break
  # }
  
  deltas <- abs(omega - omega0)
  o <- order(deltas)
  deltas <- deltas[o]
  cumspecs <- cumsum(spec[o])
  idx <- which.max(cumspecs > threshold & deltas > 0)
  
  # idx <- which.max(diff(S) < eps & S[-1] > threshold) + 1
  # idx <- which.max(S > threshold)
  
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
      auto_threshold = auto_trend_threshold,
      delta = 1e-3,
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
  
  result_df <- data.frame(matrix(nrow = 0, ncol = 0))
  result_df[1, "mss"] <- NA
  result_df[1, "period"] <- NA
  result_df[1, "index1"] <- NA
  result_df[1, "index2"] <- NA
  
  m <- mcssa(
    x = noise_est,
    L = L,
    basis = "cos",
    fixed = fixed,
    conf.level = conf.level,
    freq.range = c(auto_trend_freq, 0.5)
  )
  
  result_df[1, "pval"] <- m$p.value
  
  freq.exclude <- list(c(0, auto_trend_freq))
  
  freq <- numeric()
  freq.bins <- list()
  n_triples <- integer()
  proj_vectors <- list()
  i <- 1
  
  groups_current <- groups
  periodics_components <- integer()
  periodics_est <- rep(0, N)
  
  # Fs <- reconstruct(ssa_obj, as.list(seq_len(max(groups)))) |> unlist() |> matrix(nrow = N)
  # pgs <- Rssa:::pgram(Fs)
  # # specs <- sweep(pgs$spec, 2, colSums(pgs$spec), FUN = "/")
  # specs <- pgs$spec
  # freqs <- pgs$freq
  # spec_x <- rowSums(specs)
  
  while (m$reject && i <= maxit) {
    freq[i] <- get_signif_freq(m, n_periodics = 1, weighted = FALSE)
    delta <- calculate_delta(N, freq[i], C_max)
    freq.bins[[i]] <- c(freq[i] - delta - 1e-9, freq[i] + delta + 1e-9)
    
    result_df[i + 1, "period"] <- 1 / freq[i]
    
    n_triples[i] <- ifelse(freq[i] == 0.5, 1, 2)
    
    # freq_mask <- (freqs >= freq.bins[[i]][1]) & (freqs <= freq.bins[[i]][2])
    # spec_x <- rowSums(specs[, groups_current, drop = FALSE])
    # 
    # powers <- apply(specs[, groups_current, drop = FALSE], 2, function(col) sum(col[freq_mask]))
    # o <- order(powers, decreasing = TRUE)
    # sum_power <- 0
    # max_sum_power <- sum(spec_x[freq_mask])
    # threshold <- 0.7
    # components <- c()
    # for (oo in o) {
    #   sum_power <- sum_power + powers[oo]
    #   components <- c(components, groups_current[oo])
    #   if (sum_power / max_sum_power > threshold)
    #     break
    # }
    # groups_current <- setdiff(groups_current, components)
    
    g <- grouping.auto.pgram(ssa_obj, groups_current, "series", freq.bins[i], auto_periodic_threshold)
    if (length(g)) {
      g <- g[[1]]
      len <- min(length(g), n_triples[i])
      components <- g[seq_len(len)]
      groups_current <- setdiff(groups_current, components)
    } else {
      components <- integer()
      freq[i] <- NA
      freq.exclude[[length(freq.exclude) + 1]] <- freq.bins[[i]]
    }

    periodic <- reconstruct(ssa_obj, list(components))[[1]]
    
    noise_est <- noise_est - periodic
    periodics_est <- periodics_est + periodic
    periodics_components <- c(periodics_components, components)
    
    mss <- mean(periodic^2)
    result_df[i + 1, "index1"] <- components[1]
    result_df[i + 1, "index2"] <- components[2]
    result_df[i + 1, "mss"] <- mss
    
    model <- as.list(arfima_whittle(noise_est, c(fixed$phi, fixed$d), freq.exclude))
    
    if ((i >= maxit) || !length(groups_current))
      break
    
    m <- mcssa(
      x = noise_est,
      L = L,
      basis = "cos",
      model = model,
      conf.level = conf.level,
      freq.range = NULL,
      freq.exclude = freq.exclude
    )
    
    result_df[i + 1, "pval"] <- m$p.value
    
    i <- i + 1
  }
  
  list(
    trend = trend_est,
    trend_components = trend_components,
    periodics = periodics_est,
    periodics_components = periodics_components,
    freq.bins = freq.bins,
    freqs = na.omit(freq),
    result_df = result_df,
    signal_rank = length(periodics_components) + length(trend_components)
  )
}


ssa_aidmc <- function(x,
                      L = (N + 1) %/% 2,
                      rank = 30,
                      conf.level1 = 0.95,
                      conf.level2 = 0.95,
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
    conf.level = conf.level1,
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
  
  if (length(groups) == 0)
    return(dec)
  
  model_obj <- dec$model_obj
  if (dec$signal_rank == 0) {
    model_obj <- ssa(x, L, svd.method = "svd", neig = rank)
  }
  
  m <- auto_mcssa(
    x_resid,
    ssa_obj = model_obj,
    L = mcssa.arguments$L,
    groups = groups,
    conf.level = conf.level2,
    auto_trend_freq = auto_trend_freq,
    auto_periodic_threshold = auto_periodic_threshold,
    C_max = C_max,
    method = method,
    maxit = maxit,
    detrend = FALSE,
    fixed = mcssa.arguments$fixed,
    ...
  )
  j_max <- max(dec$signal_rank, m$periodics_components)
  if (length(m$periodics_components)) {
    wcors <- suppressWarnings(wcor(model_obj, groups = seq_len(rank))[m$periodics_components, , drop = FALSE])
    for (i in seq_along(m$periodics_components)) {
      w <- wcors[i, ]
      j_max <- max(j_max, which(w >= 0.9) |> tail())
    }
  } else {
    return(dec)
  }
  
  indices <- setdiff(seq_len(j_max), dec$trend_indices)
  eoss <- eossa_new(
    model_obj,
    nested.groups = list(indices),
    clust_type = "distance",
    delta = 1e-4
  )
  
  x_detrended = x - dec$trend

  periodic_model <- auto_periodics_mc(
    x_detrended,
    model_obj = eoss,
    periodic_indices = indices,
    conf.level = conf.level2,
    tau_threshold = tau_threshold,
    auto_threshold = auto_periodic_threshold,
    mcssa.arguments = mcssa.arguments,
    C_max = C_max,
    p_0 = p_0,
    trace = trace,
    auto_trend_freq = auto_trend_freq,
    ...
  )
  periodic_indices <- c(periodic_model$true_two_el_indices, periodic_model$true_one_el_indices)
  list(
    trend = dec$trend,
    trend_indices = dec$trend_indices,
    true_two_el_indices = periodic_model$true_two_el_indices,
    true_one_el_indices = periodic_model$true_one_el_indices,
    periodics = periodic_model$periodics,
    result_df = periodic_model$result_df,
    model_obj = eoss,
    signal_rank = length(dec$trend_indices) + length(periodic_indices)
  )
}
