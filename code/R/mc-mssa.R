# Multiple Monte Carlo SSA
# [Golyandina N. Detection of signals by Monte Carlo singular spectrum analysis:
# multiple testing // Statistics and Its Interface. — 2023. — Vol. 16, no. 1. — P. 147–157.]

# Contribution of Egor Poteshkin:
# correction of liberal criteria is implemented
# projection vectors with parameter basis="ev" and D=1 are from Toeplitz SSA of original time series
# Toeplitz MC-MSSA (D > 1) is draft

library("Rssa")
library("pracma")
library("ggplot2")
library("matrixStats")
library("magic")
library("arfima")

source("toeplitz_mssa.R")

type = 8

### Functions for Monte Carlo SSA
# Estimate AR(1) or FI(d) parameters
est_model <- function(f, model = c("ar1", "fi")) {
  model <- match.arg(model)
  
  lmodel <- "n"
  nar <- 1
  if (identical(model, "fi")) {
    lmodel <- 'd'
    nar <- 0
  }
  
  fit <- arfima(f,
                c(nar, 0, 0),
                lmodel = lmodel,
                dmean = FALSE,
                quiet = TRUE)$modes[[1]]
  
  list(
    phi = fit$phi,
    dfrac = fit$dfrac,
    sigma2 = fit$sigma2,
    N = length(f)
  )
}

# Compute squared norms of projections
projec <- function(data, L, W, kind = c("columns", "rows")) {
  if (is.list(data)) {
    # data are given by a model
    D <- length(data)
    f <- generate(D, data)
  } else {
    # data are given by a series
    f <- data
    D <- dim(f)[2]
  }
  N <- length(f[, 1]) # assert equal length in each channel
  K <- N - L + 1
  X_res <- matrix(0, nrow = L, ncol = K * D)
  for (channel in 1:D) {
    tX <- sapply(1:L, function(i) f[i:(i + K - 1), channel])
    X_res[, (1 + (channel - 1) * K):(channel * K)] <- t(tX)
  }
  if (kind == "rows") {
    p <- X_res %*% W # Projection on rows
  }
  else {
    p <- t(X_res) %*% W # Projection on columns
  }
  colSums(p ^ 2 / N) # divide by N to weaken the dependence on t.s. length
}

# Estimate vector main frequency by ESPRIT
est_freq <- function(v) {
  s <- ssa(v)
  p <- parestimate(s, list(1:2))
  freq <- p$frequencies[[1]]
  freq
}

# Generate projection vectors corresponding to eigenvectors produced by ts
basis.ev <- function(ts, L, toeplitz.method, factor.v = FALSE) {
  D <- dim(ts)[2]
  neig <- min(L, D * (length(ts[,1]) - L + 1))
  
  if (D == 1) {
    s <- ssa(ts, L, neig, kind = "toeplitz-ssa")
  } else if (identical(toeplitz.method, "sum") || identical(toeplitz.method, "block")) {
    s <- toeplitz.mssa(ts, L, D, toeplitz.method, neig)
  } else {
    s <- ssa(ts, L, neig, kind = "mssa")
  }
  
  #freq <- apply(s$U, 2, est_freq)
  #basis <- list(freq = freq)
  
  if (factor.v) {
    s$V
  }
  else {
    s$U
  }
}

matrix.toeplitz <- function(phi, L) {
  toeplitz(phi ^ (0:(L - 1)))
}

# Generate vectors for projections corresponding to eigenvectors of theoretical autocovariance matrix 
basis.toeplitz <- function(model, L, D, fa = F) {
  if (fa) {
    # here we assume that L param represents K = N - L + 1
    toepl.array <- list()
    for (channel in 1:D) {
      toepl.array[[channel]] <- model[[channel]]$delta * matrix.toeplitz(model[[channel]]$varphi, L)
    }
    toepl <- do.call("adiag", toepl.array)
    s <- svd(toepl, nv = L)
    U <- s$v
  }
  else {
    toepl <- matrix(data = 0,
                    nrow = L,
                    ncol = L)
    for (channel in 1:D) {
      toepl <- toepl + matrix.toeplitz(model[[channel]]$varphi, L)
    }
    s <- svd(toepl, L)
    U <- s$u
  }
  
  freq <- numeric(0)
  for (i in 1:L) {
    #ss <- ssa(s$U[,i], kind = "toeplitz-ssa")
    ss <- ssa(U[, i], kind = "1d-ssa")
    #estimation of the main frequency by ESPRIT
    p <- parestimate(ss, groups = list(1:2))
    freq[i] <- p$frequencies[[1]]
  }
  list(U = U, freq = freq)
}
### end

what.reject <- function(res){
  rej <- (res$v[res$idx] < res$lower | res$v[res$idx] > res$upper) & res$idx[res$idx]
  print(res$freq[res$idx][rej==TRUE])
}

### Main functions for multiple Monte Carlo SSA
# Make multiple test
do.test <- function(x,
                    projec_vectors,
                    conf.level,
                    G,
                    two.tailed = FALSE,
                    freq.range = c(0, 0.5)) {
  freq <- apply(projec_vectors, 2, est_freq)
  idx <- freq >=  freq.range[1] & freq <= freq.range[2]
  if (!any(idx))
    stop("No vectors with given frequency range, aborting")
  
  projec_vectors <- projec_vectors[, idx, drop = FALSE]
  
  P <- replicate(G, projec(x$model, x$L, projec_vectors, x$kind))
  v <- projec(x$series, x$L, projec_vectors, x$kind)
  
  if (is.vector(P))
    P <- rbind(P)
  
  x$projec_vectors <- list(
    W = projec_vectors,
    freq = freq[idx],
    contribution = v
  )
  x$freq.range <- freq.range
  
  means <- apply(P, 1, mean)
  sds <- apply(P, 1, sd)
  
  if (!two.tailed) {
    eta <- apply(P, 2, function(p) max((p - means) / sds))
    t <- max((v - means) / sds)
  }
  else {
    eta <- apply(P, 2, function(p) max(abs(p - means) / sds))
    t <- max(abs(v - means) / sds)
  }
  
  if (!is.null(conf.level)) {
    q.upper <- quantile(eta, probs = conf.level, type = type)
    x$predint <- list()
    x$predint$upper <- means + q.upper * sds
    if (!two.tailed) {
      q.lower <- 0
      x$predint$lower <- 0
      x$reject <- as.logical(t > q.upper)
    }
    else {
      q.lower <- -q.upper
      x$predint$lower <- means + q.lower * sds
      x$reject <- as.logical(t > q.upper | t < q.lower)
    }
    x$conf.level <- conf.level
  }
  x$p.value <- 1 - ecdf(eta)(t)

  x
}

#' The wrapped function for Multiple Monte Carlo SSA
#' 
#' @param f Time series
#' @param L Window length
#' @param basis Type of vectors for projection
#' @param kind Projection on columns or rows of trajectory matrix
#' @param toeplitz.method Toeplitz MSSA decomposition method
#' @param model Noise model to be fitted (will be omitted, if model0 is specified)
#' @param model0 Exact noise model
#' @param freq.range Potential signal frequency range
#' @param G Number of surrogates
#' @param conf.level Confidence level
#' @param two.tailed If TRUE performs two-tailed test
#' @param composite If TRUE performs test with composite null hypothesis (noise + nuisance signal)
mcssa <- function(f,
                  L,
                  basis = c("ev", "t"),
                  kind = c("columns", "rows"),
                  toeplitz.method = c("no", "sum", "block"),
                  model = c("ar1", "fi"),
                  model0 = list(phi = NA, dfrac = NA, sigma2 = NA, N = NA),
                  freq.range = c(0, 0.5),
                  G = 1000,
                  conf.level = 0.8,
                  two.tailed = FALSE,
                  composite = FALSE) {
  if (is.vector(f)) {
    f <- as.matrix(f)
    if (!missing(model0)) {
      model0 <- list(model0)
    }
  } else if (composite) {
    stop("mc-ssa with nuisance signal for multivariate ts is not implemented")
  }
  
  model <- match.arg(model)
  kind <- match.arg(kind)
  toeplitz.method <- match.arg(toeplitz.method)
  
  D <- dim(f)[2]
  
  if (missing(model0)) {
    model0 <- lapply(seq_len(D), function(i)
      est_model(f[, i], model))
  }
  
  this <- list(
    series = f,
    L = L,
    model = model0,
    basis = basis,
    kind = kind
  )
  class(this) <- "mcssa"
  
  if (basis == "ev") {
    f.basis <- f
    # comment next 2 lines to project vectors of original series (another version of the nuisance algorithm)
    if (composite)
      f.basis <- f - model0$signal
    if (kind == 'fa')
      projec_vectors <- basis.ev(f.basis, L, toeplitz.method, factor.v = TRUE)
    else
      projec_vectors <- basis.ev(f.basis, L, toeplitz.method)
  }
  else {
    if (kind == 'fa')
      projec_vectors <- basis.toeplitz(estModel, N - L + 1, D, factor.v = TRUE)
    else
      projec_vectors <- basis.toeplitz(estModel, L, D)
  }
  
  this <- do.test(
    this,
    projec_vectors,
    conf.level,
    G,
    two.tailed,
    freq.range
  )
  this
}


plot.mcssa <- function(x, by.order = FALSE) {
  df <-
    data.frame(
      frequency = x$projec_vectors$freq,
      contribution = x$projec_vectors$contribution,
      lower = x$predint$lower,
      upper = x$predint$upper
    )
  df <- df |>
    mutate(reject = contribution < lower | contribution > upper)

  if (by.order) {
    df <- df |> arrange(desc(contribution))
    df$index <- 1:length(df$frequency)
    p <- ggplot(df, aes(index, contribution, color = reject))
  } else {
    p <- ggplot(df, aes(frequency, contribution, color = reject))
  }
  
  p <- p +
    geom_point() +
    geom_errorbar(aes(ymin = lower, ymax = upper), color = "blue") +
    scale_color_manual(values = c("blue", "red")) +
    theme(legend.position = "none") 
  p
}

print.mcssa <- function(x) {
  N <- length(x$series[, 1])
  D <- dim(x$series)[2]
  cat("Series length:", rep(N, D))
  cat("\nWindow length:", x$L)
  cat("\nProjection vectors: ")
  if (x$basis == "ev")
    cat("based on series (liberal test)")
  else
    cat("eigenvectors of theoretical autocovariance matrix (exact test)")
  cat("\nType of projection: on", x$kind, "of trajectory matrix")
  cat("\nNumber of projection vectors:", length(x$projec_vectors$freq))
  cat("\np-value:", x$p.value)
  if (!is.null(x$conf.level))
    cat("\nNull hypothesis is", if (!x$reject) "not", "rejected")
  invisible(x)
}


# There is implemenatation of correction liberal/conservative criteria
correction <- function(p.values, alphas = 0:1000 / 1000) {
  alphaI <- sapply(alphas, function(a) mean(p.values < a))
  alphas.fun <- approxfun(alphaI, alphas, rule = 2)
  alphas.fun
}

# CI
conf.interval <- function(p.values, alpha) {
  corr <- data.frame(upper=numeric(length(alphas)), lower=numeric(length(alphas)), alpha=numeric(length(alphas)))
  for (i in alphas_idx) 
  {
    conf <- prop.test(sum(p.values < alphas[i]), n = length(p.values))
    corr$upper[i] <- conf$conf.int[2]
    corr$lower[i] <- conf$conf.int[1]
    corr$alpha[i] <- alphas[i]
  }
  left.func <- approxfun(corr$upper, corr$alpha, rule=2)
  right.func <- approxfun(corr$lower, corr$alpha, rule=2)
  
  left <- left.func(alpha)
  right <- right.func(alpha)
  
  c(left, right)
}

###  Time series generation

# Generates sinusoidal signal with specified frequency
signal.one.channel <- function(N, omega, A = 1) {
  num <- 1:N
  if (is.null(omega))
    signal <- 0 * num
  else
    signal <- A * cos(2 * pi * num * omega + runif(1, 0, pi))
  signal
}

# Generates a time series according to the model signal + noise
generate_channel <- function(model, signal = 0) {
  # TODO: implement Davies and Harte algorithm
  xi <- arfima.sim(
    n = model$N,
    model = list(phi = model$phi, dfrac = model$dfrac),
    sigma2 = model$sigma2
  )
  xi <- xi - mean(xi)
  if (!is.null(model$signal)) # composite null hypothesis
    xi <- xi + model$signal
  f <- xi + signal
  as.vector(f)
}

# Generates a multivariate ts
generate <- function(D, model, signal = matrix(0, nrow = N, ncol = D)) {
  N <- model[[1]]$N
  res <- lapply(
    seq_len(D),
    function(i) generate_channel(model[[i]], signal[, i])
  )
  f <- matrix(unlist(res), ncol = D, nrow = N)
  f
}
