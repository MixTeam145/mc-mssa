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

# Quantile algorithm 
type <- 8

### Functions for Monte Carlo SSA

# Spectral density of ARFIMA(1, d, 0) model
spec_arfima <- function(w,
                        phi = 0,
                        d = 0,
                        sigma2 = 1) {
  sigma2 * (2 * sin(pi * w)) ^ (-2 * d)  /
    abs(1 - phi * exp(-2i * pi * w)) ^ 2
}

# Maximum likelihood estimation of ARFIMA(1, d, 0) model
arfima_mle <- function(x, fixed = NULL) {
  n <- length(x)
  
  if (is.null(fixed))
    fixed <- rep(NA, 2)
  
  mask <- is.na(fixed)
  
  objective <- function(p) {
    par <- fixed
    par[mask] <- p
    r <- tacvfARFIMA(phi = par[1], dfrac = par[2], maxlag = n - 1)
    -DLLoglikelihood(r, x)
  }
  
  init <- c(0, 0)
  lower <- c(-1, -0.5) + 1e-4
  upper <- c(1, 0.5) - 1e-4
  
  opt <- optim(
    init[mask],
    objective,
    method = "L-BFGS-B",
    lower = lower[mask],
    upper = upper[mask]
  )
  
  coef <- fixed
  coef[mask] <- opt$par
  names(coef) <- c("phi", "d")
  
  r <- tacvfARFIMA(phi = coef[1], dfrac = coef[2], maxlag = n - 1)
  error <- DLResiduals(r, x)
  
  c(coef, sigma2 = mean(error^2))
}

# Whittle estimation of ARFIMA(1, d, 0) model
arfima_whittle <- function(x, fixed = NULL) {
  n <- length(x)
  m <- (n - 1) %/% 2
  
  # Periodogram
  per <- Mod(fft(x)[2:(m + 1)]) ^ 2 / n
  freq <- 1:m / n
  
  if (is.null(fixed))
    fixed <- rep(NA, 2)
  
  mask <- is.na(fixed)
  
  # Whittle loglikelihood
  objective <- function(p) {
    par <- fixed
    par[mask] <- p
    g <- spec_arfima(freq, par[1], par[2])
    sigma2 <- mean(per / g)
    loglike <- -log(sigma2) - mean(log(g)) - 1
    - loglike
  }
  
  init <- c(0, 0)
  lower <- c(-1, -0.5) + 1e-4
  upper <- c(1, 0.5) - 1e-4
  
  opt <- optim(
    init[mask],
    objective,
    method = "L-BFGS-B",
    lower = lower[mask],
    upper = upper[mask]
  )
  
  coef <- fixed
  coef[mask] <- opt$par
  names(coef) <- c("phi", "d")
  
  c(coef, sigma2 = mean(per / spec_arfima(freq, coef[1], coef[2])))
}


# Estimate AR(1) or FI(d) parameters
est_model <- function(f, model = c("ar1", "fi")) {
  model <- match.arg(model)
  
  lmodel <- "n"
  nar <- 1
  if (identical(model, "fi")) {
    lmodel <- 'd'
    nar <- 0
  }
  
  fit <- arfima::arfima(f,
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
  f <- f - colMeans(f)
  N <- length(f[, 1]) # assert equal length in each channel
  K <- N - L + 1
  
  # X_res <- matrix(0, nrow = L, ncol = K * D)
  # for (channel in seq_len(D)) {
  #   tX <- sapply(1:L, function(i) f[i:(i + K - 1), channel])
  #   X_res[, (1 + (channel - 1) * K):(channel * K)] <- t(tX)
  # }
  
  f.fft <- fft(c(f[(N-L+1):N], f[1:(N-L)]))
  
  if (kind == "rows") {
    p <- X_res %*% W # Projection on rows
  }
  else {
    # p <- t(X_res) %*% W # Projection on columns
    p <- mvfft(W * f.fft, inverse = TRUE)[L:N, ] / N # Projection on columns
  }
  colSums(Mod(p) ^ 2 / N) # divide by N to weaken the dependence on t.s. length
}

# Estimate vector main frequency by ESPRIT
est_freq <- function(v) {
  s <- ssa(v, neig = 2)
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
  
  res <- list()
  
  if (factor.v) {
    res$W <- s$V 
  }
  else {
    res$W <- s$U
  }
  
  res
}

# Generate vectors for projections as cosine vectors
basis.cos <- function(L) {
  W <- matrix(0, nrow = L, ncol = L)
  separat <- 1 / (2 * L)
  freq <- seq(0, 0.5 - separat, separat) # Grid of frequencies
  for (i in seq_along(freq)) {
    W[, i] <- cos(2 * pi * freq[i] * 1:L)
    W[, i] <- W[, i] / Norm(W[, i])
  }
  list(W = W, freq = freq)
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
      toepl <- toepl + matrix.toeplitz(model[[channel]]$phi, L)
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
  list(W = U, freq = freq)
}
### end

what.reject <- function(x) {
  rej <- x$projec_vectors$contribution < x$predint$lower |
    x$projec_vectors$contribution > x$predint$upper
  x$projec_vectors$freq[rej]
}

### Main functions for multiple Monte Carlo SSA
# Make multiple test
do.test <- function(x,
                    projec_vectors,
                    conf.level,
                    G,
                    two.tailed = FALSE) {
  
  if (length(x$freq.range)) {
    if (!length(projec_vectors$freq))
      projec_vectors$freq <- apply(projec_vectors$W, 2, est_freq)
    idx <-
      projec_vectors$freq >=  x$freq.range[1] & projec_vectors$freq <= x$freq.range[2]
  }
  else
    idx <- seq_len(ncol(projec_vectors$W))
  if (!any(idx))
    stop("No vectors with given frequency range, aborting")
  
  projec_vectors$W <- projec_vectors$W[, idx, drop = FALSE]
  
  W.fft <- mvfft(
    rbind(
      matrix(0, nrow = dim(x$series)[1] - x$L, ncol = dim(projec_vectors$W)[2]),
      projec_vectors$W[x$L:1, ]
    )
  )
  
  P <- replicate(
    G,
    projec(x$model, x$L, W.fft, x$kind),
    simplify = FALSE
  )
  P <- do.call(cbind, P)
  v <- projec(x$series, x$L, W.fft, x$kind)
  
  x$projec_vectors <- list(
    W = projec_vectors$W,
    freq = projec_vectors$freq[idx],
    contribution = v
  )
  
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
#' @param G Number of surrogates
#' @param conf.level Confidence level
#' @param two.tailed If TRUE performs two-tailed test
#' @param est.freq If TRUE estimates the main frequency of each projection vector
#' @param freq.range Potential signal frequency range
#' @param composite If TRUE performs test with composite null hypothesis (noise + nuisance signal)
mcssa <- function(f,
                  L,
                  basis = c("ev", "t", "cos"),
                  kind = c("columns", "rows"),
                  toeplitz.method = c("no", "sum", "block"),
                  model = c("ar1", "fi"),
                  model0 = list(phi = NA, dfrac = NA, sigma2 = NA, N = NA),
                  G = 1000,
                  conf.level = 0.8,
                  two.tailed = FALSE,
                  est.freq = TRUE,
                  freq.range = c(0, 0.5),
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
    model0 <- vector("list", D)
    for (channel in seq_len(D)) {
      model0[[channel]] <- est_model(f[, channel], model)
    }
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
  else if (basis == "t") {
    if (kind == 'fa')
      projec_vectors <- basis.toeplitz(estModel, N - L + 1, D, factor.v = TRUE)
    else
      projec_vectors <- basis.toeplitz(model0, L, D)
  }
  else if (D == 1) {
    projec_vectors <- basis.cos(L)
    
  }
  else {
    stop()
  }
  
  if (est.freq)
    this$freq.range <- freq.range
  else if (!identical(freq.range, c(0, 0.5)))
    warning("est.freq is set FALSE, freq.range will be omitted")
  
  this <- do.test(
    this,
    projec_vectors,
    conf.level,
    G,
    two.tailed
  )
  this
}


plot.mcssa <- function(x, by.order = FALSE, text.size = 10, point.size = 1) {
  if (!length(x$freq.range))
    warning("The main frequencies of projection vectors missing, estimating it now")
    x$projec_vectors$freq <- apply(x$projec_vectors$W, 2, est_freq)
  
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
    geom_point(size = point.size) +
    geom_errorbar(aes(ymin = lower, ymax = upper), color = "black", linewidth = 0.2) +
    scale_color_manual(values = c("blue", "red")) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      panel.background = element_blank(),
      legend.position = "none",
      axis.title = element_text(size = text.size)
    )
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
  else if (x$basis == "t")
    cat("eigenvectors of theoretical autocovariance matrix (exact test)")
  else
    cat("cosines with frequency j / (2L) (exact test)")
  cat("\nType of projection: on", x$kind, "of trajectory matrix")
  cat("\nNumber of projection vectors:", ncol(x$projec_vectors$W))
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

# Davies-Harte algorithm to simulate Gaussian process with given autocovariance function
DH.sim <- function(n, acvf, ...) {
  # Next power of two greater than 'n'
  N <- nextn(n, 2)
  
  # Autocovariance sequence
  acvs <- acvf(maxlag = N, ...)
  
  ak <- Re(fft(c(acvs, acvs[N:2])))
  
  if (any(ak < 0))
    stop("Davies-Harte nonnegativity condition is not fulfilled")
  
  # Gaussian white noise
  eps <- rnorm(2 * N)
  
  ks <- 2:N
  
  y0 <- sqrt(ak[1]) * eps[1]
  yN <- sqrt(ak[N + 1]) * eps[2 * N]
  yk <- sqrt(0.5 * ak[ks]) *
    complex(real = eps[2 * ks - 2], imaginary = eps[2 * ks - 1])
  
  y <- c(y0, yk, yN, Conj(rev(yk)))
  
  x <- Re(fft(y, inverse = TRUE)) / sqrt(2 * N)
  
  # Truncate the resulted series
  x[1:n]
}

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
  # r <- tacvfARFIMA(phi = model$phi, dfrac = model$dfrac, sigma2 = model$sigma2, maxlag = model$N - 1)
  # xi <- ltsa::DLSimulate(model$N, r)
  xi <- DH.sim(
    model$N,
    tacvfARFIMA,
    phi = model$phi,
    dfrac = model$dfrac,
    sigma2 = model$sigma2
  )
  if (!is.null(model$signal)) # composite null hypothesis
    xi <- xi + model$signal
  f <- xi + signal
  as.vector(f)
}

# Generates a multivariate ts
generate <- function(D, model, signal = matrix(0, nrow = N, ncol = D)) {
  N <- model[[1]]$N
  res <- vector("list", D)
  for (channel in seq_len(D)) {
    res[[channel]] <- generate_channel(model[[channel]], signal[, channel])
  }
  f <- do.call(cbind, res)
  f
}
