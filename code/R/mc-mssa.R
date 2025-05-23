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

source("toeplitz_mssa.R", TRUE)

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
arfima_whittle <- function(x, fixed = NULL, freq.exclude) {
  n <- length(x)
  m <- (n - 1) %/% 2
  
  # Periodogram
  per <- Mod(fft(x)[2:(m + 1)]) ^ 2 / n
  freq <- 1:m / n
  
  if (!missing(freq.exclude)) {
    idx <- sapply(freq.exclude, function(fb) freq < fb[1] | freq > fb[2])
    idx <- which(rowSums(idx) == length(freq.exclude))
    per <- per[idx]
    freq <- freq[idx]
  }
  
  # per <- per[freq >= freq.range[1] & freq <= freq.range[2]]
  # freq <- freq[freq >= freq.range[1] & freq <= freq.range[2]]
 
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
  coef["phi"] <- max(coef["phi"], 0)
  
  c(coef, sigma2 = mean(per / spec_arfima(freq, coef[1], coef[2])))
}

imvfft <- function(x) {
  mvfft(x, inverse = TRUE) / nrow(x)
}

pad <- function(x, n, after = TRUE) {
  x_padded <- matrix(0, nrow(x) + n, ncol(x))
  if (after)
    x_padded[1:nrow(x), ] <- x
  else
    x_padded[(n + 1):nrow(x_padded), ] <- x
  x_padded
}

# Compute squared norms of projections
projec <- function(x, L, W_ft, kind = c("columns", "rows")) {
  kind <- match.arg(kind)
  N <- nrow(x)
  D <- ncol(x)

  x_ft <- mvfft(x[c((N - L + 1):N, 1:(N - L)), , drop = FALSE])

  if (kind == "rows") {
    v <- matrix(0, L, ncol(W_ft[[1]]))
    for (i in seq_len(D)) {
      p <- imvfft(W_ft[[i]] * x_ft[, i])[1:L, , drop = FALSE]
      v <- v + p
    }
    v <- colSums(Mod(v)^2)
  } else {
    v <- numeric(ncol(W_ft))
    for (i in seq_len(D)) {
      p <- imvfft(W_ft * x_ft[, i])[L:N, , drop = FALSE]
      v <- v + colSums(Mod(p)^2)
    }
  }
  v / N # divide by N to weaken the dependence on series length
}

# Estimate vector main frequency by ESPRIT
estimate_freq <- function(v) {
  s <- ssa(v, neig = 2)
  p <- parestimate(s, list(1:2))
  freq <- p$frequencies[[1]]
  freq
}

# Generate projection vectors corresponding to eigenvectors produced by ts
basis.ev <- function(ts, L, decomposition.method, vectors = c("U", "V")) {
  D <- dim(ts)[2]
  neig <- min(L, D * (length(ts[, 1]) - L + 1))
  
  if (decomposition.method == "svd") {
    if (D == 1) {
      s <- ssa(ts, L, neig, kind = "1d-ssa")
    } else {
      s <- ssa(ts, L, neig, kind = "mssa")
    }
  } else {
    if (D == 1) {
      s <- ssa(ts, L, neig, kind = "toeplitz-ssa")
    } else {
      s <- toeplitz.mssa(ts, L, neig, kind = "sum")
    }
  }
  
  res <- list()
  
  if (vectors == "V") {
    res$W <- s$V
  } else {
    res$W <- s$U
  }
  res
}

# Generate projection vectors as cosine vectors
basis.cos <- function(L) {
  W <- matrix(0, nrow = L, ncol = L)
  separat <- 1 / (2 * L)
  freq <- seq(separat, 0.5, separat) # Grid of frequencies
  for (i in seq_along(freq)) {
    W[, i] <- cos(2 * pi * freq[i] * (1:L))
    W[, i] <- W[, i] / Norm(W[, i])
  }
  list(W = W, freq = freq)
}


autocov.mat <- function(model, L) {
  r <- tacvfARFIMA(
    phi = model$phi,
    dfrac = model$d,
    sigma2 = model$sigma2,
    maxlag = L - 1
  )
  toeplitz(r)
}

# Generate projection vectors corresponding to eigenvectors of theoretical autocovariance matrix 
basis.t <- function(model, L, vectors = c("U", "V")) {
  D <- length(model)
  if (vectors == "V") {
    K <- model[[1]]$N - L + 1
    Sigma <- list()
    for (channel in 1:D) {
      Sigma[[channel]] <- autocov.mat(model[[channel]], K)
    }
    Sigma <- do.call("adiag", Sigma)
    s <- svd(Sigma, nv = K)
    W <- s$v
  } else {
    Sigma <- matrix(0, L, L)
    for (channel in 1:D) {
      Sigma <- Sigma + autocov.mat(model[[channel]], L)
    }
    s <- svd(Sigma, L)
    W <- s$u
  }
  
  list(W = W)
}
### end

what.reject <- function(x) {
  rej <- x$statistic$contribution < x$predint$lower |
    x$statistic$contribution > x$predint$upper
  x$statistic$freq[rej]
}

### Main functions for multiple Monte Carlo SSA
# Make multiple test
do.test <- function(x, G, conf.level, two.tailed, freq.range) {
  N <- x$length
  L <- x$window
  D <- x$channels
  
  if (!is.null(freq.range)) {
    if (!length(x$proj_vectors$freq))
      x$proj_vectors$freq <- apply(x$proj_vectors$W, 2, estimate_freq)
    mask <-
      x$proj_vectors$freq >=  freq.range[1] &
      x$proj_vectors$freq <= freq.range[2]
  } else {
    mask <- rep(TRUE, ncol(x$proj_vectors$W))
  }
  
  if (!any(mask))
    stop("No vectors with given frequency range, aborting")
  
  x$freq.range <- freq.range
  
  # Pre-compute FFT of the reversed vectors for the fast matrix-vector product
  W <- x$proj_vectors$W
  if (x$proj.kind == "rows") {
    K <- N - L + 1
    W_ft <- list()
    for (d in seq_len(D)) {
      idx <- (d * K):((d - 1) * K + 1)
      W_ft[[d]] <- mvfft(pad(W[idx, mask, drop = FALSE], L - 1))
    }
  } else {
    W_ft <- mvfft(pad(W[L:1, mask, drop = FALSE], N - L, after = FALSE))
  }
  
  # Simulate surrogate data and calculate projection
  P <- replicate(
    G,
    projec(generate(D, x$model, demean = TRUE), L, W_ft, x$proj.kind),
    simplify = FALSE
  )
  P <- do.call(cbind, P)
  
  # Calculate projection of input time series
  v <- projec(x$series, L, W_ft, x$proj.kind)
  
  x$statistic <- list(
    freq = x$proj_vectors$freq[mask],
    contribution = v
  )
  
  means <- rowMeans(P)
  sds <- rowSds(P)
  
  if (!two.tailed) {
    eta <- apply(P, 2, function(p)
      max((p - means) / sds))
    t <- max((v - means) / sds)
  } else {
    eta <- apply(P, 2, function(p)
      max(abs(p - means) / sds))
    t <- max(abs(v - means) / sds)
  }
  
  if (!is.null(conf.level)) {
    q.upper <- quantile(eta, probs = conf.level, type = type)
    x$predint <- list()
    x$predint$upper <- means + q.upper * sds
    if (!two.tailed) {
      q.lower <- 0
      x$predint$lower <- 0
      x$reject <- t > q.upper
    } else {
      q.lower <- -q.upper
      x$predint$lower <- means + q.lower * sds
      x$reject <- t > q.upper | t < q.lower
    }
    x$conf.level <- conf.level
  }
  x$p.value <- 1 - ecdf(eta)(t)
  
  x
}

#' The wrapped function for Multiple Monte Carlo SSA
#' 
#' @param x Time series
#' @param L Window length
#' @param basis Type of vectors for projection
#' @param proj.kind Projection on columns or rows of trajectory matrix
#' @param decomposition.method Decomposition method (for SSA)
#' @param model A list of noise model parameters
#' @param fixed A list of parameters to be fixed (if `model` is missing)
#' @param G Number of surrogates
#' @param conf.level Confidence level
#' @param two.tailed If `TRUE` performs two-tailed test
#' @param freq.range Potential signal frequency range
#' @param composite If `TRUE` performs test with composite null hypothesis (noise + nuisance signal)
mcssa <- function(x,
                  L,
                  basis = c("ev", "t", "cos"),
                  proj.kind = c("columns", "rows"),
                  decomposition.method = c("toeplitz", "svd"),
                  model,
                  fixed = list(phi = NA, d = NA),
                  G = 1000,
                  conf.level = 0.8,
                  two.tailed = FALSE,
                  freq.range = c(0, 0.5),
                  composite = FALSE) {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
    if (!missing(model))
      model <- list(model)
  } else if (composite) {
    stop("mc-ssa with nuisance signal for multivariate ts is not implemented")
  }
  
  x <- scale(x, scale = FALSE)
  
  proj.kind <- match.arg(proj.kind)
  decomposition.method <- match.arg(decomposition.method)
  
  N <- dim(x)[1]
  D <- dim(x)[2]
  
  if (missing(model)) {
    model <- vector("list", D)
    freq.exclude <- list(c(0, freq.range[1]), c(freq.range[2], 0.5))
    for (channel in seq_len(D)) {
      model[[channel]] <- as.list(
        arfima_whittle(x[, channel], c(fixed$phi, fixed$d), freq.exclude)
      )
      model[[channel]]$N <- N 
    }
  }
  
  call <- match.call()
  
  this <- list(
    series = x,
    length = N,
    channels = D,
    window = L,
    model = model,
    basis = basis,
    proj.kind = proj.kind,
    call = call
  )
  class(this) <- "mcssa"
  
  if (proj.kind == "rows") {
    vectors <- "V"
  } else {
    vectors <- "U"
  }
    
  if (basis == "ev") {
    x.basis <- x
    # comment next 2 lines to project vectors of original series (another version of the nuisance algorithm)
    if (composite)
      x.basis <- x - model$signal
    proj_vectors <- basis.ev(x.basis, L, decomposition.method, vectors)
  } else if (basis == "t") {
    proj_vectors <- basis.t(model, L, vectors)
  } else if (D == 1) {
    proj_vectors <- basis.cos(L)
  } else {
    stop()
  }
  
  this$proj_vectors <- proj_vectors
  
  this <- do.test(
    this,
    G,
    conf.level,
    two.tailed,
    freq.range
  )
  this
}


plot.mcssa <- function(x, by.order = FALSE, text.size = 10, point.size = 1) {
  if (!length(x$freq.range)) {
    warning("The main frequencies of projection vectors are missing, estimating them now")
    x$statistic$freq <- apply(x$proj_vectors$W, 2, estimate_freq)
  }
  df <-
    data.frame(
      frequency = x$statistic$freq,
      contribution = x$statistic$contribution,
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
  cat("\nCall:\n", deparse(m1$call), "\n\n", sep = "")
  cat("Series length:", paste(rep(x$length, x$channels), collapse = ", "))
  cat(",\tWindow length:", x$window)
  cat("\n\nProjection vectors: ")
  if (x$basis == "ev")
    cat("based on series (liberal test)")
  else if (x$basis == "t")
    cat("eigenvectors of theoretical autocovariance matrix (exact test)")
  else
    cat("cosines with frequencies j / (2L) (exact test)")
  cat("\nType of projection: on", x$proj.kind, "of trajectory matrix")
  cat("\nNumber of projection vectors:", length(x$statistic$contribution))
  cat("\n\nP-value:", x$p.value)
  if (!is.null(x$conf.level))
    cat("\nNull hypothesis is", if (!x$reject) "not", "rejected")
  invisible(x)
}

# Correction of liberal/conservative criteria
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

# Generates a time series according to the model signal + noise
generate_channel <- function(model, signal = 0) {
  # r <- tacvfARFIMA(phi = model$phi, dfrac = model$dfrac, sigma2 = model$sigma2, maxlag = model$N - 1)
  # xi <- ltsa::DLSimulate(model$N, r)
  xi <- DH.sim(
    model$N,
    tacvfARFIMA,
    phi = model$phi,
    dfrac = model$d,
    sigma2 = model$sigma2
  )
  if (!is.null(model$signal)) # composite null hypothesis
    xi <- xi + model$signal
  f <- xi + signal
  as.vector(f)
}

# Generates a multivariate ts
generate <- function(D, model, signal = matrix(0, nrow = N, ncol = D), demean = FALSE) {
  N <- model[[1]]$N
  res <- vector("list", D)
  for (channel in seq_len(D)) {
    res[[channel]] <- generate_channel(model[[channel]], signal[, channel])
  }
  f <- do.call(cbind, res)
  if (demean) {
    f <- sweep(f, 2, colMeans(f))
  }
  f
}
