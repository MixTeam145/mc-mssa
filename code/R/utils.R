library("arfima")


# Calculate spectral density of an ARFIMA(1, d, 0) model
spec_arfima <- function(w,
                        phi = 0,
                        d = 0,
                        sigma2 = 1) {
  sigma2 * (2 * sin(pi * w)) ^ (-2 * d)  /
    abs(1 - phi * exp(-2i * pi * w)) ^ 2
}


# Calculate signal-to-noise ratio
snr <- function(s, ...) {
  n <- length(s)
  per <- Mod(fft(s))^2 / n
  freq <- 0:(n - 1) / n
  mean(per / spec_arfima(freq, ...))
}


# Fit an ARFIMA(1, d, 0) model using maximum likelihood estimation
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


# Fit an ARFIMA(1, d, 0) model using Whittle estimation
arfima_whittle <- function(x, fixed = NULL, freq.exclude) {
  n <- length(x)
  m <- (n - 1) %/% 2
  
  # Periodogram
  spec <- Mod(fft(x)[2:(m + 1)]) ^ 2 / n
  freq <- 1:m / n
  
  if (!missing(freq.exclude)) {
    idx <- sapply(freq.exclude, function(fb) freq < fb[1] | freq > fb[2])
    idx <- which(rowSums(idx) == length(freq.exclude))
    spec <- spec[idx]
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
    sigma2 <- mean(spec / g)
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
  
  c(coef, sigma2 = mean(spec / spec_arfima(freq, coef[1], coef[2])))
}


# Simulate a Gaussian process with given autocovariance function using Davies-Harte algorithm
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


# Simulate from an ARFIMA model
arfima.sim <- function(model, N) {
  x <- DH.sim(
    N,
    tacvfARFIMA,
    phi = model$phi,
    dfrac = model$d,
    sigma2 = model$sigma2
  )
  
  if (!is.null(model$signal)) # composite null hypothesis
    x <- x + model$signal
  
  as.ts(x)
}


# Generate a multivariate ts
generate <- function(model, N, D = 1, demean = FALSE) {
  if (D == 1 | (D > 1 & !is.list(model[[1]])))
    model <- replicate(D, model, simplify = FALSE)
  
  x <- ts(matrix(nrow = N, ncol = D))
  for (channel in seq_len(D))
    x[, channel] <- arfima.sim(model[[channel]], N)
  
  if (demean)
    x <- sweep(x, 2, colMeans(x))
  
  x
}


# Generate exponentially modulated harmonic
em_harmonic <- function(N, omega, C = 0, A = 1, phi = 0) {
  t <- 1:N
  alpha <- C / N
  signal <- A * exp(alpha * t) * cos(2 * pi * t * omega + phi)
  as.ts(signal)
}
