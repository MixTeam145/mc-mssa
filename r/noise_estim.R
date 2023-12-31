library(pracma)

mu2 <- function(x, N) {
  -1/N + 2/N^2 * ((N-x^N)/(1-x) - x * (1-x^(N-1)) / (1-x)^2)
}

autocov.est <- function(f, lag) {
  res <- 0
  N <- length(f)
  d <- f - mean(f)
  for (i in 1:(N-lag)) {
    res <- res + d[i] * d[i+lag]
  }
  res / (N - lag)
}

est.model.as <- function(f) {
  N <- length(f)
  c0 <- autocov.est(f,0)
  c1 <- autocov.est(f,1)
  start <- c1/c0
  est.varphi <- newtonRaphson(
    function(x) return((x-mu2(x,N))/(1-mu2(x,N))-start),
    x0 = start
    )$root
  est.delta <- sqrt(c0/(1-mu2(est.varphi,N))*(1-est.varphi^2))
  est.model <- list(varphi = est.varphi, delta = est.delta, N = M)
  est.model
}

est.model.extract <- function(f, L = (N + 1) %/% 2, alpha) {
  N <- length(f)
  m <- MonteCarloSSA(f = f,
                     L = L,
                     D = 1,
                     basis = "ev",
                     kind = "ev",
                     G = G,
                     level.conf = 1 - alphas.corrected.arima(alpha),
                     composite = "none")
  if (m$reject) {
    suspect <- which(m$v > m$upper)
    s <- ssa(f, L = 50, D = 1, kind = "toeplitz-ssa")
    r <- reconstruct(s, groups = list(signal = suspect))
    estModel <- est.model.arima(residuals(r))
  }
  else
    estModel <- est.model.arima(f)
  estModel
}

freq.est <- function(u) {
  s <- ssa(u, kind = "1d-ssa")
  p <- parestimate(s, groups = list(1:2))
  p$frequencies[[1]]
}
