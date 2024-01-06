source("mcmssa_utils.R")

est.model.extract <- function(f, L = (N + 1) %/% 2, alpha) {
  N <- length(f)
  m <- MonteCarloSSA(f = f,
                     L = L,
                     D = 1,
                     basis = "ev",
                     kind = "ev",
                     G = G,
                     level.conf = 1 - alphas_corrected_arima(alpha),
                     composite = "none")
  if (m$reject) {
    suspect <- which(m$v > m$upper)
    s <- ssa(f, L, kind = "toeplitz-ssa")
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
