# Rssa (ttps://CRAN.R-project.org/package=Rssa): only Basic MSSA is implemented
# Contribution of Egor Poteshkin: Toeplitz MSSA in two versions is implemented
# "Sum" version of Toeplitz MSSA is preferable

library("Rssa")
library("dplyr")
library("pracma")

# Calculate toeplitz cross-covariance matrix
Lcov <- function(f1, f2, K) {
  f1 <- as.vector(f1)
  f2 <- as.vector(f2)
  N <- length(f1)
  c1 <- sapply(0:(K - 1), function(i)
    sum(f1 * lag(f2, i), na.rm = TRUE) / (N - i))
  c2 <- sapply(0:(K - 1), function(i)
    sum(lag(f1, i) * f2, na.rm = TRUE) / (N - i))
  Toeplitz(c1, c2)
}

toeplitz.mssa <- function(ts, # time series
                          L, # window length
                          neig = NULL, # number of desired non-zero eigenvalues
                          kind = c("sum", "block") # decomposition method
                          ) {
  N <- dim(ts)[1] # assert equal length in each channel
  D <- dim(ts)[2]
  K <- N - L + 1
  
  kind <- match.arg(kind)
  
  if (is.null(neig))
    neig <- min(L, D * K)
  
  this <- list("F" = ts, N = N, L = L, K = K, D = D, kind = kind)
  
  traj.mat <- new.hbhmat(ts, L = c(L, 1))
  if (identical(kind, "sum")) {
    toepl.mat <- matrix(0, nrow = L, ncol = L)
    for (i in 1:D)
      toepl.mat <- toepl.mat + Lcov(ts[, i], ts[, i], L)
  } else {
    toepl.mat <- rbind()
    for (i in 1:D) {
      mat <- cbind()
      for (j in 1:D)
        mat <- cbind(mat, Lcov(ts[, i], ts[, j], K))
      toepl.mat <- rbind(toepl.mat, mat)
    }
    traj.mat <- t(traj.mat)
  }
  
  S <- eigen(toepl.mat, symmetric = TRUE)
  U <- S$vectors
  Z <- crossprod(traj.mat, U)
  sigma <- apply(Z, 2, function(x)
    sqrt(sum(x ^ 2)))
  V <- sweep(Z, 2, sigma, FUN = "/")
  
  o <- 1:min(neig, L, D * K)
  sigma <- sigma[o]
  U <- U[, o, drop = FALSE]
  V <- V[, o, drop = FALSE]
  if (identical(kind, "sum")) {
    this$U <- U
    this$V <- V
  }
  else {
    this$U <- V
    this$V <- U
  }
  this$sigma <- sigma
  this
}

diag.avg <- function(x, group) {
  X <- matrix(0, nrow = x$L, ncol = x$D * x$K)
  for (i in group) {
    X <- X + x$sigma[i] * x$U[, i] %*% t(x$V[, i])
  }
  series <- matrix(0, nrow = x$N, ncol = x$D)
  for (i in 1:x$D) {
    M <- X[, ((i-1)*x$K + 1):(i*x$K)]
    if (x$K < x$L)
      M <- t(M)
    series[, i] <- sapply(0:(x$N-1), function(k) mean(M[row(M) + col(M) == k + 2]))
  }
  ts(series)
}

toeplitz.reconstruct <- function(x, # result of toeplitz.mssa
                                 groups # desired grouping
                                 ) {
  residuals <- x$F
  out <- list()
  for (i in seq_along(groups)) {
    out[[names(groups)[i]]] <- diag.avg(x, groups[[i]])
    residuals <- residuals - out[[names(groups)[i]]]
  }
  out[['Residuals']] <- residuals
  out
}