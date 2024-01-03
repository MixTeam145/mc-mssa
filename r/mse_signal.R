library(Rssa)
source("toeplitz_mssa.R")

N <- 71
sigma <- 5
Ls <- c(12, 24, 36, 48, 60)
signal1 <- 30 * cos(2*pi * (1:N) / 12)
signal2 <- 30 * cos(2*pi * (1:N) / 12 + pi / 4)
signal <- cbind(signal1, signal2)
R <- 10000
n <- 2

ssa.errors <- function(Ls) {
  f <- signal1 + rnorm(N, sd = sigma)
  err.rec <- numeric(length(Ls));
  names(err.rec) <- Ls
  for (l in seq_along(Ls)) {
    L <- Ls[l]
    s <- ssa(f, L = L, kind = "1d-ssa")
    rec <- reconstruct(s, groups = list(1:n))[[1]]
    err.rec[l] <- mean((rec - signal1)^2)
  }
  err.rec
}

toeplssa.errors <- function(Ls) {
  f <- signal1 + rnorm(N, sd = sigma)
  err.rec <- numeric(length(Ls));
  names(err.rec) <- Ls
  for (l in seq_along(Ls)) {
    L <- Ls[l]
    s <- ssa(f, L = L, kind = "toeplitz-ssa")
    rec <- reconstruct(s, groups = list(1:n))[[1]]
    err.rec[l] <- mean((rec - signal1)^2)
  }
  err.rec
}

mssa.errors <- function(Ls) {
  f1 <- signal1 + rnorm(N, sd = sigma)
  f2 <- signal2 + rnorm(N, sd = sigma)
  f <- cbind(f1, f2)
  err.rec <- numeric(length(Ls))
  names(err.rec) <- Ls
  for (l in seq_along(Ls)) {
    L <- Ls[l]
    s <- ssa(f, L = L, kind = "mssa")
    rec <- reconstruct(s, groups = list(1:n))[[1]]
    err.rec[l] <- mean((rec - signal)^2)
  }
  err.rec
}

toeplSum.errors <- function(Ls) {
  f1 <- signal1 + rnorm(N, sd = sigma)
  f2 <- signal2 + rnorm(N, sd = sigma)
  f <- cbind(f1, f2)
  err.rec <- numeric(length(Ls))
  names(err.rec) <- Ls
  for (l in seq_along(Ls)) {
    L <- Ls[l]
    s <- toeplitz.mssa(f, L = L, D = 2, method = "sum")
    rec <- toeplitz.reconstruct(s, groups = list(signal = 1:n))[[1]]
    err.rec[l] <- mean((rec - signal)^2)
  }
  err.rec
} 

toeplBlock.errors <- function(Ls) {
  f1 <- signal1 + rnorm(N, sd = sigma)
  f2 <- signal2 + rnorm(N, sd = sigma)
  f <- cbind(f1, f2)
  err.rec <- numeric(length(Ls))
  names(err.rec) <- Ls
  for (l in seq_along(Ls)) {
    L <- Ls[l]
    s <- toeplitz.mssa(f, L = L, D = 2, method = "block")
    rec <- toeplitz.reconstruct(s, groups = list(signal = 1:n))[[1]]
    err.rec[l] <- mean((rec - signal)^2)
  }
  err.rec
} 

mres.ssa <- replicate(R, ssa.errors(Ls))
rowMeans(mres.ssa)


mres.toeplssa <- replicate(R, toeplssa.errors(Ls))
rowMeans(mres.toeplssa)


mres.mssa <- replicate(R, mssa.errors(Ls))
rowMeans(mres.mssa)


mres.toeplSum <- replicate(R, toeplSum.errors(Ls))
rowMeans(mres.toeplSum)


mres.toeplBlock <- replicate(R, toeplBlock.errors(Ls))
rowMeans(mres.toeplBlock)


signal1 <- 1.2 * (1:N)
signal2 <- 0.8 * (1:N)
signal <- cbind(signal1, signal2)
n <- 2

mres.ssa <- replicate(R, ssa.errors(Ls))
rowMeans(mres.ssa)


mres.toeplssa <- replicate(R, toeplssa.errors(Ls))
rowMeans(mres.toeplssa)


mres.mssa <- replicate(R, mssa.errors(Ls))
rowMeans(mres.mssa)


mres.toeplSum <- replicate(R, toeplSum.errors(Ls))
rowMeans(mres.toeplSum)


mres.toeplBlock <- replicate(R, toeplBlock.errors(Ls))
rowMeans(mres.toeplBlock)
