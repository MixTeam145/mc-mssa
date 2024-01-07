library(Rssa)
library(xtable)
source("toeplitz_mssa.R")

ssa.errors <- function(Ls) {
  f <- signal1 + rnorm(N, sd = sigma)
  err.rec <- numeric(length(Ls));
  names(err.rec) <- Ls
  for (l in seq_along(Ls)) {
    L <- Ls[l]
    s <- ssa(f, L = L, kind = "1d-ssa")
    rec <- reconstruct(s, groups = list(1:r1))[[1]]
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
    rec <- reconstruct(s, groups = list(1:r1))[[1]]
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
    rec <- reconstruct(s, groups = list(1:r2))[[1]]
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
    rec <- toeplitz.reconstruct(s, groups = list(signal = 1:r2))[[1]]
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
    rec <- toeplitz.reconstruct(s, groups = list(signal = 1:r2))[[1]]
    err.rec[l] <- mean((rec - signal)^2)
  }
  err.rec
} 

# Initual setup
N <- 71
sigma <- 5
Ls <- c(12, 24, 36, 48, 60)
R <- 1000

# 1. Same frequency case
signal1 <- 30 * cos(2*pi * (1:N) / 12)
signal2 <- 20 * cos(2*pi * (1:N) / 12)
signal <- cbind(signal1, signal2)
r1 <- 2
r2 <- 2

set.seed(42)
mres.ssa1 <- replicate(R, ssa.errors(Ls))
set.seed(42)
mres.toeplssa1 <- replicate(R, toeplssa.errors(Ls))
set.seed(42)
mres.mssa1 <- replicate(R, mssa.errors(Ls))
set.seed(42)
mres.toeplSum1 <- replicate(R, toeplSum.errors(Ls))
set.seed(42)
mres.toeplBlock1 <- replicate(R, toeplBlock.errors(Ls))

table1 <- matrix(nrow = 5, ncol = 5)
table1[1,] <- rowMeans(mres.ssa1)
table1[2,] <- rowMeans(mres.toeplssa1)
table1[3,] <- rowMeans(mres.mssa1)
table1[4,] <- rowMeans(mres.toeplSum1)
table1[5,] <- rowMeans(mres.toeplBlock1)
rownames(table1) <- c("ssa", "toeplitz ssa", "mssa", "toeplitz mssa, sum", "toeplitz mssa, block")
colnames(table1) <- c("12", "24", "36", "48", "60")
xtable(table1)


# 2. Different frequency case
signal1 <- 30 * cos(2*pi * (1:N) / 12)
signal2 <- 20 * cos(2*pi * (1:N) / 8)
signal <- cbind(signal1, signal2)
r1 <- 2
r2 <- 2

set.seed(42)
mres.ssa2 <- replicate(R, ssa.errors(Ls))
set.seed(42)
mres.toeplssa2 <- replicate(R, toeplssa.errors(Ls))
set.seed(42)
mres.mssa2 <- replicate(R, mssa.errors(Ls))
set.seed(42)
mres.toeplSum2 <- replicate(R, toeplSum.errors(Ls))
set.seed(42)
mres.toeplBlock2 <- replicate(R, toeplBlock.errors(Ls))

table2 <- matrix(nrow = 5, ncol = 5)
table2[1,] <- rowMeans(mres.ssa2)
table2[2,] <- rowMeans(mres.toeplssa2)
table2[3,] <- rowMeans(mres.mssa2)
table2[4,] <- rowMeans(mres.toeplSum2)
table2[5,] <- rowMeans(mres.toeplBlock2)
rownames(table2) <- c("ssa", "toeplitz ssa", "mssa", "toeplitz mssa, sum", "toeplitz mssa, block")
colnames(table2) <- c("12", "24", "36", "48", "60")
xtable(table2)


# 3. Nonstationary case
signal1 <- 1.2 * (1:N)
signal2 <- 0.8 * (1:N)
signal <- cbind(signal1, signal2)
r1 <- 2
r2 <- 2

set.seed(42)
mres.ssa3 <- replicate(R, ssa.errors(Ls))
set.seed(42)
mres.toeplssa3 <- replicate(R, toeplssa.errors(Ls))
set.seed(42)
mres.mssa3 <- replicate(R, mssa.errors(Ls))
set.seed(42)
mres.toeplSum3 <- replicate(R, toeplSum.errors(Ls))
set.seed(42)
mres.toeplBlock3 <- replicate(R, toeplBlock.errors(Ls))

table3 <- matrix(nrow = 5, ncol = 5)
table3[1,] <- rowMeans(mres.ssa3)
table3[2,] <- rowMeans(mres.toeplssa3)
table3[3,] <- rowMeans(mres.mssa3)
table3[4,] <- rowMeans(mres.toeplSum3)
table3[5,] <- rowMeans(mres.toeplBlock3)
rownames(table3) <- c("ssa", "toeplitz ssa", "mssa", "toeplitz mssa, sum", "toeplitz mssa, block")
colnames(table3) <- c("12", "24", "36", "48", "60")
xtable(table3)
