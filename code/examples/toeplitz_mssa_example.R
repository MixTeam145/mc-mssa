# Rssa (ttps://CRAN.R-project.org/package=Rssa): only Basic MSSA is implemented
# Contribution of Egor Poteshkin: Toeplitz MSSA in two versions is implemented
# "Sum" version of Toeplitz MSSA is preferable

source("R/toeplitz_mssa.R")

N <- 100 # time series length
L <- 70 # window length
D <- 2 # number of channels

signal1 <- 2 * cos(2 * pi * (1:N) / 10)
signal2 <- 3 * cos(2* pi * (1:N) / 10)
signal <- cbind(signal1, signal2)

set.seed(5)
f <- signal + cbind(rnorm(N), rnorm(N))

plot.ts(f)

s <- toeplitz.mssa(f,
                   L = L,
                   D = D,
                   method = "sum", # or "block"
                   neig = 10 # desired non-zero eigenvalues
                   )

plot(s$sigma, type = "b") # plot eigenvalues

i <- 1
plot(s$U[, i], type = "l") # plot i-th eigenvector

j <- 2
plot(s$U[, i], s$U[, j], type = "l") # plot 2-d eigenvectors

r <- toeplitz.reconstruct(s, groups = list(season = 1:2))

plot.ts(r$season, main = "Reconstructed series")

plot.ts(r$Residuals, main = "Residuals")
