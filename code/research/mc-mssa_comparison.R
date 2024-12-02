library(foreach)
library(doParallel)
library(doRNG)
library(gcplyr)
source("toeplitz_mssa.R")
source("mcmssa_utils.R")

N <- 100
Ls <- c(10, 20, 50, 80, 90)
L_idx <- 1:length(Ls)
D <- 2
G <- 1000
M <- 1000
delta <- 1

# phi = 0.7
varphi <- 0.7
model <- list(list(varphi = varphi,
                   delta = delta,
                   N = N),
              list(varphi = varphi,
                   delta = delta,
                   N = N))

# H_0
signal <- replicate(D, 0, simplify=F)

p.values_noise.me1block.ev <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "ev",
        toeplitz.kind = "block",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_noise.me1block.ev[[idx]] <- result
}
stopCluster(cluster)

p.values_noise.me1block.fa <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "fa",
        toeplitz.kind = "block",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_noise.me1block.fa[[idx]] <- result
}
stopCluster(cluster)

p.values_noise.me1sum.ev <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "ev",
        toeplitz.kind = "sum",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_noise.me1sum.ev[[idx]] <- result
}
stopCluster(cluster)

p.values_noise.me1sum.fa <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "fa",
        toeplitz.kind = "sum",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_noise.me1sum.fa[[idx]] <- result
}
stopCluster(cluster)

p.values_noise.me1.ev <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "ev",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_noise.me1.ev[[idx]] <- result
}
stopCluster(cluster)

p.values_noise.me1.fa <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "fa",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_noise.me1.fa[[idx]] <- result
}
stopCluster(cluster)

# H_1, omega = 0.075
omega <- 0.075
signal <- list(signal.one.channel(model[[1]]$N, omega), 
               signal.one.channel(model[[1]]$N, omega))

p.values_signal.me1block.ev_omega0075 <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "ev",
        toeplitz.kind = "block",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_signal.me1block.ev_omega0075[[idx]] <- result
}
stopCluster(cluster)

p.values_signal.me1block.fa_omega0075 <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "fa",
        toeplitz.kind = "block",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_signal.me1block.fa_omega0075[[idx]] <- result
}
stopCluster(cluster)

p.values_signal.me1sum.ev_omega0075 <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "ev",
        toeplitz.kind = "sum",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_signal.me1sum.ev_omega0075[[idx]] <- result
}
stopCluster(cluster)

p.values_signal.me1sum.fa_omega0075 <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "fa",
        toeplitz.kind = "sum",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_signal.me1sum.fa_omega0075[[idx]] <- result
}
stopCluster(cluster)

# H_1, omega = 0.225
omega <- 0.225
signal <- list(signal.one.channel(model[[1]]$N, omega, A = 0.4), 
               signal.one.channel(model[[1]]$N, omega, A = 0.4))

p.values_signal.me1block.ev_omega0225 <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "ev",
        toeplitz.kind = "block",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_signal.me1block.ev_omega0225[[idx]] <- result
}
stopCluster(cluster)

p.values_signal.me1block.fa_omega0225 <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "fa",
        toeplitz.kind = "block",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_signal.me1block.fa_omega0225[[idx]] <- result
}
stopCluster(cluster)

p.values_signal.me1sum.ev_omega0225 <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "ev",
        toeplitz.kind = "sum",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_signal.me1sum.ev_omega0225[[idx]] <- result
}
stopCluster(cluster)

p.values_signal.me1sum.fa_omega0225 <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "fa",
        toeplitz.kind = "sum",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_signal.me1sum.fa_omega0225[[idx]] <- result
}
stopCluster(cluster)

p.values_signal.me1.ev_omega0225 <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "ev",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_signal.me1.ev_omega0225[[idx]] <- result
}
stopCluster(cluster)

p.values_signal.me1.fa_omega0225 <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "fa",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_signal.me1.fa_omega0225[[idx]] <- result
}
stopCluster(cluster)


# phi=0.3
varphi <- 0.3
model <- list(list(varphi = varphi,
                   delta = delta,
                   N = N),
              list(varphi = varphi,
                   delta = delta,
                   N = N))

# H_0
signal <- replicate(D, 0, simplify=F)

p.values_noise.me1block.ev_phi3 <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "ev",
        toeplitz.kind = "block",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_noise.me1block.ev_phi3[[idx]] <- result
}
stopCluster(cluster)

p.values_noise.me1block.fa_phi3 <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "fa",
        toeplitz.kind = "block",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_noise.me1block.fa_phi3[[idx]] <- result
}
stopCluster(cluster)

p.values_noise.me1sum.ev_phi3 <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "ev",
        toeplitz.kind = "sum",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_noise.me1sum.ev_phi3[[idx]] <- result
}
stopCluster(cluster)

p.values_noise.me1sum.fa_phi3 <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "fa",
        toeplitz.kind = "sum",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_noise.me1sum.fa_phi3[[idx]] <- result
}
stopCluster(cluster)

p.values_noise.me1.ev_phi3 <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "ev",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_noise.me1.ev_phi3[[idx]] <- result
}
stopCluster(cluster)

p.values_noise.me1.fa_phi3 <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "fa",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_noise.me1.fa_phi3[[idx]] <- result
}
stopCluster(cluster)

# H_1, omega = 0.075
omega <- 0.075
signal <- list(signal.one.channel(model[[1]]$N, omega, A = 0.5), 
               signal.one.channel(model[[1]]$N, omega, A = 0.5))
p.values_signal.me1block.ev_phi3_omega0075 <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "ev",
        toeplitz.kind = "block",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_signal.me1block.ev_phi3_omega0075[[idx]] <- result
}
stopCluster(cluster)

p.values_signal.me1block.fa_phi3_omega0075 <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "fa",
        toeplitz.kind = "block",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_signal.me1block.fa_phi3_omega0075[[idx]] <- result
}
stopCluster(cluster)

p.values_signal.me1sum.ev_phi3_omega0075 <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "ev",
        toeplitz.kind = "sum",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_signal.me1sum.ev_phi3_omega0075[[idx]] <- result
}
stopCluster(cluster)

p.values_signal.me1sum.fa_phi3_omega0075 <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "fa",
        toeplitz.kind = "sum",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_signal.me1sum.fa_phi3_omega0075[[idx]] <- result
}
stopCluster(cluster)

p.values_signal.me1.ev_phi3_omega0075 <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "ev",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_signal.me1.ev_phi3_omega0075[[idx]] <- result
}
stopCluster(cluster)

p.values_signal.me1.fa_phi3_omega0075 <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "fa",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_signal.me1.fa_phi3_omega0075[[idx]] <- result
}
stopCluster(cluster)


alphas <- 0:1000/1000
alphas_idx <- 1:length(alphas)
clrs <- c('black', 'red', 'green', 'orange', 'purple')
lwds <- c(2, 1, 1, 1, 1)


alphaI.me1.ev <- lapply(p.values_noise.me1.ev, function(pvals) sapply(alphas, function(a) mean(pvals < a)))
alphaI.me1.fa <- lapply(p.values_noise.me1.fa, function(pvals) sapply(alphas, function(a) mean(pvals < a)))

beta.me1.ev_omega0075 <- lapply(p.values_signal.me1.ev_omega0075, function(pvals) sapply(alphas, function(a) mean(pvals < a)))
beta.me1.fa_omega0075 <- lapply(p.values_signal.me1.fa_omega0075, function(pvals) sapply(alphas, function(a) mean(pvals < a)))


pdf("../tex/img/roc_mssa_ev.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue', xlab = 'type I error', ylab = 'power')
for (l in L_idx)
  lines(alphaI.me1.ev[[l]][-1], beta.me1.ev_omega0075[[l]][-1], col = clrs[l], type = 'l', lwd = lwds[l])
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()


pdf("../tex/img/type1error_mssa_fa.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue', xlab = 'significance level', ylab = 'type I error')
for (l in L_idx)
  lines(alphas, alphaI.me1.ev[[l]], col = clrs[l], lwd = lwds[l])
lines(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue')
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()


pdf("../tex/img/roc_mssa_fa.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue', xlab = 'type I error', ylab = 'power')
for (l in L_idx)
  lines(alphaI.me1.fa[[l]][-1], beta.me1.fa_omega0075[[l]][-1], col = clrs[l], type = 'l', lwd = lwds[l])
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()


pdf("../tex/img/type1error_mssa_fa.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue', xlab = 'significance level', ylab = 'type I error')
for (l in L_idx)
  lines(alphas, alphaI.me1.fa[[l]], col = clrs[l], lwd = lwds[l])
lines(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue')
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()


roc.me1sum.fa <- list()
alpha_1.me1sum.fa <- list()
beta.me1sum.fa <- list()
for (l in L_idx)
{
  roc.me1sum.fa[[as.character(Ls[l])]] <- data.frame(fpr=numeric(length(alphas)), tpr=numeric(length(alphas)), alpha=numeric(length(alphas)))
  alpha_1.me1sum.fa[[as.character(Ls[l])]] <- data.frame(alpha=numeric(length(alphas)), fpr=numeric(length(alphas)))
  beta.me1sum.fa[[as.character(Ls[l])]] <- data.frame(alpha=numeric(length(alphas)), tpr=numeric(length(alphas)))
  
  for (i in alphas_idx) 
  {
    alpha <- alphas[i]
    
    roc.me1sum.fa[[as.character(Ls[l])]]$fpr[i] <- sum(p.values_noise.me1sum.fa[[l]] < alpha)/M
    alpha_1.me1sum.fa[[as.character(Ls[l])]]$fpr[i] <- sum(p.values_noise.me1sum.fa[[l]] < alpha)/M
    roc.me1sum.fa[[as.character(Ls[l])]]$tpr[i] <- sum(p.values_signal.me1sum.fa_omega0075[[l]] < alpha)/M
    beta.me1sum.fa[[as.character(Ls[l])]]$tpr[i] <- sum(p.values_signal.me1sum.fa_omega0075[[l]] < alpha)/M
    roc.me1sum.fa[[as.character(Ls[l])]]$alpha[i] <- alpha
    alpha_1.me1sum.fa[[as.character(Ls[l])]]$alpha[i] <- alpha
    beta.me1sum.fa[[as.character(Ls[l])]]$alpha[i] <- alpha
  }
}

pdf("../tex/img/roc_sum_fa_omega0075.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue', xlab = 'type I error', ylab = 'power')
for (l in L_idx)
  lines(roc.me1sum.fa[[l]]$fpr[-1], roc.me1sum.fa[[l]]$tpr[-1], col = clrs[l], type = 'l', lwd = lwds[l])
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()

pdf("../tex/img/type1error_sum_fa.pdf", width = 6, height = 3.5, bg = "white")
plot(alpha_1.me1sum.fa[[1]]$alpha, alpha_1.me1sum.fa[[1]]$fpr, type = 'l', lwd = lwds[1], xlab = 'significance level', ylab = 'type I error')
for (l in L_idx[-1])
  lines(alpha_1.me1sum.fa[[l]]$alpha, alpha_1.me1sum.fa[[l]]$fpr, col = clrs[l], lwd = lwds[l])
lines(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue')
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()

roc.me1sum.ev <- list()
alpha_1.me1sum.ev <- list()
beta.me1sum.ev <- list()
for (l in L_idx)
{
  roc.me1sum.ev[[as.character(Ls[l])]] <- data.frame(fpr=numeric(length(alphas)), tpr=numeric(length(alphas)), alpha=numeric(length(alphas)))
  alpha_1.me1sum.ev[[as.character(Ls[l])]] <- data.frame(alpha=numeric(length(alphas)), fpr=numeric(length(alphas)))
  beta.me1sum.ev[[as.character(Ls[l])]] <- data.frame(alpha=numeric(length(alphas)), tpr=numeric(length(alphas)))
  
  for (i in alphas_idx) 
  {
    alpha <- alphas[i]
    
    roc.me1sum.ev[[as.character(Ls[l])]]$fpr[i] <- sum(p.values_noise.me1sum.ev[[l]] < alpha)/M
    alpha_1.me1sum.ev[[as.character(Ls[l])]]$fpr[i] <- sum(p.values_noise.me1sum.ev[[l]] < alpha)/M
    roc.me1sum.ev[[as.character(Ls[l])]]$tpr[i] <- sum(p.values_signal.me1sum.ev_omega0075[[l]] < alpha)/M
    beta.me1sum.ev[[as.character(Ls[l])]]$tpr[i] <- sum(p.values_signal.me1sum.ev_omega0075[[l]] < alpha)/M
    roc.me1sum.ev[[as.character(Ls[l])]]$alpha[i] <- alpha
    alpha_1.me1sum.ev[[as.character(Ls[l])]]$alpha[i] <- alpha
    beta.me1sum.ev[[as.character(Ls[l])]]$alpha[i] <- alpha
  }
}

pdf("../tex/img/roc_sum_ev_omega0075.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue', xlab = 'type I error', ylab = 'power')
for (l in L_idx)
  lines(roc.me1sum.ev[[l]]$fpr[-1], roc.me1sum.ev[[l]]$tpr[-1], col = clrs[l], type = 'l', lwd = lwds[l])
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()

pdf("../tex/img/type1error_sum_ev.pdf", width = 6, height = 3.5, bg = "white")
plot(alpha_1.me1sum.ev[[1]]$alpha, alpha_1.me1sum.ev[[1]]$fpr, type = 'l', lwd = lwds[1], xlab = 'significance level', ylab = 'type I error')
for (l in L_idx[-1])
  lines(alpha_1.me1sum.ev[[l]]$alpha, alpha_1.me1sum.ev[[l]]$fpr, col = clrs[l], lwd = lwds[l])
lines(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue')
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()


roc.me1block.ev <- list()
alpha_1.me1block.ev <- list()
beta.me1block.ev <- list()
for (l in L_idx)
{
  roc.me1block.ev[[as.character(Ls[l])]] <- data.frame(fpr=numeric(length(alphas)), tpr=numeric(length(alphas)), alpha=numeric(length(alphas)))
  alpha_1.me1block.ev[[as.character(Ls[l])]] <- data.frame(alpha=numeric(length(alphas)), fpr=numeric(length(alphas)))
  beta.me1block.ev[[as.character(Ls[l])]] <- data.frame(alpha=numeric(length(alphas)), tpr=numeric(length(alphas)))
  
  for (i in alphas_idx) 
  {
    alpha <- alphas[i]
    
    roc.me1block.ev[[as.character(Ls[l])]]$fpr[i] <- sum(p.values_noise.me1block.ev[[l]] < alpha)/M
    alpha_1.me1block.ev[[as.character(Ls[l])]]$fpr[i] <- sum(p.values_noise.me1block.ev[[l]] < alpha)/M
    roc.me1block.ev[[as.character(Ls[l])]]$tpr[i] <- sum(p.values_signal.me1block.ev_omega0075[[l]] < alpha)/M
    beta.me1block.ev[[as.character(Ls[l])]]$tpr[i] <- sum(p.values_signal.me1block.ev_omega0075[[l]] < alpha)/M
    roc.me1block.ev[[as.character(Ls[l])]]$alpha[i] <- alpha
    alpha_1.me1block.ev[[as.character(Ls[l])]]$alpha[i] <- alpha
    beta.me1block.ev[[as.character(Ls[l])]]$alpha[i] <- alpha
  }
}

pdf("../tex/img/roc_block_ev_omega0075.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue', xlab = 'type I error', ylab = 'power')
for (l in L_idx)
  lines(roc.me1block.ev[[l]]$fpr[-1], roc.me1block.ev[[l]]$tpr[-1], col = clrs[l], type = 'l', lwd = lwds[l])
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()

pdf("../tex/img/type1error_block_ev.pdf", width = 6, height = 3.5, bg = "white")
plot(alpha_1.me1block.ev[[1]]$alpha, alpha_1.me1block.ev[[1]]$fpr, type = 'l', lwd = lwds[1], xlab = 'significance level', ylab = 'type I error')
for (l in L_idx[-1])
  lines(alpha_1.me1block.ev[[l]]$alpha, alpha_1.me1block.ev[[l]]$fpr, col = clrs[l], lwd = lwds[l])
lines(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue')
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()


roc.me1block.fa <- list()
alpha_1.me1block.fa <- list()
beta.me1block.fa <- list()
for (l in L_idx)
{
  roc.me1block.fa[[as.character(Ls[l])]] <- data.frame(fpr=numeric(length(alphas)), tpr=numeric(length(alphas)), alpha=numeric(length(alphas)))
  alpha_1.me1block.fa[[as.character(Ls[l])]] <- data.frame(alpha=numeric(length(alphas)), fpr=numeric(length(alphas)))
  beta.me1block.fa[[as.character(Ls[l])]] <- data.frame(alpha=numeric(length(alphas)), tpr=numeric(length(alphas)))
  
  for (i in alphas_idx) 
  {
    alpha <- alphas[i]
    
    roc.me1block.fa[[as.character(Ls[l])]]$fpr[i] <- sum(p.values_noise.me1block.fa[[l]] < alpha)/M
    alpha_1.me1block.fa[[as.character(Ls[l])]]$fpr[i] <- sum(p.values_noise.me1block.fa[[l]] < alpha)/M
    roc.me1block.fa[[as.character(Ls[l])]]$tpr[i] <- sum(p.values_signal.me1block.fa_omega0075[[l]] < alpha)/M
    beta.me1block.fa[[as.character(Ls[l])]]$tpr[i] <- sum(p.values_signal.me1block.fa_omega0075[[l]] < alpha)/M
    roc.me1block.fa[[as.character(Ls[l])]]$alpha[i] <- alpha
    alpha_1.me1block.fa[[as.character(Ls[l])]]$alpha[i] <- alpha
    beta.me1block.fa[[as.character(Ls[l])]]$alpha[i] <- alpha
  }
}

pdf("../tex/img/roc_block_fa_omega0075.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue', xlab = 'type I error', ylab = 'power')
for (l in L_idx)
  lines(roc.me1block.fa[[l]]$fpr[-1], roc.me1block.fa[[l]]$tpr[-1], col = clrs[l], type = 'l', lwd = lwds[l])
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()

pdf("../tex/img/type1error_block_fa.pdf", width = 6, height = 3.5, bg = "white")
plot(alpha_1.me1block.fa[[1]]$alpha, alpha_1.me1block.fa[[1]]$fpr, type = 'l', lwd = lwds[1], xlab = 'significance level', ylab = 'type I error')
for (l in L_idx[-1])
  lines(alpha_1.me1block.fa[[l]]$alpha, alpha_1.me1block.fa[[l]]$fpr, col = clrs[l], lwd = lwds[l])
lines(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue')
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()



# MSSA
signal <- replicate(D, 0, simplify=F)

p.values_noise.me1.ev_phi3 <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "ev",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_noise.me1.ev_phi3[[idx]] <- result
}
stopCluster(cluster)

p.values_noise.me1.fa_phi3 <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "fa",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_noise.me1.fa_phi3[[idx]] <- result
}
stopCluster(cluster)

omega <- 0.075
signal <- list(signal.one.channel(model[[1]]$N, omega, A = 0.5), 
               signal.one.channel(model[[1]]$N, omega, A = 0.5))
p.values_signal.me1.ev_phi3_omega0075 <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "ev",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_signal.me1.ev_phi3_omega0075[[idx]] <- result
}
stopCluster(cluster)

p.values_signal.me1.fa_phi3_omega0075 <- list()
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
for (idx in L_idx) {
  result <- foreach (
    i = 1:M,
    .combine = 'c',
    .export = c('Norm', 'rowQuantiles', 'Toeplitz', 'lag'),
    .packages = "Rssa",
    .options.snow = opts
  ) %dopar% {
    f <- generate(model, signal, D)
    res <-
      MonteCarloSSA(
        f = f,
        L = Ls[idx],
        model = model,
        basis = "ev",
        kind = "fa",
        D = D,
        G = G,
        level.conf = NULL
      )
    res$p.value
  }
  p.values_signal.me1.fa_phi3_omega0075[[idx]] <- result
}
stopCluster(cluster)

alphaI.me1sum.ev_phi3 <- lapply(p.values_noise.me1sum.ev_phi3, function(pvals) sapply(alphas, function(a) mean(pvals < a))) 
beta.me1sum.ev_phi3_omega0075 <- lapply(p.values_signal.me1sum.ev_phi3_omega0075, function(pvals) sapply(alphas, function(a) mean(pvals < a)))

alphaI.me1sum.fa_phi3 <- lapply(p.values_noise.me1sum.fa_phi3, function(pvals) sapply(alphas, function(a) mean(pvals < a))) 
beta.me1sum.fa_phi3_omega0075 <- lapply(p.values_signal.me1sum.fa_phi3_omega0075, function(pvals) sapply(alphas, function(a) mean(pvals < a)))

alphaI.me1block.ev_phi3 <- lapply(p.values_noise.me1block.ev_phi3, function(pvals) sapply(alphas, function(a) mean(pvals < a))) 
beta.me1block.ev_phi3_omega0075 <- lapply(p.values_signal.me1block.ev_phi3_omega0075, function(pvals) sapply(alphas, function(a) mean(pvals < a)))

alphaI.me1block.fa_phi3 <- lapply(p.values_noise.me1block.fa_phi3, function(pvals) sapply(alphas, function(a) mean(pvals < a))) 
beta.me1block.fa_phi3_omega0075 <- lapply(p.values_signal.me1block.fa_phi3_omega0075, function(pvals) sapply(alphas, function(a) mean(pvals < a)))


pdf("../tex/img/roc_sum_ev_phi3_omega0075.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue', xlab = 'type I error', ylab = 'power')
for (l in L_idx)
  lines(alphaI.me1sum.ev_phi3[[l]][-1], beta.me1sum.ev_phi3_omega0075[[l]][-1], col = clrs[l], type = 'l', lwd = lwds[l])
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()


pdf("../tex/img/type1error_sum_ev_phi3.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue', xlab = 'significance level', ylab = 'type I error')
for (l in L_idx)
  lines(alphas, alphaI.me1sum.ev_phi3[[l]], col = clrs[l], lwd = lwds[l])
lines(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue')
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()


pdf("../tex/img/roc_sum_fa_phi3_omega0075.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue', xlab = 'type I error', ylab = 'power')
for (l in L_idx)
  lines(alphaI.me1sum.fa_phi3[[l]][-1], beta.me1sum.fa_phi3_omega0075[[l]][-1], col = clrs[l], type = 'l', lwd = lwds[l])
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()


pdf("../tex/img/type1error_sum_fa_phi3.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue', xlab = 'significance level', ylab = 'type I error')
for (l in L_idx)
  lines(alphas, alphaI.me1sum.fa_phi3[[l]], col = clrs[l], lwd = lwds[l])
lines(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue')
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()


pdf("../tex/img/roc_block_ev_phi3_omega0075.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue', xlab = 'type I error', ylab = 'power')
for (l in L_idx)
  lines(alphaI.me1block.ev_phi3[[l]][-1], beta.me1block.ev_phi3_omega0075[[l]][-1], col = clrs[l], type = 'l', lwd = lwds[l])
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()


pdf("../tex/img/type1error_block_ev_phi3.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue', xlab = 'significance level', ylab = 'type I error')
for (l in L_idx)
  lines(alphas, alphaI.me1block.ev_phi3[[l]], col = clrs[l], lwd = lwds[l])
lines(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue')
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()


pdf("../tex/img/roc_block_fa_phi3_omega0075.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue', xlab = 'type I error', ylab = 'power')
for (l in L_idx)
  lines(alphaI.me1block.fa_phi3[[l]][-1], beta.me1block.fa_phi3_omega0075[[l]][-1], col = clrs[l], type = 'l', lwd = lwds[l])
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()


pdf("../tex/img/type1error_block_fa_phi3.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue', xlab = 'significance level', ylab = 'type I error')
for (l in L_idx)
  lines(alphas, alphaI.me1block.fa_phi3[[l]], col = clrs[l], lwd = lwds[l])
lines(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue')
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()


alphaI.me1sum.ev <- lapply(p.values_noise.me1sum.ev, function(pvals) sapply(alphas, function(a) mean(pvals < a)))
beta.me1sum.ev_omega0225 <- lapply(p.values_signal.me1sum.ev_omega0225, function(pvals) sapply(alphas, function(a) mean(pvals < a)))

alphaI.me1sum.fa <- lapply(p.values_noise.me1sum.fa, function(pvals) sapply(alphas, function(a) mean(pvals < a)))
beta.me1sum.fa_omega0225 <- lapply(p.values_signal.me1sum.fa_omega0225, function(pvals) sapply(alphas, function(a) mean(pvals < a)))

alphaI.me1block.ev <- lapply(p.values_noise.me1block.ev, function(pvals) sapply(alphas, function(a) mean(pvals < a)))
beta.me1block.ev_omega0225 <- lapply(p.values_signal.me1block.ev_omega0225, function(pvals) sapply(alphas, function(a) mean(pvals < a)))

alphaI.me1block.fa <- lapply(p.values_noise.me1block.fa, function(pvals) sapply(alphas, function(a) mean(pvals < a)))
beta.me1block.fa_omega0225 <- lapply(p.values_signal.me1block.fa_omega0225, function(pvals) sapply(alphas, function(a) mean(pvals < a)))

alphaI.me1.ev <- lapply(p.values_noise.me1.ev, function(pvals) sapply(alphas, function(a) mean(pvals < a)))
beta.me1.ev_omega0225 <- lapply(p.values_signal.me1.ev_omega0225, function(pvals) sapply(alphas, function(a) mean(pvals < a)))

alphaI.me1.fa <- lapply(p.values_noise.me1.fa, function(pvals) sapply(alphas, function(a) mean(pvals < a)))
beta.me1.fa_omega0225 <- lapply(p.values_signal.me1.fa_omega0225, function(pvals) sapply(alphas, function(a) mean(pvals < a)))


plot(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue', xlab = 'type I error', ylab = 'power', cex.lab=1.8, cex.axis=1.8, cex.sub=1.8)
for (l in L_idx)
  lines(alphaI.me1.ev[[l]][-1], beta.me1.ev_omega0225[[l]][-1], col = clrs[l], type = 'l', lwd = lwds[l])
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds, cex = 2)



pdf("../tex/img/roc_block_ev_omega0225.pdf", width = 15, height = 5, bg = "white")
plot(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue', xlab = 'type I error', ylab = 'power', cex.lab=1.8, cex.axis=1.8, cex.sub=1.8)
for (l in L_idx)
  lines(alphaI.me1block.ev[[l]][-1], beta.me1block.ev_omega0225[[l]][-1], col = clrs[l], type = 'l', lwd = lwds[l])
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds, cex = 2)
dev.off()

pdf("../tex/img/roc_block_fa_omega0225.pdf", width = 15, height = 5, bg = "white")
plot(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue', xlab = 'type I error', ylab = 'power', cex.lab=1.8, cex.axis=1.8, cex.sub=1.8)
for (l in L_idx)
  lines(alphaI.me1block.fa[[l]][-1], beta.me1block.fa_omega0225[[l]][-1], col = clrs[l], type = 'l', lwd = lwds[l])
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds, cex = 2)
dev.off()

pdf("../tex/img/roc_sum_ev_omega0225.pdf", width = 15, height = 5, bg = "white")
plot(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue', xlab = 'type I error', ylab = 'power', cex.lab=1.8, cex.axis=1.8, cex.sub=1.8)
for (l in L_idx)
  lines(alphaI.me1sum.ev[[l]][-1], beta.me1sum.ev_omega0225[[l]][-1], col = clrs[l], type = 'l', lwd = lwds[l])
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds, cex = 2)
dev.off()

pdf("../tex/img/roc_sum_fa_omega0225.pdf", width = 15, height = 5, bg = "white")
plot(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue', xlab = 'type I error', ylab = 'power', cex.lab=1.8, cex.axis=1.8, cex.sub=1.8)
for (l in L_idx)
  lines(alphaI.me1sum.fa[[l]][-1], beta.me1sum.fa_omega0225[[l]][-1], col = clrs[l], type = 'l', lwd = lwds[l])
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds, cex = 2)
dev.off()



alphaI.me1.ev_phi3 <- lapply(p.values_noise.me1.ev_phi3, function(pvals) sapply(alphas, function(a) mean(pvals < a)))
beta.me1.ev_phi3 <- lapply(p.values_signal.me1.ev_phi3_omega0075, function(pvals) sapply(alphas, function(a) mean(pvals < a)))

alphaI.me1.fa_phi3 <- lapply(p.values_noise.me1.fa_phi3, function(pvals) sapply(alphas, function(a) mean(pvals < a)))
beta.me1.fa_phi3 <- lapply(p.values_signal.me1.fa_phi3_omega0075, function(pvals) sapply(alphas, function(a) mean(pvals < a)))

plot(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue', xlab = 'type I error', ylab = 'power', cex.lab=1.8, cex.axis=1.8, cex.sub=1.8)
for (l in L_idx)
  lines(alphaI.me1.fa_phi3[[l]][-1], beta.me1.fa_phi3[[l]][-1], col = clrs[l], type = 'l', lwd = lwds[l])
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)


beta.me1.ev_omega0225 <- lapply(p.values_signal.me1.ev_omega0225, function(pvals) sapply(alphas, function(a) mean(pvals < a)))
beta.me1.fa_omega0225 <- lapply(p.values_signal.me1.fa_omega0225, function(pvals) sapply(alphas, function(a) mean(pvals < a)))


plot(c(0, 1), c(0, 1), type='l', lty = 2, col = 'blue', xlab = 'type I error', ylab = 'power', cex.lab=1.8, cex.axis=1.8, cex.sub=1.8)
for (l in L_idx)
  lines(alphaI.me1.fa_phi3[[l]][-1], beta.me1.fa_phi3[[l]][-1], col = clrs[l], type = 'l', lwd = lwds[l])
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds, cex = 2)
