source("mcmssa_utils.R")
source("noise_estim.R")
source("auto_trend_ssa.R")
library(foreach)
library(doSNOW)
library(parallel)
library(doRNG)


N <- 100
Ls <- c(10, 20, 50, 80, 90)
D <- 1
G <- 1000
M <- 1000
pb <- txtProgressBar(max = M, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Пример 1 (phi=0.7, omega=0.075)
varphi <- 0.7
delta <- 1
model <- list(varphi = varphi,
              delta = delta,
              N = N)
omega.unknown <- 0.075
signal.unknown <- signal.one.channel(N, omega.unknown, A = 1)

# Случай, когда известный сигнал - периодика
omega.known <- 0.25
signal.known <- signal.one.channel(N, omega.known, A = 3)

# Сигнал и параметры шума известны точно
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
p.values_h0sin <- list()
for (k in seq_along(Ls)) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known)
      model.signal <- model
      model.signal$signal <- signal.known
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h0sin[[k]] <- result
}
stopCluster(cluster)


cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
p.values_h1sin <- list()
for (k in 1:5) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known + signal.unknown)
      model.signal <- model
      model.signal$signal <- signal.known
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h1sin[[k]] <- result
}
stopCluster(cluster)

# Параметры шума неизвестны, сигнал известен точно
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
p.values_h0sin_est_noise <- list()
for (k in 1:5) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known)
      model.signal <- est.model.arima(f - signal.known)
      model.signal$signal <- signal.known
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h0sin_est_noise[[k]] <- result
}
stopCluster(cluster)


cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
p.values_h1sin_est_noise <- list()
for (k in 1:5) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known + signal.unknown)
      model.signal <- est.model.arima(f - signal.known)
      model.signal$signal <- signal.known
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h1sin_est_noise[[k]] <- result
}
stopCluster(cluster)

# Неизвестен сигнал и параметры шума
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
p.values_h0sin_est_noise_signal <- list()
for (k in 1:5) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known)
      s <- ssa(f, kind = "toeplitz-ssa")
      freq <- numeric()
      for (i in 1:nu(s)) {
        ss <- ssa(s$U[,i], kind = "1d-ssa")
        p <- parestimate(ss, groups = list(1:2))
        freq[i] <- p$frequencies[[1]]
      }
      components <- order(abs(omega.known - freq))[1:2]
      r <- reconstruct(s, groups = list(components))
      model.signal <- est.model.arima(resid(r))
      model.signal$signal <- r[[1]]
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h0sin_est_noise_signal[[k]] <- result
}
stopCluster(cluster)


cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
p.values_h1sin_est_noise_signal <- list()
for (k in 1:5) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known + signal.unknown)
      s <- ssa(f, kind = "toeplitz-ssa")
      freq <- numeric()
      for (i in 1:nu(s)) {
        ss <- ssa(s$U[,i], kind = "1d-ssa")
        p <- parestimate(ss, groups = list(1:2))
        freq[i] <- p$frequencies[[1]]
      }
      components <- order(abs(omega.known - freq))[1:2]
      r <- reconstruct(s, groups = list(components))
      model.signal <- est.model.arima(resid(r))
      model.signal$signal <- r[[1]]
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h1sin_est_noise_signal[[k]] <- result
}
stopCluster(cluster)


# Случай, когда известный сигнал - тренд
signal.known <- 0.2 * exp(0.05 * (1:N))

# Сигнал и параметры шума известны точно
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
p.values_h0trend <- list()
for (k in seq_along(Ls)) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct', 'grouping.auto'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known)
      model.signal <- model
      model.signal$signal <- signal.known
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h0trend[[k]] <- result
}
stopCluster(cluster)


cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
p.values_h1trend <- list()
for (k in seq_along(Ls)) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct', 'grouping.auto'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known + signal.unknown)
      model.signal <- model
      model.signal$signal <- signal.known
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h1trend[[k]] <- result
}
stopCluster(cluster)


# Параметры шума неизвестны, сигнал известен точно
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
p.values_h0trend_est_noise <- list()
for (k in seq_along(Ls)) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct', 'grouping.auto'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known)
      model.signal <- est.model.arima(f - signal.known)
      model.signal$signal <- signal.known
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h0trend_est_noise[[k]] <- result
}
stopCluster(cluster)


cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
p.values_h1trend_est_noise  <- list()
for (k in seq_along(Ls)) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct', 'grouping.auto'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known + signal.unknown)
      model.signal <- est.model.arima(f - signal.known)
      model.signal$signal <- signal.known
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h1trend_est_noise[[k]] <- result
}
stopCluster(cluster)


# Неизвестен сигнал и параметры шума
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
p.values_h0trend_est_noise_signal <- list()
for (k in seq_along(Ls)) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct', 'grouping.auto'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known)
      s <- ssa(f)
      auto <- auto_trend_model(s, method = "basic.auto", signal_rank = 1)
      model.signal <- est.model.arima(resid(auto))
      model.signal$signal <- auto$trend
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h0trend_est_noise_signal[[k]] <- result
}
stopCluster(cluster)


cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
p.values_h1trend_est_noise_signal <- list()
for (k in seq_along(Ls)) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct', 'grouping.auto'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known + signal.unknown)
      s <- ssa(f)
      auto <- auto_trend_model(s, method = "basic.auto", signal_rank = 1)
      model.signal <- est.model.arima(resid(auto))
      model.signal$signal <- auto$trend
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h1trend_est_noise_signal[[k]] <- result
}
stopCluster(cluster)

alphas <- c(seq(0, 0.025, 0.0001), seq(0.025, 1, 0.001))
alphas_idx <- seq_along(alphas)
clrs <- c('black', 'red', 'green', 'orange', 'purple')
lwds <- c(2, 1, 1, 1, 1)

alphaI_sin_new <- lapply(p.values_h0sin_new, function(p) sapply(alphas, function(a) sum(p < a) / M))
beta_sin <- lapply(p.values_h1sin, function(p) sapply(alphas, function(a) sum(p < a) / M))

alphaI_sin_est_noise_new <- lapply(p.values_h0sin_est_noise_new, function(p) sapply(alphas, function(a) sum(p < a) / M))
beta_sin_est_noise <- lapply(p.values_h1sin_est_noise, function(p) sapply(alphas, function(a) sum(p < a) / M))

alphaI_sin_est_noise_signal_new <- lapply(p.values_h0sin_est_noise_signal_new, function(p) sapply(alphas, function(a) sum(p < a) / M))
beta_sin_est_noise_signal <- lapply(p.values_h1sin_est_noise_signal, function(p) sapply(alphas, function(a) sum(p < a) / M))


alphaI_trend <- lapply(p.values_h0trend, function(p) sapply(alphas, function(a) sum(p < a) / M))
beta_trend <- lapply(p.values_h1trend, function(p) sapply(alphas, function(a) sum(p < a) / M))

alphaI_trend_est_noise <- lapply(p.values_h0trend_est_noise, function(p) sapply(alphas, function(a) sum(p < a) / M))
beta_trend_est_noise <- lapply(p.values_h1trend_est_noise, function(p) sapply(alphas, function(a) sum(p < a) / M))

alphaI_trend_est_noise_signal <- lapply(p.values_h0trend_est_noise_signal, function(p) sapply(alphas, function(a) sum(p < a) / M))
beta_trend_est_noise_signal <- lapply(p.values_h1trend_est_noise_signal, function(p) sapply(alphas, function(a) sum(p < a) / M))


pdf("../tex/img/type1error_sin.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0,1),c(0,1), type="l", col = "blue", lty = 2, xlab = 'significance level', ylab = 'type I error')
for (i in 1:5)
  lines(alphas, alphaI_sin[[i]], lwd = lwds[i], col = clrs[i])
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()

pdf("../tex/img/roc_sin.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0,1),c(0,1), type="l", col="blue", lty = 2, xlab = 'type I error', ylab = 'power')
for (i in seq_along(Ls)) {
  lines(alphaI_sin[[i]][-1], beta_sin[[i]][-1], lwd = lwds[i], col = clrs[i])
}
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()


pdf("../tex/img/type1error_sin_est_noise.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0,1),c(0,1), type="l", col = "blue", lty = 2, xlab = 'significance level', ylab = 'type I error')
for (i in 1:5)
  lines(alphas, alphaI_sin_est_noise[[i]], lwd = lwds[i], col = clrs[i])
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()

pdf("../tex/img/roc_sin_est_noise.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0,1),c(0,1), type="l", col="blue", lty = 2, xlab = 'type I error', ylab = 'power')
for (i in seq_along(Ls)) {
  lines(alphaI_sin_est_noise[[i]][-1], beta_sin_est_noise[[i]][-1], lwd = lwds[i], col = clrs[i])
}
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()


pdf("../tex/img/type1error_sin_est_noise_signal.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0,1),c(0,1), type="l", col = "blue", lty = 2, xlab = 'significance level', ylab = 'type I error')
for (i in 1:5)
  lines(alphas, alphaI_sin_est_noise_signal[[i]], lwd = lwds[i], col = clrs[i])
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()

pdf("../tex/img/roc_sin_est_noise_signal.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0,1),c(0,1), type="l", col="blue", lty = 2, xlab = 'type I error', ylab = 'power')
for (i in seq_along(Ls)) {
  lines(alphaI_sin_est_noise_signal[[i]][-1], beta_sin_est_noise_signal[[i]][-1], lwd = lwds[i], col = clrs[i])
}
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()


pdf("../tex/img/type1error_trend.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0,1),c(0,1), type="l", col = "blue", lty = 2, xlab = 'significance level', ylab = 'type I error')
for (i in 1:5)
  lines(alphas, alphaI_trend[[i]], lwd = lwds[i], col = clrs[i])
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()

pdf("../tex/img/roc_trend.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0,1),c(0,1), type="l", col="blue", lty = 2, xlab = 'type I error', ylab = 'power')
for (i in seq_along(Ls)) {
  lines(alphaI_trend[[i]][-1], beta_trend[[i]][-1], lwd = lwds[i], col = clrs[i])
}
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()


pdf("../tex/img/type1error_trend_est_noise.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0,1),c(0,1), type="l", col = "blue", lty = 2, xlab = 'significance level', ylab = 'type I error')
for (i in 1:5)
  lines(alphas, alphaI_trend_est_noise[[i]], lwd = lwds[i], col = clrs[i])
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()

pdf("../tex/img/roc_trend_est_noise.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0,1),c(0,1), type="l", col="blue", lty = 2, xlab = 'type I error', ylab = 'power')
for (i in seq_along(Ls)) {
  lines(alphaI_trend_est_noise[[i]][-1], beta_trend_est_noise[[i]][-1], lwd = lwds[i], col = clrs[i])
}
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()

pdf("../tex/img/type1error_trend_est_noise_signal.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0,1),c(0,1), type="l", col = "blue", lty = 2, xlab = 'significance level', ylab = 'type I error')
for (i in 1:5)
  lines(alphas, alphaI_trend_est_noise_signal[[i]], lwd = lwds[i], col = clrs[i])
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()

pdf("../tex/img/roc_trend_est_noise_signal.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0,1),c(0,1), type="l", col="blue", lty = 2, xlab = 'type I error', ylab = 'power')
for (i in seq_along(Ls)) {
  lines(alphaI_trend_est_noise_signal[[i]][-1], beta_trend_est_noise_signal[[i]][-1], lwd = lwds[i], col = clrs[i])
}
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()

legend <- c("Exact", "Est (noise)", "Est (noise + signal)")
pdf("type1error_sin_copm.pdf", width = 15, height = 5, bg = "white")
plot(c(0,1),c(0,1), type="l", col="blue", lty = 2, xlab = 'significance level', ylab = 'type I error', cex.lab=1.8, cex.axis=1.8, cex.sub=1.8)
lines(alphas, alphaI_sin[[5]], col = "red")
lines(alphas, alphaI_sin_est_noise[[5]], col = "purple")
lines(alphas, alphaI_sin_est_noise_signal[[5]], col = "green")
legend(x = "bottomright", legend, col = c("red", "purple", "green"), lty = 1, lwd = 1, cex = 2)
dev.off()

legend <- c("Exact", "Est (noise)", "Est (noise + signal)")
pdf("../tex/img/roc_sin_copm.pdf", width = 15, height = 5, bg = "white")
plot(c(0,1),c(0,1), type="l", col="blue", lty = 2, xlab = 'type I error', ylab = 'power', cex.lab=1.8, cex.axis=1.8, cex.sub=1.8)
lines(alphaI_sin[[5]], beta_sin[[5]], col = "red")
lines(alphaI_sin_est_noise[[5]], beta_sin_est_noise[[5]], col = "purple")
lines(alphaI_sin_est_noise_signal[[5]], beta_sin_est_noise_signal[[5]], col = "green")
legend(x = "bottomright", legend, col = c("red", "purple", "green"), lty = 1, lwd = 1, cex = 2)
dev.off()

legend <- c("Exact", "Est (noise)", "Est (noise + signal)")
pdf("type1error_trend_copm.pdf", width = 15, height = 5, bg = "white")
plot(c(0,1),c(0,1), type="l", col="blue", lty = 2, xlab = 'significance level', ylab = 'type I error', cex.lab=1.8, cex.axis=1.8, cex.sub=1.8)
lines(alphas, alphaI_trend[[5]], col = "red")
lines(alphas, alphaI_trend_est_noise[[5]], col = "purple")
lines(alphas, alphaI_trend_est_noise_signal[[5]], col = "green")
legend(x = "bottomright", legend, col = c("red", "purple", "green"), lty = 1, lwd = 1, cex = 2)
dev.off()

legend <- c("Exact", "Est (noise)", "Est (noise + signal)")
pdf("../tex/img/roc_trend_copm.pdf", width = 15, height = 5, bg = "white")
plot(c(0,1),c(0,1), type="l", col="blue", lty = 2, xlab = 'type I error', ylab = 'power', cex.lab=1.8, cex.axis=1.8, cex.sub=1.8)
lines(alphaI_trend[[5]], beta_trend[[5]], col = "red")
lines(alphaI_trend_est_noise[[5]], beta_trend_est_noise[[5]], col = "purple")
lines(alphaI_trend_est_noise_signal[[5]], beta_trend_est_noise_signal[[5]], col = "green")
legend(x = "bottomright", legend, col = c("red", "purple", "green"), lty = 1, lwd = 1, cex = 2)
dev.off()

# Пример 2 (phi=0.3, omega=0.075)
varphi <- 0.3
delta <- 1
D <- 1
G <- 1000
M <- 1000
model <- list(varphi = varphi,
              delta = delta,
              N = N)

omega.unknown <- 0.075
signal.unknown <- 0.7 * cos(2 * pi * (1:N) * omega.unknown)

# Случай, когда известный сигнал - периодика
omega.known <- 0.25
signal.known <- signal.one.channel(N, omega.known, A = 3)

# Сигнал и параметры шума известны точно
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
p.values_h0sin_phi3 <- list()
for (k in seq_along(Ls)) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known)
      model.signal <- model
      model.signal$signal <- signal.known
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h0sin_phi3[[k]] <- result
}
stopCluster(cluster)


cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
p.values_h1sin_phi3 <- list()
for (k in 1:5) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known + signal.unknown)
      model.signal <- model
      model.signal$signal <- signal.known
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h1sin_phi3[[k]] <- result
}
stopCluster(cluster)


# Параметры шума неизвестны, сигнал известен точно
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
p.values_h0sin_phi3est <- list()
for (k in 1:5) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known)
      model.signal <- est.model.arima(f - signal.known)
      model.signal$signal <- signal.known
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h0sin_phi3est[[k]] <- result
}
stopCluster(cluster)


cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
p.values_h1sin_phi3est <- list()
for (k in 1:5) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known + signal.unknown)
      model.signal <- est.model.arima(f - signal.known)
      model.signal$signal <- signal.known
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h1sin_phi3est[[k]] <- result
}
stopCluster(cluster)


# Неизвестен сигнал и параметры шума
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
p.values_h0sin_phi3est_signal <- list()
for (k in 1:5) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known)
      s <- ssa(f, kind = "toeplitz-ssa")
      freq <- numeric()
      for (i in 1:nu(s)) {
        ss <- ssa(s$U[,i], kind = "1d-ssa")
        p <- parestimate(ss, groups = list(1:2))
        freq[i] <- p$frequencies[[1]]
      }
      components <- order(abs(omega.known - freq))[1:2]
      r <- reconstruct(s, groups = list(components))
      model.signal <- est.model.arima(resid(r))
      model.signal$signal <- r[[1]]
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h0sin_phi3est_signal[[k]] <- result
}
stopCluster(cluster)


cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
p.values_h1sin_phi3est_signal <- list()
for (k in 1:5) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known + signal.unknown)
      s <- ssa(f, kind = "toeplitz-ssa")
      freq <- numeric()
      for (i in 1:nu(s)) {
        ss <- ssa(s$U[,i], kind = "1d-ssa")
        p <- parestimate(ss, groups = list(1:2))
        freq[i] <- p$frequencies[[1]]
      }
      components <- order(abs(omega.known - freq))[1:2]
      r <- reconstruct(s, groups = list(components))
      model.signal <- est.model.arima(resid(r))
      model.signal$signal <- r[[1]]
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h1sin_phi3est_signal[[k]] <- result
}
stopCluster(cluster)


# Случай, когда известный сигнал - тренд
signal.known <- 0.2 * exp(0.05 * (1:N))

# Сигнал и параметры шума известны точно
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
p.values_h0trend_phi3 <- list()
for (k in seq_along(Ls)) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct', 'grouping.auto'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known)
      model.signal <- model
      model.signal$signal <- signal.known
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h0trend_phi3[[k]] <- result
}
stopCluster(cluster)


cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
p.values_h1trend_phi3 <- list()
for (k in seq_along(Ls)) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct', 'grouping.auto'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known + signal.unknown)
      model.signal <- model
      model.signal$signal <- signal.known
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h1trend_phi3[[k]] <- result
}
stopCluster(cluster)


# Параметры шума неизвестны, сигнал известен точно
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
p.values_h0trend_phi3est <- list()
for (k in seq_along(Ls)) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct', 'grouping.auto'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known)
      model.signal <- est.model.arima(f - signal.known)
      model.signal$signal <- signal.known
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h0trend_phi3est[[k]] <- result
}
stopCluster(cluster)


cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
p.values_h1trend_phi3est  <- list()
for (k in seq_along(Ls)) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct', 'grouping.auto'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known + signal.unknown)
      model.signal <- est.model.arima(f - signal.known)
      model.signal$signal <- signal.known
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h1trend_phi3est[[k]] <- result
}
stopCluster(cluster)


# Неизвестен сигнал и параметры шума
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
p.values_h0trend_phi3est_signal <- list()
for (k in seq_along(Ls)) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct', 'grouping.auto'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known)
      s <- ssa(f)
      #auto <- auto_trend_model(s, method = "basic.auto", signal_rank = 1)
      r <- reconstruct(s, list(trend = 1))
      model.signal <- est.model.arima(resid(r))
      model.signal$signal <- r$trend
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h0trend_phi3est_signal[[k]] <- result
}
stopCluster(cluster)


cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
p.values_h1trend_phi3est_signal <- list()
for (k in seq_along(Ls)) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct', 'grouping.auto'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known + signal.unknown)
      s <- ssa(f)
      #auto <- auto_trend_model(s, method = "basic.auto", signal_rank = 1)
      r <- reconstruct(s, list(trend = 1))
      model.signal <- est.model.arima(resid(r))
      model.signal$signal <- r$trend
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h1trend_phi3est_signal[[k]] <- result
}
stopCluster(cluster)


alphaI_sin_phi3 <- lapply(p.values_h0sin_phi3, function(p) sapply(alphas, function(a) sum(p < a) / M))
beta_sin_phi3 <- lapply(p.values_h1sin_phi3, function(p) sapply(alphas, function(a) sum(p < a) / M))

alphaI_sin_phi3est <- lapply(p.values_h0sin_phi3est, function(p) sapply(alphas, function(a) sum(p < a) / M))
beta_sin_phi3est <- lapply(p.values_h1sin_phi3est, function(p) sapply(alphas, function(a) sum(p < a) / M))

alphaI_sin_phi3est_signal <- lapply(p.values_h0sin_phi3est_signal, function(p) sapply(alphas, function(a) sum(p < a) / M))
beta_sin_phi3est_signal <- lapply(p.values_h1sin_phi3est_signal, function(p) sapply(alphas, function(a) sum(p < a) / M))


alphaI_trend_phi3 <- lapply(p.values_h0trend_phi3, function(p) sapply(alphas, function(a) sum(p < a) / M))
beta_trend_phi3 <- lapply(p.values_h1trend_phi3, function(p) sapply(alphas, function(a) sum(p < a) / M))

alphaI_trend_phi3est <- lapply(p.values_h0trend_phi3est, function(p) sapply(alphas, function(a) sum(p < a) / M))
beta_trend_phi3est <- lapply(p.values_h1trend_phi3est, function(p) sapply(alphas, function(a) sum(p < a) / M))

alphaI_trend_phi3est_signal <- lapply(p.values_h0trend_phi3est_signal, function(p) sapply(alphas, function(a) sum(p < a) / M))
beta_trend_phi3est_signal <- lapply(p.values_h1trend_phi3est_signal, function(p) sapply(alphas, function(a) sum(p < a) / M))


pdf("../tex/img/type1error_sin_phi3est_signal.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0,1),c(0,1), type="l", col = "blue", lty = 2, xlab = 'significance level', ylab = 'type I error')
for (i in 1:5)
  lines(alphas, alphaI_sin_phi3est_signal[[i]], lwd = lwds[i], col = clrs[i])
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()

pdf("../tex/img/roc_sin_phi3est_signal.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0,1),c(0,1), type="l", col="blue", lty = 2, xlab = 'type I error', ylab = 'power')
for (i in seq_along(Ls)) {
  lines(alphaI_sin_phi3est_signal[[i]][-1], beta_sin_phi3est_signal[[i]][-1], lwd = lwds[i], col = clrs[i])
}
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()


pdf("../tex/img/type1error_trend_phi3est_signal.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0,1),c(0,1), type="l", col = "blue", lty = 2, xlab = 'significance level', ylab = 'type I error')
for (i in 1:5)
  lines(alphas, alphaI_trend_phi3est_signal[[i]], lwd = lwds[i], col = clrs[i])
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()

pdf("../tex/img/roc_trend_phi3est_signal.pdf", width = 6, height = 3.5, bg = "white")
plot(c(0,1),c(0,1), type="l", col="blue", lty = 2, xlab = 'type I error', ylab = 'power')
for (i in seq_along(Ls)) {
  lines(alphaI_trend_phi3est_signal[[i]][-1], beta_trend_phi3est_signal[[i]][-1], lwd = lwds[i], col = clrs[i])
}
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
dev.off()


# Пример 3 (phi=0.7, omega=0.225)
varphi <- 0.7
delta <- 1
model <- list(varphi = varphi,
              delta = delta,
              N = N)
omega.unknown <- 0.225
signal.unknown <- 0.4 * cos(2 * pi * (1:N) * omega.unknown)

# Случай, когда известный сигнал - периодика
omega.known <- 0.25
signal.known <- signal.one.channel(N, omega.known, A = 3)

# Сигнал и параметры шума известны точно
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
p.values_h1sin_phi7_omega0225 <- list()
for (k in 1:5) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known + signal.unknown)
      model.signal <- model
      model.signal$signal <- signal.known
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h1sin_phi7_omega0225[[k]] <- result
}
stopCluster(cluster)


# Параметры шума неизвестны, сигнал известен точно
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
p.values_h1sin_phi7est_omega0225 <- list()
for (k in 1:5) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known + signal.unknown)
      model.signal <- est.model.arima(f - signal.known)
      model.signal$signal <- signal.known
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h1sin_phi7est_omega0225[[k]] <- result
}
stopCluster(cluster)


# Неизвестен сигнал и параметры шума
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
p.values_h1sin_phi7est_signal_omega0225 <- list()
for (k in 1:5) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known + signal.unknown)
      s <- ssa(f, kind = "toeplitz-ssa")
      freq <- numeric()
      for (i in 1:nu(s)) {
        ss <- ssa(s$U[,i], kind = "1d-ssa")
        p <- parestimate(ss, groups = list(1:2))
        freq[i] <- p$frequencies[[1]]
      }
      components <- order(abs(omega.known - freq))[1:2]
      r <- reconstruct(s, groups = list(components))
      model.signal <- est.model.arima(resid(r))
      model.signal$signal <- r[[1]]
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h1sin_phi7est_signal_omega0225[[k]] <- result
}
stopCluster(cluster)


# Случай, когда известный сигнал - тренд
signal.known <- 0.2 * exp(0.05 * (1:N))

# Сигнал и параметры шума известны точно
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
p.values_h1trend_phi7_omega0225 <- list()
for (k in 1:5) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known + signal.unknown)
      model.signal <- model
      model.signal$signal <- signal.known
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h1trend_phi7_omega0225[[k]] <- result
}
stopCluster(cluster)


# Параметры шума неизвестны, сигнал известен точно
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
registerDoRNG(seed = 1, once = FALSE)
p.values_h1trend_phi7est_omega0225 <- list()
for (k in 1:5) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known + signal.unknown)
      model.signal <- est.model.arima(f - signal.known)
      model.signal$signal <- signal.known
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h1trend_phi7est_omega0225[[k]] <- result
}
stopCluster(cluster)


# Неизвестен сигнал и параметры шума
cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)
p.values_h1trend_phi7est_signal_omega0225 <- list()
for (k in 1:5) {
  result <- foreach (i=1:M,
                     .export=c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'reconstruct'),
                     .combine='cbind',
                     .options.snow=opts) %dopar%
    {
      f <- one.channel.ts(model, signal.known + signal.unknown)
      s <- ssa(f)
      #auto <- auto_trend_model(s, method = "basic.auto", signal_rank = 1)
      r <- reconstruct(s, list(trend = 1))
      model.signal <- est.model.arima(resid(r))
      model.signal$signal <- r$trend
      m <- MonteCarloSSA(f = f,
                         model = model.signal,
                         L = Ls[k],
                         D = 1,
                         basis = "ev",
                         kind = "ev",
                         G = G,
                         level.conf = NULL,
                         composite = TRUE)
      m$p.value
    }
  p.values_h1trend_phi7est_signal_omega0225[[k]] <- result
}
stopCluster(cluster)


beta_sin_phi7_omega0225 <- lapply(p.values_h1sin_phi7_omega0225, function(p) sapply(alphas, function(a) sum(p < a) / M))
beta_sin_phi7est_omega0225 <- lapply(p.values_h1sin_phi7est_omega0225, function(p) sapply(alphas, function(a) sum(p < a) / M))
beta_sin_phi7est_signal_omega0225 <- lapply(p.values_h1sin_phi7est_signal_omega0225, function(p) sapply(alphas, function(a) sum(p < a) / M))

beta_trend_phi7_omega0225 <- lapply(p.values_h1trend_phi7_omega0225, function(p) sapply(alphas, function(a) sum(p < a) / M))
beta_trend_phi7est_omega0225 <- lapply(p.values_h1trend_phi7est_omega0225, function(p) sapply(alphas, function(a) sum(p < a) / M))
beta_trend_phi7est_signal_omega0225 <- lapply(p.values_h1trend_phi7est_signal_omega0225, function(p) sapply(alphas, function(a) sum(p < a) / M))


pdf("../tex/img/roc_sin_phi7est_signal_omega0225.pdf", width = 15, height = 5, bg = "white")
plot(c(0,1),c(0,1), type="l", col="blue", lty = 2, xlab = 'type I error', ylab = 'power', cex.lab=1.8, cex.axis=1.8, cex.sub=1.8)
for (i in seq_along(Ls)) {
  lines(alphaI_sin_est_noise_signal_new[[i]][-1], beta_sin_phi7est_signal_omega0225_new[[i]][-1], lwd = lwds[i], col = clrs[i])
}
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds, cex = 2)
dev.off()

pdf("../tex/img/roc_trend_phi7est_signal_omega0225.pdf", width = 15, height = 5, bg = "white")
plot(c(0,1),c(0,1), type="l", col="blue", lty = 2, xlab = 'type I error', ylab = 'power', cex.lab=1.8, cex.axis=1.8, cex.sub=1.8)
for (i in seq_along(Ls)) {
  lines(alphaI_trend_est_noise_signal[[i]][-1], beta_trend_phi7est_signal_omega0225[[i]][-1], lwd = lwds[i], col = clrs[i])
}
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds, cex = 2)
dev.off()
