library("Rssa")
library("pracma")
library("ggplot2")
library("matrixStats")
library("magic")
library("stats")
library(LSTS)
# source("toeplitz_mssa.R")
#source("auto_trend_ssa.R")

type = 8

#generates a series according to the model signal + AR(1)
one.channel.ts <- function(model, signal) {
  if (model$varphi == 0)
    xi = rnorm(model$N, sd = model$delta)
  else
    xi <- arima.sim(
      n = model$N,
      list(ar = model$varphi),
      sd = model$delta,
      n.start = 1,
      start.innov =
        rnorm(1, sd = model$delta / sqrt(1 - model$varphi^2))
    )
  if (!is.null(model$signal))
    xi <- xi + model$signal
  f <- xi + signal
  as.vector(f)
}

# Estimation of AR parameters, without signal extraction
# est.model.arima <-  function(f) {
#   param <- list()
#   ar.par <- arima(
#     f,
#     order = c(1, 0, 0),
#     include.mean = FALSE,
#     method = "ML"
#   )
#   param$varphi <- coef(ar.par)
#   if (sqrt(ar.par$var.coef) > abs(ar.par$coef))
#     param$varphi <- 0
#   param$delta <- sqrt(ar.par$sigma2)
#   #print(param)
#   estModel <-
#     list(varphi = param$varphi,
#          delta = param$delta,
#          N = length(f))
#   estModel
# }

est.model.arima <-  function(f, 
                             freq.range=c(0, 0.5),
                             method="mss") {
  if (method == "pg"){
    spec_values <- spectrum(f, taper=0, detrend=FALSE)$spec
    sigma2 <- mean(spec_values[ceiling(freq.range[1] * length(f)):length(spec_values)])
  } else if (method == "mss") {
    sigma2 <- mean(f ^ 2)
  }
  delta <- sqrt(sigma2)
  #print(param)
  estModel <-
    list(varphi = 0,
         delta = delta,
         N = length(f))
  estModel
}
###end

est.model.arima.red <- function(f,
                                freq.range=c(0, 0.5),
                                mode="MLE",
                                opt_init_phi=0.5){
  # opt_fun <- function(params, f, omega0) {
  #   fixed_freq <- omega0 * 2 * pi
  #   phi <- params[1]
  #   sigma2 <- params[2]
  #   
  #   if (abs(phi) >= 1 || sigma2 <= 0) {
  #     return(1e10)
  #   }
  #   
  #   n <- length(f)
  #   freq <- seq(0, pi, length.out = floor(n/2) + 1)
  #   S_y <- sigma2 / (1 + phi^2 - 2 * phi * cos(freq))
  #   S_y_filtered <- ifelse(freq >= fixed_freq, S_y, 0)
  #   I <- abs(fft(f))^2 / n
  #   I <- I[1:(floor(n/2) + 1)]
  #   positive_indices <- which(S_y_filtered > 0)
  #   I <- I[positive_indices]
  #   S_y_filtered <- S_y_filtered[positive_indices]
  #   log_likelihood <- sum(log(S_y_filtered) + I / S_y_filtered)
  #   
  #   return(log_likelihood)
  # }
  
  lsts_periodogram <- function(y, include.taper = FALSE) {
    series <- y - mean(y, na.rm = TRUE)
    N <- sum(is.na(series))
    series[is.na(series)] <- 0
    n <- length(series)
    if (include.taper == TRUE) {
      a <- 0.5 * (1 - cos(2 * pi * (0:(n - 1)) / n))
      series <- series * a
    }
    aux <- Mod(fft(series))^2
    m <- n / 2
    periodogram <- (aux[2:(m + 1)]) / (2 * pi * (n - N))
    if (include.taper == TRUE) {
      periodogram <- (aux[2:(m + 1)]) / (3 * pi * (n - N) / 4)
    }
    lambda <- (2 * pi * (1:m)) / n
    return(list(periodogram = periodogram, lambda = lambda))
  }
  
  opt_fun_fixed <- function(x, series, order = c(p = 1, q = 0), omega_0=-1,
                            ar.order = NULL, ma.order = NULL, sd.order = NULL,
                            d.order = NULL, include.d = FALSE, N = NULL,
                            S = NULL, include.taper = TRUE) {
    y <- series
    T. <- length(y) # Длина исходного ряда
    
    if (is.null(N)) {
      N <- trunc(T.^0.8) # Длина окна для подсчета периодограмм
    }
    if (is.null(S)) {
      S <- trunc(0.2 * N) # Сдвиг между блоками
    }
    
    M <- trunc((T. - N) / S + 1) # Количество блоков
    p <- c(0, 0)
    
    if (length(x) != sum(p + 1)) {
      stop("error in the number of parameters")
    } else {
      lik <- 0
      for (j in 1:M) {
        u <- (N / 2 + S * (j - 1)) / T.
        aux <- lsts_periodogram(y[(1 + S * (j - 1)):(N + S * (j - 1))], include.taper = TRUE)
        I <- aux$periodogram
        
        phi <- numeric()
        k <- 1
        phi[1] <- x[1]
        
        theta <- numeric()
        
        d <- 0
        sigma <- x[2]
        
        f <- spectral.density(ar = phi, ma = theta, d = d, sd = sigma, lambda = aux$lambda) # Здесь lambda --- это частоты 2 pi k / n, k = 1, ... , n / 2
        mask <- aux$lambda > 2 * pi * omega_0
        f <- f[mask]
        I <- I[mask]
        
        lik <- sum(log(f) + I / f) / N + lik
      }
      
      lik <- lik / M
      
      if (is.na(sigma) | sigma <= 0) {
        lik <- 1e10
      }
      
      lik
    }
  }
  loglik <- function(x, omega_0) opt_fun_fixed(c(x[1],x[2]), f, order = c(1,0), omega_0=omega_0)

  if (mode == "MLE"){
    init_params <- c(phi = opt_init_phi, sigma2 = var(f))
    
    result <- optim(
      par = init_params,
      fn = loglik,
      omega_0 = freq.range[1],
      method = "L-BFGS-B",
      lower = c(-0.999, 1e-6),
      upper = c(0.999, Inf)
    )
    
    estimated_phi <- result$par[1]
    estimated_sigma2 <- result$par[2] ^ 2
    
    estModel <-
      list(varphi = estimated_phi,
           delta = sqrt(estimated_sigma2),
           N = length(f))
  } else if (mode == "arima"){
    param <- list()
      ar.par <- arima(
        f,
        order = c(1, 0, 0),
        include.mean = FALSE,
        method = "ML"
      )
      param$varphi <- coef(ar.par)
      if (sqrt(ar.par$var.coef) > abs(ar.par$coef))
        param$varphi <- 0
      param$delta <- sqrt(ar.par$sigma2)
      estModel <-
        list(varphi = param$varphi,
             delta = param$delta,
             N = length(f))
  } else if (mode == "real"){
    estModel <-
      list(varphi = 0.7,
           delta = 1,
           N = length(f))
  }
  estModel
}

###Functions for Monte Carlo SSA
# Computes squared norms of projections to column vectors of U
projec <- function(data, L, D, U, kind=c("ev", "fa")) {
  if (is.list(data)) {
    # data are given by a model
    f <- generate(data, replicate(D, 0, simplify=F), D)
  } else {
    # data are given by a series
    f <- data
  }
  if (D == 1)
    N <- length(f)
  else
    N <- length(f[,1]) # assert equal length in each channel
  K <- N - L + 1
  X_res <- matrix(0, nrow = L, ncol = K*D)
  if (D == 1) {
    X_res <- t(sapply(1:L, function(i) f[i:(i+K-1)]))
  }
  else
    for (channel in 1:D) {
      tX <- sapply(1:L, function(i) f[i:(i+K-1),channel])
      X_res[,(1+(channel - 1)*K):(channel*K)] <- t(tX)
    }
  if (kind=='fa') {
    W <- X_res %*% U #Projection
  }
  else {
    W <- t(X_res) %*% U #Projection
  }
  colSums(W ^ 2 / N) #divide by N to weaken the dependence on t.s. length
}

# Generate vectors for projections corresponding to eigenvectors produced by t.s.
basis.ev <- function(ts, L, factor.v = T, toeplitz.kind = "no") {
  D <- dim(ts)[2]
  neig <- min(L, D * (length(ts[,1]) - L + 1))
  if (D == 1)
      s <- ssa(ts, L = L, neig = neig, kind = "toeplitz-ssa")
  else
    if (toeplitz.kind == "sum" || toeplitz.kind == "block")
      s <- toeplitz.mssa(ts, L = L, D = D, method = toeplitz.kind)
    else
      s <- ssa(ts, L = L, svd.method = "svd", neig = neig, kind = "mssa")
  freq <- numeric(0)
  for (i in 1:nu(s)){
    #ss <- ssa(s$U[,i], kind = "toeplitz-ssa")
    ss <- ssa(s$U[,i], kind = "1d-ssa")
    #estimation of the main frequency by ESPRIT
    p <- parestimate(ss, groups = list(1:2))
    freq[i] <- p$frequencies[[1]]
  }
  if (factor.v) {
    return(list(U = s$V, freq = freq)) # [,1:min(nu(s), neig)]
  }
  else {
    return(list(U = s$U, freq = freq)) # [,1:min(nu(s), neig)]
  }
}


# Generate vectors for projections as sine and cosine vectors
basis.sin0 <- function(L, freq.range=c(0, 0.5)) {
  numb <- L
  U <- matrix(0, nrow = L, ncol = numb)
  idx <- 1:L
  
  separat <- (freq.range[2] - freq.range[1]) / (L-0.5)
  from <- freq.range[1] + separat / 2
  freq <-
    seq(from = from,
        to = 0.5,
        by = separat) # Grid of frequencies

  for (i in 1:length(freq)) {
    U[, i] <- cos(2 * pi * freq[i] * idx)
    U[, i] <- U[, i] / Norm(U[, i])
  }  
  list(U = U, freq = freq)
}

basis.sin1 <- function(L) {
  numb <- L
  U <- matrix(0, nrow = L, ncol = numb)
  idx <- 1:L
  separat <- 1 / (2 * L)
  from <- 1 / (2 * L)
  freq <-
    seq(from = separat,
        to = 0.5,
        by = separat) # Grid of frequencies
  for (i in 1:length(freq)) {
    U[, i] <- sin(2 * pi * freq[i] * idx)
    U[, i] <- U[, i] / Norm(U[, i])
  }
  list(U = U, freq = freq)
}

basis.toeplitz.sin <- function(phi, L, omega=0.2, Amp=1) {
  toepl <- matrix(0, L, L)
  numb <- L
  for (i in 1:L) {
    for (j in 1:L) {
      toepl[i, j] <- phi ^ (abs(i - j))+ 0.5*Amp*cos(2*pi*omega*abs(i - j))
    }
  }
  s <- svd(toepl, L)
  
  U <- s$u
  freq <- numeric(0)
  for (i in 1:L) {
    #ss <- ssa(s$U[,i], kind = "toeplitz-ssa")
    ss <- ssa(U[, i], kind = "1d-ssa")
    #estimation of the main frequency by ESPRIT
    p <- parestimate(ss, groups = list(1:2))
    freq[i] <- p$frequencies[[1]]
  }
  list(U = U, freq = freq)
}
###end

# Generate vectors for projections corresponding to eigenvectors teoretical matrix 
matrix.toeplitz <- function(phi, L) {
  toepl <- matrix(0, L, L)
  numb <- L
  for (i in 1:L) {
    for (j in 1:L) {
      toepl[i, j] <- phi ^ (abs(i - j))
    }
  }
  toepl
}

basis.toeplitz.m <- function(phis, L, D, fa=F) {
  if (fa) {
    # here we assume that L param represents K = N - L + 1
    toepl.array <- list() 
    for (channel in 1:D) {
      toepl.array[[channel]] <- matrix.toeplitz(phis[[channel]][[1]], L)
    }
    toepl <- do.call("adiag", toepl.array)
    s <- svd(toepl, nv=L)
    U <- s$v
  }
  else {
    toepl <- matrix(data=0, nrow=L, ncol=L)
    for (channel in 1:D) {
      toepl <- toepl + matrix.toeplitz(phis[[channel]][[1]], L)
    }
    s <- svd(toepl, L)
    U <- s$u
  }
  
  freq <- numeric(0)
  for (i in 1:L) {
    #ss <- ssa(s$U[,i], kind = "toeplitz-ssa")
    ss <- ssa(U[, i], kind = "1d-ssa")
    #estimation of the main frequency by ESPRIT
    p <- parestimate(ss, groups = list(1:2))
    freq[i] <- p$frequencies[[1]]
  }
  list(U = U, freq = freq)
}


###end

what.reject <- function(res){
  rej <- (res$v[res$idx] < res$lower | res$v[res$idx] > res$upper) & res$idx[res$idx]
  print(res$freq[res$idx][rej==TRUE])
}

###Main functions for multiple Monte Carlo SSA
# Make multiple test
do.ci <-
  function(f,
           plan, 
           kind=c("ev", "fa"),
           model,
           level.conf,
           L,
           G,
           D,
           two.tailed = FALSE,
           composite,
           transf = function(x) {
             return(x)
           },
           inv.transf = function(x) {
             return(x)
           },
           weights = 1) {
    P <- replicate(G, projec(data = model, L = L, D = D, U = plan$U, kind=kind))
    v <- projec(data = f, L = L, D = D, U = plan$U, kind=kind)
    
    idx <- plan$freq >=  plan$range[1] & plan$freq <= plan$range[2]
    if (!(TRUE %in% idx))
      warning("no vectors with given frequency range")
    X <- transf(P[idx, , drop = FALSE])
    x <- transf(v[idx, drop = FALSE])
    
    if (is.vector(X))
      dim(X) <- c(length(X), 1)
    
    res <- list()
    
    res$freq <- plan$freq[idx]
    
    ci <- list()
    ci$means <- apply(X, 1, mean)
    
    ci$sds <- apply(X, 1, sd)
    
    #print(c("sd", ci$sds))
    if (weights[1] == "equal")
      weights <- ci$sds
    ci$sds[weights == 0] <- 1000 #something large
    ci$sds[weights != 0] <- ci$sds[weights != 0] / weights[weights != 0]
    #print(c("weight", weights))
    #print(c("sd", ci$sds))
    stats.max <-
      apply(X, 2, function(vv)
        max((vv - ci$means) / ci$sds))
    stats.max.abs <-
      apply(X, 2, function(vv)
        max(abs(vv - ci$means) / ci$sds))
    #print(stats.max)
    if (!is.null(level.conf))
    {
      if (two.tailed == FALSE) {
        ci$q.upper <- quantile(stats.max, probs = level.conf, type = type)
        ci$q.lower <- 0
      } else{
        ci$q.upper <-
          quantile(stats.max.abs, probs = level.conf, type = type)
        ci$q.lower <- -ci$q.upper
      }
      
      stat.v.max <- max((x - ci$means) / ci$sds)
      if (two.tailed == TRUE)
        stat.v.max <- max(abs(x - ci$means) / ci$sds)
      if (two.tailed == TRUE)
        res$reject <-
        as.logical(stat.v.max > ci$q.upper | stat.v.max < ci$q.lower)
      if (two.tailed == FALSE)
        res$reject <- as.logical(stat.v.max > ci$q.upper)
      res$freq.max <- NA
      #res$freq <- plan$freq
      if (res$reject == TRUE)
        res$freq.max <-
        plan$freq[idx][which.max((x - ci$means) / ci$sds)]
      res$freq.max
      res$upper <- inv.transf(ci$means + ci$q.upper * ci$sds)
      res$lower <- 0
      if (two.tailed == TRUE)
        res$lower <- inv.transf(ci$means + ci$q.lower * ci$sds)
        # res$lower <- inv.transf(max(0, ci$means + ci$q.lower * ci$sds))
      res$plan <- plan
      res$v <- x
      res$f <- f
      res$idx <- idx
      res
    }
    else 
    {
      if (two.tailed == FALSE)
      {
        stat.v.max <- max((x - ci$means) / ci$sds)
        F_ <- ecdf(stats.max)
        res$p.value <- 1 - F_(stat.v.max)
        res
      }
      else {
        stat.v.max <- max(abs(x - ci$means) / ci$sds)
        F_ <- ecdf(stats.max)
        res$p.value <- 1 - F_(stat.v.max)
        res
      }
    }
  }

#make single test
do.ci.single <-
  function(f,
           plan,
           model,
           level.conf,
           L,
           G,
           two.tailed = FALSE,
           transf = function(x) {
             return(x)
           },
           inv.transf = function(x) {
             return(x)
           }) {
    P <- replicate(G, projec(data = model, L = L, U = plan$U))
    v <- projec(data = f, L = L, U = plan$U)
    
    idx <- plan$freq >=  plan$range[1] & plan$freq <= plan$range[2]
    if (!(TRUE %in% idx))
      warning("no vectors with given frequency range")
    X <- transf(P[idx, , drop = FALSE])
    x <- transf(v[idx, drop = FALSE])
    #print(plan$freq[idx])
    
    if (is.vector(X))
      dim(X) <- c(length(X), 1)
    
    res <- list()
    ci <- list()
    if (!is.null(level.conf))
    {
      if (two.tailed == FALSE) {
        ci$q.upper <- rowQuantiles(X, probs = level.conf, type = type)
        ci$q.lower <- 0
      } else{
        ci$q.upper <- rowQuantiles(X, probs = (1 + level.conf) / 2, type = type)
        ci$q.lower <-
          rowQuantiles(X, probs = (1 - level.conf) / 2, type = type)
      }
      
      #stat.v.max <- max(x-ci$q.upper)
      #if(two.tailed == FALSE) stat.v.min <- 0 else stat.v.min <- min(x)
      if (two.tailed == FALSE)
        res$reject <- TRUE %in% as.logical(x > ci$q.upper)
      if (two.tailed == TRUE)
        res$reject <- TRUE %in% as.logical(x > ci$q.upper | x < ci$q.lower)
      res$freq.max <- NA
      res$freq <- plan$freq
      if (res$reject)
        res$freq.max <-
        plan$freq[idx][which.max((x - ci$q.upper) / ci$q.upper)]
      res$freq.max
      res$upper <- inv.transf(ci$q.upper)
      res$lower <- 0
      if (two.tailed == TRUE)
        res$lower <- inv.transf(ci$q.lower)
      res$plan <- plan
      res$v <- v
      res$f <- f
      res$idx <- idx
      res
    }
    else
    {
      if (two.tailed == FALSE) {
        p.values <- numeric(length(x))
        for(i.inner in 1:length(x))
        {
          F_ <- ecdf(X[i.inner,])
          p.values[i.inner] <- 1 - F_(x[i.inner])
        }
        res$p.value <- min(p.values)
        res
      }
      else {
        p.values <- numeric(length(x))
        for(i.inner in 1:length(x))
        {
           F_ <- ecdf(X[i.inner,])
           p.values[i.inner] <- (1 - F_(x[i.inner])) * 2
           if (p.values[i.inner] > 1) {
              p.values[i.inner] <- F_(x[i.inner]) * 2
           }
        }
        res$p.value <- min(p.values)
        res
      }
    }
  }

#plot by dominant frequency
plot.ci <- function(res, log_ = FALSE) {
  df <-
    data.frame(
      frequency = res$freq,
      contribution = res$v,
      ci.lower = res$lower,
      ci.upper = res$upper
    )
  df_reject <- df |> filter(contribution < ci.lower | contribution > ci.upper)
  p <-
    ggplot(df, aes(frequency, contribution)) +
    geom_point() +
    geom_point(data = df_reject, aes(x=frequency, y=contribution), color = "red") +
    geom_errorbar(aes(ymin = res$lower, ymax = res$upper), color = "blue") +
    theme_bw()
  if (log_)
    p <- p + scale_y_sqrt()
  p
  # v <- res$v
  # freq <- res$freq
  # idx <- res$idx
  # v <- res$v
  # sp <- spec.ar(res$f, order = 1, plot = FALSE)
  # if (is.null(lim))
  
  #  lim = c(0, max(res$upper, v))
  # plot(
  #   sp$spec ~ sp$freq,
  #   type = "l",
  #   ylim = lim,
  #   xlab = "frequency",
  #   ylab = "contribution"
  # ) #spectral density of AR(1)
  # if(log_)
  # {
  #   print('log')
  #   log.lower <- log(res$lower)
  #   log.upper <- log(res$upper)
  #   log.lower[res$lower <= 0] <- -6
  # 
  #   plot(
  #     log(v[idx]) ~ freq[idx],
  #     ylim = c(min(log.lower, v), max(log.upper, v)),
  #     xlab = "frequency",
  #     ylab = "log(contribution)"
  #   )
  #   
  #   segments(freq[idx], log.lower, freq[idx], log.upper, col = "red")
  # segments(freq[idx], -5, freq[idx], log(res$upper), col = "red")
  # }
  #prediction intervals
  # lines(v[idx] ~ freq[idx], type = "p") #Squared projection norm for the original time series
  # else 
  # {
    # 
    # plot(
    #   v ~ freq,
    #   ylim = c(min(res$lower, v), max(res$upper, v)),
    #   xlab = "frequency",
    #   ylab = "contribution"
    # )
    # 
    # segments(freq, res$lower, freq, res$upper, col = "red")
  # }
}


#plot by dominant contribution (the projection norm)
plot.ci.by.order <- function(res, lim = NULL) {
  v <- res$v
  freq <- res$freq
  num <- 1:length(res$freq)
  idx <- res$idx
  v <- res$v
  idx.order <- idx[order(v[idx], decreasing = TRUE)]
  v.order <- v[order(v[idx], decreasing = TRUE)]
  lower.order <- 0
  if (length(res$lower) > 1)
    lower.order <- res$lower[order(v[num[idx]], decreasing = TRUE)]
  upper.order <- res$upper[order(v[idx], decreasing = TRUE)]
  num.order <- num[order(v[idx], decreasing = TRUE)]
  #sp <- spec.ar(res$f, order = 1, plot = FALSE)
  #lim <-c(min(res$lower,v[idx]), max(res$upper,v[idx]))
  if (is.null(lim))
    lim = c(0, max(res$upper, v))
  #plot(sp$spec*(N) ~ sp$freq, type = "l", ylim = lim) #spectral density of AR(1)
  plot(
    v.order[idx] ~ num[idx],
    type = "p",
    ylim = lim,
    xlab = "number by order",
    ylab = "contribution"
  ) #Squared projection norm for the original time series
  segments(num[idx], lower.order, num[idx], upper.order, col = "red")
  #prediction intervals
}

# The wrapped function for Multiple Monte Carlo SSA
MonteCarloSSA <-
  function(f,
           L,
           D,
           noise_type="white",
           basis = c("ev", "sin0", "sin1", "t"),
           kind=c("ev", "fa"),
           toeplitz.kind = "no",
           model = NULL,
           freq.range = c(0, 0.5),
           G = 1000,
           level.conf = 0.8,
           two.tailed = FALSE,
           weights = 1,
           composite = FALSE,
           sigma_est_method="mss",
           mode="MLE",
           opt_init_phi=0.3) {
    f <- as.matrix(f)
    if (is.null(model))
    {
      estModel <- list()
      if (D == 1){
        if (noise_type == "white")
          estModel <- est.model.arima(f, freq.range=freq.range, method=sigma_est_method)
        else if (noise_type == "red")
          estModel <- est.model.arima.red(f, freq.range=freq.range, mode=mode, opt_init_phi=opt_init_phi)
      }
      else {
        for (channel in 1:D){
          if (noise_type == "white")
            estModel[[channel]] <- est.model.arima(f[,channel], freq.range=freq.range, method=sigma_est_method)
          else if (noise_type == "red")
            estModel[[channel]] <- est.model.arima(f[,channel], freq.range=freq.range, mode=mode, opt_init_phi=opt_init_phi)
        }
      }
    }
    else
      estModel <- model
    if (basis == "ev") {
      #f.basis <- f
      #if (composite) 
        #f.basis <- f - model$signal
      if (kind == 'fa')
        basis <- basis.ev(f, L, factor.v=T, toeplitz.kind = toeplitz.kind)
      else
        basis <- basis.ev(f, L, factor.v=F, toeplitz.kind = toeplitz.kind)
    }
    else if (basis == "sin0") {
      basis <- basis.sin0(L, freq.range = freq.range)
    }
    else if (basis == "sin1"){
      basis <- basis.sin1(L)
    }
    else {
      if (kind == 'fa')
        basis <- basis.toeplitz.m(estModel, estModel[[1]]$N - L + 1, D, fa=T)
      else
        basis <- basis.toeplitz.m(estModel, L, D, fa=F)
    }
    # basis <- basis.toeplitz(estModel, L*D)
    
    plan <- list(U = basis$U,
                 freq = basis$freq,
                 range = freq.range)
    res <-
      do.ci(
        f,
        plan = plan,
        kind = kind,
        model = estModel,
        level.conf = level.conf,
        L = L,
        G = G,
        D = D,
        two.tailed = two.tailed,
        weights = weights,
        composite = composite
      )
    #,transf = function(x){return(log(x))}, inv.transf = function(x){return(exp(x))})
    res
  }

MonteCarloSSA.single <-
  function(f,
           L,
           basis = c("ev", "sin0", "sin1", "t"),
           model = NULL,
           freq.range = c(0, 0.5),
           G = 1000,
           level.conf = 0.8,
           two.tailed = FALSE) {
    if (is.null(model))
      estModel <- est.model.arima(f)
    else
      estModel <- model
    if (basis == "ev")
      basis <-
        basis.ev(f, L)
    else if (basis == "sin")
      basis <- basis.sin(L)
    else
      basis <- basis.toeplitz(estModel[[1]], L)
    plan <- list(U = basis$U,
                 freq = basis$freq,
                 range = freq.range)
    
    res <-
      do.ci.single(
        f,
        plan = plan,
        model = estModel,
        level.conf = level.conf,
        L = L,
        G = G,
        two.tailed = two.tailed
      )
    res
  }

# There is implemenatation of correction radical/conservative criteria
correction <- function(p.values) {
  corr <- data.frame(alpha=numeric(length(alphas)), fpr=numeric(length(alphas)))
  for (i in alphas_idx) 
  {
    alpha.corr <- alphas[i]
    
    corr$fpr[i] <- sum(p.values < alpha.corr)/M
    corr$alpha[i] <- alpha.corr
  }
  
  alphas.fun <- approxfun(corr$fpr, corr$alpha, rule=2)
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

# TS Generation

# generate sinusoidal signal with specified frequency
signal.one.channel <- function(N, omega, A=1) {
  num <- 1:N
  if (is.null(omega))
    y.clear <- 0*num
  else
    # y.clear <- amplitude(gamma, x) * sin(2*pi*(1:N)/x+ runif(1,0,pi))
    y.clear <- A * cos(2 * pi * num * omega + runif(1,0,pi))
  y.clear
}

# Generate multichanel t.s.
generate <- function(model, signal, D) {
  # res<-replicate(D, one.channel.ts(model, 0))
  res <- list()
  if (D == 1)
    res <- one.channel.ts(model, as.numeric(signal))
  else {
    for (channel in 1:D)
    {
      res[[channel]] <- one.channel.ts(model[[channel]], signal[[channel]])
    }
    res <- matrix(unlist(res), ncol = D, nrow = model[[1]]$N)
  }
  # unlist(res, use.names=FALSE)
  # unlist(res, recursive =F, use.names=FALSE)
  #matrix(unlist(res), ncol = D, nrow = model[[1]]$N)
  # as.vector(unlist(res, recursive =F, use.names=FALSE))
  res
}
