# Multiple Monte Carlo SSA
# [Golyandina N. Detection of signals by Monte Carlo singular spectrum analysis:
# multiple testing // Statistics and Its Interface. — 2023. — Vol. 16, no. 1. — P. 147–157.]

# Contribution of Egor Poteshkin:
# correction of liberal criteria is implemented
# projection vectors with parameter basis="ev" and D=1 are from Toeplitz SSA of original time series
# Toeplitz MC-MSSA (D > 1) is draft

library("Rssa")
library("pracma")
library("ggplot2")
library("matrixStats")
library("magic")

source("toeplitz_mssa.R")

type = 8

# Estimation of AR(1) parameters, without signal extraction
est.model.arima <-  function(f) {
  param <- list()
  ar.par <- arima(
    f,
    order = c(1, 0, 0),
    include.mean = FALSE,
    method = "CSS-ML"
  )
  param$varphi <- coef(ar.par)
  if (sqrt(ar.par$var.coef) > abs(ar.par$coef))
    param$varphi <- 0
  param$delta <- sqrt(ar.par$sigma2)
  estModel <-
    list(varphi = param$varphi,
         delta = param$delta,
         N = length(f))
  estModel
}
###end

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
  N <- length(f[,1]) # assert equal length in each channel
  K <- N - L + 1
  X_res <- matrix(0, nrow = L, ncol = K*D)
  for (channel in 1:D) {
    tX <- sapply(1:L, function(i) f[i:(i + K - 1), channel])
    X_res[, (1 + (channel - 1) * K):(channel * K)] <- t(tX)
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
basis.ev <- function(ts, L, factor.v = T, toeplitz.kind) {
  D <- dim(ts)[2]
  neig <- min(L, D * (length(ts[,1]) - L + 1))
  if (D == 1)
    s <- ssa(ts, L = L, neig = neig, kind = "toeplitz-ssa")
  else
    if (identical(toeplitz.kind, "sum") || identical(toeplitz.kind, "block"))
      s <- toeplitz.mssa(ts, L = L, D = D, method = toeplitz.kind, neig =  neig)
  else
    s <- ssa(ts, L = L, neig = neig, kind = "mssa")
  freq <- numeric(0)
  for (i in 1:nu(s)){
    #ss <- ssa(s$U[,i], kind = "toeplitz-ssa")
    ss <- ssa(s$U[,i], kind = "1d-ssa")
    #estimation of the main frequency by ESPRIT
    p <- parestimate(ss, groups = list(1:2))
    freq[i] <- p$frequencies[[1]]
  }
  if (factor.v) {
    return(list(U = s$V, freq = freq))
  }
  else {
    return(list(U = s$U, freq = freq))
  }
}

# Generate vectors for projections corresponding to eigenvectors teoretical matrix 
matrix.toeplitz <- function(phi, L) {
  toeplitz(phi ^ (0:(L - 1)))
}

basis.toeplitz <- function(model, L, D, fa = F) {
  if (fa) {
    # here we assume that L param represents K = N - L + 1
    toepl.array <- list()
    for (channel in 1:D) {
      toepl.array[[channel]] <- model[[channel]]$delta * matrix.toeplitz(model[[channel]]$varphi, L)
    }
    toepl <- do.call("adiag", toepl.array)
    s <- svd(toepl, nv = L)
    U <- s$v
  }
  else {
    toepl <- matrix(data = 0,
                    nrow = L,
                    ncol = L)
    for (channel in 1:D) {
      toepl <- toepl + matrix.toeplitz(model[[channel]]$varphi, L)
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
    
    if (weights[1] == "equal")
      weights <- ci$sds
    ci$sds[weights == 0] <- 1000 # something large
    ci$sds[weights != 0] <- ci$sds[weights != 0] / weights[weights != 0]
    
    stats.max <-
      apply(X, 2, function(vv)
        max((vv - ci$means) / ci$sds))
    stats.max.abs <-
      apply(X, 2, function(vv)
        max(abs(vv - ci$means) / ci$sds))
    
    if (!is.null(level.conf)) {
      if (two.tailed == FALSE) {
        ci$q.upper <- quantile(stats.max, probs = level.conf, type = type)
        ci$q.lower <- 0
      } else {
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
      if (res$reject == TRUE)
        res$freq.max <-
        plan$freq[idx][which.max((x - ci$means) / ci$sds)]
      
      res$upper <- inv.transf(ci$means + ci$q.upper * ci$sds)
      res$lower <- 0
      if (two.tailed == TRUE)
        res$lower <- inv.transf(ci$means + ci$q.lower * ci$sds)
      
      res$plan <- plan
      res$v <- x
      res$f <- f
      res$idx <- idx
    }
    if (two.tailed == FALSE) {
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

#plot by dominant frequency
plot.ci <- function(res, log_ = FALSE) {
  df <-
    data.frame(
      frequency = res$freq,
      contribution = res$v,
      ci.lower = res$lower,
      ci.upper = res$upper
    )
  df_reject <- df |> filter(contribution < ci.lower |
                              contribution > ci.upper)
  p <-
    ggplot(df, aes(frequency, contribution)) +
    geom_point() +
    geom_point(data = df_reject,
               aes(x = frequency, y = contribution),
               color = "red") +
    geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), color = "blue") +
    theme_bw()
  if (log_)
    p <- p + scale_y_log10()
  p
}

#plot by dominant contribution (the projection norm)
plot.ci.by.order <- function(res, log_ = FALSE) {
  df <-
    data.frame(
      contribution = res$v,
      ci.lower = res$lower,
      ci.upper = res$upper
    ) |> arrange(desc(contribution))
  
  df$num <- 1:length(res$freq)
  
  df_reject <- df |> filter(contribution < ci.lower |
                              contribution > ci.upper)
  
  p <-
    ggplot(df, aes(num, contribution)) +
    geom_point() +
    geom_point(data = df_reject,
               aes(x = num, y = contribution),
               color = "red") +
    geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), color = "blue") +
    theme_bw()
  if (log_)
    p <- p + scale_y_log10()
  p
}

# The wrapped function for Multiple Monte Carlo SSA
MonteCarloSSA <-
  function(f, # time series
           L, # window length
           D = 1, # number of channels
           basis = c("ev", "t"), # vectors for projection
           kind = c("ev", "fa"), # left or right vectors (if D=1 "ev" and "fa" are equal up to replacement L->N-L+1)
           toeplitz.kind = c("no", "sum", "block"), # for toeplitz mc-mssa
           model = NULL, # AR(1) model
           freq.range = c(0, 0.5), #
           G = 1000, # number of surrogates
           level.conf = 0.8,
           two.tailed = FALSE,
           weights = 1,
           composite = FALSE # mc-ssa with nuisance signal
  ) {
    
    if (D == 1) {
      f <- as.matrix(f)
      if (!is.null(model))
        model <- list(model)
    }
    
    kind <- match.arg(kind)
    
    if (composite & D != 1)
      stop("mc-ssa with nuisance signal for multivariate ts is not implemented")
    
    if (is.null(model)) {
      estModel <- list()
      for (channel in 1:D)
        estModel[[channel]] <- est.model.arima(f[,channel])
    }
    else
      estModel <- model
    if (basis == "ev") {
      f.basis <- f
      # comment next 2 lines to project vectors of original series (another version of the nuisance algorithm)
      if (composite) 
        f.basis <- f - estModel[[1]]$signal
      if (kind == 'fa')
        basis <- basis.ev(f.basis, L, factor.v=T, toeplitz.kind = toeplitz.kind)
      else
        basis <- basis.ev(f.basis, L, factor.v=F, toeplitz.kind = toeplitz.kind)
    }
    else {
      if (kind == 'fa')
        basis <- basis.toeplitz(estModel, estModel[[1]]$N - L + 1, D, fa=T)
      else
        basis <- basis.toeplitz(estModel, L, D, fa=F)
    }
    
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
    res
  }

# There is implemenatation of correction liberal/conservative criteria
correction <- function(p.values, alphas = 0:1000 / 1000) {
  alphaI <- sapply(alphas, function(a) mean(p.values < a))
  alphas.fun <- approxfun(alphaI, alphas, rule=2)
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

# Generates sinusoidal signal with specified frequency
signal.one.channel <- function(N, omega, A = 1) {
  num <- 1:N
  if (is.null(omega))
    signal <- 0 * num
  else
    signal <- A * cos(2 * pi * num * omega + runif(1, 0, pi))
  signal
}

# Generates a series according to the model signal + AR(1)
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
  if (!is.null(model$signal)) # composite null hypothesis
    xi <- xi + model$signal
  f <- xi + signal
  as.vector(f)
}

# Generates a multichanel ts
generate <- function(model, signal, D) {
  res <- list()
  for (channel in 1:D)
    res[[channel]] <- one.channel.ts(model[[channel]], signal[[channel]])
  
  res <- matrix(unlist(res), ncol = D, nrow = model[[1]]$N)
  res
}
