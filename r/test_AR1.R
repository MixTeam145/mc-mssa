library(hwwntest)
library(wavethresh)
library(sarima)
#library(normwhn.test)

test.ar1.wavelet <- function(ts0, phi) {
  # ar.par <- arima(
  #   ts0,
  #   order = c(1, 0, 0),
  #   include.mean = FALSE,
  #   method = "CSS-ML"
  # ) # does not estimate parameters for nonstationary t.s. with trend
  #phi.est <- coef(ar.par) #test is conservative for estimated phi
  phi.est <- phi #test is exact for true phi
  Sigma <- toeplitz(phi.est ^ (0:(length(ts0) - 1)))
  s.svd <- eigen(Sigma)
  Sigma.root.inv <- s.svd$vectors %*% diag(sqrt(1 / s.svd$values)) %*% t(s.svd$vectors)
  ts.pre <- t(t(ts0) %*% Sigma.root.inv)
  genwwn.test(ts.pre, filter.number = 10)$p.value #this works
  #bartlettB.test(ts.pre)$p.value #this works
  
  #normwhn.test(ts.pre)$p.value
  #x.acf <- autocorrelations(ts.pre)
  #whiteNoiseTest(x.acf, h0 = "iid", x = x, method = "LjungBox")$test[3]
}

test.ar1 <- function(ts0) {
  # ar.par <- arima(
  #   ts0,
  #   order = c(1, 0, 0),
  #   include.mean = FALSE,
  #   method = "CSS-ML"
  # ) #does not estimate parameters for nonstationary t.s. with trend
  # phi.est <- coef(ar.par) #test is conservative for estimated phi
  phi.est <- phi #test is exact for true phi
  Sigma <- toeplitz(phi.est ^ (0:(length(ts0) - 1)))
  s.svd <- eigen(Sigma)
  Sigma.root.inv <- s.svd$vectors %*% diag(sqrt(1 / s.svd$values)) %*% t(s.svd$vectors)
  ts.pre <- t(t(ts0) %*% Sigma.root.inv)
  Box.test(ts.pre)$p.value
}
