library(hwwntest)
library(wavethresh)
library(sarima)
#library(normwhn.test)

test.ar1.wavelet <- function(ts0){
  ar.par <- arima(
    ts0,
    order = c(1, 0, 0),
    include.mean = FALSE,
    method = "CSS-ML"
  ) #does not estimate parameters for nonstationary t.s. with trend
  #phi.est <- coef(ar.par) #test is conservative for estimated phi
  phi.est <- phi #test is exact for true phi
  Sigma <- toeplitz(phi.est^(0:(length(ts0)-1)))
  s.svd <- eigen(Sigma)
  Sigma.root.inv <- s.svd$vectors %*% diag(sqrt(1/s.svd$values)) %*% t(s.svd$vectors)
  ts.pre <- t(t(ts0) %*% Sigma.root.inv) 
  genwwn.test(ts.pre, filter.number = 10)$p.value #this works
  #bartlettB.test(ts.pre)$p.value #this works
  
  #normwhn.test(ts.pre)$p.value
  #x.acf <- autocorrelations(ts.pre)
  #whiteNoiseTest(x.acf, h0 = "iid", x = x, method = "LjungBox")$test[3]
}

phi <- 0.7
N <- 128
M <- 1000
for (A in seq(0,2,0.5)){
  p <- numeric(0)
  set.seed(5)
  for(i in 1:M){
    ts.ar <- arima.sim(n = N, model = list(ar=c(phi)))
    #ts0 <- ts.ar + A*(1:N) 
    ts0 <- ts.ar + A*sin(2*pi*(1:N) * 0.075)
    p <- c(p, test.ar1.wavelet(ts0))
  }
  print(sum(p < 0.5) / M)
  
  #print(c(A,mean(p)), quote = FALSE)
}


phi <- 0.7
N <- 128
M <- 1000
#A <- 1
#model <- list(varphi = phi, delta = 1, N = N)

p.value.noise <- c()
set.seed(5)
for (i in 1:M) {
  #ts.ar <- one.channel.ts(model, 0)
  ts.ar <- arima.sim(n = N, model = list(ar = c(phi)))
  #ts0 <- ts.ar + A*(1:N)
  ts0 <- ts.ar
  p.value.noise <- c(p.value.noise, test.ar1.wavelet(ts0))
}

#print(sum(p.value.noise < 0.2) / M)

p.value.signal <- c()
set.seed(5)
for (i in 1:M) {
  #ts.ar <- one.channel.ts(model, 0)
  ts.ar <- arima.sim(n = N, model = list(ar = c(phi)))
  #ts0 <- ts.ar + A*(1:N)
  ts0 <- ts.ar + A * cos(2 * pi * (1:N) * omega)
  p.value.signal <- c(p.value.signal, test.ar1.wavelet(ts0))
}

#print(sum(p.value.signal < 0.2) / M)

alphas <- 0:1000 / 1000
alphaI.wavelet <- sapply(alphas, function(a) sum(p.value.noise < a) / M)
beta.wavelet <- sapply(alphas, function(a) sum(p.value.signal < a) / M)

pdf("alphaI_wavelet.pdf", width = 15, height = 8)
plot(c(0,1),c(0,1), type="l", col = "blue", lty = 2, xlab = 'significance level', ylab = 'type I error', cex.lab=1.8, cex.axis=1.8, cex.sub=1.8)
lines(alphas, alphaI.wavelet, col = "orange", lwd = 2) 
lines(alphas, alphaI[[1]], col = "purple", lwd = 2)
legend(x="bottomright", c("MC-SSA (L=120)", "Whitening + genwwn.test"), col = c("purple", "orange"), lty = 1, lwd = c(2, 2), cex = 2)
dev.off()

pdf("roc_wavelet_omega01.pdf", width = 15, height = 8)
plot(c(0,1),c(0,1), type="l", col = "blue", lty = 2, xlab = 'type I error', ylab = 'power', cex.lab=1.8, cex.axis=1.8, cex.sub=1.8)
lines(alphaI.wavelet, beta.wavelet, lwd = 2, col = "orange")
lines(alphaI[[1]][-1], beta[[1]][-1], col = "purple", lwd = 2)
legend(x="bottomright", c("MC-SSA (L=120)", "Whitening + genwwn.test"), col = c("purple", "orange"), lty = 1, lwd = c(2, 2), cex = 2)
dev.off()
