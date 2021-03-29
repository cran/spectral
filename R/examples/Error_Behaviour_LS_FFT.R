### numerischer Fehler -- Test für LS verfahren
#
#
#
#

require(spectral)
dx <- 1e-3
x <- seq(0,1,by = dx)
# x <- sample(x = x,size = 1e2,replace = F)
x <- sort(x)
n <- length(x)

w <- 2*pi*20.25 # frequency
v <- 0.1        # noise standard deviation

res <- NULL

e <- 4/pi*qnorm(1 - 1e-2) * v * sqrt(2/length(x))

set.seed(Sys.time())

for(i in 1:1e5)
{
  n <- rnorm(length(x),0,v)
  p <- rnorm(1) # random phase
  y <- cos(w*x + p) + n


  l <- spec.lomb(x = x, y = y, f = rep(w / (2*pi),2))
  dA <- 1 - l$A[1]

  res <- rbind(res,c(e=e,dA=dA))
  cat(i,"\r")

}
res <- as.data.frame(res)
res$e <- e
res$cond <- e < abs(res$dA)

plot(res$dA,ylim=range(res$dA)*2,col=res$cond+1)
abline(h = e*c(-1,1),lty = 2)

print(sum(res$cond))


#### FFT Spectrum mit Lomb rekonstruieren ####
#
#
#
#
############################################
x0 <- seq(0,1,by=1e-3)

# 10% Lücken setzen
x <- x0[x < 0.1 | x > 0.2]

y <- sin(pi*x)*(cos(2*pi*10*x) + rnorm(length(x),sd = 0.5))

# FFT spectrum generieren
df <- 1/diff(range(x))
fmax <- length(x)/diff(range(x))

f <- seq(0,fmax-df,by = df)
LS <- spec.lomb(x = x, y = y, f=f)

LS$A[LS$p > 0.0000001] <- 0
LS$A[LS$f > fmax] <- 0

plot(LS)
ft <- LS$A*exp(1i*(2*pi*LS$f + LS$phi))

yfs <- fft(ft,inverse = T)

plot(x,y,"l")
points(seq(0,1,length.out = length(yfs)),0.5*Re(yfs))

