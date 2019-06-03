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

w <- 2*pi*20.25
v <- 0.1
windows()

res <- NULL

e <- 4/5*qnorm(1-1e-3) * v * 1/(sqrt(length(x))) # *(1-(w*diff(range(x)))^2/(6*length(x))))


for(i in 1:1e4)
{
  n <- rnorm(length(x),0,v)
  p <- rnorm(1)
  y <- cos(w*x+p) + n
  plot(x,y)
  curve(cos(w*x+p),add=T,n = length(x))

  SOT <- (w)*x

  # tau <- atan2(sum(sin(2 * SOT)), sum(cos(2 * SOT))) / 2

  # 0.5*atan2(sum(sin(2*w*x)),sum(cos(2*w*x)))
  wt <- c(1,min(diff(x))/diff(x))

  ### für gleichmäßig abgetastete Punkte ist das OK
  tau <- (max(w*x) - min(w*x))/2
  tau <- tau - (tau %/% (pi/2) ) * pi/2

  ###
  # tau <- atan2((sin(w*max(x)) - sin(w*min(x)))/n, (cos(w*max(x)) - cos(w*min(x)))/n)

  arg <- SOT - tau


  cs <- cos(arg)
  ss <- sin(arg)

  R <- sum(y * cs)
  I <- sum(y * ss)

  C <- sum(cs ^ 2)
  S <- sum(ss ^ 2)

  A <- sqrt((R/C)^2 + (I/S)^2)

  dA <- abs(1-A)

  res <- rbind(res,c(e=e,dA=dA))
  print(i)

}
res <- as.data.frame(res)
res$e <- e
res$cond <- e*5/4<res$dA

plot(res$dA,ylim=range(res$dA)*2,col=res$cond+1)
abline(h = 5/4*e,lty = 2)

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

