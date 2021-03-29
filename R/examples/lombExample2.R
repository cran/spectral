## under sampling with jitter ##
require(spectral)


f0 <- 20

x <- sin(2*pi*1.1*seq(0,1,by = 1/20)[-1])
# x <- runif(length(x),0,max(x))
# x <- x + runif(length(x),sd = min(diff(x))/4)
x <- sort(x)
x <- x - min(x); x <- x/max(x)


fy <- function(x) cos(2*pi*f0*x)
y <- fy(x) + rnorm(length(x),sd = 0.05)



f <- seq(0,500,by = 1)

lmb <- spec.lomb(x = x,y = y, f = f,mode = "generalized")
ft <- spec.fft(x = x, y = y)

Aft <- lmb$A * (cos(2*pi*lmb$f + lmb$phi) + 1i * sin(2*pi*lmb$f + lmb$phi))
Aft[lmb$p > 0.005] <- 0

yft <- fft(Aft,inverse = T)

par(mfrow = c(2,1))
plot(x,y,"b")
curve(fy,add = T, col = "darkgrey",n = 1000)
lines(seq(0,max(x),length.out= length(yft)),base::Re(yft),lty=3,lwd=2)


plot(ft,type = "p",xlim = 50*c(-1,1),ylim = c(0,1))
lines(lmb$f,lmb$A/2)
abline(v = f0, lty = 3)
abline(h = 0.5)
legend("topright",c("FT","LOMB"),pch=c(1,NA),lty = c(NA,1))
