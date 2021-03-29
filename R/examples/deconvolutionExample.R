### Deconvolution ###
#
#
# we define a test function with gaps and noise
# we show:
# - the aliased Fourier spectrum and for comparison Lomb Spectrum
# - the corrected spectrum
#

## definition of the sampling series
x <- seq(0,pi/2,by = 5e-3)
n <- rnorm(length(x),sd = 0.1)

## definition of the test function
## with 2 frequencies
yf <- function(x)
{
  cos(2*pi*5.123*x) +
  cos(2*pi*17*x)
}

y <- yf(x)
y <- y - mean(y)

## define strongly correlated gaps
i <- NULL

i <- c(i,which(sin(2*pi / 0.3 * x) - 0.5 > 0))
i <- c(i,which(sin(2*pi / 0.04 * x + 1.123) - 0.5 > 0))
i <- sort(unique(i))


xs <- x
ys <- yf(xs) + n # add some noise
ys[i] <- NA

## for comparison we calculate a Lomb-Spectrum
LT <- spec.lomb(x = xs,y = ys
                ,f = seq(0,250,by = 0.02)
                ,mode = "generalized"
)

WS <- deconvolve(x = xs, y = ys,SNR.enable = 1,SNR.level = 1)
FT <- spec.fft(x = xs, y = ys,center = FALSE)
FTS <- spec.fft(x = xs, y = is.na(ys),center = FALSE)
LTS <- spec.lomb(x = xs, y = is.na(ys),f = seq(0,250,by = 0.02))

### results ###
#
# - signal spectrum (solid) dominant peaks at around f0 = {5, 17}
# - (minor) alias peaks (grey line, FFT dots) at f0 +/- fs
# - sampling spectrum (dashed) with fs = {3.3, 25} (dominant modes)
# - deconvolved spectrum (solid black) rejects the aliases and sampling
#
#

### time series

par(mfrow = c(1,1),mar = c(4,4,3,0.3))
curve(yf,0,max(x), col = "grey",n = 1000
      ,xlim = c(0,max(x)),ylim = c(-2,3)
      ,xlab = "Time", ylab = "y(t)"
      ,main = "Fragmented Time Series"
      )
points(xs,ys)
points(xs[is.na(ys)],yf(xs[is.na(ys)]),pch = 16,cex = 0.5)

legend("topright",c("y(t)","y(tn) + n(tn)","NA's")
       ,lty = c(1,NA,NA)
       ,lwd = c(1,NA,NA)
       ,pch = c(NA,1,16)
       ,col = c("darkgrey","black","black")
       ,bg = "white"
       ,cex = 0.8
)

## plot spectra
par(mfrow = c(1,1),mar = c(4,4,3,0.3))
with(FT,plot(fx,PSD,type="p",log = "x"
     # ,col="grey"
     ,xlim = c(1,100),ylim = c(1e-2,0.75)
     ,xlab = "f", ylab = "PSD"
     ,pch = 1
     ,lwd = 1
     ,main = "Spectra"
))
with(LT,lines(f,PSD,col = "grey",lwd = 4))
with(WS,lines(f,S, lwd = 2, col = "black"))
with(LTS,lines(f,PSD,lty = 2))
abline(h = c(1,0.5),lty = 3)
legend("topright",c("Fourier","Lomb","Decon.","Sampling")
       ,lty = c(NA,1,1,2)
       ,lwd = c(2,2,2,2)
       ,pch = c(1,NA,NA,NA)
       ,col = c("black","grey","black","black")
       ,bg = "white"
       ,cex = 0.8
       ,ncol = 2
)
