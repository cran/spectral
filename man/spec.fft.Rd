% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/specFFT.R
\name{spec.fft}
\alias{spec.fft}
\title{1D/2D/nD (multivariate) spectrum of the Fourier transform}
\usage{
spec.fft(y = NULL, x = NULL, z = NULL, center = T)
}
\arguments{
\item{y}{1D data vector, y coordinate of a 2D matrix, nD (even 2D) array
or object of class \code{fft}}

\item{x}{x-coordinate of the data in \code{y} or \code{z}. If \code{y} is an array, \code{x} must be a named list \code{x = list(x = ..., y = ...)}.}

\item{z}{optional 2D matrix}

\item{center}{logical vector, indicating which axis to center in frequency space}
}
\value{
An object of the type \code{fft} is returned. It contains the
spectrum \code{A}, with "reasonable" frequency vectors along each ordinate. \code{psd} represents
the standardized power spectral density, [0,1]. The false alarm probability (FAP)
\code{p} is given similar to the Lomb-Scargle method, see \link{spec.lomb}.
}
\description{
This function calculates the Fourier spectrum and power spectral density
of a given data object. The dimension of the array can be of arbitary size
e. g. 3D or 4D.
}
\details{
The function returns an user friendly object, which contains as much frequency
vectors as ordinates of the array. \code{spec.fft} provides the
ability to center the spectrum along multiple axis. The amplitude output is already
normalized to the sample count and the frequencies are given in terms of
\eqn{1/\Delta x}-units.
}
\section{Missing Values}{


Given a regualar grid \eqn{x_i = \delta x \cdot i} there might be missing values
marked with \code{NA}, which are treated by the function as 0's.
This "zero-padding" leads to a loss of signal energy being
roughly proportional to the number of missing values.
The correction factor is then \eqn{(1 - Nna/N)} as long as \eqn{Nna / N < 0.2}.
If the locations of missing values are randomly
distributed the implemented procedure workes quite robust. If correalted
gaps are present, the proposed correction is faulty and
scales wrong. This is because a convolution of the incomplete
sampling vector with the the signal takes place. An aliasing effect
takes place distorting the spectral content.

To be compatible with the underlying Fourier transform, the amplitudes
are not affected by this rescaling.
Only the power spectral density (PSD) is corrected in terms of the energy
content, which is experimental for the moment.
}

\examples{
# 1D Example with two frequencies
#################################

x <- seq(0, 1, length.out = 1e3)
y <- sin(4 * 2 * pi * x) + 0.5 * sin(20 * 2 * pi * x)
FT <- spec.fft(y, x)
par(mfrow = c(2, 1))
plot(x, y, type = "l", main = "Signal")
plot(
  FT,
  ylab = "Amplitude",
  xlab = "Frequency",
  type = "l",
  xlim = c(-30, 30),
  main = "Spectrum"
)
summary(FT)

# 2D example with a propagating wave
####################################

x <- seq(0, 1, length.out = 50)
y <- seq(0, 1, length.out = 50)

# calculate the data
m <- matrix(0, length(x), length(y))
for (i in 1:length(x))
  for (j in 1:length(y))
    m[i, j] <- sin(4 * 2 * pi * x[i] + 10 * 2 * pi * y[j])

# calculate the spectrum
FT <- spec.fft(x = x, y = y, z = m)

# plot
par(mfrow = c(2, 1))
rasterImage2(x = x,
           y = y,
           z = m,
           main = "Propagating Wave")
plot(
  FT,
  main = "2D Spectrum",
  palette = "wb"
  ,
  xlim = c(-20, 20),
  ylim = c(-20, 20),
  zlim = c(0, 0.51)
  ,
  xlab = "fx",
  ylab = "fy",
  zlab = "A",
  ndz = 3,
  z.adj = c(0, 0.5)
  ,
  z.cex = 1
)
summary(FT)

# 3D example with a propagating wave
####################################

# sampling vector
x <- list(x = seq(0,2,by = 0.1)[-1]
          ,y = seq(0,1, by = 0.1)[-1]
          ,z = seq(0,1, by = 0.1)[-1]
)

# initializing array
m <- array(data = 0,dim = sapply(x, length))

for(i in 1:length(x$x))
  for(j in 1:length(x$y))
    for(k in 1:length(x$z))
      m[i,j,k] <- cos(2*pi*(1*x$x[i] + 2*x$y[j] + 2*x$z[k])) + sin(2*pi*(1.5*x$x[i]))^2

FT <- spec.fft(x = x, y = m, center = c(TRUE,TRUE,FALSE))

par(mfrow = c(2,2))
# plotting m = 0
rasterImage2( x = FT$fx
              ,y = FT$fy
              ,z = abs(FT$A[,,1])
              ,zlim = c(0,0.5)
              ,main="m = 0"
              )

# plotting m = 1
rasterImage2( x = FT$fx
              ,y = FT$fy
              ,z = abs(FT$A[,,2])
              ,zlim = c(0,0.5)
              ,main="m = 1"
)

# plotting m = 2
rasterImage2( x = FT$fx
              ,y = FT$fy
              ,z = abs(FT$A[,,3])
              ,zlim = c(0,0.5)
              ,main="m = 2"
)
rasterImage2( x = FT$fx
              ,y = FT$fy
              ,z = abs(FT$A[,,4])
              ,zlim = c(0,0.5)
              ,main="m = 3"
)

summary(FT)


# calculating the derivative with the help of FFT
################################################
#
# Remember, a signal has to be band limited.
# !!! You must use a window function !!!
#

# preparing the data
x <- seq(-2, 2, length.out = 1e4)
dx <- mean(diff(x))
y <- win.tukey(x) * (-x ^ 3 + 3 * x)

# calcualting spectrum
FT <- spec.fft(y = y, center = TRUE)
# calculating the first derivative
FT$A <- FT$A * 2 * pi * 1i * FT$fx
# back transform
dm <- spec.fft(FT)

# plot
par(mfrow=c(1,1))
plot(
  x,
  c(0, diff(y) / dx),
  type = "l",
  col = "grey",
  lty = 2,
  ylim = c(-4, 3)
)
# add some points to the line for the numerical result
points(approx(x, Re(dm$y) / dx, n = 100))
# analytical result
curve(-3 * x ^ 2 + 3,
      add = TRUE,
      lty = 3,
      n = length(x))

legend(
  "topright",
  c("analytic", "numeric", "spectral"),
  title = "diff",
  lty = c(3, 2, NA),
  pch = c(NA, NA, 1),
  col=c("black","grey","black")
)
title(expression(d / dx ~ (-x ^ 3 + 3 * x)))
}
\seealso{
\link{plot.fft}
}
\concept{Fourier Transform}
