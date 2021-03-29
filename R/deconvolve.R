#' Deconvolve Sampling Spectrum for Equidistant Sampling
#'
#' The function removes the probable alias
#' peaks in the power spectral density. These projections originate
#' from correlated gaps, missing values and interactions with noise.
#' The function should be considered as *experimental* but with didactic background.
#'
#' In the special case of a non complete equidistant grid
#' containing the data and missing values (\code{NA}), this function
#' performs the deconvolution of \code{Y = fft(y)} from the
#' sampling spectrum of the aquisition series \code{x}. The data
#' is assumed to exist on a equidistant grid with missing values and
#' gaps.
#'
#' Given a one dimensional vector \code{y} of data this function reverses
#' the spectral convolution of \eqn{Y = S * X + N}, if * describes the convolution
#' operation and \eqn{Y = F(y)} denotes the discrete Fourier transform
#' via the operator \eqn{F(.)}. If, the sampling series \code{x} is considered to be purely deterministic,
#' which should be the case for captured data, or the distortions
#' (missing values, gaps) are *correlated* (see example), then there exists an analytic inversion of
#' the convolution. Given the general definition of power spectral density
#' \eqn{|Y|^2 = |S * X + N|^2} the challenge is to prove \eqn{|S * X + N|^2 ~ |S|^2 * |X|^2 + |N|^2}.
#' Here \eqn{N} describes a stochastic term of gaussian noise. This issue is
#' solved in correlation space, where convolution becomes a multiplitation. The
#' auto correlation function (acf) of \eqn{y} is given by \eqn{Ry = F(|Y|^2)}.
#' As a remark, IF we consider the special case of
#' equispaced sampling, modeled by the Diraq distribution \eqn{\delta}(x),
#' it is easy to show that the correlation function of a product
#' is the product of individual correaltaion functions, \eqn{F(|S*X|^2) = F(|S|^2) . F(|X|^2)}.
#'
#' The aim is, to approximate \eqn{S} as the "true" spectrum. To the cost
#' of the phase information, the result is the standardized power spectral
#' density. The spectral noise term \eqn{F(N)} is approximated by a theshold in
#' Fourier space. Here \code{SNR.level} sets the factor of \code{mean(fft(y))} below
#' which noise level is assumed. Above this value, the signal should be present.
#' As a parameter to play with, \code{SNR.enable} enables or disables the noise term.
#' This parameter was introduced to be consistent with present approaches,
#' not considering the presence of noise.
#'
#'
#' @param x sampling instances
#' @param y values
#' @param SNR.enable binary value, include or exclude the noise
#' @param SNR.level theshold in the sense of a multiple of mean()
#'                  noise level
#'
#' @return list of frequency \code{f} and spectral density function \code{S}
#' @concept deconvolution
#' @concept convolution
#' @concept power spectral density
#' @example R/examples/deconvolutionExample.R
#' @export
deconvolve <- function(x,y, SNR.enable = T, SNR.level = 1)
{
  N <- length(y)
  M <- sum(is.na(y))

  f <- seq(0,1/min(diff(x)),length.out = N)

  # Zero padding of missing values
  tmp <- y
  tmp[is.na(y)] <- 0


  # Amplitude Spectrum
  Y <- fft(tmp) / N

  Sy <- Re(Y * Conj(Y))
  Ry <- Re(fft(Sy,inverse = T))

  ## estimate noise spectrum
  Sn <- Sy

  # find all signal amplitudes to cut them
  cond <- Sn > mean(Sy)*SNR.level

  # mast the probable signal amplitudes
  Sn[cond] <- abs(rnorm(sum(cond),sd = sqrt(2)*mean(Sy)))
  Rn <- Re(fft(Sn,inverse = T))

  # define sampling function as series of binaries
  w <- rep(1,N)
  w[is.na(y)] <- 0

  # amplitude spectrum of sampling function
  W <- fft(w) / N

  Sw <- Re(W * Conj(W))
  Rw <- Re(fft(Sw,inverse = T))

  # first bin Rw[1] = 1
  # because it includes energy content of all possible sampling positions,
  # it must be 1 in the normalized case, other whise N^2.
  # the later forces Y = fft(y) and Sy <- Re(Y * Conj(Y)) / N
  Rw[1] <- 1
  # estimation of SNR^-1
  iSNR <- Rn / Ry
  iSNR[1] <- 0 # inverses SNR!!!

  iSNR[is.na(iSNR)] <- 1

  # estimation of original spectrum
  # normalization to N only makes sense if it is normalized too, W <- fft(w) / N
  # remind to check/change the rest of the definition for the difference between
  # normalized (Fourier Series) and Fourier Transform
  Sp <- 2*abs( fft( ( 1 - SNR.enable*iSNR ) / Rw * Ry )) / N

  return(list(f = f, S = Sp) )
}
