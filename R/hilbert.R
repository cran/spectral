#' The Hilbert transformation
#'
#' The Hilbert transform is a phase shifter, which represents the complex complement
#' to a real vauled signal. It is calculated in the complex frequency space of the
#' signal by using the Fourier transform. Finally, calculating \eqn{f = y + i*H(y)}
#' gives the analytic signal, with a one sided spectrum. (See \code{\link{analyticFunction}})
#'
#' @usage H(x)
#'
#' @param x real valued time series
#' @return A numeric real valued vector is returned
#' @export
H <- function(x)
{
  x <- as.numeric(x)
  # first calculate the normalized FFT
  X <- fft(x) / length(x)

  # then we need a virtual spatial vector which is symmetric with respect to
  # f = 0. The signum function will do that. The advantage is, that we need not
  # take care of the odd-/evenness of the length of our dataset
  xf <- 0:(length(X) - 1)
  xf <- xf - mean(xf)

  # because the negative Frequencies are located in the upper half of the
  # FFT-data vector it is nesccesary to use "-sign". This will mirror the relation
  # The "-0.5" effect is that the Nyquist frequency, in case of odd data set lenghts,
  # is not rejected.
  Xh <- -1i * X * -sign(xf - 0.5)
  return(fft(Xh, inverse = T))
}
