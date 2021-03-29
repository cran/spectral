#' Analytic function
#'
#' In general a causal real valued signal in time has negative
#' frequencies, when a Fourier transform is applied. To overcome this,
#' a complex complement can be calculated to compensate the negative
#' frequency spectrum. The result is called analytic signal or analytic
#' function, which provides a one sided spectrum.
#'
#' An analytic function \eqn{xa} is composed of the real valued signal
#' representation \eqn{y} and its Hilber transform \eqn{H(y)} as the
#' complex complement \deqn{xa(t) = x(t)+i H(x(t))}.
#' In consequence, the analytic function has a one sided spectrum,
#' which is more natural. Calculating the discrete Fourier transform
#' of such a signal will give a complex vector, which is only non zero
#' until the half of the length. Every component higher than the half
#' of the sampling frequency is zero. Still, the analytic signal
#' and its spectrum are a unique representation of the original signal
#' \eqn{x(t)}. The new properties enables us to do certain filtering
#' and calculations more efficient in the spectral space compared to the
#' standard FFT approach. Some examples are:
#'
#' \describe{
#'  \item{Filtering}{because the spectrum is one sided, the user must
#'  only modifiy values in the lower half of the vector. This strongly
#'  reduces mistakes in indexing.
#'  See \code{\link{filter.fft}}}
#'
#'  \item{Envelope functions}{Since the Hilbert transform is a perfect phase shifter
#'  by pi/2, the envelope of a band limited signal can be calculated.
#'  See \code{\link{envelope}}}
#'
#'  \item{Calculations}{Deriving and integrating on band limited discrete data becomes
#'  possible, without taking the symmetry of the discrete Fourier transform into
#'  account. The secound example of the \code{\link{spec.fft}} function calculates
#'  the derivative as well, but plays with a centered spectrum and its corresponding
#'  "true" negative frequencies}
#' }
#' A slightly different approach on the analytic signal can be found in R. Hoffmann
#' "Signalanalyse und -erkennung" (Chap. 6.1.2). Here the signal \eqn{x(t)} is split
#' into the even and odd part. According to Marko (1985) and Fritzsche (1995)
#' this two parts can be composed to the analytic signal, which lead to the
#' definition with the Hilbert transform above.
#'
#' @references
#'  R. Hoffmann, Signalanalyse und -erkennung: eine Einfuehrung fuer
#'  Informationstechniker, Berlin; Heidelberg: Springer, 1998.
#' @references
#'  H. Marko, Systemtheorie: Methoden und Anwendungen fuer ein- und mehrdimensionale
#'  Systeme. 3. Aufl., Berlin: Springer, 1995.
#' @references
#'  G. Fritzsche, Signale und Funktionaltransformationen - Informationselektronik.
#'  Berlin: VEB Verlag Technik, 1985
#'
#'
#' @param x real valued data vector
#' @return Complex valued analytic function
#' @export
analyticFunction <- function(x)
{
  # normalized FFT-spectrum
  X <- fft(x) / length(x)

  # virtual spatial vector is symmetrical to 0.
  # simplifies signum-function, because data
  # need not to be of a whole-numbered length
  xf <- 0:(length(X) - 1)
  xf <- xf - mean(xf)

  X <- X * (1 - sign(xf - 0.5))
  # correct DC value
  X[1] <- 0.5 * X[1]

  return(fft(X, inverse = T))
}

# # The extended approach to calculate the one-sided
# # spectrum of a signal is the the  decomposition in
# # even and odd functions, which is documented in
# # 'Signalanalyse und -Erkennung' by R. Hoffmann (1998)
# #
# # first decompose the signal x in even and odd parts
# xe <- 1 / 2 * (x + x[length(x):1])
# xo <- 1 / 2 * (x - x[length(x):1])
# # calculate the ffts
# XRe <- fft(xe) / length(xe)
# XIo <- fft(xo) / length(xe)
# # Now the spectrum must be mirrored at fc/2,
# # which corresponds to the half of the sampling frequency.
# # Doing so we must invoke kind of a sign-function to change
# # signs. In the end this equals to the Hilbert transform calculated
# # in the spectral domain.
# XRo <- XRe * c(rep(1, length(XRe) / 2), rep(-1, length(XRe) / 2))
# XIe <- XIo * c(rep(1, length(XIo) / 2), rep(-1, length(XIo) / 2))
# X <- XRe + XRo + XIe + XIo
# return(fft(X, inverse = T))
# # Of course the shorter way to an one-sided spectrum is the pure
# # Hilbert transform:
# # x+1i*H(x)
# # which is way faster too.
