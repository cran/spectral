#' Calculates the envelope of a band limited signal
#'
#' The envelope of an amplitude modulated signal can be calculated by using the
#' Hilbert-transform \eqn{H(y)} of the signal or the analytic signal.
#'
#' An amplitude  modulated function \eqn{y(x) = A(x) * cos(\omega * x)} can be
#' demodulated as follows:
#'
#' \deqn{A(x)^2 = y(x)^2 + H(y(x))^2}
#'
#' If the signal is not band limited, strange things can happen. See the ripple
#' at the edges in the example below. Pay attention, that the envelope is always
#' the real part of the returned value.
#'
#' @param y numeric vector of the signal
#' @return real valued envelope function of the signal
#' @example R/examples/envelopeExample.R
#' @export
envelope<-function(y)
{
  if(is.complex(y))
  {
    y <- base::Re(y)
    warning("Used only real part!")
  }
  return(sqrt(y^2+H(y)^2))
}
