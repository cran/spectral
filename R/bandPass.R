#' Simple bandpass function
#'
#' This function represents a simple weighting function for spectral filtering.
#'
#' The band pass is represented troughout a polynom in the form
#' \deqn{w = 1 - a * (f-fc)^n} with the degree \eqn{n}. The parameter \eqn{fc}
#' controlls the center frequency and \eqn{a} scales the required band width \code{BW}.
#' outside the band width the result is forced to zero.
#'
#' @param f vector of frequencies
#' @param fc center frequency
#' @param BW bandwidth, with \code{w[f < (fc - BW) | f > (fc+BW)] == 0}
#' @param n degree of the polynom, \code{n} can be real, e.g. \code{n = 2.5}
#' @return This function returns a weight vector, which is to apply to the frequency
#' vector \code{f} in a top level function
BP <- function(f, fc, BW, n = 3)
{
  # Find roots
  a <- 1 / abs(BW ^ n)
  w <- 1 - abs(a * abs(f - fc) ^ n)
  w[w < 0] <- 0
  return(w)
}
