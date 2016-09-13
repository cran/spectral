#' Estimate the local frequencies
#'
#' A \code{waterfall}-diagramm displays the local frequency in dependence of
#' or spatial vector. One can then locate an event in time or space.
#'
#' Each frequency is evaluated by calculating the demodulation. This is equivalent
#' to the envelope function of the bandpass filtered signal. The frequency of
#' interest defines automatically the center frequency of the applied bandpass
#' with the bandwidth \eqn{BW}:
#' \deqn{BW = f0 / 4, BW < 4df -> BW = 4df, BW > width * df -> BW = width * df}
#' The minimal frequency is \eqn{df} and \eqn{f0} denotes the center
#' frequency of the bandpass.
#' With increasing frequency the bandwidth becomes wider, which lead to a variable
#' resolution in space and frequency. This is comparable to the wavelet transform,
#' which scales the wavelet according to the frequency.
#' However, the necessary bandwidth is changed by frequency to take the
#' uncertainty principle into account. Slow oscillating events are measured precisely
#' in frequency and fast changing processes can be determined more exact in space.
#' This means for a signal with steady
#' increasing frequency the \code{waterfall} function will produce a diagonally
#' stripe. See the examples below.
#'
#' @param y numeric real valued data vector
#' @param x numeric real valued spatial vector. (time or space)
#' @param nf steepness of the bandpass filter, degree of the polynomial.
#' @param width normalized (to \eqn{df}) maximum width of the bandpass.
#'
#' @return a special \code{fft}-object is returned. It has mode "waterfall" and
#'         \code{x} and \code{fx} present, so it is only plotable.
#' @example R/examples/waterfallExample.r
#' @export
waterfall <-
  function(y = stop("y value is missing"),
           x = NULL,
           nf = 3,
           width = 10)
  {
    if (!is.vector(y))
      stop("y must be a vector")
    if (!is.null(x) & !is.vector(x))
      stop("x must be a vector")
    if (is.null(x))
      x <- seq(0, 1, length.out = length(y))

    bandwidth <- function(f0)
    {
      # calculates the bandwidth 'w' for fast_envelope
      # limits for the return value w are [4*df, width*df]
      w <- f0 / 4
      w <- ifelse(w < 4 * df, 4 * df, ifelse(w > width * df, width * df, w))
      return(w)
    }

    fast_envelope <- function(y, x, fc, BW, nf)
    {
      Y <- BP(Y.f, fc, BW, nf) * fft(y) / n * sY.f
      hk <- base::Mod(fft(Y + 1i * Y, inverse = T) / sqrt(2))
      return(hk)
    }

    n  <- length(y)
    df <- 1 / diff(range(x))  # x length
    fs <- (n - 1) * df  # sampling rate

    # get frequency vector up to nyquist frequency
    f <- seq(0, fs / 2, length.out = as.integer(n / 2))

    Y.f <- seq(0, (n - 1) * df, length.out = n)
    sY.f <- (1 - sign(Y.f - mean(Y.f)))
    sY.f[1] <- 1

    # calculate the waterfall data
    WF <-
      sapply(f, function(f0)
        return(fast_envelope(y, x, f0, bandwidth(f0), nf)))

    result <-
      list(
        x = x,
        fx = f,
        A = WF,
        mode = "waterfall",
        center = F,
        inverse = F
      )
    class(result) <- "fft"

    return(result)
  }
