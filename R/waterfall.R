#' Estimate the local frequencies
#'
#' A \code{waterfall}-diagramm displays the local frequency in dependence of
#' or spatial vector. One can then locate an event in time or space.
#'
#' Each frequency is evaluated by calculating the amplitude demodulation, which
#' is equivalent to the envelope function of the band pass filtered signal.
#' The frequency of interest defines automatically the center frequency \eqn{fc} of the
#' applied band pass with the bandwidth \eqn{BW}:
#' \deqn{BW = fc / width, BW < width -> BW = width, BW > width -> BW = fc / width}
#' The frequency is normalized so the minimal frequency is \eqn{1}.
#' With increasing frequency the bandwidth becomes wider, which lead to a variable
#' resolution in space and frequency. This is comparable to the wavelet
#' (or Gabor) transform,
#' which scales the wavelet (window) according to the frequency.
#' However, the necessary bandwidth is changed by frequency to take the
#' uncertainty principle into account. Slow oscillating events are measured precisely
#' in frequency and fast changing processes can be determined more exact in space.
#' This means for a signal with steady
#' increasing frequency the \code{waterfall} function will produce a diagonally
#' stripe. See the examples below.
#'
#' @section Missing values:
#'
#' Given a regualar grid \eqn{x_i = \delta x \cdot i} there might be missing values
#' marked with \code{NA}, which are treated by the function as 0's.
#' This "zero-padding" leads to a loss of signal energy being
#' roughly proportional to the number of missing values.
#' The correction factor is then \eqn{(1 - Nna/N)} as long as \eqn{Nna / N < 0.2}.
#' As long as the locations of missing values are randomly
#' distributed the implemented procedure workes quite robust. If, in any case,
#' the distribution becomes correlated the proposed correction is faulty and
#' projects the wrong energies.
#'
#' The amplitudes and PSD values are compensated to show up an estimate of the
#' "correct" value. Therefore this method is experimental
#'
#' @param y numeric real valued data vector
#' @param x numeric real valued spatial vector. (time or space)
#' @param nf steepness of the bandpass filter, degree of the polynomial.
#' @param type type of weightening function: "poly", "sinc", "bi-cubic","gauss", can be abbreviated
#' @param width normalized maximum "inverse" width of the bandpass \eqn{bw = fc/width}.
#'
#' @return a special \code{fft}-object is returned. It has mode "waterfall" and
#'         \code{x} and \code{fx} present, so it is only plotable.
#' @example R/examples/waterfallExample.r
#' @import pbapply
#' @export
waterfall <-
  function(y = stop("y value is missing"),
           x = NULL,
           nf = 3,
           type = "b",
           width = 7)
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
      # limits for the return value w are [width, f0/width]
      w <- f0 / width
      # w <- ifelse(w < 4, 4, ifelse(w > width, width, w))
      w <- ifelse(w < width, width, w)
      return(w)
    }

    fast_envelope <- function(Y, fc, BW, nf, Y.f, type)
    {
      Y_bp <- BP(Y.f, fc, BW, nf, type) * Y
      hk <- base::Mod(fft(Y_bp + 1i * Y_bp, inverse = TRUE) / sqrt(2) )
      return(hk)
    }

    N  <- length(y)
    Nna <- sum(is.na(y))


    df <- 1  / diff(range(x))  # x length
    fs <- (N - 1) * df  # sampling rate

    sdy <- sd(y,na.rm=T)

    y[is.na(y)] <- 0 # zero padding

    # get frequency vector up to normalized nyquist frequency
    f <- seq(0, (N - 1) / 2, length.out = as.integer(N / 2))

    Y.f <- seq(0, (N - 1), length.out = N)
    sY.f <- (1 - sign(Y.f - mean(Y.f)))
    sY.f[1] <- 1

    # # generate progressbar
    # pb <- progress_bar$new(
    #   format = "[:bar] :current / :total | :eta"
    #   ,total = length(f)
    #   ,force = TRUE
    #
    # )

    # making up progressbar
    pboptions(style = 3, char = "=")
    # calculate the waterfall data
    WF <-
      pbsapply(f, function(f0,Y,nf,Y.f,type)
      {
        # pb$tick()
        return(fast_envelope(Y, f0, bandwidth(f0), nf, Y.f,type))
      }
        ,Y = fft(y) / length(y) * sY.f
        ,nf = nf
        ,Y.f = Y.f
        ,type = type
      )

    WF <- WF / (1 - Nna / N) # compensate missing values

    PSD <- N / (N - 1) / ( 2 * (sdy)^2 ) * ( WF )^2


    result <-
      list(
        x = x,
        fx = f * df,
        A = WF,
        PSD = PSD,
        mode = "waterfall",
        center = FALSE,
        inverse = FALSE
      )
    class(result) <- "fft"

    return(result)
  }
