#' Lomb-Scargle Periodigram
#'
#' The Lomb-Scargle periodigram represents an
#' statistical estimator for the amplitude and
#' phase for a given frequency.
#'
#' A given time series does not need to be evenly sampled. This means a time series mainly
#' consists of data pairs \code{x} and \code{y}, which store the data and the sampling
#' position (e.g. in time). Additionally, this method enables the user to analyse
#' the data with respect to a given frequency vector, which can be artifical dense.
#' \describe{
#'  \item{\code{ofac}}{If the user does not provide a corresponding frequency vector, the
#'  \code{ofac} parameter causes the function to estimate \deqn{nf = ofac*length(x)/2}
#'  equidistant frequencies.}
#'  \item{\code{p}-value}{The \code{p}-value gives the probability, wheater the
#'  estimated amplitude is NOT significant. However, if \code{p} tends to zero the
#'  amplidutde is significant. The user must decide which maximum value is acceptable,
#'  until an amplitude is not valid.}
#' }
#'
#' @param y data vector
#' @param x sampling vector
#' @param f frequeny vector
#' @param ofac in case \code{f=NULL} this value controlls the amount of frequency
#' oversampling.
#'
#' @return The \code{spec.lomb} function returns an object of the type \code{lomb},
#' which is a \code{list} containg the following parameters:
#' \describe{
#'  \item{\code{A}}{A vector with amplitude spectrum}
#'  \item{\code{f}}{corresponding frequency vector}
#'  \item{\code{phi}}{phase vector}
#'  \item{\code{x,y}}{original data}
#'  \item{p}{p-value as statistical measure}
#' }
#' @references
#' A. Mathias, F. Grond, R. Guardans, D. Seese, M. Canela, H. H. Diebner, und G. Baiocchi,
#' "Algorithms for spectral analysis of irregularly sampled time series",
#' Journal of Statistical Software, 11(2), pp. 1--30, 2004.
#' @references
#'  J. D. Scargle, "Studies in astronomical time series analysis. II - Statistical
#'  aspects of spectral analysis of unevenly spaced data", The Astrophysical Journal,
#'  263, pp. 835--853, 1982.
#'
#' @example R/examples/lombExample.R
#' @seealso \code{\link{filter.lomb}}
#' @export
spec.lomb <-
  function(y = stop("Missing y-Value"),
           x = stop("Missing x-Value"),
           f = NULL,
           ofac = 4)
  {
    f_not_defined  <-  is.null(f)
    fzero  <-  F
    if (length(x) != length(y))
      stop("x-y lengths differ")

    # delete "missing values"
    c   <-  is.finite(x) & is.finite(y)
    x_org  <-  x
    y_org  <-  y
    x  <-  x[c]
    y  <-  y[c]

    mean_y <-  mean(y)
    var_y  <-  var(y)
    y_  <-  y - mean_y

    # make x symmetric around 0
    x_ave  <-  (x[1] + tail(x, 1)) / 2
    x   <-  x - x_ave
    nx   <-  length(x)

    # minimal frequency
    f0 <- 1 / (2 * diff(range(x)) * 4)

    # mean sampling frequency
    fsm <- nx / (diff(range(x)))

    # set up frequency vector acording to Nyquists sampling theorem
    # band limitation is set to f[sampling]/2
    if (f_not_defined)
    {
      f <- seq(f0, 0.5 * fsm, length.out = (nx * ofac - 1))
      addMean <- 1
    }

    if (min(f) == 0)
    {
      # remove f == 0 and remember
      f <- f[-which(f == 0)]
      fzero = T
    }

    nf <- length(f)
    omega  <-  2 * pi * f

    # do lomb-scargle
    result <- lmb(x, y_, omega)

    # in case of f is not defined or f contains zero:
    # manually add DC value
    if (f_not_defined | fzero)
    {
      f <- c(0, f)
      A <- c(mean_y, result$amp)
      phi <- c(0, result$phase)
    } else{
      f <- f
      A <- result$amp
      phi <- result$phi
    }

    # M equals to the degrees of freedom OR the independent frequencies
    M <- nx / 2

    p <- M * exp(-(A ^ 2 * nx / (4 * var_y)))
    # to reduce approximation errors for small p
    p[is.na(p)] <- 1
    # large p can be calculated as follows
    cond <- p > 0.01
    p[cond] <- 1 - (1 - exp(-(A[cond] ^ 2 * nx / (4 * var_y)))) ^ M

    res <- list(
      f = f,
      A = A,
      phi = phi,
      x = x_org,
      y = y_org,
      p = p
    )
    class(res) <- "lomb"

    return(res)
  }
