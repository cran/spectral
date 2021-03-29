#' Filter in the frequency domain
#'
#' This function provides a method to band pass filter in the frequency domain.
#'
#' A signal \eqn{y} is meant to be equaly spaced and causal, which means it starts
#' at \eqn{t=0}. For times \eqn{y < 0} the signal is not defined. The filtering
#' itself takes place with the analytic function of \eqn{y} which provides an one sided
#' spectrum. Applying the Fourier transform, all properties of \eqn{y} will be
#' preserved.
#'
#' The band pass is represented throughout a function in the form of four different types, i.e.
#' "polynom", "sin(x)/x", "bi-cubic", "gauss". A detailed
#' description about these types can be found in \code{\link{BP}}.
#'
#' Setting \code{fc = 0} one can achieve a low pass filter.
#'
#' @param y numeric data vector
#' @param x optional x-coordinate
#' @param fc center frequency of the bandpass
#' @param BW bandwith of the bandpass
#' @param n parameter to control the stiffness of the bandpass
#' @param type type of weightening function: "poly", "sinc", "bi-cubic","gauss", can be abbreviated
#'
#' @examples
#' ## noisy signal with amplitude modulation
#' x <- seq(0,1, length.out=500)
#'
#' # original data
#' y_org <- (1+sin(2*2*pi*x))*sin(20*2*pi*x)
#'
#' # overlay some noise
#' y_noise <- y_org+rnorm(length(x),sd=0.2)
#'
#' # filter the noisy data
#' y_filt <- filter.fft(y_noise,x,fc=20,BW=4,n=50)
#'
#' # plot results
#' plot(x,y_noise,type="l",lwd=1,col="darkgrey",lty=2,ylab="y",main="Spectral filtering")
#' lines(x,y_org,lwd=5,col="grey")
#' lines(x,y_filt)
#' legend("topright",c("org","noisy","filtered"),col=c("grey","darkgrey","black")
#'         ,lty=c(1,2,1),lwd=c(5,1,1))
#' @export
filter.fft <-
  function(y = stop("y-value is missing"),x = NULL,fc = 0,BW = 0,n = 3,type = "poly")
  {
    if (!is.vector(y))
      stop("y must be a vector")
    if (!is.null(x) & !is.vector(x))
      stop("x must be a vector")

    # calculate analytical function
    y.ana <- analyticFunction(y)
    # calculate spectrum
    FT <- spec.fft(y = y.ana,x = x,center = F)
    # calculate weights for filtering
    w <- BP(FT$fx,fc,BW,n, type)

    # do the filtering
    FT$A <- w * FT$A

    # back transform
    Y <- spec.fft(FT)

    return(Y$y)
  }
