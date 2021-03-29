#' Simple bandpass function
#'
#' This function represents a simple weightening procedure for spectral filtering accoring
#' to the type (\code{"poly", "sinc", "bi-cubic", "gauss"}) provided.
#'
#' The band pass is represented troughout a function in the form of four different types:
#'
#' 1. polynominial function \deqn{w = 1 - |((f - fc) / BW)|^n } with the degree \eqn{n}. The parameter \eqn{fc}
#' controlls the center frequency and desired band width \code{BW}. Outside the band width
#' \deqn{|f - fc| > BW} the result is forced to zero. With \code{n =  1.6} a quasi sinc-filter
#' without side bands can be constructed. A quasi rectangular window can be gained by setting
#' \code{n > 5}.
#'
#' 2. sinc function corresponds to a rectangular observation window in time domain
#' with \deqn{\Delta T ~ 1/BW}. It values ALL frequencies according to the si(x) function.
#' Calculation speed might be reduced.
#'
#' 3. bi-cubic encounters 2nd order interpolation kernel, providing a quasi rectangular observation window.
#'
#' 4. exponential Gauss curve. Here the band width is defined as the value of 90% damping.
#'
#'
#' @param f vector of frequencies
#' @param fc center frequency
#' @param BW bandwidth, with \code{w[ abs(f - fc) > BW ] == min}
#' @param n degree of the polynom, \code{n} can be real, e.g. \code{n = 1.6 as sinc alike}
#' @param type Type of weightening function: "poly", "sinc", "bi-cubic", "exp", can be abbreviated
#'
#' @return This function returns a weight vector [0..1], which is to apply to the frequency
#' vector \code{f} in a top level function
#' @export
#' @examples
#'
#' f <- seq(-50,50,by = 1e-2)
#' fc <- 0.3
#' BW <- 0.75
#'
#'
#' par(mfrow = c(2,1))
#'
#' curve(BP(x,fc = fc, BW = BW, type = "p"), -2,2, ylim = c(-0.2,1)
#'       ,main = "Filter weights"
#'       ,xlab = "fx",ylab = "w"
#' )
#' curve(BP(x,fc = fc, BW = BW, type = "s"), add = TRUE, lty = 2)
#' curve(BP(x,fc = fc, BW = BW, type = "b"), add = TRUE, lty = 3)
#' curve(BP(x,fc = fc, BW = BW, type = "g"), add = TRUE, lty = 4)
#'
#' abline(v = c(fc,fc+BW,fc-BW), lty = 3, col = "grey")
#'
#' # the corresponding Fourier-Transforms
#'
#' ty <- c("p","s","b","g")
#' A0 <- integrate(BP,fc = fc, BW = BW, type = "s",lower = -2,upper = 2)$value
#'
#' plot(NA,NA,xlab = "x", ylab = "|A|"
#'      ,main = "corresponding convolution kernels"
#'      ,xlim = 2*c(-1,1),ylim = c(0, sqrt(2)*A0/(length(f)*BW*min(diff(f))) )
#' )
#' for(i in 1:length(ty))
#' {
#'   FT <- spec.fft(y = BP(f,fc,BW,type = ty[i]))
#'   lines(FT$fx * length(FT$fx) / diff(range(f)),Mod(FT$A),lty = i)
#'
#' }
BP <- function(f, fc, BW, n = 3, type = "poly")
{
  type_org <- type
  w <- 0 * f

  if(BW == 0)
  {
    return(w)
  }

  type <- substr(tolower(type),1,1)

  # polynominial
  if(type == "p")
  {
    w <- 1 - abs( ((f - fc) / BW) )^n
    w[w < 0] <- 0
    return(w)
  }

  # sinc type
  if(type == "s")
  {
    x <- f / BW
    x <- x - fc / BW
    w <- sin(pi *  x)/(pi *  x)
    w[x == 0] <- 1
    return(w)
  }

  # bi-cubic
  if(type == "b")
  {
    x <- f / BW
    x <- x - fc/BW

    c <- abs(x) < 1
    w[c] <- 3/2 * abs(x[c])^3 - 5/2 * abs(x[c])^2 + 1

    c <- abs(x) >= 1 & abs(x) < 2
    w[c] <- -0.5 * abs(x[c])^3 + 5/2 * abs(x[c])^2 - 4*abs(x[c]) + 2

    return(w)
  }

  # exp.- Gauss
  if(type == "g")
  {
    w <- exp(-2.3025*((f - fc)/BW)^2)
    return(w)
  }

  stop("'",type_org, "' is not a propper band pass type")
}
