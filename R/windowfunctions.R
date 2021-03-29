#' Windowfunctions
#'
#' Some typical windowfunctions are defined below:
#'
#' \describe{
#'    \item{\code{win.cos()}}{cosine window}
#'    \item{\code{win.tukey()}}{Tukey window}
#'    \item{\code{win.hann()}}{Hann window}
#'    \item{\code{win.nutt()}}{Nutt window}
#' }
#'
#' A window function weights a given dataset in a way, that the new data set is
#' coerced to be periodic. This method reduces the leakage effects of the
#' discrete Fourier transform.
#'
#' @return All window functions return a wighting vector with the same length
#' as the provided data vector.
#'
#' @examples
#' y <- 1:100
#' y_cos <- y * win.cos(y)
#' y_tuk <- y * win.tukey(y)
#' y_han <- y * win.hann(y)
#'
#' # Plot the original data
#' plot(y,main="Effect of window functions")
#' legend("topleft",c("original","cos","tukey","han"),pch=c(1,16,17,18))
#' points(y_cos,pch=16)
#' points(y_tuk,pch=17)
#' points(y_han,pch=18)
#'
#'
#' @name Windowfunctions
NULL

#' Cosine window function
#'
#' This window function returns a vector of weights with means of a cosine window
#'
#' @param n data vector to be windowed
#' @seealso \code{\link{Windowfunctions}}
#' @export
win.cos <- function(n)
{
  M <- length(n)
  n <- 0:(length(n) - 1)
  w <- sin(pi * n / (M))
  return(w)
}

#' Tukey window function
#'
#' This window function returns a vector of weights with means of a
#' Tukey-window. In contrast to a cosine window this function is more steep
#' at the beginning and the end. And it is 1 in the middle.
#'
#' @param n data vector to be windowed
#' @param a width of the rising and falling edge as ratio of the total data length
#' @seealso \code{\link{Windowfunctions}}
#' @export
win.tukey <- function(n, a = 0.5)
{
  M <- length(n) - 1
  n <- 0:(length(n) - 1)
  w <- rep(0, M)
  c <- n < a * (M) / 2
  w[c] <- 0.5 * (1 + cos(pi * (2 * n[c] / (a * (
    M
  )) - 1)))
  c <- a * (M) / 2 <= n & n <= (M) * (1 - a / 2)
  w[c] <- 1
  c <- (M) * (1 - a / 2) <= n
  w[c] <- 0.5 * (1 + cos(pi * (2 * n[c] / (a * (
    M
  )) - 2 / a + 1)))
  return(w)
}

#' Hanning window function
#'
#' This window function returns a vector of weights with means of a
#' generlized Hann-window.
#'
#' @param n data vector to be windowed
#' @param a order of the window, default \code{a = 2}
#' @seealso \code{\link{Windowfunctions}}
#' @export
win.hann <- function(n, a = 2)
{
  M <- length(n) - 1
  n <- 0:(M)

  # w <- factorial(a)/factorial(2*a-1) * (1 - cos(2*pi*n / M))^a / 2^(a+1)
  # w <- 0.5*(1 - cos(2*pi*n/(M-1)))
  w <- sin(pi*n/M)^a

  # w <- 2^a*factorial(a)^2 / factorial(2*a)*(1 + cos(2*pi*n/M))^a
  return(w)
}

#' Nuttall window function
#'
#' This window function returns a vector of weights with means of a
#' Nuttall-window.
#'
#' This window function provides a continuous first derivative everywhere,
#' like the Hann window. Adopted from the idea of Hann this window consists of
#' up to 5 trigonometric polynominial terms, i.e.
#'
#' \deqn{w_{n} = a_1 - a_2  \cos(2\pi n/M) + a_3  \cos(4\pi n/M) - a_4 \cos(6\pi n/M)
#'               + a_5 \cos(8\pi n/M) }
#'
#' Different sets of coefficients:
#'
#' \tabular{ll}{
#'   \strong{Nuttall(Default)} \tab  \code{c(0.355768, 0.487396, 0.144232, 0.012604,0)} \cr
#'   \strong{Blackman-Nuttall} \tab  \code{c(0.3635819, 0.4891775, 0.1365995, 0.0106411,0)} \cr
#'   \strong{Blackman-Harris}  \tab  \code{c(0.35875, 0.48829, 0.14128, 0.01168,0)} \cr
#'   \strong{Flat-Top}         \tab  \code{c(0.211557895, 0.41663158, 0.277263158, 0.083578947, 0.006947368)} \cr
#' }
#'
#'
#'
#' @param n data vector to be windowed
#' @param a coefficients default \code{a = c(0.355768, 0.487396, 0.144232, 0.012604,0)}
#' @seealso \code{\link{Windowfunctions}}
#' @export
win.nutt <- function(n, a = c(0.355768, 0.487396, 0.144232, 0.012604,0))
{
  M <- length(n) - 1
  n <- 0:(M)

  w <- a[1] - a[2] * cos(2*pi*n/M) + a[3] * cos(4*pi*n/M) - a[4] * cos(6*pi*n/M) + a[5] * cos(8*pi*n/M)
  return(w)
}
