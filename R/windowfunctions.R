#' Windowfunctions
#'
#' Some typical windowfunctions are defined below:
#'
#' \describe{
#'    \item{\code{win.cos()}}{cosine window}
#'    \item{\code{win.tukey()}}{the tukey window}
#' }
#'
#' A window function weights a given dataset in a way, that the new data set is
#' coerced to be periodic. This method reduces the leakage effects of the Fourier
#' transform.
#'
#' @return All window functions return a wighting vector with the same length
#' as the provided data vector.
#'
#' @examples
#' y <- 1:100
#' y_cos <- y * win.cos(y)
#' y_tuk <- y * win.tukey(y)
#'
#' # Plot the original data
#' plot(y,main="Effect of window functions")
#' legend("topleft",c("original","cos","tukey"),pch=c(16,17))
#' points(y_cos,pch=16)
#' points(y_tuk,pch=17)
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
