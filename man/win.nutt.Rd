% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/windowfunctions.R
\name{win.nutt}
\alias{win.nutt}
\title{Nuttall window function}
\usage{
win.nutt(n, a = c(0.355768, 0.487396, 0.144232, 0.012604, 0))
}
\arguments{
\item{n}{data vector to be windowed}

\item{a}{coefficients default \code{a = c(0.355768, 0.487396, 0.144232, 0.012604,0)}}
}
\description{
This window function returns a vector of weights with means of a
Nuttall-window.
}
\details{
This window function provides a continuous first derivative everywhere,
like the Hann window. Adopted from the idea of Hann this window consists of
up to 5 trigonometric polynominial terms, i.e.

\deqn{w_{n} = a_1 - a_2  \cos(2\pi n/M) + a_3  \cos(4\pi n/M) - a_4 \cos(6\pi n/M)
              + a_5 \cos(8\pi n/M) }

Different sets of coefficients:

\tabular{ll}{
  \strong{Nuttall(Default)} \tab  \code{c(0.355768, 0.487396, 0.144232, 0.012604,0)} \cr
  \strong{Blackman-Nuttall} \tab  \code{c(0.3635819, 0.4891775, 0.1365995, 0.0106411,0)} \cr
  \strong{Blackman-Harris}  \tab  \code{c(0.35875, 0.48829, 0.14128, 0.01168,0)} \cr
  \strong{Flat-Top}         \tab  \code{c(0.211557895, 0.41663158, 0.277263158, 0.083578947, 0.006947368)} \cr
}
}
\seealso{
\code{\link{Windowfunctions}}
}
