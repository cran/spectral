% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyticFunction.R
\name{analyticFunction}
\alias{analyticFunction}
\title{Analytic function}
\usage{
analyticFunction(x)
}
\arguments{
\item{x}{real valued data vector}
}
\value{
Complex valued analytic function
}
\description{
In general a causal real valued signal in time has negative
frequencies, when a Fourier transform is applied. To overcome this,
a complex complement can be calculated to compensate the negative
frequency spectrum. The result is called analytic signal or analytic
function, which provides a one sided spectrum.
}
\details{
An analytic function \eqn{xa} is composed of the real valued signal
representation \eqn{y} and its Hilber transform \eqn{H(y)} as the
complex complement \deqn{xa(t) = x(t)+i H(x(t))}.
In consequence, the analytic function has a one sided spectrum,
which is more natural. Calculating the discrete Fourier transform
of such a signal will give a complex vector, which is only non zero
until the half of the length. Every component higher than the half
of the sampling frequency is zero. Still, the analytic signal
and its spectrum are a unique representation of the original signal
\eqn{x(t)}. The new properties enables us to do certain filtering
and calculations more efficient in the spectral space compared to the
standard FFT approach. Some examples are:

\describe{
 \item{Filtering}{because the spectrum is one sided, the user must
 only modifiy values in the lower half of the vector. This strongly
 reduces mistakes in indexing.
 See \code{\link{filter.fft}}}

 \item{Envelope functions}{Since the Hilbert transform is a perfect phase shifter
 by pi/2, the envelope of a band limited signal can be calculated.
 See \code{\link{envelope}}}

 \item{Calculations}{Deriving and integrating on band limited discrete data becomes
 possible, without taking the symmetry of the discrete Fourier transform into
 account. The secound example of the \code{\link{spec.fft}} function calculates
 the derivative as well, but plays with a centered spectrum and its corresponding
 "true" negative frequencies}
}
A slightly different approach on the analytic signal can be found in R. Hoffmann
"Signalanalyse und -erkennung" (Chap. 6.1.2). Here the signal \eqn{x(t)} is split
into the even and odd part. According to Marko (1985) and Fritzsche (1995)
this two parts can be composed to the analytic signal, which lead to the
definition with the Hilbert transform above.
}
\references{
R. Hoffmann, Signalanalyse und -erkennung: eine Einfuehrung fuer
 Informationstechniker, Berlin; Heidelberg: Springer, 1998.

H. Marko, Systemtheorie: Methoden und Anwendungen fuer ein- und mehrdimensionale
 Systeme. 3. Aufl., Berlin: Springer, 1995.

G. Fritzsche, Signale und Funktionaltransformationen - Informationselektronik.
 Berlin: VEB Verlag Technik, 1985
}
