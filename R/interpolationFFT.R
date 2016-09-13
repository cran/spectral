#' interpolates data using the Fourier back transform
#'
#' There are two way to interpolate data from a given spectrum.
#' Frist, one can do zero padding to cover \code{n} new data points. Or, secound
#' the complex amplitude with the associated frequency is taken and evaluated
#' at given points \code{xout}. Doing that for all frequencies and amplitudes
#' will give the interpolation.
#'
#' @param y numeric data vector to be interpolated
#' @param x numeric data vector with reference points
#' @param n number of new points
#' @param xout a vector new points
#' @return A list with a \code{x} and \code{y} component is returned. The \code{e99}
#' value evaluates the error of the interpolation with respect to linear approximation
#' with the \code{approx()} function.
interpolate.fft<-function(y,x=NULL,n=NULL,xout=NULL)
{
  yout=NULL
  if(is.null(n) & is.null(xout)) stop("Missing values for n or xout")
  # Wenn n definiert -> oversampling
  if(!is.null(n))
  {
    ### calculate the offset of x ###
    # because the FFT contains the
    # mean part in the first bin,
    # we have to compensate that in
    # the spatial vector

    o<-(n/length(x)-1)*diff(range(x))/n
    xout<-seq(min(x),max(x)+o,length.out=n)

    # calculate the spectrum
    y<-analyticFunction(y)
    Y<-fft(y)/length(y)

    # doing zero padding
    if(n >= length(y))
      YOUT <- c(Y,rep(0,n-length(y)))
    if(n < length(y))
      YOUT <- Y[1:n]

    # reconstruct new data vector
    yout <- base::Re(fft(YOUT,inverse=T))

  }

  # in case the new sampling pints xout are defined
  # we have to do the inverse transform manually
  if(!is.null(xout) & is.null(n))
  {
    o <- min(x)
    FT <- spec.fft(y,x-o)
    # select the ranges for error approximation
    n <- length(xout)
    y <- y[which(x>min(xout) & x<max(xout))]
    x <- x[which(x>min(xout) & x<max(xout))]

    yout <- 0*xout
    for(i in 1:length(FT$A))
      yout <- yout-FT$A[i]*exp(1i*2*pi*FT$fx[i]*(xout-o))

    yout <- base::Re(yout)
  }

  y <- base::Re(y)
  # calculate the error with respect to linear aproximation
  # use the larger data set for linear interpolation
  if(n <= length(x))
    e<-abs(1-abs( yout / approx(x,y,xout=xout,rule=2)$y) )
  if(n > length(x))
    e<-abs(1-abs( y / approx(xout,yout,xout=x,rule=2)$y) )

  e[which(e>1)]<-1
  emax<-mean(e)+qt(1-0.05/2,length(e)-1)/sqrt(length(e))*sd(e)

  return(list(x=xout,y=yout,e99=emax))
}
