#' 1D/2D/nD (multivariate) spectrum of the Fourier transform
#'
#' This function calculates the Fourier spectrum and power spectral density
#' of a given data object. The dimension of the array can be of arbitary size
#' e. g. 3D or 4D.
#'
#' The function returns an user friendly object, which contains as much frequency
#' vectors as ordinates of the array. \code{spec.fft} provides the
#' ability to center the spectrum along multiple axis. The amplitude output is already
#' normalized to the sample count and the frequencies are given in terms of
#' \eqn{1/\Delta x}-units.
#'
#' @section Missing Values:
#'
#' Given a regualar grid \eqn{x_i = \delta x \cdot i} there might be missing values
#' marked with \code{NA}, which are treated by the function as 0's.
#' This "zero-padding" leads to a loss of signal energy being
#' roughly proportional to the number of missing values.
#' The correction factor is then \eqn{(1 - Nna/N)} as long as \eqn{Nna / N < 0.2}.
#' If the locations of missing values are randomly
#' distributed the implemented procedure workes quite robust. If correalted
#' gaps are present, the proposed correction is faulty and
#' scales wrong. This is because a convolution of the incomplete
#' sampling vector with the the signal takes place. An aliasing effect
#' takes place distorting the spectral content.
#'
#' To be compatible with the underlying Fourier transform, the amplitudes
#' are not affected by this rescaling.
#' Only the power spectral density (PSD) is corrected in terms of the energy
#' content, which is experimental for the moment.
#'
#' @param y 1D data vector, y coordinate of a 2D matrix, nD (even 2D) array
#'          or object of class \code{fft}
#' @param x x-coordinate of the data in \code{y} or \code{z}. If \code{y} is an array, \code{x} must be a named list \code{x = list(x = ..., y = ...)}.
#' @param z optional 2D matrix
#' @param center logical vector, indicating which axis to center in frequency space
#'
#' @return An object of the type \code{fft} is returned. It contains the
#' spectrum \code{A}, with "reasonable" frequency vectors along each ordinate. \code{psd} represents
#' the standardized power spectral density, [0,1]. The false alarm probability (FAP)
#' \code{p} is given similar to the Lomb-Scargle method, see \link{spec.lomb}.
#'
#' @seealso \link{plot.fft}
#' @example R/examples/specFFTExample.r
#'
#' @concept Fourier Transform
#'
#' @export
spec.fft <- function(y = NULL, # vector or nD-array
                     x = NULL, # list of ordinates
                     z = NULL, # optional in 2D case
                     center = T # single or vector
                    )
{
  mode <- ""
  res <- NULL # stores result
  inverse = F


  # If an object of class "fft" is given, all variables are going
  # to be assigned properly
  if (class(y) == "fft")
  {
    if (y$mode == "waterfall")
      stop("Found time dependend analysis. Cannot work with that.")

    # x,y = spatial vectors, get defined if necessary

    # looking for frequency vectors
    nf <- grep(pattern = "f",x = substr(names(y),1,1))

    # reconstructing ordinates
    x <- list()
    for(i in nf)
    {
      Tmax <- 1/min(diff(y[[i]]))
      x[[length(x) + 1]] <- seq(0,Tmax - 1/diff(range(y[[i]])),length.out = length(y[[i]]))
    }
    names(x) <- substr(names(y)[nf],2,nchar(names(y)[nf]))

    center <- y$center
    mode <- y$mode
    inverse <- T
    y <- y$A
  }

  ### Cases of input data
  #
  # 1D: x,y = vector z = NULL
  #       y = vector, x,z = NULL
  #
  # 2D: x,y = vector z = matrix
  #       y = matrix
  #
  # nD: x = List of ordinates, y = array, z = NUll
  #
  #
  if(mode(inverse) != "logical" | length(inverse) != 1)
    stop("Wrong datatype or length of inverse.")

  if(!inverse)
  {
    # 1D case
    if ( is.vector(y) & is.null(z) )
    {
      if(!is.vector(x))
      {
        if(!is.null(x))
          stop("Wrong input data on x")
        if(is.null(x))
          x <- 1:length(y)
      }
      x <- list(x = x)
      mode = "1D"
    }

    # 2D case
    if( is.matrix(y) & !is.list(x))
    {
      if(!is.null(x) | !is.null(z))
        stop("Wrong input data!")

      z <- y
      y <- NULL
    }

    if( is.matrix(z))
    {
      if(is.null(x))
        x <- 1:dim(z)[1]
      if(is.null(y))
        y <- 1:dim(z)[2]

      if(length(x) != dim(z)[1] | length(y) != dim(z)[2])
        stop("x,y dimensions do not fit to z!")

      x <- list(x = x, y = y)
      y <- z
      z <- NULL
    }

    # nD case
    if(is.array(y))
    {
      if(!is.null(x))
        if(length(dim(y)) != length(x))
          stop("nD array Dimensions do not fit to length of x")

      if(!is.null(z))
        warning("Ignoring z input!")

      # generate x ordinates for nD array
      if(is.null(x))
      {
        x <- list()
        for(i in dim(y))
          x[[length(x) + 1]] <- 1:i
        names(x) <- paste(1:length(x))
      }

      mode <- ifelse(length(dim(y)) == 2,"2D","nD")
    }

    # correct center variable
    if(mode(center) == "numeric")
    {
      tmp <- rep(F, length(dim(y)))
      tmp[center] <- T
      center <- tmp
    }
    if(mode(center) != "logical")
      stop("Wrong data type of center. A logical vector is assumed.")
    if(length(center) == 1)
      center <- rep(center,length(dim(as.array(y))))

    # end of variable checks

    ### calculate frequencies
    f <- lapply(x,function(x)
    {
      nx <- length(x)
      seq(0, (nx - 1) / (min(diff(x)) * nx), length.out = nx)
    }
    )
    names(f) <- paste("f",names(x),sep="")

    ### center axis

    for(j in which(center))
    {
      # centered FFT via binary modulation with fs/2
      # in nD case after Gonzalez & Wintz (1977) Digital Image Processing p.53
      # centering remains a modulation as well
      if(mode != "1D")
      {
        s <- paste(paste(rep(" ,",length(dim(y))-1),collapse = "")," ",sep = "")
        substr(s,2*(j-1) + 1,2*(j-1) + 1) <- "i"
        s <- paste("y[",s,"] <- y[",s,"] * (-1)^(i)")
        for(i in 1:dim(y)[j])
          eval(parse(text = s))
      }
      else
      {
        y <- y * (-1) ^ (1:length(y))
      }

      ### recalculate frequencies
      nf <- length(f[[j]])
      Ts <- min(diff(x[[j]]))
      f[[j]] <- seq(-(nf) / (2 * Ts * nf), (nf - 2) / (2 * Ts * nf), length.out = nf)
    }
  }

  ### Doing FFT ###
  N <- length(y)

  tmp <- y
  tmp[is.na(tmp)] <- 0

  Y <- fft(tmp, inverse = inverse)

  ### correct inverse or non inverse mode ###
  if(!inverse)
  {
    Nna <- sum(is.na(y))
    sdy <- sd(y,na.rm=T)

    res <- f
    res$A <- Y / N
    ### add standardized PSD ###
    res$PSD <- N / (N - 1) / ( 2 * sdy^2 ) * ( 2 * abs(res$A) / (1 - Nna / N) )^2

    ### add false Alarm Probability ###
    M <- -6.362 + 1.193 * N + 0.00098 * N^2
    tmp <- (1 - res$PSD)^((N - 3) / 2)
    res$p <- M * tmp
    res$p[is.na(res$p)] <- 1
    if(any(res$p > 0.001))
      res$p[res$p > 0.001] <- 1 - (1 - tmp[res$p > 0.001])^M

    res$mode <- mode
    res$center <- center

    class(res) <- "fft"
  }
  if(inverse)
  {
    # un-modulate incase of centering
    for(j in which(center))
    {
      # centered FFT via binary modulation with fs/2
      # in nD case after Gonzalez & Wintz (1977) Digital Image Processing p.53
      # centering remains a modulation as well
      if(mode != "1D")
      {
        s <- paste(paste(rep(" ,",length(dim(y))-1),collapse = "")," ",sep = "")
        substr(s,2*(j-1) + 1,2*(j-1) + 1) <- "i"
        s <- paste("Y[",s,"] <- Y[",s,"] * (-1)^(i)")
        for(i in 1:dim(y)[j])
          eval(parse(text = s))
      }
      else
      {
        Y <- Y * (-1) ^ (1:length(Y))
      }
    }

    if(mode == "1D")
    {
      res <- x
      res$y <- Y
    }
    if(mode == "2D")
    {
      res <- x
      res$z <- Y
    }
    else
    {
      res$x <- x
      res$y <- Y
    }


  }

  return(res)
}
