#' Lomb-Scargle Periodigram
#'
#' The Lomb-Scargle periodigram represents a statistical estimator for the
#' amplitude and phase at a given frequency. This function takes also multivariate
#' (n-dimensional) input data.
#'
#' Since the given time series does not need to be evenly sampled, the data
#' mainly consists of data pairs \code{x1, x2, x3, ...} (sampling points) and (one)
#' corresponding value \code{y}, which stores the realisation/measurement data.
#' As can be seen from the data definition above multivariate (n-dimensional)
#' input data is allowed and properly processed.
#'
#' Two different methods are implemented: the standard Lomb-Scargle method with
#'
#' \eqn{y(t) = a * cos(\omega (t - \tau)) + b * sin(\omega (t - \tau))}
#'
#' as model function and the generalized Lomb-Scargle (after Zechmeister 2009)
#' method with
#'
#' \eqn{y(t) = a * cos(\omega t) + b * sin(\omega t) + c}
#'
#' as model function, which investigates a floating average parameter \eqn{c}
#' as well.
#'
#' Both methods can be supplied by an artifical dense frequency vector \code{f}.
#' In conjunction with the resulting phase information the user might be able to
#' build a "Fourier" spectrum to reconstruct or interpolate the timeseries in equally
#' spaced sampling. Remind the band limitation which must be fulfilled for this.
#'
#' \describe{
#'  \item{f}{The frequencies should be stored in a 1D vector or -- in case of
#'   multivariate analysis -- in a \code{data.frame} structure to preserve variable names}
#'  \item{\code{ofac}}{If the user does not provide a corresponding frequency
#'  vector, the \code{ofac} parameter causes the function to estimate
#'  \deqn{nf = ofac*length(x)/2} equidistant frequencies.}
#'  \item{\code{p}-value}{The \code{p}-value (aka false alarm probability FAP)
#'  gives the probability, wheater the estimated amplitude is NOT significant.
#'  However, if \code{p} tends to zero the
#'  amplidutde is significant. The user must decide which maximum value is acceptable,
#'  until an amplitude is not valid.}
#' }
#'
#' If missing values \code{NA} or \code{NaN} appear in any column, the whole row
#' is excluded from calculation.
#'
#' In general the function calculates everything in a vectorized manner, which
#' speeds up the procedure. If the memory requirement is more than \code{maxMem},
#' the calculation is split into chunks which fit in the memory size.
#'
#'
#'
#' @param f optional frequency vector / data frame. If not supplied \code{f} is calculated.
#' @param x sampling vector or data frame \code{data.frame(x1, x2, x3, ...)}
#' @param y input data vector or data frame \code{data.frame(x1, x2, ..., val)}
#' @param ofac in case \code{f=NULL} this value controlls the amount of frequency
#' oversampling.
#' @param w weights for data. It must be a 1D vector.
#' @param mode \code{"normal"} calculates the normal Lomb-Scargle periodogram;
#' \code{"generalized"} calculates the generalized Lomb-Scargle periodogram including
#' floating average and weights.
#' @param maxMem sets the amount of memory (in MB) to utilize, as a rough approximate.
#'
#' @return The \code{spec.lomb} function returns an object of the class \code{lomb},
#' which is a \code{list} containg the following parameters:
#' \describe{
#'  \item{\code{A}}{A vector with amplitude spectrum}
#'  \item{\code{f}}{corresponding frequency vector}
#'  \item{\code{phi}}{phase vector}
#'  \item{\code{PSD}}{power spectral density normalized to the sample variance}
#'  \item{\code{floatAvg}}{floating average value only in case of
#'                          \code{mode == "generalized"}}
#'  \item{\code{x,y}}{original data}
#'  \item{p}{p-value as statistical measure}
#' }
#'
#' @references
#' A. Mathias, F. Grond, R. Guardans, D. Seese, M. Canela, H. H. Diebner,
#' and G. Baiocchi, "Algorithms for spectral analysis of irregularly sampled
#' time series", Journal of Statistical Software, 11(2), pp. 1--30, 2004.
#' @references
#' J. D. Scargle, "Studies in astronomical time series analysis. II - Statistical
#' aspects of spectral analysis of unevenly spaced data", The Astrophysical Journal,
#' 263, pp. 835--853, 1982.
#' @references
#' M. Zechmeister and M. Kurster, "The generalised Lomb-Scargle periodogram.
#' A new formalism for the floating-mean and Keplerian periodograms",
#' Astronomy & Astrophysics, 496(2), pp. 577--584, 2009.
#'
#' @concept Lomb Scargle
#' @concept Periodogram
#' @example R/examples/lombExample.R
#' @seealso \code{\link{filter.lomb}}
#' @export
#' @import lattice
spec.lomb <- function(x = NULL, y=stop("Missing y-Value"),f = NULL, ofac = 1
                      ,w = NULL, mode = "normal", maxMem = 100)
{
  f_not_defined <- is.null(f)
  fzero <- F

  dat <- NULL

  # put data into dataframe
  if(!is.data.frame(y) & (is.null(x) || is.vector(x)))
  {
    # space normalized periodogram
    if(is.null(x))
      x <- 1:length(y)

    if(!is.vector(x))
      stop("x is not an atomic vector")

    if(length(x)!=length(y)) stop("x-y lengths differ")

    # coerce data to data frame with last (rightmost column) as values
    dat <- data.frame(x = x, val = y)
  }
  if(!is.data.frame(y) & (is.data.frame(x)))
  {
    dat <- as.data.frame(cbind(x,y))
  }
  if(is.data.frame(y))
  {
    if(!is.null(x))
      warning("ignoring x-values, because y is data.frame including depending variables")
    dat <- y
  }

  names(dat)[dim(dat)[2]] <- "val"

  ### checking for weigtheing vector w
  if(is.null(w) | length(w) == 1)
  {
    w <- rep(1,dim(dat)[1])
  }

  x_org <- x
  y_org <- y

  # delete "missing values"
  c <- apply(dat,1,function(x) !any(is.na(x)))
  w <- w[c]
  dat <- dat[c , ]

  # force sum(w) == 1
  w <- w/sum(w)

  mean_val <- mean(dat$val)
  var_val  <- var(dat[,dim(dat)[2]])
  nt <- dim(dat)[1]

  ### calculating frequencies
  if(!is.null(f))
  {
    if( !is.data.frame(f) & length(f) == dim(dat)[2]-1 )
    {
      dim(f) <- c(1,length(f))
      f <- as.data.frame(f)
      names(f) <- paste("f",1:dim(f)[2],sep="")
      warning("coerced frequency to 1xN data.frame")
    }

    if( !( (is.vector(f) & is.vector(x) & is.vector(y)) |
           (is.data.frame(f) & is.data.frame(y)) )
      )
      stop("frequency vector does not fit to data dimensions")

    if(ifelse(is.vector(f), length(f),dim(f)[1]) == (dim(dat)[2]-1))
    {
      f <- t(f)
      warning("transposed frequency vector")
    }

  }
  else # calculating the frequencies
  {
    fzero <- T
    f <- apply(as.matrix(dat[,-dim(dat)[2]]),2,function(x)
    {
      # minimal frequency #
      f0 <- 1/(2*diff(range(x)))

      # Sampling frequency #
      fsm <- length(x)^(1/(dim(dat)[2]-1)) / (diff(range(x)))

      ### calculation of frequencies ###
      fmin <- ifelse(dim(dat)[2] == 2, f0, -0.5*fsm + f0)
      return( seq(fmin,0.5*fsm - f0,by = f0/ofac) )
    }
    )
    f <- as.data.frame(f)
    f <- expand.grid(f)
    names(f) <- ifelse(length(names(f)) == 1, "f", paste("f",1:dim(f)[2],sep=""))
  }

  # Lomb rechnen
  f <- as.matrix(f)
  dat <- as.matrix(dat)

  res <- NULL

  # estimate the memory usage and divide the analysis in appropriate fractions
  memSize <- 3 * max(dim(f)) * max(dim(dat)) * 8 / 1024^2
  nFrac <- ceiling(memSize / maxMem)
  fInd <- split(1:dim(f)[1], rep(1:nFrac, each = ceiling(dim(f)[1]/nFrac)
                                 ,length.out = dim(f)[1])
                )

  if(mode == "generalized")
  {
    Y <- sum(w*dat[,dim(dat)[2]])
    hYY <- sum(w*dat[,dim(dat)[2]]^2)
  }

  if(mode == "normal")
  {
    dat[,dim(dat)[2]] <- dat[,dim(dat)[2]] - mean_val
  }

  res <- lapply(fInd,function(ind)
  {
    if(mode == "generalized")
    {
      return( gLmb(f[ind,],dat = dat,w = w,Y = Y,hYY = hYY) )
    }

    if(mode == "normal")
    {
      return( lmb(f[ind,],dat = dat,var_val = var_val) )
    }
  })

  res <- do.call(rbind,res)
  res <- cbind(f,res)


  # column of a data frame is a list, so rbind() works better with a list
  if(fzero)
  {
    res <- rbind(as.list(numeric(dim(res)[2])),res)
  }

  if(min(abs(rowSums(f))) == 0) # if we find f == 0 exatly
  {
    # correct f == 0 value
    res$A[which(rowSums(f) == 0)] <- mean_val
    res$phi[which(rowSums(f) == 0)] <- NA
    res$PSD[which(rowSums(f) == 0)] <- 0
    res$p[which(rowSums(f) == 0)] <- NA
  }

  # False alarm Probability
  #
  # based on the normalized PSD in [0,1] after Zechmeister 2009
  #
  # Prob(p > p0) = (1 - p0)^((N-3)/2)
  #
  Prob <- (1 - res$PSD)^((nt-3)/2)

  # Degrees of freedom
  # after Press (1989,1992)
  # M <- nt/2

  # after Horn and Baliunas (1986)
  M <- -6.362 + 1.193*nt + 0.00098*(nt)^2

  # minimize calculation errors at small p
  res$p <- M*Prob
  res$p[is.na(res$p)] <- 1
  # 'large' p calculate like...
  res$p[res$p > 0.001] <- 1 - ( 1 - Prob[res$p > 0.001])^M


  res <- as.list(res)
  res$x <- x_org
  res$y <- y_org

  class(res) <- "lomb"

  return(as.list(res))
}

#' Lomb-Scargle estimation function
#'
#' calculates the standard Lomb-Scargle estimation. The calculation is vectorized
#' to enhance calculation speed.
#'
#'
#' @param f frequency
#' @param dat spatial vector including locations and values
#' @param var_val variance of the data
lmb <- function(f, dat, var_val)
{
  ### Input definitions ###
  #
  # f = c(f1, f2, f3, ...)
  # dat = data.frame(x1, x2, x3, ..., val)
  #
  #########################

  val <- dat[, dim(dat)[2]]
  ### Calculating phase (e.g. omega * time) ###

  # if f is a single vector
  if(is.null(dim(f)) & dim(dat)[2] > 2)
  {
    dim(f) <- c(1,length(f))
  }

    # matrix product takes sum over frequencies, if multivariate
  SOT <- 2 * pi * tcrossprod(dat[, -dim(dat)[2]], f) # takes the matrix product
  dSOT <- dim(SOT)

  tau <- atan2(  .colSums( sin(2 * SOT), dSOT[1], dSOT[2] )
               , .colSums( cos(2 * SOT), dSOT[1], dSOT[2] )
               ) / 2

  ### berechne Argument der Summen ###
  arg <- SOT - matrix(tau, nrow = dSOT[1], ncol = dSOT[2],byrow = T)

  cs <- cos(arg)
  ss <- sin(arg)
  rm(SOT,arg)

  ### calculate individual parts ###
  R <- .colSums(dat[,dim(dat)[2]] * cs, dSOT[1], dSOT[2])
  I <- .colSums(dat[,dim(dat)[2]] * ss, dSOT[1], dSOT[2])

  C <- .colSums(cs ^ 2, dSOT[1], dSOT[2])
  S <- .colSums(ss ^ 2, dSOT[1], dSOT[2])

  # plausibility test for cs = 1 and ss = 1
  # prevent diff by 0 like error
  C[which(C < 1e-20)] <- 1
  S[which(S < 1e-20)] <- 1

  l <- sqrt(R ^ 2 + I ^ 2)
  phi  <- -tau - atan2((I / l), (R / l))
  A  <- sqrt((R / C) ^ 2 + (I / S) ^ 2)

  # normalized Periodogram (Hocke 1998, Press 1992 and Zechmeister 2009)
  PSD <- dim(dat)[1]/ (4*var_val) * A^2 * 2 / (dim(dat)[1] - 1)

  return(data.frame(A = A, phi = phi, PSD = PSD))
}

#' generalized Lomb-Scargle estimation function
#'
#' calculates the generalized Lomb-Scargle estimation after Zechmeister et al. (2009)
#'
#' This method is based on the generalized approach
#'
#' \eqn{y(t) = a*cos(w*t) + b*sin(w*t) + c}
#'
#' which contains the floating average value \eqn{c} of the model function above.
#' The calculation is vectorized to enhance calculation
#' speed.
#'
#'
#' @param f frequency
#' @param dat spatial vector including locations and values
#' @param w vector of weights
#' @param Y weighted sum of values
#' @param hYY weighted sum of squared values
gLmb <- function(f, dat, w, Y, hYY)
{
  ### Input definitions ###
  #
  # f = c(f1, f2, f3, ...)
  # dat = data.frame(x1, x2, x3, ..., val)
  #
  #########################

  # val <- dat[, dim(dat)[2]]

  # if f is a single vector
  if(is.null(dim(f)) & dim(dat)[2] > 2)
  {
    dim(f) <- c(1,length(f))
  }
  ### Calculating phase (e.g. omega * time) ###

  # matrix product takes sum over frequencies, if multivariate
  SOT <- 2 * pi * tcrossprod(dat[, -dim(dat)[2]],f) # takes the matrix product
  dSOT <- dim(SOT)

  ## put the data into a Matrix
  cSOT <- cos(SOT)
  sSOT <- sin(SOT)
  rm(SOT)

  ms1 <- w * cSOT               # C
  ms2 <- w * sSOT               # S

  s1 <- .colSums(ms1,dSOT[1],dSOT[2])
  s2 <- .colSums(ms2,dSOT[1],dSOT[2])
  s3 <- .colSums( dat[, dim(dat)[2]] * ms1 ,dSOT[1],dSOT[2])  # YC
  s4 <- .colSums( dat[, dim(dat)[2]] * ms2 ,dSOT[1],dSOT[2])  # YS
  s5 <- .colSums( ms1 * cSOT ,dSOT[1],dSOT[2]) # CC
  s6 <- .colSums( ms1 * sSOT ,dSOT[1],dSOT[2]) # CS
#  s7 <- .colSums( ms2 * sSOT ,dSOT[1],dSOT[2]) # SS, ersetzt durch SS = 1 - CC
  rm(cSOT,sSOT,ms1,ms2)

  # Calculationg the coefficients
  YY <- hYY - Y^2
  YC <- s3 - Y*s1
  YS <- s4 - Y*s2
  CC <- s5 - s1^2
  CS <- s6 - s1*s2
  # SS <- s7 - s2^2 # 1-s[5] = SS = 1-CC
  SS <- 1 - s5 - s2^2 # 1-s[5] = SS = 1-CC

  D <- CC*SS - CS^2

  a <- ( YC * SS - YS * CS ) / D
  b <- ( YS * CC - YC * CS ) / D
  c <- Y - a*s1 - b*s2

  phi  <- -atan2(b, a)
  A  <- sqrt(a ^ 2 + b ^ 2)
  PSD <- 1/(YY*D)*(SS*YC^2 + CC*YS^2 - 2*CS*YC*YS)

  return(data.frame(A = A, phi = phi, PSD = PSD,floatAvg = c))
}
