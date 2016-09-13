#' Filter and reconstruction of data analysed via spec.lomb
#'
#' Given an object of class \code{lomb}, this function allows the
#' reconstruction of the input signal using (a) a frequency selection
#' of single or multiple frequency (ranges), and/or (b) the most
#' significant peaks in the periodogram.
#'
#' To properly reconstruct the signal out of the calculated
#' \code{lomb}-object, three different methods are available, which are
#' controlled by the \code{filt}-argument.
#' \enumerate{
#' \item If \code{filt=NULL}, the most significant values in the (dense) spectrum
#' are used.
#'
#' \item If \code{filt=c(f1, .., fn)}, the given frequencies are used. The corresponding
#' phase is approximated.
#'
#' \item If \code{class(filt)=="matrix"}, each row of the 2 x n matrix defines a
#' frequency range. With in each range the "significant" frequencies are selected for
#' reconstruction.
#' }
#' Prior to the reconstruction the \code{filter.lomb}-function calculates the
#' most significant amplitudes and corresponding phases. As a measure to select
#' the "correct" frequencies, the \code{threshold} argument can be adjusted.
#' The corresponding phases of the underlying sine/cosine-waves are estimated by
#' one of the four following methods.
#' \enumerate{
#' \item \code{phase=="nextnb"}... use the phase of the bin of nearest neighbour.
#' \item \code{phase=="lin"}... linear interpolation between the two closest bins.
#' \item \code{phase=="lockin"}... principle of lock-in amplification, also known as
#' quadrature-demodulation technique.
#' \item \code{phase=="fit"}... non-linear least squares fit with \code{stats::nls}
#' }
#'
#' @param l lomb object
#' @param newx vector of new values at which the restored function is to be evaluated
#' @param threshold statistical threshold in terms of a standard deaviation of
#' the amplidudes. It determines which frequencies are used. Lower values give
#' more frequencies.
#' @param filt vector or matix of frequencies (ranges) in which to select the frequencies
#' @param phase set the method to determine the phase at a given frequency (moegliche werte???)
#' @return This function returns a list which contains the reconstruction according to the
#' \code{lomb}-object and \code{newx} for the given data \code{x} and \code{y}. The returned
#' object contains the following:
#' \describe{
#' \item{\code{x,y}}{reconstructed signal}
#' \item{\code{f,A,phi}}{used parameters from the \code{lomb}-object}
#' \item{\code{p}}{corresponding significance values}
#' }
#'
#' @import stats
#' @importFrom utils tail
#' @export
filter.lomb <-
  function(l = stop("No Lomb-Data"),
           newx = NULL,
           threshold = 6,
           filt = NULL,
           phase = "nextnb")
  {
    # test variables and validation
    if (!is.null(filt))
    {
      if (!is.vector(filt) &
          !is.matrix(filt))
        stop("Vector or Matrix for filt expected")
      if (is.matrix(filt))
        if ((dim(filt)[2]) > 2)
          stop("Filtermatrix should be a 2 x n matrix")
    }
    if (!is.finite(threshold) |
        length(threshold) > 1)
      stop("Threshold is supposed to be a single Number")
    if (threshold > 6)
      warning("Observig more than 6 sigma!")
    if (!(phase %in% c("nextnb", "lin", "lockin", "fit")))
      stop("Wrong phase argument. Use 'nextnb', 'lin', 'lockin' or 'fit'.")
    # test the lomb object for completeness
    if (sum(unlist(lapply(c("f", "A", "p", "x", "y"), function(x, l)
      is.null(l[[x]]), l))) > 0)
      stop("Missing data in Lomb-list [f, A, phi, p, x, y]")

    # Reconstruct time vector of signal
    if (is.null(newx))
      newx <-
        seq(min(l$x, na.rm = T), max(l$x, na.rm = T), length.out = length(l$f))
    tave <- (newx[1] + tail(newx, 1)) / 2
    newx <-
      newx - tave # for the reconstruction newx must be symmetrical around 0

    ### EITHER: filter frequencies and determine maxima of the peaks ###
    f <- l$f[l$f > 0]
    A <- l$A[l$f > 0]
    p <- l$p[l$f > 0]
    phi <- l$phi[l$f > 0]


    ##### OR: find significant frequencies and filter them ####
    # prerequisite: we need enough frequencies calculated by spec.lomb for valid results
    if (length(f) > 10 * length(filt))
    {
      ### find all peaks using local maxima ###
      a1 <- c(NA, diff(A)) > 0
      a2 <- c(diff(A), NA) < 0
      a3 <- (a1 & a2)

      Apk <- A[a3]
      Apk <- Apk[!is.na(Apk)]
      fpk <- f[a3]
      fpk <- fpk[!is.na(fpk)]
      ppk <- phi[a3]
      ppk <- ppk[!is.na(ppk)]

      ## sort the amplitudes ##
      index <- rev(order(Apk))
      sApk <- Apk[index]
      sfpk <- fpk[index]

      # if there is more than one peak:
      # 1.) eliminate side peaks
      # 2.) test against noise
      if (length(index) > 1)
      {
        ## calculate variance of noise ##
        vsApk <-
          sd(sApk[round(0.05 * length(sApk)):length(sApk)], na.rm = T)
        n <- threshold
        # for small frequencies noise is small as well (uncertainty relation)
        fmin <-
          abs(3 / diff(range(l$x, na.rm = TRUE)))
        for (i in 1:(length(Apk)))
        {
          if (Apk[i] < n * (1 - 0.5 * (fmin > fpk[i])) * vsApk)
            Apk[i] <- 0
        }
        ch <- which(Apk > 0)

        ## False Friends -- eliminate side peaks ##
        tmp <- NULL
        for (i in ch)
        {
          falsefriend <- FALSE
          for (j in ch)
            if (i != j)
              if ((abs(fpk[i] - fpk[j]) < 0.1 * fpk[i] &
                   Apk[i] < 0.6 * Apk[j]))
                falsefriend <- TRUE
              if (!falsefriend)
                tmp <- c(tmp, i)
        }
        index <- tmp
      }
      # select amplitudes #
      A <- Apk[index]
      # select frquencies #
      f <- fpk[index]

      p <- p[a3]
      p <- p[!is.na(p)][index]

      phi <- ppk[index]
    }

    ### determine phase for signal reconstruction ###

    ### next neighbour ###
    if (phase == "nextnb")
    {
      phi <- NULL
      for (i in f)
        # get phase for the frequency bin closest to desired frequency
        phi <- c(phi, l$phi[abs(l$f - i) == min(abs(l$f - i))])
    }


    ### Phase linear approximation ###
    if (phase == "lin")
      phi <- stats::approx(l$f, l$phi, xout = f)$y

    ### quadrature demodulation: phase estimation ###
    if (phase == "lockin")
    {
      phi <- NULL
      for (i in 1:length(f))
      {
        ### correct the length of x ###
        # increases the precision of the estimation
        k <- pi / 4
        npi_max <-
          as.integer((2 * pi * f[i]) * max(l$x) / k) # how many times f will fit into data
        npi_min <- as.integer((2 * pi * f[i]) * min(l$x) / k)

        index <- which(l$x < (npi_max * (k) / (2 * pi * f[i]))
                       & l$x > (npi_min * (k) / (2 * pi * f[i]))) ## calculate index
        if (length(index) < 4)
          index <- 1:length(x)

        # calculate COS for lock-in
        arg1 <-
          2 / A[i] * mean(cos(2 * pi * f[i] * l$x[index]) * l$y[index], na.rm = T)
        if (abs(arg1) > 1)
          arg1 <- sign(arg1)

        # calculate SIN for lock-in
        arg2 <-
          2 / A[i] * mean(sin(2 * pi * f[i] * l$x[index]) * l$y[index], na.rm = T)
        if (abs(arg2) > 1)
          arg2 <- sign(arg2)

        N <- sqrt(arg1 ^ 2 + arg2 ^ 2)

        phi <- c(phi, -atan2(arg2 / N, arg1 / N))
      }
    }

    ### fit the phase to the signal ###
    if (phase == "fit")
    {
      phi <- NULL
      for (i in 1:length(f))
      {
        # next neighbour as starting value for optimization
        tmp_phi <- l$phi[which.min(abs(l$f - f[i]))] + pi / 2

        # non-linear optimization
        nfit <-
          stats::nls(l$y ~ A[i] * cos(2 * pi * f[i] * l$x + P), start = list(P =
                                                                               tmp_phi))
        phi <- c(phi, coef(nfit)[1])
      }
    }


    # frequency selection according to filt
    if (!is.null(filt))
    {
      tmp <- NULL
      if (is.vector(filt))
        # reconstruct next neighbour to filt
      {
        filt <- sort(filt)
        for (i in filt)
          tmp <- c(tmp, f[abs(f - i) == min(abs(f - i))])

      }

      if (is.matrix(filt))
      {
        filt <- matrix(filt[order(filt[, 1]), ], ncol = 2) # sort using column 1
        for (i in 1:dim(filt)[1])
          tmp <- c(tmp, f[f > filt[i, 1] & f < filt[i, 2]])
      }

      A <- A[f %in% tmp]
      phi <- phi[f %in% tmp]
      p <- p[f %in% tmp]
      f <- f[f %in% tmp]
    }

    # do the reconstruction
    rel <- NULL
    for (x in newx)
    {
      rel <- c(rel, sum(A * (cos(2 * pi * f * x + phi))))
    }
    # correct offset
    rel <- rel + l$A[l$f == 0]
    return(list(
      x = newx + tave,
      y = rel,
      f = c(0, f),
      A = c(l$A[l$f == 0], A),
      phi = c(0, phi),
      p = c(0, p)
    ))
  }
