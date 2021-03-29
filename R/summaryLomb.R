#' Summarize Lomb objects
#'
#' The function summarizes properties from the Lomb object.
#'
#' The false alarm probability threshold \code{p0} value will
#' adjust the number of peaks.
#'
#' The \code{effectiveBandWidth} describes the coverage of processed
#' frequencies by the \code{spec.lomb} function. If the ratio to
#' \code{averageSampling} is almost 2, then the Nyquist criterion can
#' be assumed to be fullfilled. If the ratio is much less than 2
#' then only a fraction of information is analysed.
#'
#' The \code{minFreqStep} is an estimate of the minimum frequency
#' step determined from the Lomb-Object.
#'
#' Average sampling is calculated from the median distance between
#' two spatial points.
#'
#' The possible frequency resolution originates also from the spatial
#' (temporal) input data by \code{1/(diff(range(x)))}, if \code{x} is
#' the spatial (temporal) coordinate.
#'
#' @param object \code{lomb} object
#' @param p0 False Alarm Probability threshold, default 1\%
#' @param ... not used
#'
#' @return a list of significant values of the spectral analysis
#' @examples
#' # see spec.lomb() example
#' @export
summary.lomb <- function(object,p0 = 0.01,...)
{
  out <- list()
  # determine the dimension of data
  out$nD <- which(names(object) == "A") - 1

  tmp <- NULL
  if(is.null(object$x))
  {
    tmp <- lapply(object$y[,1:(dim(object$y)[2] - 1)]
                   ,function(x) length(unique(x))
                  )
  }
  else
    tmp <- length(unique(object$x))
  out$N <- as.numeric(tmp)

  # determine average sampling resolution
  out$resolution <- as.data.frame(lapply(object[1:out$nD],function(x) diff(range(x))))
  rownames(out$resolution) <- "effectiveBandWidth"

  out$resolution <- rbind(out$resolution, minFreqStep = lapply(object[1:out$nD],function(x) min(abs(diff(unique(x))))))

  tmp <- NULL
  if(is.null(object$x))
  {
    tmp <- (lapply(object$y[,1:(dim(object$y)[2] - 1)]
                                            ,function(x) 1/median(diff(sort(unique(x))),na.rm = T)))
  }
  else
    tmp <- 1/median(diff(object$x),na.rm = T)

  out$resolution <- rbind(out$resolution, averageSampling = as.numeric(tmp))

  tmp <- NULL
  if(is.null(object$x))
  {
    tmp <- (lapply(object$y[,1:(dim(object$y)[2] - 1)]
                   ,function(x) 1/diff(range(unique(x),na.rm = T))))
  }
  else
    tmp <- 1/diff(range(object$x,na.rm = T))

  out$resolution <- rbind(out$resolution, frequencyResolution = as.numeric(tmp))

  # find the maximum amplitude in the spectrum
  # HINT: use linear interpolation to project on a regular grid
  # HINT: Project PSD on the cordinate axis. && all axis to get the maxima

  i <- which.max(object$PSD)

  # out$maxAmp <- as.data.frame(object[1:(out$nD + 4)])[i,]
  out$maxAmp <- as.data.frame(object[1:( which(names(object) == "p") )])[i,]
  rownames(out$maxAmp) <- ""

  tmp <- data.frame(i = 1:length(object$PSD), PSD = object$PSD, max = 0 * object$PSD)

  for(i in 1:out$nD)
  {
    j <- order(object[[i]])
    tmp2 <- split(tmp, object[[i]][j])

    tmp3 <- lapply(tmp2, function(x)
    {
      i <- which.max(x$PSD)
      return(x[i,])
    })
    tmp3 <- do.call(rbind,tmp3)

    i.max <- amax(tmp3$PSD)
    tmp$max[ tmp3$i[i.max] ] <- tmp$max[ tmp3$i[i.max] ] + 1
  }
  # tmp[tmp < out$nD] <- 0
  tmp <- tmp[tmp$max >= out$nD - 1,]

  out$Peaks <- data.frame(object[1:( which(names(object) == "p") )],row.names = NULL)[tmp$i,]
  out$Peaks <- out$Peaks[out$Peaks$p < p0,]
  out$Peaks <- out$Peaks[order(out$Peaks$p),]
  rownames(out$Peaks) <- NULL

  # exclude peaks
  i <- 0
  while(i < dim(out$Peaks)[1])
  {
    i <- i+1
    tmp <- apply(out$Peaks,1,function(x)
    {
      sqrt(sum(( (out$Peaks[i,1:out$nD] - x[1:out$nD]) / out$resolution["averageSampling",] * out$N )^2))
    })

    j <- which(tmp < 1)[-1]
    if(length(j) > 0)
      out$Peaks <- out$Peaks[-j,]
  }

  ## print summary
  cat("Summary of Lomb-Object\n----------------------\n")
  cat("Dimension:",out$nD,"\n")
  cat("# of Freq:", length(object[[1]]),"\n")

  print(out$resolution)

  cat("\nMost Dominant Frequency:\n")
  print(out$maxAmp)

  cat("\nMost Dominant Peaks:\n")
  print(out$Peaks)

  message("Hint: The dominant frequencies are not optimized and therefore only an estimate!")

  return(invisible(out))

}
