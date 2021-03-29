#' Summarize FFT objects
#'
#' The function summarizes properties from the \code{class(fft)} object.
#'
#' The false alarm probability threshold \code{p0} value can be
#' changed to modify the amount of significant peaks.
#'
#'
#' @param object \code{lomb} object
#' @param p0 False Alarm Probability (FAP) threshold, default 1\%
#' @param ... not used
#'
#' @return a list of significant values of the spectral analysis
#' @examples
#' # see spec.fft() example
#'
#' @export
summary.fft <- function(object, p0 = 0.01,...)
{
  out <- list()
  # determine the dimension of data
  out$nD <- which(names(object) == "A") - 1

  out$N <- length(object$A)

  # determine average sampling resolution
  out$resolution <- as.data.frame(lapply(object[1:out$nD],function(x) diff(range(x))))
  rownames(out$resolution) <- "SamplingFrequency"

  out$resolution <- rbind(out$resolution, minFreqStep = lapply(object[1:out$nD],function(x) min(abs(diff(unique(x))))))

  # correct sampling rate / effective Bandwidth
  out$resolution[1,] <- colSums(out$resolution)

      tmp <- lapply(object[1:out$nD], function(x) 1:length(x) )
      tmp <- expand.grid(tmp)# data.frame(i = 1:length(object$PSD), PSD = object$PSD, max = 0 * object$PSD)
    tmp$A <- apply(tmp,1,function(n) abs(object$A)[array(n,dim = c(1,out$nD))])
  tmp$PSD <- apply(tmp,1,function(n) object$PSD[array(n,dim = c(1,out$nD))])
    tmp$p <- apply(tmp,1,function(n) object$p[array(n,dim = c(1,out$nD))])
    tmp$i <- 1:length(object$PSD)
  tmp$max <- 0 * tmp$PSD

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
  tmp <- tmp[tmp$max >= out$nD - 1,]

  for(i in 1:out$nD)
  {
    tmp[,i] <- object[[i]][tmp[,i]]
  }

  out$Peaks <- tmp[,-(dim(tmp)[2] - 0:1)]
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
  cat("Summary of FFT-Object\n----------------------\n")
  cat("Dimension:",out$nD,"\n")
  cat("# of Freq:", length(object$A),"\n")

  print(out$resolution)

  cat("\nMost Dominant Peaks:\n")
  print(out$Peaks)

  message("Hint: The dominant frequencies are not optimized and therefore only an estimate!")

  return(invisible(out))
}
