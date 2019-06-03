#' plot method for Lomb-Scargle periodograms
#'
#' This method plots a standard Lomb-Scargle periodogram, which contains the
#' normalized power spectra \code{PSD} and the corresponding false alarm
#' probability \code{p}. For more details refer to Zechmeister et al. (2009).
#'
#' The \code{plot.lomb} function is a wrapper function for R's standard scatter \code{plot}
#' To switch off certain properties, simply overwrite the parameter. For example
#' \code{log = ""} will reset the plot axis back to non-log scale.
#'
#' @param x object of class \code{lomb}
#' @param main setting the title of the plot
#' @param FAPcol color of the FAP line
#' @param FAPlwd line width of the FAP line
#' @param FAPlim limits to the FAP
#' @param FAPlab label of the right vertical axis
#' @param FAPlty line type for the FAP graph
#' @param legend.pos position of the legend
#' @param legend.cex cex value for the legend
#' @param legend.on logical, wheater to draw a legend or not
#' @param legend.text legend text
#' @param legend.lwd line width
#' @param legend.lty line type
#' @param legend.col color vector of the legend elements
#' @param ... further parameters to the plot function
#'
#' @examples
#' # See spec.lomb
#'
#' @references
#' M. Zechmeister and M. Kurster, "The generalised Lomb-Scargle periodogram.
#' A new formalism for the floating-mean and Keplerian periodograms",
#' Astronomy & Astrophysics, 496(2), pp. 577--584, 2009.
#'
#' @seealso \code{\link{spec.lomb}}
#' @inheritParams graphics::plot.default
#' @inheritParams graphics::par
#' @import graphics
#' @export
plot.lomb <-
  function(x
           ,FAPcol = 1,FAPlwd = 1,FAPlty = "dashed" ,FAPlim = c(1,1e-3),FAPlab = "FAP"
           ,legend.pos = "topleft",legend.cex = 1,legend.on = T
           ,legend.text = c("Spectrum","False Alarm Probability")
           ,legend.lwd = NULL,legend.lty = NULL,legend.col = NULL
           ,xlab = "Frequency" ,ylab = "Normalized PSD",main="",...)
  {
    ## checking common plot parameters and set defauls values
    params <- list(...)

    # restore each plot parameter from par
    lty <- lwd <- col <- NULL

    for(n in c("lwd","lty","col"))
      if(!(n %in% names(params)))
        assign(n,par()[[n]])
    else
      assign(n,params[[n]])

    i <- which(x$f == 0)
    # set up defaults
    if(!("xlim" %in% names(params)))
    {
      assign("xlim",NULL)
      if(length(i) == 0)
        xlim = range(x$f)
      if (length(i) > 0)
        xlim = range(x$f[-i])
    }
    else
      assign("xlim",params$xlim)

    if(!("ylim" %in% names(params)))
      assign("ylim",c(0,max(x$PSD)))
    else
      assign("ylim",params$ylim)

    if(!("type" %in% names(params)))
      assign("type","l")
    else
      assign("type",params$type)

    if(!("log" %in% names(params)))
      assign("log","x")
    else
      assign("log",params$log)

    if (is.null(FAPlim))
      FAPlim <- c(1,min(x$p))

    # deleting x == 0 values in log-scale plot
    if(log == "x" & length(which(x$f == 0)) > 0)
    {
      i <- which(x$f == 0)
      x$f <- x$f[-i]
      x$A <- x$A[-i]
      x$p <- x$p[-i]
      x$PSD <- x$PSD[-i]
    }

    if(is.null(legend.lwd))
      legend.lwd <- c(lwd,FAPlwd)
    if(is.null(legend.lty))
      legend.lty <- c(lty,FAPlty)
    if(is.null(legend.col))
      legend.col <- c(col,FAPcol)


    par(mar = par()$mar * c(1,1,1,0) + c(0, 0, 0, par()$mgp[1] + 1))
    args <- list(
      x = x$f, y = x$PSD
      ,type = type
      ,log = log
      ,xlab = xlab ,ylab = ylab
      ,xlim = xlim
      ,ylim = ylim
      ,main = main
    )
    args <- append(args,params[!(names(params) %in% names(args))])
    do.call(plot,args)

    ### Plot the False Alarm Probability

    op <- par(no.readonly = TRUE)
    par(new = T)
    plot(
      x = x$f, y = x$p
      ,log = paste0(log,"y")
      ,type = "l"
      ,ylim = FAPlim
      ,xlim = xlim
      ,xaxt = "n",yaxt = "n"
      ,xlab = "",ylab = ""
      ,lty = FAPlty
      ,col = FAPcol
      ,lwd = FAPlwd
    )
    grid()
    axis(4);mtext(FAPlab,side = 4,line = par("mgp")[1])

    if (legend.on)
    {
      args <- list(
        legend.pos
        ,legend.text #,c("Spectrum","False Alarm Probability")
        ,lty = legend.lty
        ,lwd = legend.lwd
        ,col = legend.col
        ,cex = legend.cex
      )
      args <- append(args,params[(names(params) %in% names(formals(legend)))])
      do.call(legend,args)
    }
    par(op)

  }
