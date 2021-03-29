#' Plot \code{fft}-objects
#'
#' This is a wrapper function to plot \code{fft}-class objects.
#'
#' @param x Object of the class \code{fft}
#' @param ... further arguments to the plot functions
#' @import rasterImage
#' @seealso \code{\link{spec.fft}}
#' @examples
#' # See spec.fft
#' @export
plot.fft<-function(x,...)
{
  if(x$mode=="1D")
    plot(x$fx,base::Mod(x$PSD),...)
  if(x$mode == "2D")
    rasterImage::rasterImage2(x=x$fx,y=x$fy,z=base::Mod(x$A),...)
  if(x$mode == "waterfall")
    rasterImage::rasterImage2(x=x$x,y=x$fx,z=base::Mod(x$PSD)
                              ,zlab = "PSD",zlim = c(0,1),z.adj = c(0,0.5)
                              ,...)
}
