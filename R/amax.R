#' Local Maxima
#'
#' Determines all local maxima from a real valued vector.
#'
#' The purpose is to detect all local maxima in a real valued 1D vector.
#' If the first element \code{x[1]} is the global maximum, it is ignored,
#' because there is no information about the previous element. If there
#' is a plateau, the first edge is detected.
#'
#' @param x numeric vector
#'
#' @return returns the indicies of local maxima. If \code{x[1] = max}, then
#'         it is ignored.
#' @export
#'
#' @examples
#'
#' a <- c(1,2,3,2,1,5,5,4)
#' amax(a) # 3, 6
amax<-function(x)
{
  a1 <- c(0,x,0)
  a2 <- c(x,0,0)
  a3 <- c(0,0,x)
  e <- which((a1 >= a2 & a1 > a3)[2:(length(x))])
  if(!is.na(e[1] == 1))
    if(e[1]==1)
      e <- e[-1]
  if(length(e) == 0) e <- NaN
  return (e)
}
