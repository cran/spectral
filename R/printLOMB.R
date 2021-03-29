#' Lomb-Plotting Function
#'
#' It calls the summary function.
#'
#' @param x lomb object
#' @param ... not used
#'
#' @return This function returns nothing
#' @export
#'
#' @examples
#' # see summary.lomb() function
print.lomb <- function(x,...)
{
  summary(x)
}
