#' Setting up multithread BLAS library
#'
#' The number of used cores is set to \code{RhpcBLASctl::get_num_cores()} on the
#' attach event of the package.
#'
#' The algebra system of R relies on a BLAS library
#' which can be set to use many threads / cores. This
#' feature is considered as experimental since there are many differences across
#' the operating systems R is running on. If there is an issue and there is a need
#' to run R in multi-thread mode, consider to install
#' a different optimized version of BLAS.
#' If necessary, the number of cores required can also be changed manually
#' by calling \code{blas_set_num_threads(nCores)} and
#' \code{omp_set_num_threads(nCores)}.
#'
#'
#' This function is invoked automatically
#'
#' @param libname  a character string giving the library directory where the package defining the namespace was found.
#' @param pkgname  a character string giving the name of the package.
#'
#' @rdname onAttach
#' @import RhpcBLASctl
.onAttach <- function(libname,pkgname)
{
  # require(RhpcBLASctl)

  nCores <- get_num_cores()
  packageStartupMessage(paste("Detecting",nCores,"cores"))

  blas_set_num_threads(nCores)
  omp_set_num_threads(nCores)

  # detect cache size
  # .... to do
}

#' Reset multithread BLAS to default
#'
#' This function is invoked automatically
#'
#' @param libpath a character string giving the complete path to the package.
#'
#' @rdname onDetach
#' @import RhpcBLASctl
.onDetach <- function(libpath)
{
  blas_set_num_threads(1)
  omp_set_num_threads(1)
}

