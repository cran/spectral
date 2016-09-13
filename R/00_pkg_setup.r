### Compiler settings for Package comilation on windows

# Sys.setenv("PKG_LIBS"="-lgomp -lpthread")
# Sys.setenv("PKG_CXXFLAGS"="-fopenmp")

#' @useDynLib spectral
#' @importFrom Rcpp sourceCpp
NULL

### Hoock's ###
.onUnload <- function (libpath) {
  library.dynam.unload("spectral", libpath)
}

