#' 1D/2D spectrum of the Fourier transform
#'
#' This function calculates the Fourier spectrum of a given 1D or 2D data object.
#' It returns a user friendly object, which contains one or two frequency vectors
#' to map the complex amplitudes (vector or matirx) to the corresponding
#' frequencies. The output is already normalized and the frequencies can be seen
#' in terms of \eqn{1/\Delta x}-units.
#'
#' @param y 1D data vector, y coordinate of a 2D matrix or object of class \code{fft}
#' @param x x-coordinate of the data in y vector or z matrix
#' @param z optional 2D matrix
#' @param center logical parameter, if the spectrum should be centered or not
#' @param inverse logical parameter, if the back transformation should be performed
#'
#' @return An object of the type \code{fft} is returned. This contains the
#' original dataset and the corresponding spectrum, with "reasonable" frequency vectors.
#'
#' @seealso \link{plot.fft}
#' @example R/examples/specFFTExample.r
#' @export
spec.fft <- function(y = NULL,
                     x = NULL,
                     z = NULL,
                     center = T,
                     inverse = F)
{
  mode <- "1D"
  dimZ <- NULL
  fx <- NULL
  fy <- NULL

  # If an object of class "fft" is given, all variables are going
  # to be assigned properly
  if (class(y) == "fft")
  {
    if (y$mode == "waterfall")
      stop("Found time dependend analysis. Cannot work with that.")
    # x,y = spatial vectors, get defined if necessary
    if (!is.null(y$x) | !is.null(y$y))
    {
      mode   <- y$mode
      center <- y$center
      z <- y$z
      x <- y$x
      y <- y$y
    }
    # fx, fy frequency vectors are defined -> inverse
    if (!is.null(y$fx) | !is.null(y$fy))
    {
      mode    <- y$mode
      center  <- y$center
      inverse <- T

      if (mode == "1D")
      {
        x <- y$fx
        y <- y$A

      }
      if (mode == "2D")
      {
        z <- y$A
        x <- y$fx
        y <- y$fy
      }
    }
  }

  # If y is a matrix, correct the assignment of the variables.
  if (is.matrix(y))
  {
    z <- y
    y <- NULL
    x <- NULL
  }

  # check variables
  if (!is.null(z))
  {
    if (!is.matrix(z))
      stop("z is not a matrix")
    mode <- "2D"
    dimZ <- dim(z)
  }

  if (is.null(y))
  {
    ifelse(mode == "2D", y <- 1:dimZ[2], stop("y-value is missing"))
    if (inverse)
      y <- seq(0, 1, length.out = length(y))
  }
  if (is.null(x))
  {
    if (mode == "2D")
      x <- 1:dimZ[1]
    if (mode == "1D")
      x <- 1:length(y)
    if (inverse)
      x <- seq(0, 1, length.out = length(y))
  }

  nx <- length(x)
  ny <- length(y)

  if (mode == "2D")
  {
    if ((dimZ[1] == ny) & (dimZ[2] == nx) & (dimZ[1] != dimZ[2]))
    {
      z <- t(z)
      dimZ <- dim(z)
      warning("!!! I transposed z !!!")
    }
    if ((dimZ[1] != nx) &
        (dimZ[2] != ny))
      stop("x,y dimensions do not fit to z")
  }

  if ((mode == "1D") & (nx != ny))
    stop("x,y lengths differ")

  # end of variable checks

  if (center)
  {
    if (mode == "1D")
    {
      # centered FFT via binary modulation with fs/2
      if (!inverse)
        y <- y * (-1) ^ (1:ny)

      Ts <- min(diff(x))
      fx <- seq(-(nx) / (2 * Ts * nx), (nx - 2) / (2 * Ts * nx), length.out = nx)

      A <- fft(y, inverse = inverse) * 1 / ny
    }

    if (mode == "2D")
    {
      # centered spectrum: Gonzalez & Wintz (1977) Digital Image Processing p.53
      # is a modulation as well
      if (!inverse)
        z <- z * (-1) ^ (row(z) + col(z))

      Ts  <- min(diff(x))
      fx  <- seq(-(nx) / (2 * Ts * nx), (nx - 2) / (2 * Ts * nx), length.out = nx)

      Ty  <- min(diff(y))
      fy  <- seq(-(ny) / (2 * Ty * ny), (ny - 2) / (2 * Ty * ny), length.out = ny)

      A <- fft(z, inverse = inverse) / length(z)
    }
  }
  if (!center)
  {
    if (mode == "1D")
    {
      fx <- seq(0, (nx - 1) / (min(diff(x)) * nx), length.out = nx)
      A  <- fft(y, inverse = inverse) / ny
    }

    if (mode == "2D")
    {
      fx <- seq(0, (nx - 1) / (min(diff(x)) * nx), length.out = nx)
      fy <- seq(0, (ny - 1) / (min(diff(y)) * ny), length.out = ny)

      A <- fft(z, inverse = inverse) / length(z)
    }
  }

  if (!inverse)
  {
    if (mode == "1D")
      res <- list(
        fx = fx,
        A = A,
        mode = mode,
        center = center
      )
    if (mode == "2D")
      res <- list(
        fx = fx,
        fy = fy,
        A = A,
        mode = mode,
        center = center
      )
  }

  if (inverse)
  {
    x <- seq(0, (nx - 1) / diff(range(x)), length.out = nx)

    if (mode == "1D")
    {
      y <- A * ny # denormalization
      if (center)
        y <- y * (-1) ^ (1:ny)
      res <- list(
        x = x,
        y = y,
        mode = mode,
        center = center
      )
    }

    if (mode == "2D")
    {
      y <- seq(0, (ny - 1) / diff(range(y)), length.out = ny)
      A <-
        A * length(A) # denormalization
      if (center)
        A <- A * (-1) ^ (row(A) + col(A))
      res <- list(
        x = x,
        y = y,
        z = A,
        mode = mode,
        center = center
      )
    }

  }
  class(res) <- "fft"
  return(res)
}
