#'
#'   adaptive.R
#'
#'  Adaptive kernel smoothing
#'

densityAdaptiveKernel <- function(X, ...) {
  UseMethod("densityAdaptiveKernel")
}

