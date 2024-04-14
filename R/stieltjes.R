#'    stieltjes
#'
#'    Stieltjes integration
#'
#'   Copyright (c) 2000-2023 Adrian Baddeley, Rolf Turner and Ege Rubak
#'   GNU Public Licence (>= 2.0)
#' 
#'    $Revision: 1.2 $ $Date: 2023/11/05 02:02:53 $


stieltjes <- function(f, M, ...) {
  ## stieltjes integral of f(x) dM(x)
  stopifnot(is.function(f))
  StieltjesCalc(M, f, ...)
}

StieltjesCalc <- function(M, f, ...) {
  UseMethod("StieltjesCalc")
}

StieltjesCalc.stepfun <- function(M, f, ...) {
  stopifnot(inherits(M, "stepfun"))
  envM <- environment(M)
  #' jump locations
  x <- get("x", envir=envM)
  #' values of integrand
  fx <- f(x, ...)
  #' jump amounts
  xx <- c(-Inf, (x[-1L] + x[-length(x)])/2, Inf)
  dM <- diff(M(xx))
  #' integrate f(x) dM(x)
  f.dM <- fx * dM
  result <- sum(f.dM[is.finite(f.dM)])
  return(list(result))
}

