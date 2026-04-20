#'
#' farebro.R
#'
#'   Farebrother's algorithm for the distribution of a linear combination
#'   of chi-squared random variables
#'
#'   R.W. Farebrother (1984)
#'   Algorithm AS 204:
#'   The distribution of a positive linear combination of \chi^2
#'   random variables. Applied Statistics 33 #3 (1984) 332--339
#'
#'   C implementation by Adrian Baddeley 2026
#'
#'   $Revision: 1.4 $ $Date: 2026/04/20 07:46:34 $

farebro <- function(lambda, mult=1, delta=0, x, ...,
                    maxit=1e4, eps = 1e-4, warn=TRUE) {
  stopifnot(all(lambda > 0))
  stopifnot(all(mult > 0))
  stopifnot(all(delta >= 0))
  df <- data.frame(lambda=lambda, mult=mult, delta=delta)
  tol <- log(.Machine$double.xmin)/2
  nx <- length(x)
  if(nx > .Machine$integer.max) stop("Cannot handle long vectors")
  if(has.zeroes <- (min(x) <= 0)) {
    xorig <- x
    positive <- (x > 0)
    x <- x[positive]
    nx <- length(x)
  }
  if(nx == 0) {
    result <- data.frame(x=numeric(0), p=numeric(0), d=numeric(0))
  } else {
    z <- .C(SK_farebro,
            lambda = as.double(df$lambda),
            mult   = as.integer(df$mult),
            delta  = as.double(df$delta),
            n      = as.integer(nrow(df)),
            x      = as.double(x),
            nx     = as.integer(nx),
            mode   = as.double(-1.0),
            maxit  = as.integer(maxit),
            eps    = as.double(eps),
            tol    = as.double(tol),
            ifault = as.integer(integer(nx)),
            density = as.double(numeric(nx)),
            probability = as.double(numeric(nx)),
            PACKAGE="spatstat.univar")
    ifault <- z$ifault
    if(warn && any(ifault != 0)) {
      pre <- "Farebrother algorithm returned error code:"
      for(falsch in unique(setdiff(ifault, 0))) 
        splat(pre, paren(falsch), .FarebrotherError(falsch))
    }
    result <- data.frame(x=x, p=z$probability, d=z$density, ifault=ifault)
  }
  if(has.zeroes) {
    fullresult <- data.frame(x=xorig, p=NA, d=NA, ifault=0)
    fullresult[positive,] <- result
    result <- fullresult
  }
  return(result)
}

.FarebrotherError <- function(i) {
  if(i < 0) {
    paste(.FBE[["negative"]], "for i = ", -i)
  } else if(!is.na(k <- match(as.character(i), names(.FBE)))) {
    .FBE[[k]]
  } else {
    paste("unrecognised error code", i)
  }
}

.FBE <- c(
  negative = paste("one or more of the constraints",
                   "lambda[i] > 0, mult[i] > 0 and delta[i] > 0",
                   "is not satisfied"),
  "0" = "normal exit",
  "1" = "non-fatal underflow of a0",
  "2" = paste("one or more of the constraints x > 0, maxit > 0 and eps > 0",
              "is not satisfied"),
  "3" = "current estimate of probability is less than -1",
  "4" = paste("required accuracy could not be attained",
            "within the specified maximum number of iterations"),
  "5" = "cumulative probability does not lie between 0 and 1",
  "6" = "probability density is negative"
)

.FBE[["9"]] <- paste(.FBE[[4]], "AND", .FBE[[5]])
.FBE[["10"]] <- paste(.FBE[[4]], "AND", .FBE[[6]])
