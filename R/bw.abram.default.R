#'
#'   bw.abram.default.R
#'
#'   Abramson bandwidths for numeric data
#'
#'  Copyright (c) 2020-2023 Adrian Baddeley, Tilman Davies and Martin Hazelton
#'  GNU Public Licence (>= 2.0)

bw.abram <- function(X, h0, ...) {
  UseMethod("bw.abram")
}

bw.abram.default <- function(X, h0,
                             ..., 
                             at=c("data", "grid"), 
                             pilot=NULL, hp=h0, trim=5,
                             smoother=density.default){
  check.nvector(X)
  at <- match.arg(at)

  if(missing(h0) || is.null(h0)) {
    h0 <- bw.nrd0(X)
  } else {
    check.1.real(h0)
    stopifnot(h0 > 0)
  }

  check.1.real(trim)
  stopifnot(trim > 0)

  if(is.numeric(pilot) && at == "grid")
    stop("'pilot' must be a function or density estimate, when at='grid'")
  
  if(is.null(pilot)) {
    ## compute pilot density by smoothing X
    if(!missing(smoother)) {
      if(is.character(smoother)) {
        smoother <- get(smoother, mode="function")
      } else stopifnot(is.function(smoother))
    }
    pilot <- smoother(X, hp, ...)
    xx <- pilot$x
    px <- pilot$y
  } else if(is.numeric(pilot)) {
    check.nvector(pilot, length(X))
  } else if(inherits(pilot, "density")) {
    xx <- pilot$x
    px <- pilot$y
  } else if(is.function(pilot)) {
    stuff <- resolve.defaults(list(...),
                              list(from=min(X), to=max(X), n=512))
    xx <- with(stuff, seq(from, to, length.out=n))
    px <- pilot(xx)
  } else {
    stop("'pilot' should be a numeric vector, a function, or a density object")
  }
  
  #' evaluate pilot at data points
  if(is.numeric(pilot)) {
    pilot.X <- pilot
  } else if(is.function(pilot)) {
    pilot.X <- pilot(X)
  } else if(inherits(pilot, "density")) {
    pilot.X <- approx(xx, px, X, rule=2)$y
  } else {
    ## unreachable, but..
    stop("Format of 'pilot' is not recognised")
  }

  if(!all(is.finite(pilot.X)))
    stop("Evaluation of pilot density at X yielded NA, NaN or infinite values")
  
  #' clip bad values
  if(min(pilot.X) <= 0)
    pilot.X[pilot.X<=0] <- min(pilot.X[pilot.X>0])
  
  #' geometric mean re-scaler (Silverman, 1986; ch 5).  
  gamma <- exp(mean(log(pilot.X)))^(-0.5)

  #' evaluate Abramson bandwidths
  switch(at,
         data = {
           bw <- h0 * pmin((pilot.X^(-0.5))/gamma,trim)           
         },
         grid = {
           if(min(px) <= 0)
             px[px<=0] <- min(px[px>0])
           bb <- h0 * pmin((px^(-0.5))/gamma, trim)
           bw <- approxfun(xx, bb, rule=2)
         })

  return(bw)
}           
