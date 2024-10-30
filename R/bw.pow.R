#'
#'   bw.pow.R
#'
#'  Adaptive bandwidths proportional to x^POW
#'
#'  Copyright (c) 2021-2023 Tilman Davies, Martin Hazelton and Adrian Baddeley
#'  GNU Public Licence (>= 2.0)
#'

bw.pow <- function(X, h0, POW=0.75, trim=5, ...){
  check.nvector(X)
  
  if(missing(h0) || is.null(h0)) {
    h0 <- bw.nrd0(X)
  } else {
    check.1.real(h0)
    stopifnot(h0 > 0)
  }
  
  check.1.real(trim)
  stopifnot(trim > 0)
  
  # POW <- 0.75
  
  #' geometric mean of distances for scaling purposes
  gamma <- exp(mean(log(X^POW)))
  
  #' compute variable bandwidths to be proportional to distance
  bw <- h0 * pmin(X^POW/gamma,trim)
  
  return(bw)
}
