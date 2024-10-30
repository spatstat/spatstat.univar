#'
#'    access to internal functions, for debugging
#'
#'   kermom()  should agree with kernel.moment()
#' 
#'  Copyright (c) 2023 Adrian Baddeley, Tilman Davies and Martin Hazelton
#'  GNU Public Licence (>= 2.0)

kermom <- function(m, r, kernel="gaussian", mean = 0,
                   sd=1/kernel.factor(kernel)) {
  kernel <- match.kernel(kernel)
  kerncode <- switch(kernel,
                     gaussian=1,
                     rectangular=2,
                     triangular=3,
                     epanechnikov=4,
                     biweight=5,
                     cosine=6,
                     optcosine=7,
                     0)
  df <- data.frame(r=r, mean=mean, sd=sd)
  nr <- nrow(df)
  z <- .C(SK_kermom,
          nx = as.integer(nr),
          x = as.double(df$r),
          mean = as.double(df$mean),
          sd = as.double(df$sd),
          m = as.integer(m),
          kerncode = as.integer(kerncode),
          y = as.double(numeric(nr)),
          errcode = as.integer(integer(1)),
          PACKAGE="spatstat.univar")
  if(z$errcode != 0) stop(paste("Error code", z$errcode))
  return(z$y)
}


          
