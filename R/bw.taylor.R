#'
#'  bw.taylor.R
#'
#'  Bandwidth selection for univariate data
#'
#'  Copyright (c) 2024 Tilman M Davies and Adrian Baddeley
#'  GNU Public Licence (>= 2.0)
#' 
#'  $Revision: 1.4 $ $Date: 2024/09/24 03:09:56 $
#' 

bw.taylor <- local({

  bw.taylor <- function(x, ..., srange=NULL, useC=TRUE) {
    x <- as.numeric(as.vector(x))
    n <- length(x)
    if(n < 2) return(NA)
    dr <- diff(range(x))
    if(dr == 0) return(NA)
    if(is.null(srange)){
      srange <- dr * c(1/n, 2/sqrt(n))
    } else {
      check.range(srange)
      srange <- pmax(srange, sqrt(.Machine$double.eps))
    }
    if(useC) {
      z <- optimise(unibootC,interval=srange,x=x,n=n)
    } else {
      z <- optimise(uniboot,interval=srange,x=x,n=n)
    }
    result <- z$minimum
    return(result)
  }

  uniboot <- function(h,x,n,ij=FALSE){
    ijd <- outer(x,x,"-")^2
    if(!ij) diag(ijd) <- NA
    s1 <- sum(exp(-ijd/(8*h^2)),na.rm=TRUE)
    s2 <- sum(exp(-ijd/(6*h^2)),na.rm=TRUE)
    s3 <- sum(exp(-ijd/(4*h^2)),na.rm=TRUE)
    return((1/(2*n^2*h*sqrt(2*pi)))*(s1-4/sqrt(3)*s2+sqrt(2)*s3+n*sqrt(2)))
  }

  unibootC <- function(h, x, n, ij=FALSE) {
    diagok <- if(ij) 1L else 0L
    z <- .C(SK_taylorboot,
            x = as.double(x),
            n = as.integer(length(x)),
            h = as.double(h),
            diagok = as.integer(diagok),
            value = as.double(numeric(1)),
            PACKAGE="spatstat.univar")
    value <- z$value
    return(value)
  }

  bw.taylor
})

