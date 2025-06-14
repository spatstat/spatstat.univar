#'
#'     weightedStats.R
#'
#'   weighted versions of hist, var, median, quantile
#'
#'  $Revision: 1.11 $  $Date: 2024/05/06 00:50:47 $
#'


#'
#'    whist      weighted histogram
#'

whist <- function(x, breaks, weights=NULL, method=c("C", "interpreted")) {
    N <- length(breaks)
    if(length(x) == 0) 
      h <- numeric(N+1)
    else {
      # classify data into histogram cells (breaks need not span range of data)
      cell <- findInterval(x, breaks, rightmost.closed=TRUE)
      # values of 'cell' range from 0 to N.
      nb <- N + 1L
      if(is.null(weights)) {
        ## histogram
        h <- tabulate(cell+1L, nbins=nb)
      } else {
        ##  weighted histogram
        method <- match.arg(method)
        switch(method,
               interpreted = {
                 cell <- factor(cell, levels=0:N)
                 h <- unlist(lapply(split(weights, cell), sum, na.rm=TRUE))
               },
               C = {
                 h <- .Call(SK_Cwhist,
                            as.integer(cell),
                            as.double(weights),
                            as.integer(nb),
                            PACKAGE="spatstat.univar")
               })
      }
    }
    h <- as.numeric(h)
    y <- h[2:N]
    attr(y, "low") <- h[1]
    attr(y, "high") <- h[N+1]
    return(y)
}

#' wrapper for computing weighted variance of a vector
#' Note: this includes a factor 1 - sum(v^2) in the denominator
#' where v = w/sum(w). See help(cov.wt)

weighted.var <- function(x, w, na.rm=TRUE) {
  bad <- is.na(w) | is.na(x)
  if(any(bad)) {
    if(!na.rm) return(NA_real_)
    ok <- !bad
    x <- x[ok]
    w <- w[ok]
  }
  cov.wt(matrix(x, ncol=1),w)$cov[]
}

#' weighted median

weighted.median <- function(x, w, na.rm=TRUE, type=4, collapse=TRUE) {
  unname(weighted.quantile(x, probs=0.5, w=w, na.rm=na.rm, type=type, collapse=collapse))
}

#' weighted quantile

weighted.quantile <- function(x, w, probs=seq(0,1,0.25), na.rm=TRUE, type=4, collapse=TRUE) {
  x <- as.numeric(as.vector(x))
  w <- as.numeric(as.vector(w))
  if(length(x) == 0)
    stop("No data given")
  stopifnot(length(x) == length(w))
  if(is.na(m <- match(type, c(1,2,4))))
    stop("Argument 'type' must equal 1, 2 or 4", call.=FALSE)
  type <- c(1,2,4)[m]
  if(anyNA(x) || anyNA(w)) {
    ok <- !(is.na(x) | is.na(w))
    x <- x[ok]
    w <- w[ok]
  }
  if(length(x) == 0)
    stop("At least one non-NA value is required")
  stopifnot(all(w >= 0))
  if(all(w == 0)) stop("All weights are zero", call.=FALSE)
  #'
  oo <- order(x)
  x <- x[oo]
  w <- w[oo]
  Fx <- cumsum(w)/sum(w)
  #' 
  if(collapse && anyDuplicated(x)) {
    dup <- rev(duplicated(rev(x)))
    x <- x[!dup]
    Fx <- Fx[!dup]
  }
  #'
  nx <- length(x)
  if(nx > 1) {
    result <- switch(as.character(type),
                     "1" = approx(Fx, x, xout=probs, ties="ordered", rule=2,
                                  method="constant", f=1)$y,
                     "2" = {
                       j <- approx(Fx, 1:nx, xout=probs, ties="ordered",
                                   rule=2, method="constant", f=0)$y
                       j <- as.integer(j)
                       g <- probs - Fx[j]
                       jplus1 <- pmin(j+1, nx)
                       ifelse(g == 0, (x[j]+x[jplus1])/2, x[jplus1])
                     },
                     "4" = approx(Fx, x, xout=probs, ties="ordered", rule=2,
                                  method="linear")$y)
  } else {
    result <- rep.int(x, length(probs))
  }
  names(result) <- paste0(format(100 * probs, trim = TRUE), "%")
  return(result)
}

