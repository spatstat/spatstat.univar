#'
#'     weightedStats.R
#'
#'   weighted versions of hist, var, median, quantile
#'
#'   whist
#'   weighted.var
#'   weighted.median
#'   weighted.quantile
#' 
#'  $Revision: 1.24 $  $Date: 2025/07/13 01:09:56 $
#'
#' --------------------------------------------------------
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
  if(missing(w) || is.null(w))
    w <- rep(1, length(x))
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

weighted.median <- function(x, w, na.rm=TRUE, type=2, collapse=FALSE) {
  if(missing(w) || is.null(w))
    w <- rep(1, length(x))
  unname(weighted.quantile(x, probs=0.5, w=w,
                           na.rm=na.rm, type=type, collapse=collapse))
}

#' weighted quantile

weighted.quantile <- function(x, w, probs=seq(0,1,0.25),
                              na.rm=TRUE, type=4, collapse=FALSE) {
  x <- as.numeric(as.vector(x))
  if(missing(w) || is.null(w)) {
    w <- rep(1, length(x))
  } else {
    w <- as.numeric(as.vector(w))
    stopifnot(length(x) == length(w))
  }
  if(length(x) == 0) {
    warning("No data values given; quantiles returned as NA", call.=FALSE)
    return(rep(NA_real_, length(probs)))
  }
  if(anyNA(x) || anyNA(w)) {
    if(na.rm) {
      ok <- !(is.na(x) | is.na(w))
      x <- x[ok]
      w <- w[ok]
      if(length(x) == 0) {
        warning("All data values are NA; quantiles returned as NA", call.=FALSE)
        return(rep(NA_real_, length(probs)))
      } 
    } else {
      warning("Data contain NA values; quantiles are NA", call.=FALSE)      
      return(rep(NA_real_, length(probs)))
    }
  }
  rs <- range(sign(w))
  if(rs[1L] <= 0) {
    #' some non-positive weights
    if(rs[1L] < 0) stop("Some weights are negative", call.=FALSE) 
    if(rs[2L] == 0) {
      warning("All weights are zero; quantiles are NA", call.=FALSE)
      return(rep(NA_real_, length(probs)))
    }
  }
  #' type of quantile
  type <- as.integer(type)
  supportedtypes <- 1:5
  documentedtypes <- 1:4
  if(is.na(m <- match(type, supportedtypes)))
    stop(paste("Argument", sQuote("type"),
               "must equal", commasep(documentedtypes, "or")),
         call.=FALSE)
  type <- supportedtypes[m]
  #' warn about experimental code
  experimental <- c(3, 5)
  if(!is.na(match(type, experimental)))
    warning(paste("Implementation of weighted quantile type", type,
                  "is experimental and may be changed"),
            call.=FALSE)
  #'
  oo <- order(x)
  x <- x[oo]
  w <- w[oo]
  Fx <- cumsum(w)/sum(w)
  #' 
  if(collapse && anyDuplicated(x)) {
    last <- rev(!duplicated(rev(x)))
    x <- x[last]
    Fx <- Fx[last]
  }
  #'
  nx <- length(x)
  if(nx > 1) {
    result <- switch(type,
    {
      #' 1
      approx(Fx, x, xout=probs, ties="ordered",
             rule=2, method="constant",
             f=1)$y
    },
    {
      #' 2
      j <- approx(Fx, 1:nx, xout=probs, ties="ordered",
                  rule=2, method="constant",
                  f=0)$y
      j <- as.integer(j)
      ## j is position immediately to left (or j=1)
      g <- probs - Fx[j]
      jplus1 <- pmin(j+1, nx)
      ifelse(g < 0, x[j],
                    ifelse(g == 0,
                          (x[j]+x[jplus1])/2,
                          x[jplus1]))
    },
    {
      #' 3
      FxLag <- Fx - diff(c(0, Fx))/2
      j <- approx(FxLag, 1:nx, xout=probs, ties="ordered",
                  rule=2, method="constant",
                  f=0)$y
      j <- as.integer(j)
      jplus1 <- pmin(j+1, nx)
      ##
      giszero <- !is.na(match(probs, FxLag))
      choice <- ifelse(giszero & (j %% 2 == 0), j, jplus1)
      x[choice]
    },
    {
      #' 4
      approx(Fx, x, xout=probs, ties="ordered", rule=2,
             method="linear")$y
    },
    {
      #' 5
      FxLag <- Fx - diff(c(0, Fx))/2
      approx(FxLag, x, xout=probs, ties="ordered", rule=2,
             method="linear")$y
    })
  } else {
    result <- rep.int(x, length(probs))
  }
  names(result) <- paste0(format(100 * probs, trim = TRUE), "%")
  return(result)
}

