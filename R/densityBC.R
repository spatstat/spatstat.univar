#'
#'
#' densityBC.R
#'
#'  Kernel smoothing with boundary correction at zero
#'
#' Copyright (c) 2024 Adrian Baddeley, Tilman Davies and Martin Hazelton
#' GNU Public Licence (>= 2.0)
#' 

densityBC <- function(x, kernel="epanechnikov", bw=NULL, 
                    ...,
                    h=NULL,
                    adjust=1,
                    weights = rep(1, length(x))/length(x),
                    from=0, to=max(x), n=256,
                    zerocor=c("none", "weighted", "convolution",
                              "reflection", "bdrykern", "JonesFoster"),
                    fast=FALSE, internal=list()) {

  xname <- short.deparse(substitute(x))
  trap.extra.arguments(..., .Context="In densityBC()")

  zerocor <- match.arg(zerocor)
  ker <- match.kernel(kernel)

  ## internal option to suppress construction of $call 
  thecall <- if(isFALSE(internal$addcall)) NULL else match.call() 
  DEBUG <- isTRUE(internal$debug)

  if(DEBUG) {
    splat("\tdensityBC,", zerocor, ", h =", format(h), ", bw =", format(bw))
    splat("\t\trange(x) = ", prange(range(x)))
    started <- proc.time()
  }
  

  ## ...........   validate arguments  ..............................

  x <- as.vector(x)
  stopifnot(is.numeric(x))
  nx <- length(x)

  stopifnot(is.numeric(weights))
  if(length(weights) == 1) {
    weights <- rep.int(weights, length(x))
  } else stopifnot(length(weights) == length(x))

  if(!missing(adjust)) {
    check.1.real(adjust)
    stopifnot(adjust > 0)
  }

  if(is.null(from)) from <- min(0, x)
  if(is.null(to)) to <- max(0, x)
  
  stopifnot(is.numeric(from) && length(from) == 1)
  stopifnot(is.numeric(to) && length(to) == 1)
  stopifnot(from < to)
  stopifnot(is.numeric(n) && length(n) == 1)
  stopifnot(n >= 2)
  nr <- as.integer(n)

  
  ## ............... determine bandwidth ...........................

  if(!is.null(h) && !is.null(bw))
    stop("Only one of the arguments h and bw should be given")

  if(is.null(h) && is.null(bw))
    bw <- "nrd0"
  
  if(!is.null(h)) {
    if(is.character(h)) {
      bw <- h
      h <- NULL
    } else check.bandwidth(h, "the user-supplied bandwidth 'h'")
  }

  if(!is.null(bw)) {
    if(is.function(bw)) {
      ## bandwidth selection function
      bw <- bw(x)
      check.bandwidth(bw, "the bandwidth returned by bw(x)")
    } else if(is.character(bw)) {
      ## bandwidth selection rule -- copied from density.default
      if (nx < 2) 
        stop("need at least 2 points to select a bandwidth automatically")
      bwdescrip <- paste("the bandwidth returned by rule", sQuote(bw))
      bw <- switch(tolower(bw),
                   nrd0 = bw.nrd0(x),
                   nrd = bw.nrd(x), 
                   ucv = bw.ucv(x),
                   bcv = bw.bcv(x),
                   sj = ,
                   'sj-ste' = bw.SJ(x, method = "ste"),
                   'sj-dpi' = bw.SJ(x, method = "dpi"), 
                   stop(paste("unknown bandwidth rule", sQuote(bw)),
                        call.=FALSE))
      check.bandwidth(bw, bwdescrip)
    } else check.bandwidth(bw, "the user-supplied bandwidth 'bw'")
  }

  ## if bw given, determine h (or vice versa)
  cker <- kernel.factor(ker)
  if(!is.null(bw)) {
    h <- bw * cker
  } else {
    bw <- h/cker
  }

  ## finally, adjust to actual values
  h  <- h * adjust
  bw <- bw * adjust

  ## .............. initialise function table ......................
  r <- seq(from, to, length.out=nr)
  
  if(nx == 0) {
    f <- rep(0, nr)
    result <- list(x=r, y=f, bw=bw,
                   n=nr, call=thecall, data.name=xname, has.na=FALSE)
    class(result) <- c("density", class(result))
    return(result)
  }

  ## ........... Jones-Foster estimate --- combination of two methods -----

  if(zerocor == "JonesFoster") {
    estconv <- densityBC(x, kernel=kernel, bw=bw, weights = weights,
                       from=from, to=to, n=n,
                       internal=list(addcall=FALSE),
                       zerocor="convolution")
    estbdry <- densityBC(x, kernel=kernel, bw=bw, weights = weights,
                       from=from, to=to, n=n,
                       internal=list(addcall=FALSE),
                       zerocor="bdrykern")
    result <- estconv
    ecy <- estconv$y
    eby <- estbdry$y
    result$y <- ifelse(ecy <= 0, 0, ecy * exp(eby/ecy - 1))
    result$call <- thecall
    result$data.name <- xname
    if(DEBUG) {
      elapsed <- proc.time() - started
      splat("\tdensityBC returned", paren(zerocor), "time taken",
            elapsed[[3]], "sec")
    }
    return(result)
  }

  ## ......... start processing x ..................................
  
  ## divide by halfwidth; henceforth the kernel has unit halfwidth
  xscal <- x/h
  rscal <- r/h

  ## ......... rejig the data to implement boundary correction at r = 0 .....

  if(zerocor != "none") {
    if(any(x < 0))
      stop("negative x values are illegal when boundary correction selected")

    ## threshold is halfwidth of support of kernel
    thresh <- if(ker == "gaussian") 3 else 1

    if(zerocor=="weighted") {
      # identify x[i] whose kernels need renormalising
      x.is.small <- (xscal <= thresh)
      # divide weights[i] by right tail of kernel 
      if(any(x.is.small)) {
        mass <- pkernel(-xscal[x.is.small], ker, sd=1/cker, lower.tail=FALSE)
        weights[x.is.small] <- ifelse(weights[x.is.small] <= 0, 0,
                                      weights[x.is.small] / mass)
      }
    }

    if(zerocor == "reflection") {
      ## reflect input data about 0
      xscal <- c(xscal, -xscal)
      weights <- rep.int(weights, 2)
      nx <- length(xscal)
    }

    if(zerocor == "bdrykern" && fast) {
      ## split data according to whether boundary kernel is operative
      ## This occurs if abs(r-x) < thresh and abs(r-0) < thresh
      ## so can only be excluded when x > 2 * thresh
      x.is.large <- (xscal > 2 * thresh)
      xlarge <- xscal[x.is.large]
      wlarge <- weights[x.is.large]

      x.is.small <- ! x.is.large
      xsmall <- xscal[x.is.small]
      wsmall <- weights[x.is.small]
    }
  } 

  ## ............. compute density estimate ...................................

  ffast <- fslow <- 0
  
  if(fast) {
    ## use FFT for some or all of the calculation

    if(zerocor != "bdrykern") {
      d <- unnormdensity(xscal, weights=weights,
                         from=rscal[1], to=rscal[nr], n=nr,
                         kernel=kernel, bw=1/cker)
    } else {
      d <- unnormdensity(xlarge, weights=wlarge,
                         from=rscal[1], to=rscal[nr], n=nr,
                         kernel=kernel, bw=1/cker)
    } 
    ffast <- d$y
  } 
  
  if(!fast || zerocor == "bdrykern") {
    ## use the bespoke C code
    kerncode <- switch(ker,
                       gaussian=1,
                       rectangular=2,
                       triangular=3,
                       epanechnikov=4,
                       biweight=5,
                       cosine=6,
                       optcosine=7,
                       0)
    if(zerocor != "bdrykern") {
      ## standard fixed bandwidth kernel estimate with appropriately rigged data
      res <- .C(SK_fcolonel,
                kerncode=as.integer(kerncode),
                nx = as.integer(nx),
                x = as.double(xscal),
                w = as.double(weights),
                nr = as.integer(nr),
                r = as.double(rscal),
                f = as.double(numeric(nr)),
                errcode = as.integer(0),
                PACKAGE="spatstat.univar")
      fslow <- res$f
      if(DEBUG) {
        ## Try older, slightly slower implementation
        res2 <- .C(SK_colonel,
                   kerncode=as.integer(kerncode),
                   nx = as.integer(nx),
                   x = as.double(xscal),
                   w = as.double(weights),
                   nr = as.integer(nr),
                   r = as.double(rscal),
                   f = as.double(numeric(nr)),
                   errcode = as.integer(0),
                   PACKAGE="spatstat.univar")
        fslow2 <- res2$f
        splat("\tCalling C functions 'colonel' and 'fcolonel';")
        splat("\t\tdiscrepancy range", prange(range(fslow - fslow2)))
      }
    } else {
      ## boundary kernel
      nu0 <- kernel.moment(0, rscal, ker)
      nu1 <- kernel.moment(1, rscal, ker)
      nu2 <- kernel.moment(2, rscal, ker)
      ## safety check
      if(length(nu0) != nr || length(nu1) != nr || length(nu2) != nr)
        stop("internal error: kernel.moment did not yield result of correct length")
      ## 
      if(fast) {
        ## most data already processed. Only need to handle data close to 0.
        xx <- xsmall
        ww <- wsmall
        nn <- length(xsmall)
      } else {
        ## handle all data
        xx <- xscal
        ww <- weights
        nn <- nx
      }
      ## go
      res <- .C(SK_fbcolonel,
                kerncode=as.integer(kerncode),
                nx = as.integer(nn),
                x = as.double(xx),
                w = as.double(ww),
                nr = as.integer(nr),
                r = as.double(rscal),
                nu0 = as.double(nu0),
                nu1 = as.double(nu1),
                nu2 = as.double(nu2),
                a = as.double(numeric(nr)),
                b = as.double(numeric(nr)),
                f = as.double(numeric(nr)),
                errcode = as.integer(0),
                PACKAGE="spatstat.univar")
      fslow <- res$f
      if(DEBUG) {
        ## Try older, slightly slower implementation
        res2 <- .C(SK_bcolonel,
                   kerncode=as.integer(kerncode),
                   nx = as.integer(nn),
                   x = as.double(xx),
                   w = as.double(ww),
                   nr = as.integer(nr),
                   r = as.double(rscal),
                   nu0 = as.double(nu0),
                   nu1 = as.double(nu1),
                   nu2 = as.double(nu2),
                   a = as.double(numeric(nr)),
                   b = as.double(numeric(nr)),
                   f = as.double(numeric(nr)),
                   errcode = as.integer(0),
                   PACKAGE="spatstat.univar")
        fslow2 <- res2$f
        splat("\tCalling C functions 'bcolonel' and 'fbcolonel';")
        splat("\t\tdiscrepancy range", prange(range(fslow - fslow2)))
      }
    }

    ## check for errors
    if(res$errcode != 0) {
      ## error codes are defined in src/interfacecodes.h
      whinge <- switch(res$errcode,
                       "illegal length of vector",
                       "unrecognised kernel",
                       "unknown error code")
      stop(paste("Internal error in C function call:", whinge))
    }
  }

  ## combine contributions
  
  f <- ffast + fslow

  ## ...... post-process -------------------------------

  
  ## deal with effect of rescaling
  f <- f/h

  ## correct density estimates

  if(zerocor == "convolution") {
    r.is.small <- (rscal <= thresh)
    if(any(r.is.small)) {
        mass <- pkernel(-rscal[r.is.small], ker, sd=1/cker, lower.tail=FALSE)
        f[r.is.small] <- ifelse(f[r.is.small] <= 0, 0,
                                f[r.is.small] / mass)
      }
  }

  ## wrap up
  result <- list(x=r, y=f, bw=bw, n=nr,
                 call=thecall, data.name=xname, has.na=FALSE)
  class(result) <- c("density", class(result))
  if(DEBUG) {
    elapsed <- proc.time() - started
    splat("\tdensityBC returned", paren(zerocor), "time taken",
          elapsed[[3]], "sec")
  }
  return(result)
}
