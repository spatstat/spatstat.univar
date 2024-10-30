#'
#'   adaptive.R
#'
#'  Adaptive kernel smoothing
#'

densityAdaptiveKernel <- function(X, ...) {
  UseMethod("densityAdaptiveKernel")
}

densityAdaptiveKernel.default <-
  function(X, bw, ...,  weights=NULL,
           zerocor=c("none", "weighted", "convolution",
                     "reflection", "bdrykern", "JonesFoster"),
           at=c("grid", "data"),
           ngroups=Inf, fast=TRUE) {
  check.nvector(X)
  at <- match.arg(at)
  zerocor <- match.arg(zerocor)

  if(at == "grid")
    Xname <- short.deparse(substitute(X))

  nX <- length(X)
  if(nX == 0) 
    switch(at,
           data = return(numeric(nX)),
           grid = return(unnormdensity(X, ..., weights=0)))

  switch(zerocor,
         none = {
           from.default <- min(X)
         },
         {
           if(min(X) < 0)
             stop(paste("For boundary correction on the positive axis,",
                        "X must contain only nonnegative values"),
                  call.=FALSE)
           from.default <- 0
         })

  if(!missing(ngroups)) {
    if(is.null(ngroups)) {
      ngroups <- max(1L, floor(sqrt(nX)))
    } else {
      check.1.real(ngroups)
      if(is.finite(ngroups)) check.1.integer(ngroups)
      if(ngroups > nX) ngroups <- Inf
    }
  }

  if(at == "data" && ngroups == Inf) {
    warning("This case is not yet implemented; setting ngroups=nX")
    ngroups <- nX
  }
  if(weighted <- !is.null(weights)) {
    check.nvector(weights, nX, oneok=TRUE, vname="weights")
    if(length(weights) == 1) weights <- rep(weights, nX)
  } else weights <- rep(1/nX,nX)

  ## determine bandwidth for each data point
  if(missing(bw) || is.null(bw)) {
    ## default is Abramson rule
    bw <- bw.abram.default(X, at="data", ...)
  } else if(is.character(bw) && length(bw) == 1) {
    switch(bw,
           pow = ,
           bw.pow = {
             bw <- bw.pow(X,...)
           },
           abram = ,
           bw.abram = {
             bw <- bw.abram.default(X, at="data", ...)
           },
           stop(paste("Unrecognised bandwidth rule", sQuote(bw)), call.=FALSE))
  } else if(is.numeric(bw)) {
    check.nvector(bw, nX, oneok=TRUE, vname="bw")
    if(length(bw) == 1) bw <- rep(bw, nX)
  } else if(is.function(bw)) {
    ## adaptive bandwidth, function applied to data X
    bw <- bw(X)
    if(!is.numeric(bw))
      stop("The function bw() did not return a numeric vector", call.=FALSE)
    if(anyNA(bw))
      stop("Some computed bandwidths were NA", call.=FALSE)
  } else stop(paste("Argument 'bw' should be a numeric vector of bandwidths,",
                    "a function to compute bandwidths,",
                    "or a character string specifying the bandwidth rule"),
              call.=FALSE)

  if(zerocor == "JonesFoster") {
    #' Jones-Foster estimate involves two corrected estimates
    resC <- densityAdaptiveKernel.default(X, bw, ..., weights=weights,
                                          zerocor="convolution", at=at,
                                          ngroups=ngroups, fast=fast)
    resB <- densityAdaptiveKernel.default(X, bw, ..., weights=weights,
                                          zerocor="bdrykern", at=at,
                                          ngroups=ngroups, fast=fast)
    switch(at,
           data = {
             result <- resC * exp(resB/resC - 1)
           },
           grid = {
             result <- resC
             result$y <- resC$y * exp(resB$y/resC$y - 1)
             result$call <- match.call()
           })
    return(result)
  }

  if(ngroups == Inf && at == "grid") {
    ## use brute force C code
    stuff <- resolve.defaults(list(...),
                              list(from=from.default, to=max(X), n=512,
                                   kernel="gaussian"))
    kernel <- with(stuff, match.kernel(kernel))
    r <- with(stuff, seq(from, to, length.out=n))
    nr <- length(r)
    kerncode <- switch(kernel,
                       gaussian=1,
                       rectangular=2,
                       triangular=3,
                       epanechnikov=4,
                       biweight=5,
                       cosine=6,
                       optcosine=7,
                       0)
    ## go
    switch(zerocor,
           none = {
             res <- .C(SK_adaptiveKDE,
                       kerncode=as.integer(kerncode),
                       nx = as.integer(nX),
                       x = as.double(X),
                       sd = as.double(bw),
                       w = as.double(weights),
                       nr = as.integer(nr),
                       r = as.double(r),
                       f = as.double(numeric(nr)),
                       errcode = as.integer(0),
                       PACKAGE="spatstat.univar")
           },
           weighted = {
             res <- .C(SK_adaptKDEweight,
                       kerncode=as.integer(kerncode),
                       nx = as.integer(nX),
                       x = as.double(X),
                       sd = as.double(bw),
                       w = as.double(weights),
                       nr = as.integer(nr),
                       r = as.double(r),
                       f = as.double(numeric(nr)),
                       errcode = as.integer(0),
                       PACKAGE="spatstat.univar")
           },
           reflection = {
             res <- .C(SK_adaptKDEreflect,
                       kerncode=as.integer(kerncode),
                       nx = as.integer(nX),
                       x = as.double(X),
                       sd = as.double(bw),
                       w = as.double(weights),
                       nr = as.integer(nr),
                       r = as.double(r),
                       f = as.double(numeric(nr)),
                       errcode = as.integer(0),
                       PACKAGE="spatstat.univar")
           },
           convolution = {
             res <- .C(SK_adaptKDEconvol,
                       kerncode=as.integer(kerncode),
                       nx = as.integer(nX),
                       x = as.double(X),
                       sd = as.double(bw),
                       w = as.double(weights),
                       nr = as.integer(nr),
                       r = as.double(r),
                       f = as.double(numeric(nr)),
                       errcode = as.integer(0),
                       PACKAGE="spatstat.univar")
           },
           bdrykern = {
             res <- .C(SK_adaptKDEbdry,
                       kerncode=as.integer(kerncode),
                       nx = as.integer(nX),
                       x = as.double(X),
                       sd = as.double(bw),
                       w = as.double(weights),
                       nr = as.integer(nr),
                       r = as.double(r),
                       f = as.double(numeric(nr)),
                       errcode = as.integer(0),
                       PACKAGE="spatstat.univar")
           },
           stop("Internal error: zerocor not recognised")
           )
    if(res$errcode != 0)
      stop(paste("Internal error code", res$errcode))
    result <- structure(list(x=res$r,
                             y=res$f,
                             bw=exp(mean(log(bw))),
                             bwvalues=bw,
                             n=length(X),
                             call=match.call(),
                             data.name=Xname,
                             has.na=FALSE),
                        class=c("adaptivedensity", "density"))
  } else {
    #' divide bandwidths into groups and use FFT 
    if(ngroups == nX) {
      ## every data point is a separate group
      groupid <- 1:nX
      qmid <- bw
    } else {
      ## usual case
      p <- seq(0,1,length=ngroups+1)
      qbands <- quantile(bw, p)
      groupid <- findInterval(bw,qbands,all.inside=TRUE)
      #' map to middle of group
      pmid <- (p[-1] + p[-length(p)])/2
      qmid   <- quantile(bw, pmid)
    }

    group <- factor(groupid, levels=1:ngroups)
    Y <- split(X, group)
    W <- split(weights, group) 

    densityargs <- resolve.defaults(list(...),
                                    list(from=from.default, to=max(X), n=512,
                                         zerocor=zerocor, fast=fast))

    Z <- mapply(densityBC,
                x=Y,
                bw=as.list(qmid),
                weights=W,
                MoreArgs=densityargs,
                SIMPLIFY=FALSE)

    switch(at,
           data = {
             ftot <- numeric(nX)
             for(i in seq_along(Z)) {
               densi <- Z[[i]]
               ## interpolate
               fi <- approx(densi$x, densi$y, X)$y
               ## accumulate
               ftot <- ftot + fi
             }
             result <- ftot
           },
           grid = {
             fvalues <- sapply(Z, getElement, name="y")
             ftot <- if(length(Z) == 1) as.numeric(fvalues) else
                     .rowSums(fvalues, nrow(fvalues), ncol(fvalues))
             xvalues <- Z[[1]]$x
             result <- structure(list(x=xvalues, y=ftot,
                                      bw=exp(mean(log(bw))),
                                      bwvalues=bw,
                                      n=length(X),
                                      call=match.call(),
                                      data.name=Xname,
                                      has.na=FALSE),
                                 class=c("adaptivedensity", "density"))
           })
  }
  return(result)
}


plot.adaptivedensity <- function(x, ..., xlab) {
  if(missing(xlab)) 
    xlab <- paste("N =", x$n, "  Bandwidths ",
                  prange(formatC(range(x$bwvalues))))
  class(x) <- "density"
  plot(x, ..., xlab=xlab)
}

