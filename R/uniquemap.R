#'
#'   uniquemap.R
#'
#'   Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2019
#'   Licence: GNU Public Licence >= 2
#'
#'   $Revision: 1.20 $  $Date: 2024/06/26 08:50:53 $

uniquemap <- function(x) { UseMethod("uniquemap") }

uniquemap.default <- function(x) {
  result <- seqn <- seq_along(x)
  if(length(x) <= 1) return(result)
  if(is.atomic(x) && (is.factor(x) || (is.vector(x) && is.numeric(x)))) {
    if(is.factor(x)) x <- as.integer(x)
    o <- order(x, seqn)
    isfirst <- c(TRUE, (diff(x[o]) != 0))
    omap <- cumsum(isfirst)
    result <- seqn
    result[o] <- o[isfirst][omap]
    return(result)
  }
  dup <- duplicated(x)
  ux <- x[!dup]
  mapdup <- match(x[dup], ux)
  result[dup] <- which(!dup)[mapdup]
  return(result)
}

uniquemap.matrix <- function(x) {
  n <- nrow(x)
  result <- seqn <- seq_len(n)
  if(n <= 1)
    return(result)
  #' faster algorithms for special cases
  nc <- ncol(x)
  if(nc == 1L) return(uniquemap(x[,1]))
  if(is.numeric(x)) {
    if(nc == 2L) {
      oo <- order(x[,1], x[,2], seqn)
      xx <- x[oo, , drop=FALSE]
      isfirst <- c(TRUE, (diff(xx[,1]) != 0) | (diff(xx[,2]) != 0))
    } else {
      ## y <- asplit(x, 2) would require R 3.6.0
      y <- split(as.vector(x), factor(as.vector(col(x)), levels=1:nc))
      oo <- do.call(order, append(unname(y), list(seqn)))
      xx <- x[oo, , drop=FALSE]
      isfirst <- c(TRUE, matrowany(apply(xx, 2, diff) != 0))
    }
    uniqueids <- seqn[oo][isfirst]
    lastunique <- cumsum(isfirst)
    result[oo] <- uniqueids[lastunique]
    return(result)
  }
  #' non-numeric matrix e.g. character
  if(!anyDuplicated(x))
    return(result)
  dup <- duplicated(x)
  uni <- which(!dup)
  for(j in which(dup)) {
    for(i in uni[uni < j]) {
      if(IdenticalRowPair(i, j, x)) {
        result[j] <- i
        break
      }
    }
  }
  return(result)
}

uniquemap.data.frame <- function(x) {
  n <- nrow(x)
  result <- seqn <- seq_len(n)
  if(n <= 1)
    return(result)
  #' faster algorithms for special cases
  nc <- ncol(x)
  if(nc == 1L) return(uniquemap(x[,1]))
  if(all(sapply(x, is.numeric))) {
    if(nc == 2L) {
      oo <- order(x[,1], x[,2], seqn)
      xx <- x[oo, , drop=FALSE]
      isfirst <- c(TRUE, (diff(xx[,1]) != 0) | (diff(xx[,2]) != 0))
    } else {
      oo <- do.call(order, append(unname(as.list(x)), list(seqn)))
      xx <- x[oo, , drop=FALSE]
      isfirst <- c(TRUE, matrowany(apply(xx, 2, diff) != 0))
    }
    uniqueids <- seqn[oo][isfirst]
    lastunique <- cumsum(isfirst)
    result[oo] <- uniqueids[lastunique]
    return(result)
  }
  #' general case
  if(!anyDuplicated(x))
    return(result)
  dup <- duplicated(x)
  uni <- which(!dup)
  for(j in which(dup)) {
    for(i in uni[uni < j]) {
      if(IdenticalRowPair(i, j, x)) {
        result[j] <- i
        break
      }
    }
  }
  return(result)
}

## utility to check whether two rows are identical

IdenticalRowPair <- function(i,j, a, b=a) {
  #' i and j are row indices (single integers)
  ai <- a[i,]
  bj <- b[j,]
  row.names(ai) <- row.names(bj) <- NULL
  identical(ai, bj)
}

## vectorised

IdenticalRows <- function(i, j, a, b=a) {
  #' i and j are row index vectors of equal length
  #' result[k] = identical( a[i[k],]  , b[j[k],] )
  Mo <- if(missing(b)) list(a=a) else list(a=a, b=b)
  mapply(IdenticalRowPair, i=i, j=j, MoreArgs=Mo,
         SIMPLIFY=TRUE, USE.NAMES=FALSE)
}

