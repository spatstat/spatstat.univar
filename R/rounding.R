#'   detect rounding of data
#'
#'   $Revision: 1.4 $ $Date: 2024/06/14 05:20:33 $
#'

rounding <- function(x) {
  UseMethod("rounding")
}

rounding.default <- function(x) {
  # works for numeric, complex, matrix etc
  if(all(x == 0))
    return(0)
  if(isTRUE(all.equal(x, round(x)))) { 
    # integers: go up
    k <- 0
    smallk <- -log10(.Machine$double.xmax)
    repeat {
      if(k < smallk || !all(x == round(x, k-1)))
        return(k)
      k <- k-1
    }
  } else {
    # not integers: go down
    k <- 1
    bigk <- -log10(.Machine$double.eps)
    repeat {
      if(k > bigk || all(x == round(x, k)))
        return(k)
      k <- k+1
    }
  }
}


## least significant digit in decimal representation
lastdigit <- function(x) {
  x <- abs(as.numeric(x))
  z <- (x * 10^sapply(x, rounding.default)) %% 10
  return(z)
}

## most significant digit in decimal representation
firstdigit <- function(x) {
  x <- abs(as.numeric(x))
  z <- trunc(x/10^(floor(log10(ifelse(x == 0, 1, x)))))
  return(z)
}

## number of digits in decimal representation
ndigits <- function(x) {
  x <- abs(as.numeric(x))
  z <- pmax(1, ceiling(log10(ifelse(x == 0, 1, x * 10^sapply(x, rounding.default)))))
  return(z)
}

