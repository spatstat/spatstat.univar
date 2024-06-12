#'   detect rounding of data
#'
#'   $Revision: 1.1 $ $Date: 2024/06/12 06:23:36 $
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
      if(k < smallk || !isTRUE(all.equal(x, round(x, k-1))))
        return(k)
      k <- k-1
    }
  } else {
    # not integers: go down
    k <- 1
    bigk <- -log10(.Machine$double.eps)
    repeat {
      if(k > bigk || isTRUE(all.equal(x, round(x, k))))
        return(k)
      k <- k+1
    }
  }
}
