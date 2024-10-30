#'
#' kernelsBC.R
#'
#' Boundary-corrected kernels on the positive half-line
#'
#' Copyright (c) 2008-2023 Adrian Baddeley, Tilman Davies and Martin Hazelton
#' GNU Public Licence (>= 2.0)

dkernelBC <- function(x,  mean, sd=1, kernel="gaussian",
                      zerocor=c("none", "weighted", "convolution",
                              "reflection", "bdrykern")) {
  kernel <- match.kernel(kernel)
  zerocor <- match.arg(zerocor)
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(mean))
  stopifnot(is.numeric(sd) && length(sd) == 1 && sd > 0)
  #' compute uncorrected density
  fx <- dkernel(x, mean=mean, sd=sd, kernel=kernel)
  #' now adjust it
  switch(zerocor,
         none = {},
         weighted = {
           ## divide by mass of original kernel on positive half line
           fx <- ifelse(fx <= 0, 0,
                        fx/(1 - pkernel(0, mean=mean, sd=sd, kernel=kernel)))
         },
         convolution = {
           ## divide by mass of kernel at query point
           fx <- ifelse(fx <= 0, 0,
                        fx/(1 - pkernel(0, x, sd=sd, kernel=kernel)))
         },
         reflection = {
           ## add density of reflected kernel
           fx <- fx + dkernel(x, mean= -mean, sd=sd, kernel=kernel)
         },
         bdrykern = {
           h <- sd * kernel.factor(kernel)
           p <- x/h
           u <- (x-mean)/h
           nu0x <- kernel.moment(0, p, kernel=kernel) # default is template
           nu1x <- kernel.moment(1, p, kernel=kernel)
           nu2x <- kernel.moment(2, p, kernel=kernel)
           denomx <- nu0x * nu2x - nu1x^2
           ax <- nu2x/denomx
           bx <- nu1x/denomx
           fx <- ifelse(fx <= 0, 0,
                        fx * (ax - bx * u))
         })
  return(fx)
}

