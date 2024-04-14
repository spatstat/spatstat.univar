#'
#'   Header for all (concatenated) test files
#'
#'   Require spatstat.univar
#'   Obtain environment variable controlling tests.
#'
#'   $Revision: 1.5 $ $Date: 2020/04/30 05:31:37 $

require(spatstat.univar)
FULLTEST <- (nchar(Sys.getenv("SPATSTAT_TEST", unset="")) > 0)
ALWAYS   <- TRUE
cat(paste("--------- Executing",
          if(FULLTEST) "** ALL **" else "**RESTRICTED** subset of",
          "test code -----------\n"))
# tests/weightedstats.R
# $Revision: 1.2 $ $Date: 2023/11/05 01:40:53 $

local({
  if(ALWAYS) { # depends on hardware
    ## whist()
    ## check agreement between C and interpreted code for whist()
    set.seed(98123)
    x <- runif(1000)
    w <- sample(1:5, 1000, replace=TRUE)
    b <- seq(0,1,length=101)
    aC <- whist(x,b,w, method="C")
    aR <- whist(x,b,w, method="interpreted")
    if(!all(aC == aR))
      stop("Algorithms for whist disagree")
  }
  if(FULLTEST) {
    ## cases of 'unnormdensity()'
    x <- rnorm(20) 
    d0 <- unnormdensity(x, weights=rep(0, 20))
    dneg <- unnormdensity(x, weights=c(-runif(19), 0))
  }

})

#'
#'   tests/parzen.R
#'
#'   Tests of the Parzen-Rosenblatt estimator
#'       (fixed bandwidth, no boundary correction)
#'
#'  $Revision: 1.1 $ $Date: 2023/10/22 02:39:49 $

local({
  if(FULLTEST) {
    #' code in kernels.R
    kernames <- c("gaussian", "rectangular", "triangular",
                  "epanechnikov", "biweight", "cosine", "optcosine")
    X <- rnorm(20)
    U <- runif(20)
    for(ker in kernames) {
      dX <- dkernel(X, ker)
      fX <- pkernel(X, ker)
      qU <- qkernel(U, ker)
      m0 <- kernel.moment(0, 0, ker)
      m1 <- kernel.moment(1, 0, ker)
      m2 <- kernel.moment(2, 0, ker)
      m3 <- kernel.moment(3, 0, ker)
      fa <- kernel.factor(ker)
      sq <- kernel.squint(ker)
    }
  }
})

local({
  if(ALWAYS) {
    ## unnormdensity
    x <- rnorm(20) 
    d0 <- unnormdensity(x, weights=rep(0, 20))
    dneg <- unnormdensity(x, weights=c(-runif(19), 0))
  }
})
#
# tests/NAinCov.R
#
# Testing the response to the presence of NA's in covariates
#
# $Revision: 1.9 $ $Date: 2023/11/05 01:45:36 $

if(FULLTEST) {
local({
  #' quantile.ewcdf
  f <- ewcdf(runif(100), runif(100))
  qf <- quantile(f, probs=c(0.1, NA, 0.8))
  #' quantile.density
  f <- density(runif(100))
  qf <- quantile(f, probs=c(0.1, NA, 0.8))
})
}
