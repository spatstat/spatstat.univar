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
# $Revision: 1.10 $ $Date: 2024/09/30 23:13:54 $

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

#' tests/direct.R
#'
#' Check output of densityBC() by comparing with
#' kernel estimates computed directly in R code.
#'
#' $Revision: 1.2 $ $Date: 2023/10/22 02:39:35 $

INTERACTIVE <- FALSE

#' ...... direct implementation ........................
#'
#' Biweight boundary kernel (for estimation at point r)

bdry.bwt <- function(x,r,h=1){
		u <- x/h
		k <- (15/(16*h))*(1-u^2)^2*(u^2<1)
		p <- r/h
		p[p>1] <- 1
		a0 <- (3*p^5 - 10*p^3 + 15*p + 8)/16
		a1 <- (5*p^6 - 15*p^4 + 15*p^2 -5)/32
		a2 <- (15*p^7 - 42*p^5 + 35*p^3 + 8)/112
		bk <- (a2-a1*u)*k*(u<p)/(a0*a2-a1^2)
		bk
}

#' Standard biweight kernel

bwt <- function(x,r,h=1){
  # r is ignored
		u <- x/h
		k <- (15/(16*h))*(1-u^2)^2*(u^2<1)
		k
}

#' Kernel estimate 

kernR <- function(x, h, kname="bwt", from, to, nr=200) {
  ker <- get(kname, mode="function")
  fhat <- numeric(0)
  rvalues <- seq(from, to, length=nr)
  for (r in rvalues) 
    fhat <- c(fhat,mean(ker(r-x, r, h=h)))
  return(list(x=rvalues, y=fhat))
}

#' ...................... RUN EXAMPLE .......................

sim.dat <- rexp(500)

from <- 0.01
to <- 5
nr <- 200

if(INTERACTIVE) opa <- par(ask=TRUE)

cat("--- Fixed bandwidth ---\n")

fhatR <- kernR(sim.dat, 0.4, "bwt", from=from, to=to, nr=nr)
fhatC <- densityBC(sim.dat, "biweight", h=0.4, from=from, to=to, n=nr)

stopifnot(length(fhatR$x) == length(fhatC$x))
stopifnot(all(fhatR$x == fhatC$x))
rvalues <- fhatR$x
ftrue <- dexp(rvalues)

if(INTERACTIVE) {
  plot(c(from, to), c(0, 1.1), type="n", xlab="r", ylab="density",
       main="Fixed bandwidth")
  lines(rvalues, fhatR$y, col=1)
  lines(rvalues, fhatC$x, col=2)
  lines(rvalues, ftrue, col=3)
}

cat("Range of absolute discrepancies between estimates (fixed bandwidth):\n")
print(range(fhatC$y - fhatR$y))

cat("\n--- Boundary kernel ---\n")

fhat.R <- kernR(sim.dat, 0.4, "bdry.bwt", from=from, to=to, nr=nr)
fhat.C <- densityBC(sim.dat, "biweight", h=0.4, from=from, to=to, n=nr,
                  zerocor="bdrykern")

stopifnot(length(fhat.R$x) == length(fhat.C$x))
stopifnot(all(fhat.R$x == fhat.C$x))

rvalues <- fhat.R$x
ftrue <- dexp(rvalues)

if(INTERACTIVE) {
  plot(c(from, to), c(0, 1.1), type="n", xlab="r", ylab="density",
       main="Boundary kernel")
  lines(rvalues, fhat.R$y, col=1)
  lines(rvalues, fhat.C$y, col=2)
  lines(rvalues, ftrue, col=3)
}

cat("Range of absolute discrepancies between estimates (boundary kernel):\n")
print(range(fhat.C$y - fhat.R$y))

rel.discrep <- (fhat.C$y - fhat.R$y)/ftrue
rel.discrep <- rel.discrep[ftrue > 0.1]
cat("Range of relative discrepancies between estimates (boundary kernel):\n")
print(range(rel.discrep))
if(max(abs(rel.discrep)) > 0.01)
  stop("Relative discrepancies between C and R code exceed 1 percent")

if(INTERACTIVE) par(opa)

#' tests/kermom.R
#'
#' Check R function kernel.moment() against C function kermom ()
#'
#' $Revision: 1.2 $ $Date: 2024/10/31 10:36:27 $
#'

moo <- 1
sdee <- 0.5
xx <- seq(moo - 4 * sdee, moo + 4 * sdee, length=512)

kernames <- c("gaussian", "rectangular", "triangular",
              "epanechnikov", "biweight", "cosine", "optcosine")

eps <- sqrt(.Machine$double.eps)

for(m in 0:2) {
  cat("Incomplete moment of order", m, fill=TRUE)
  for(ker in kernames) {
    Rvalues <- kernel.moment(m, xx, ker, mean=moo, sd=sdee)
    Cvalues <- kermom(m, xx, ker, mean=moo, sd=sdee)
    discrep <- max(abs(Rvalues-Cvalues))
    if(discrep > eps) 
      stop(paste("kernel.moment and kermom disagree",
                 "for m =", m, "for kernel", sQuote(ker),
                 "\n\tDiscrepancy", discrep))
    cat("Kernel", sQuote(ker), "\tdiscrepancy", discrep, fill=TRUE)
  }
}
