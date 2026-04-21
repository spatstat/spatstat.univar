#'
#'    socs.R
#'
#'  Code for the distribution of a weighted sum of chi^2_1
#'
#'  Based on Wood's approximation (F, gamma or inverse-gamma)
#'  and Farebrother's algorithm (truncated infinite series)
#' 
#'   R.W. Farebrother (1984)
#'   Algorithm AS 204:
#'   The distribution of a positive linear combination of \chi^2
#'   random variables. Applied Statistics 33 #3 (1984) 332--339
#' 
#'   A.T.A. Wood (1989)
#'   An F approximation to the distribution of a
#'   linear combination of chi-squared variables
#'   Communications in Statistics - Simulation and Computation
#'   18:4, 1439-1456
#'
#' 
#'  Copyright (c) Adrian Baddeley 2026
#'  GNU Public Licence (>= 2.0)
#'
#' $Revision: 1.21 $ $Date: 2026/04/21 03:44:59 $

dsocs <- function(x, lambda, log=FALSE,
                  method=c("Wood", "Farebrother")) {
  method <- match.arg(method)
  if(method == "Farebrother") {
    a <- farebro(lambda, x=x, warn=FALSE)
    d <- a$d
    if(log) d <- log(ifelse(is.finite(d) & d > 0, d, 1))
    if(any(bad <- (a$ifault == 0))) {
      ## fall back on Wood's approximation
      d[bad] <- dsocs(x[bad], lambda, log=log, method="Wood")
    }
    return(d)
  }
  with(woodCalc(lambda), {
    switch(case,
           I = {
             ## F approximation: X ~ eta * F(df1, df2)
             y <- df(x/eta, df1, df2, log=log)
             y <- if(log) y - log(eta) else y/eta
           },
           II = ,
           IV = {
             ## gamma approximation: X ~ Gamma(alpha, beta)
             y <- dgamma(x, alpha, beta, log=log)
           },
           III = {
             ## inverse gamma approximation: 1/X ~ Gamma(alpha, beta)
             logy <- numeric(length(x))
             ok <- (x > 0)
             logy[ok] <- dgamma(1/x[ok], alpha, beta, log=TRUE) - 2 * log(x[ok])
             logy[!ok] <- -Inf
             y <- if(log) logy else exp(logy)
           })
    return(y)
  })
}

psocs <- function(q, lambda, lower.tail=TRUE, log.p=FALSE,
                  method=c("Wood", "Farebrother")) {
  method <- match.arg(method)
  if(method == "Farebrother") {
    a <- farebro(lambda, x=q, warn=FALSE)
    p <- a$p
    if(!lower.tail) p <- 1 - p
    if(log.p) p <- log(ifelse(is.finite(p) & p > 0, p, 1))
    if(any(bad <- (a$ifault != 0))) {
      ## fall back on Wood's approximation
      p[bad] <- psocs(q[bad], lambda, lower.tail=lower.tail, log.p=log.p,
                      method="Wood")
    }
    return(p)
  }
  with(woodCalc(lambda), {
    switch(case,
           I = {
             ## F approximation: X ~ eta * F(df1, df2)
             p <- pf(q/eta, df1, df2, lower.tail=lower.tail, log.p=log.p)
           },
           II = ,
           IV = {
             ## gamma approximation: X ~ Gamma(alpha, beta)
             p <- pgamma(q, alpha, beta, lower.tail=lower.tail, log.p=log.p)
           },
           III = {
             ## inverse gamma approximation: 1/X ~ Gamma(alpha, beta)
             p <- pgamma(1/q, alpha, beta, lower.tail=!lower.tail, log.p=log.p)
           })
    return(p)
  })
}

qsocs <- function(p, lambda, lower.tail=TRUE, log.p=FALSE,
                  method=c("Wood", "Farebrother")) {
  method <- match.arg(method)
  ## First calculate using Wood's approximation
  qwood <- with(woodCalc(lambda), {
    switch(case,
           I = {
             ## F approximation: X ~ eta * F(df1, df2)
             q <- eta * qf(p, df1, df2,
                           lower.tail=lower.tail, log.p=log.p)
           },
           II = ,
           IV = {
             ## gamma approximation: X ~ Gamma(alpha, beta)
             q <- qgamma(p, alpha, beta,
                         lower.tail=lower.tail, log.p=log.p)
           },
           III = {
             ## inverse gamma approximation: 1/X ~ Gamma(alpha, beta)
             q <- 1/qgamma(p, alpha, beta,
                           lower.tail= !lower.tail, log.p=log.p)
           }
           )
           q
  })
  switch(method,
         Wood = {
           q <- qwood
         },
         Farebrother = {
           q <- numeric(length(p))
           ## handle p = 0 or 1
           pzero <- if(lower.tail) 0 else 1
           pInf  <- if(lower.tail) 1 else 0
           if(log.p) {
             pzero <- log(pzero)
             pInf <- log(pInf)
           }
           if(any(trivial <- (p == pzero | p == pInf))) {
             q[ p == pzero ] <- 0
             q[ p == pInf  ] <- Inf
           }
           ## use Wood approximation to determine search interval
           qlower <- 0.5 * qwood
           qupper <- 1.5 * qwood
           ## search
           for(i in which(!trivial)) {
             q[i] <- uniroot(function(x, ..., prob) { prob - psocs(x, ...) },
                             interval = c(qlower[i], qupper[i]),
                             prob = p[i],
                             lambda = lambda,
                             lower.tail = lower.tail,
                             log.p = log.p,
                             tol = 1e-6,
                             method="Farebrother")$root
           }
         })
  return(q)
}

rsocs <- function(n, lambda, approx=FALSE) {
  approx <- isTRUE(approx)
  if(!approx) {
    ## exact
    m <- length(lambda)
    x <- matrix(rnorm(n * m)^2, n, m) %*% lambda
    return(x)
  }
  with(woodCalc(lambda), {
    switch(case,
           I = {
             ## F approximation: eta * F(df1, df2)
             x <- eta * rf(n, df1, df2)
           },
           II = ,
           IV = {
             ## gamma approximation: X ~ Gamma(alpha, beta)
             x <- rgamma(n, alpha, beta)
           },
           III = {
             ## inverse gamma approximation: 1/X ~ Gamma(alpha, beta)
             x <- 1/rgamma(n, alpha, beta)
           })
    return(x)
  })
}

woodCalc <- function(lambda, zero=sqrt(.Machine$double.eps), verbose=FALSE) {
  #' cumulants 2^(m-1) (m-1)! sum(lambda^m)
  kappa1 <- sum(lambda)
  kappa2 <- 2 * sum(lambda^2)
  kappa3 <- 8 * sum(lambda^3)
  #' Wood (1989)
  tau1 <- 4 * kappa2^2 * kappa1 + kappa3 * (kappa2 - kappa1^2)
  tau2 <- kappa3 * kappa1 - 2 * kappa2^2
  if(tau1 > zero && tau2 > zero) {
    ## Case I: F approximation
    alpha1 <- 2 * kappa1 * (kappa1 * kappa3 + kappa1^2 * kappa2 - kappa2^2)/tau1
    beta <- tau1/tau2
    alpha2 <- 3 + 2 * kappa2 * (kappa2 + kappa1^2)/tau2
    df1 <- 2 * alpha1
    df2 <- 2 * alpha2
    eta <- beta * df1/df2
    result <- list(case="I", df1=df1, beta=beta, df2=df2, eta=eta)
  } else if(abs(tau2) <= zero) {
    ## Case II: Gamma (Satterthwaite-Welsh) approximation
    alphaSW <- kappa1^2/kappa2
    betaSW <- kappa1/kappa2
    result <- list(case="II", alpha=alphaSW, beta=betaSW)
  } else if(tau1 <= zero) {
    ## Case III: Inverse gamma approximation
    alphastar <- 2 + kappa1^2/kappa2^2
    betastar <- kappa1^3/kappa2 + kappa1
    result <- list(case="III", alpha=alphastar, beta=betastar)
  } else {
    ## Fallback: gamma approx
    alphaSW <- kappa1^2/kappa2
    betaSW <- kappa1/kappa2
    result <- list(case="IV", alpha=alphaSW, beta=betaSW)
  }
  if(verbose) 
    splat("Case", result$case, ":", 
          paste(names(result)[-1], "=", signif(as.numeric(result[-1]), 4),
                collapse=", "))
  return(result)
}

#' Create an 'htest' object which represents the result of a test
#' with test statistic 'statistic' referred to the distribution of
#' a weighted sum of chi^2 variables with coefficients 'lambda'.
#' 
#' Any additional arguments '...' will be added to the object (name=value)
#' 
#' Arguments handled by print.htest include
#'        parameter
#'        alternative         'two.sided', 'less', 'greater'
#'        null.value          numeric value or vector
#'        conf.int            num[2] with attribute 'conf.level'
#'        estimate            numeric value or vector
#'
#' Common additional arguments include 'observed', 'expected', 'residuals'

socsTest <- function(statistic, lambda, data.name="x", ...,
                     approxmethod=c("Wood", "Farebrother"),
                     testmethod=NULL) {
  approxmethod <- match.arg(approxmethod)
  stopifnot(length(statistic) == 1)
  p <- psocs(unname(statistic), lambda, lower.tail=FALSE, log.p=FALSE,
             method=approxmethod)
  switch(approxmethod,
         Farebrother = {
           calcname <- "Farebrother algorithm AS204"
         },
         Wood = {
           v <- woodCalc(lambda)
           calcname <- switch(v$case,
                              I = "Wood's F approximation",
                              II = ,
                              IV = "Satterthwaite-Welsh gamma approximation",
                              III = "Inverse-gamma approximation",
                              NULL)
         })
  result <-
    structure(list(
      statistic = statistic,
      p.value   = p,
      method    = c(testmethod, calcname),
      data.name = data.name,
      ...
    ),
    class="htest")
  return(result)
}
