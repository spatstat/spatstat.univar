#'
#'    socs.R
#'
#'  Code for approximating the distribution of a weighted sum of chi^2_1
#'
#'  Based on
#'    A.T.A. Wood (1989)
#'    An F approximation to the distribution of a
#'    linear combination of chi-squared variables
#'    Communications in Statistics - Simulation and Computation
#'    18:4, 1439-1456
#'
#'  Copyright (c) Adrian Baddeley 2026
#'
#' $Revision: 1.12 $ $Date: 2026/04/17 02:29:26 $

dsocs <- function(x, lambda, log=FALSE) {
  with(socsCalc(lambda), {
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

psocs <- function(q, lambda, lower.tail=TRUE, log.p=FALSE) {
  with(socsCalc(lambda), {
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

qsocs <- function(p, lambda, lower.tail=TRUE, log.p=FALSE) {
  with(socsCalc(lambda), {
    switch(case,
           I = {
             ## F approximation: X ~ eta * F(df1, df2)
             q <- eta * qf(p, df1, df2, lower.tail=lower.tail, log.p=log.p)
           },
           II = ,
           IV = {
             ## gamma approximation: X ~ Gamma(alpha, beta)
             q <- qgamma(p, alpha, beta, lower.tail=lower.tail, log.p=log.p)
           },
           III = {
             ## inverse gamma approximation: 1/X ~ Gamma(alpha, beta)
             q <- 1/qgamma(p, alpha, beta, lower.tail= !lower.tail, log.p=log.p)
           }
           )
    return(q)
  })
}

rsocs <- function(n, lambda, approx=FALSE) {
  approx <- isTRUE(approx)
  if(!approx) {
    ## exact
    m <- length(lambda)
    x <- matrix(rnorm(n * m)^2, n, m) %*% lambda
    return(x)
  }
  with(socsCalc(lambda), {
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

socsCalc <- function(lambda, zero=sqrt(.Machine$double.eps), verbose=FALSE) {
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
    ## Case II: Gamma (Satterthwaite-Welch) approximation
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

socsTest <- function(statistic, lambda, data.name="x", ..., method=NULL) {
  v <- socsCalc(lambda)
  q <- unname(statistic)
  with(v, {
    switch(case,
           I = {
             ## F approximation: X ~ eta * F(df1, df2)
             p <- pf(q/eta, df1, df2, lower.tail=FALSE, log.p=log.p)
             m <- "Wood's F approximation"
           },
           II = ,
           IV = {
             ## gamma approximation: X ~ Gamma(alpha, beta)
             p <- pgamma(q, alpha, beta, lower.tail=FALSE, log.p=log.p)
             m <- "Satterthwaite-Welch gamma approximation"
           },
           III = {
             ## inverse gamma approximation: 1/X ~ Gamma(alpha, beta)
             p <- pgamma(1/q, alpha, beta, lower.tail=TRUE, log.p=log.p)
             m <- "Inverse-gamma approximation"
           })
    result <-
      structure(list(
        statistic = statistic,
        p.value   = p,
        method    = c(method, m),
        data.name = data.name,
        ...
      ),
      class="htest")
    return(result)
    })
}
