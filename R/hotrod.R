#'
#'   hotrod.R
#'
#'  Heat kernel for a one-dimensional rod
#'
#'  Copyright (c) Greg McSwiggan and Adrian Baddeley 2017-2020
#'
#'  $Revision: 1.6 $ $Date: 2022/05/22 00:01:01 $

hotrod <- function(len, xsource, xquery, sigma, 
                   ends=c("insulated", "absorbing"),
                   nmax=20) {
  ends <- match.arg(ends)
  len <- as.numeric(len)
  xsource <- as.numeric(xsource)
  xquery <- as.numeric(xquery)
  sigma <- as.numeric(sigma)
  nmax <- as.integer(nmax)
  df <- data.frame(len=len, xsource=xsource, xquery=xquery,
                   sigma=sigma, nmax=nmax)
  n <- nrow(df)
  switch(ends,
         insulated = {
           z <- with(df,
                     .C(SK_hotrodInsul,
                        n = as.integer(n),
                        a = as.double(len),
                        x = as.double(xsource),
                        y = as.double(xquery),
                        s = as.double(sigma),
                        m = as.integer(nmax),
                        z = as.double(numeric(n)),
                        PACKAGE="spatstat.univar")$z)
         },
         absorbing = {
           z <- with(df,
                     .C(SK_hotrodAbsorb,
                        n = as.integer(n),
                        a = as.double(len),
                        x = as.double(xsource),
                        y = as.double(xquery),
                        s = as.double(sigma),
                        m = as.integer(nmax),
                        z = as.double(numeric(n)),
                        PACKAGE="spatstat.univar")$z)
         })
  return(z)
}

