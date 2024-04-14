#
#	kmrs.R
#
#	S code for Kaplan-Meier, Reduced Sample and Hanisch
#	estimates of a distribution function
#	from _histograms_ of censored data.
#
#	kaplan.meier()
#	reduced.sample()
#       km.rs()
#
#	$Revision: 3.30 $	$Date: 2023/11/04 04:46:07 $
#
#	The functions in this file produce vectors `km' and `rs'
#	where km[k] and rs[k] are estimates of F(breaks[k+1]),
#	i.e. an estimate of the c.d.f. at the RIGHT endpoint of the interval.
#

"kaplan.meier" <-
function(obs, nco, breaks, upperobs=0) {
#	obs: histogram of all observations : min(T_i,C_i)
#	nco: histogram of noncensored observations : T_i such that T_i <= C_i
# 	breaks: breakpoints (vector or 'breakpts' object, see breaks.S)
#       upperobs: number of observations beyond rightmost breakpoint
#  
        breaks <- as.breakpts(breaks)

	n <- length(obs)
	if(n != length(nco)) 
		stop("lengths of histograms do not match")
	check.hist.lengths(nco, breaks)
#
#	
#   reverse cumulative histogram of observations
	d <- revcumsum(obs) + upperobs
#
#  product integrand
	s <- ifelseXB(d > 0, 1 - nco/d, 1)
#
	km <- 1 - cumprod(s)
#  km has length n;  km[i] is an estimate of F(r) for r=breaks[i+1]
#	
	widths <- diff(breaks$val)
        lambda <- numeric(n)
        pos <- (s > 0)
        lambda[pos] <- -log(s[pos])/widths[pos]
#  lambda has length n; lambda[i] is an estimate of
#  the average of \lambda(r) over the interval (breaks[i],breaks[i+1]).
#	
	return(list(km=km, lambda=lambda))
}

"reduced.sample" <-
function(nco, cen, ncc, show=FALSE, uppercen=0)
#	nco: histogram of noncensored observations: T_i such that T_i <= C_i
#	cen: histogram of all censoring times: C_i
#	ncc: histogram of censoring times for noncensored obs:
#		C_i such that T_i <= C_i
#
#	Then nco[k] = #{i: T_i <= C_i, T_i \in I_k}
#	     cen[k] = #{i: C_i \in I_k}
#	     ncc[k] = #{i: T_i <= C_i, C_i \in I_k}.
#
#       The intervals I_k must span an interval [0,R] beginning at 0.
#       If this interval did not include all censoring times,
#       then `uppercen' must be the number of censoring times
#       that were not counted in 'cen'.
{
	n <- length(nco)
	if(n != length(cen) || n != length(ncc))
		stop("histogram lengths do not match")
#
#	denominator: reverse cumulative histogram of censoring times
#		denom(r) = #{i : C_i >= r}
#	We compute 
#		cc[k] = #{i: C_i > breaks[k]}	
#	except that > becomes >= for k=0.
#
	cc <- revcumsum(cen) + uppercen
#
#
#	numerator
#	#{i: T_i <= r <= C_i }
#	= #{i: T_i <= r, T_i <= C_i} - #{i: C_i < r, T_i <= C_i}
#	We compute
#		u[k] = #{i: T_i <= C_i, T_i <= breaks[k+1]}
#			- #{i: T_i <= C_i, C_i <= breaks[k]}
#		     = #{i: T_i <= C_i, C_i > breaks[k], T_i <= breaks[k+1]}
#	this ensures that numerator and denominator are 
#	comparable, u[k] <= cc[k] always.
#
	u <- cumsum(nco) - c(0,cumsum(ncc)[1:(n-1)])
	rs <- u/cc
#
#	Hence rs[k] = u[k]/cc[k] is an estimator of F(r) 
#	for r = breaks[k+1], i.e. for the right hand end of the interval.
#
        if(!show)
          return(rs)
        else
          return(list(rs=rs, numerator=u, denominator=cc))
}

"km.rs" <-
function(o, cc, d, breaks) {
#	o: censored lifetimes min(T_i,C_i)
#	cc: censoring times C_i
#	d: censoring indicators 1(T_i <= C_i)
#	breaks: histogram breakpoints (vector or 'breakpts' object)
#
  breaks <- as.breakpts(breaks)
  bval <- breaks$val
# compile histograms (breakpoints may not span data)
  obs <- whist( o,     breaks=bval)
  nco <- whist( o[d],  breaks=bval)
  cen <- whist( cc,    breaks=bval)
  ncc <- whist( cc[d], breaks=bval)
# number of observations exceeding largest breakpoint
  upperobs <- attr(obs, "high")
  uppercen <- attr(cen, "high")
# go
  km <- kaplan.meier(obs, nco, breaks, upperobs=upperobs)
  rs <- reduced.sample(nco, cen, ncc, uppercen=uppercen)
#
  return(list(rs=rs, km=km$km, hazard=km$lambda,
              r=breaks$r, breaks=bval))
}

"km.rs.opt" <-
function(o, cc, d, breaks, KM=TRUE, RS=TRUE) {
#	o: censored lifetimes min(T_i,C_i)
#	cc: censoring times C_i
#	d: censoring indicators 1(T_i <= C_i)
#	breaks: histogram breakpoints (vector or 'breakpts' object)
#
  breaks <- as.breakpts(breaks)
  bval <- breaks$val
  out <- list(r=breaks$r, breaks=bval)
  if(KM || RS)
    nco <- whist( o[d],  breaks=bval)
  if(KM) {
    obs <- whist( o,     breaks=bval)
    upperobs <- attr(obs, "high")
    km <- kaplan.meier(obs, nco, breaks, upperobs=upperobs)
    out <- append(list(km=km$km, hazard=km$lambda), out)
  }
  if(RS) {
    cen <- whist( cc,    breaks=bval)
    ncc <- whist( cc[d], breaks=bval)
    uppercen <- attr(cen, "high")
    rs <- reduced.sample(nco, cen, ncc, uppercen=uppercen)
    out <- append(list(rs=rs), out)
  }
  return(out)
}


