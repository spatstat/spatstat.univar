/*
  farebro.c

  Farebrother's algorithm for the distribution of a linear combination
  of chi-squared random variables

  R.W. Farebrother (1984)
  Algorithm AS 204:
  The distribution of a positive linear combination of \chi^2 random variables.
  Applied Statistics 33 #3 (1984) 332--339

  Translated to C by Adrian Baddeley 2026
  Modified to accept a vector of function arguments 'x'

  Copyright (c) Adrian Baddeley 2026
  GNU Public Licence (>= 2.0)

  $Revision: 1.7 $ $Date: 2026/04/21 02:42:25 $

*/

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <math.h>

#undef DEBUG

double centnorm(double x) {
  /* Prob( - x < Z < x ) */
  double y;
  y = 1.0 - 2.0 * pnorm(x,
			(double) 0.0, (double) 1.0,
			(int) 0, (int) 0);
  return(y);
}

void farebro(
	     /* input data */
	     double *lambda,     /* vector of coefficients */
	     int *mult,          /* vector of multiplicities */
	     double *delta,      /* vector of noncentrality parameters */
	     int *n,             /* length of the above vectors */
	     double *x,          /* (vector) argument of CDF */
	     int *nx,            /* length of vector x */
	     /* algorithm control */
	     double *mode,       /* algorithm parameter */
	     int *maxit,         /* maximum number of terms in sum */
	     double *eps,        /* desired level of accuracy */
	     double *tol,        /* underflow threshold */
	     /* outputs (vectors of length nx) */
	     int *ifault,        /* error code */
	     double *density,    /* probability density f(c) */
	     double *probability /* cumulative probability F(c) */
	     ) {
  int N, Nx, Maxit;
  double c, Mode, Eps, Tol;

  int i, j, k, m, converged, faultcode;
  int i1, m1, k1, kfixed; /* 1-based */
  double a0, a0inv, z, beta, eps2, hold, hold2, sum, sum1;
  double dans, lans, pans, prbty, dnsty;
  double *gamma, *theta, *a, *b;

#ifdef DEBUG
  int jj;
#endif  
  
  N     = *n;
  Nx    = *nx;
  Maxit = *maxit;
  Mode  = *mode;
  Eps   = *eps;
  Tol   = *tol;

  /* initialise fault codes to 'OK' */
  for(j = 0; j < Nx; j++) ifault[j] = 0;
  
  /* validate arguments */
  if(N < 1 || Maxit < 1 || Eps <= 0.0) {
    for(j = 0; j < Nx; j++) ifault[j] = 2;
    return;
  }
  for(i = 0; i < N; i++) {
    if(lambda[i] <= 0.0 || mult[i] < 1 || delta[i] < 0.0) {
      for(j = 0; j < Nx; j++)
	ifault[j] = -i-1; /* 1-based index */
      return;
    }
  }
  /* detect non-positive arguments; handle later */
  for(j = 0; j < Nx; j++) {
    if(x[j] <= 0.0) ifault[j] = 2;
  }

  /* allocate temporary storage */
  gamma = (double *) R_alloc(N, sizeof(double));
  theta = (double *) R_alloc(N, sizeof(double));
  a = (double *) R_alloc(Maxit, sizeof(double));
  b = (double *) R_alloc(Maxit, sizeof(double));

  /* preliminaries */
  beta = sum = lambda[0];
  for(i = 0; i < N; i++) {
    hold = lambda[i];
    if(beta > hold) beta = hold;
    if(sum < hold) sum = hold;
  }
  beta = (Mode > 0.0) ? (Mode * beta) : (2.0/(1.0/beta + 1.0/sum));
#ifdef DEBUG
  Rprintf("beta=%lf\n", beta);
#endif  
  k = 0;
  sum = 1.0;
  sum1 = 0.0;
  for(i = 0; i < N; i++) {
    hold = beta/lambda[i];
    gamma[i] = 1.0 - hold;
    sum *= pow(hold, (double) mult[i]);
    sum1 += delta[i];
    k += mult[i];
    theta[i] = 1.0;
  }
#ifdef DEBUG
  for(jj = 0; jj < N; jj++)
    Rprintf("gamma[%d]=%lf\n", jj, gamma[jj]);
#endif  
  kfixed = k;
  a0 = exp(0.5 * (log(sum) - sum1));
  if(a0 <= 0) {
    for(j = 0; j < Nx; j++) {
      probability[j] = 0.0;
      density[j]     = 0.0;
      ifault[j]      = 1;
    }
    return;
  }
#ifdef DEBUG
  Rprintf("a0=%lf, k=%d\n", a0, k);
#endif  

  /* start calculating */

  for(j = 0; j < Nx; j++) {
    if(ifault[j] == 0) {
      converged = 0;
      c = x[j];
      z = c/beta;
      /* reset theta[] */
      for(i = 0; i < N; i++)
	theta[i] = 1.0;
#ifdef DEBUG
  Rprintf("z=%lf\n", z);
#endif  
      /* 
	 Evaluate CDF and PDF of chi-squared with k d.f.
      */
      if(kfixed % 2 == 0) {
	i1 = 2; 
	lans = -z/2.0;
	dans = exp(lans);
	pans = 1.0 - dans;
      } else {
	i1 = 1;
	lans = -(z + log(z))/2.0 - M_LN_SQRT_PId2;
	dans = exp(lans);
	pans = centnorm(sqrt(z));
      }
#ifdef DEBUG
      Rprintf("pans=%lf, dans=%lf\n", pans, dans);
#endif  
      k1 = kfixed - 2;
      for(; i1 <= k1; i1 += 2) {
	if(lans >= Tol) {
	  dans *= z/i1;
	} else {
	  /* avoid underflow */
	  lans += log(z/i1);
	  dans = exp(lans);
	}
	pans -= dans;
      }
      /* Evaluate successive terms of expansion */
      prbty = pans;
      dnsty = dans;
      eps2 = Eps/a0;
      a0inv = 1.0/a0;
      sum = a0inv - 1.0;
#ifdef DEBUG
      Rprintf("prbty=%lf, dnsty=%lf, sum=%lf\n", prbty, dnsty, sum);
#endif  
      for(m = 0, m1 = 1; m < Maxit; m++, m1++) {
	sum1 = 0.0;
	for(i = 0; i < N; i++) {
	  hold = theta[i];
	  theta[i] = hold2 = hold * gamma[i];
	  sum1 += hold2 * mult[i] + m1 * delta[i] * (hold - hold2);
	}
#ifdef DEBUG
	Rprintf("\tTerm %d:\n", m);
	for(jj = 0; jj < N; jj++)
	  Rprintf("\t\ttheta[%d]=%lf\n", jj, theta[jj]);
#endif  
	sum1 /= 2.0;
	b[m] = sum1;
	if(m1 >= 2) {
	  for(i=m-1; i >= 0; --i) {
	    sum1 += b[i] * a[m-i-1];
	  }
	}
	sum1 /= (double) m1;
	a[m] = sum1;
#ifdef DEBUG
	for(jj = 0; jj < m; jj++) {
	  Rprintf("\t\ta[%d] = %lf, b[%d]=%lf\n",
		  jj, a[jj], jj, b[jj]);
	}
#endif
	k1 += 2;
	if(lans >= Tol) {
	  dans *= z/k1;
	} else {
	  /* underflow */
	  lans += log(z/k1);
	  dans = exp(lans);
	}
#ifdef DEBUG
	Rprintf("\tTerm %d: dans=%lf, sum1=%lf\n", m1, dans, sum1);
#endif	
	pans -= dans;
	sum -= sum1;
	dnsty += dans * sum1;
	sum1 *= pans;
	prbty += sum1;
#ifdef DEBUG
	Rprintf("\tEnd of term %d: prbty=%lf, dnsty=%lf, sum=%lf, sum1=%lf\n",
		m1, prbty, dnsty, sum, sum1);
#endif  
	if(prbty < -a0inv) {
	  ifault[j] = 3;
	  break; 
	}
	if(fabs(pans * sum) < eps2 &&
	   fabs(sum1) < eps2) {
	  converged = 1;
	  break;
	}
	/* end of loop over terms */
      }
#ifdef DEBUG
	Rprintf("\tEnd of iteration: prbty=%lf, dnsty=%lf\n",
		prbty, dnsty);
#endif	
      if(ifault[j] == 0) {
	dnsty = a0 * dnsty/(2 * beta);
	prbty *= a0;
	probability[j] = prbty;
	density[j] = dnsty; 
	/* error code? */
	faultcode = (converged == 0) ?  4 : 0;
	if(prbty < 0.0 || prbty > 1.0) {
	  faultcode += 5;
	} else if(dnsty < 0.0) {
	  faultcode += 6;
	}
	ifault[j] = faultcode;
      }
    }
  }
  return;
}
	     
 
