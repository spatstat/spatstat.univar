/* 

   colonel.c

   Kernel density estimation, fixed bandwidth,
   with weights and (optional) linear boundary correction

   Brute force algorithm:

   colonel()    kernel estimate, uncorrected
   bcolonel()   kernel estimate, linear boundary correction at zero

   Faster algorithm for equally-spaced 'r' values:

   fcolonel()    kernel estimate, uncorrected
   fbcolonel()   kernel estimate, linear boundary correction at zero

   Data assumed to be rescaled so that kernel has halfwidth = 1 

   This code does not call the self-contained kernel functions 
   defined in 'kernels.c'. Instead the kernels are coded in-line.
   This implementation is slightly faster for large datasets, 
   as it avoids repeated function calls and it makes some efficiencies. 
   However, it requires pre-scaling of the data.

   Copyright (c) 2008-2024 Adrian Baddeley, Tilman Davies and Martin Hazelton
   GNU Public Licence (>= 2.0)

*/

#include <Rmath.h>
#include "kerconstants.h"
#include "interfacecodes.h"

#define PI         M_PI
#define TWO_PI     M_2PI
#define HALF_PI    M_PI_2
#define QUARTER_PI M_PI_4
#define RECIPROCAL_SQRT_TWO_PI M_1_SQRT_2PI  

double sqrt(double), exp(double), cos(double);

#define ABS(X) (((X) < 0) ? (-(X)) : (X))
#define INSIDE(X)  (((X) >= -1) && ((X) <= 1))

#define GAUSSTHRESH ((double) 8.0)

void colonel(
	     int *kerncode,   /* integer code for kernel		    */
	     int *nx,	      /* number of data values			    */
	     double *x,	      /* vector of data values			    */
	     double *w,	      /* vector of weights for data		    */
	     int *nr,	      /* number of r values			    */
	     double *r,	      /* vector of r values (argument of density)   */
	     double *f,	      /* vector of values of density estimate f(r)  */
	     int *errcode)    /* integer code for errors                    */
{ 
  int i, j, Nx, Nr, bb;
  double xi, wi, uij, temp, kvalue, root2pi;

  /* extract arguments and validate */
  Nx = *nx;
  Nr = *nr;

  *errcode = OK;
  if(Nx < 0 || Nr <= 0) {
    *errcode = ERR_NEGATIVE_LENGTH;
    return;
  }

  /* initialise f(r) */
  for(j = 0; j < Nr; j++) 
    f[j] = 0.0;

  if(Nx == 0)
    return;

  /* go */
  if(*kerncode == GAUSSIAN) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      for(j = 0; j < Nr; j++) {
	uij = xi - r[j];
	kvalue = exp(- uij * uij/2);
	f[j] += kvalue * wi;
      }
    }
    /* divide by constant factor */
    for(j = 0; j < Nr; j++) 
      f[j] *= RECIPROCAL_SQRT_TWO_PI;
    return;
  } else if(*kerncode == RECTANGULAR) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      for(j = 0; j < Nr; j++) {
	uij = xi - r[j];
	if(ABS(uij) < 1.0)
	  f[j] += wi;
      }
    }
    /* divide by constant factor */
    for(j = 0; j < Nr; j++) 
      f[j] /= 2.0;
    return;
  } else if(*kerncode == TRIANGULAR) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      for(j = 0; j < Nr; j++) {
	uij = xi - r[j];
	kvalue = 1.0 - ABS(uij);
	if(kvalue > 0)
	  f[j] += kvalue * wi;
      }
    }
    /* no constant factor */
    return;
  } else if(*kerncode == EPANECHNIKOV) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      for(j = 0; j < Nr; j++) {
	uij = xi - r[j];
	kvalue = (1.0 - uij * uij);
	if(kvalue > 0) 
	  f[j] += kvalue * wi;
      }
    }
    /* multiply by constant factor */
    for(j = 0; j < Nr; j++) 
      f[j] *= 3.0/4.0;
    return;
  } else if(*kerncode == BIWEIGHT) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      for(j = 0; j < Nr; j++) {
	uij = xi - r[j];
	temp = (1 - uij * uij);
	if(temp > 0) {
	  kvalue = temp * temp;
	  f[j] += kvalue * wi;
	}
      }
    }
    /* multiply by constant factor */
    for(j = 0; j < Nr; j++) 
      f[j] *= 15.0/16.0;
    return;
  } else if(*kerncode == COSINE) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      for(j = 0; j < Nr; j++) {
	uij = xi - r[j];
	uij = ABS(uij);
	if(uij < 1.0) {
	  kvalue = 1.0 + cos(PI * uij);
	  f[j] += kvalue * wi;
	}
      }
    }
    /* divide by constant factor */
    for(j = 0; j < Nr; j++) 
      f[j] /= 2.0;
    return;
  } else if(*kerncode == OPTCOSINE) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      for(j = 0; j < Nr; j++) {
	uij = xi - r[j];
	uij = ABS(uij);
	if(uij < 1.0) {
	  kvalue = cos(uij * HALF_PI);
	  f[j] += kvalue * wi;
	}
      }
    }
    /* multiply by constant factor */
    for(j = 0; j < Nr; j++) 
      f[j] *= QUARTER_PI;
    return;
  } else {
    /* unrecognised kernel */
    *errcode = ERR_UNKNOWN_KERNEL;
    return;
  }
}

/* 
   ......................................................................
   fcolonel()   Faster version of colonel() for equally-spaced 'r' values
   ......................................................................
*/

void fcolonel(
	     int *kerncode,   /* integer code for kernel		    */
	     int *nx,	      /* number of data values			    */
	     double *x,	      /* vector of data values			    */
	     double *w,	      /* vector of weights for data		    */
	     int *nr,	      /* number of r values			    */
	     double *r,	      /* vector of r values (argument of density)   */
	     double *f,	      /* vector of values of density estimate f(r)  */
	     int *errcode)    /* integer code for errors                    */
{ 
  int i, j, k, Nx, Nr, bb;
  double dr, xi, wi, vij, temp, kvalue, root2pi;

  /* extract arguments and validate */
  Nx = *nx;
  Nr = *nr;

  *errcode = OK;
  if(Nx < 0 || Nr <= 0) {
    *errcode = ERR_NEGATIVE_LENGTH;
    return;
  }

  /* spacing between 'r' values */
  dr = r[Nr-1]/((double) Nr);
  
  /* initialise f(r) */
  for(j = 0; j < Nr; j++) 
    f[j] = 0.0;

  if(Nx == 0)
    return;

  /* go */
  if(*kerncode == GAUSSIAN) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      k = floor((xi - GAUSSTHRESH)/dr); /* r[k] is 8 std. dev. below x[i] */
      if(k < 0) k = 0;
      if(k < Nr) {
	for(j = k; j < Nr; j++) {
	  vij = r[j] - xi;
	  if(vij > GAUSSTHRESH) break;  /* r[j] is 8 std. dev. above x[i] */
	  kvalue = exp(- vij * vij/2);
	  f[j] += kvalue * wi;
	}
      }
    }
    /* divide by constant factor */
    for(j = 0; j < Nr; j++) 
      f[j] *= RECIPROCAL_SQRT_TWO_PI;
    return;
  } else if(*kerncode == RECTANGULAR) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      k = floor((xi - (double) 1.0)/dr); /* r[k] is 1 halfwidth below x[i] */
      if(k < 0) k = 0;
      if(k < Nr) {
	for(j = k; j < Nr; j++) {
	  vij = r[j] - xi;
	  if(vij >  1.0) break;  /* r[j] is 1 halfwidth above x[i] */
	  if(vij > -1.0)
	    f[j] += wi;
	}
      }
    }
    /* divide by constant factor */
    for(j = 0; j < Nr; j++) 
      f[j] /= 2.0;
    return;
  } else if(*kerncode == TRIANGULAR) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      k = floor((xi - (double) 1.0)/dr); /* r[k] is 1 halfwidth below x[i] */
      if(k < 0) k = 0;
      if(k < Nr) {
	for(j = k; j < Nr; j++) {
	  vij = r[j] - xi;
	  if(vij >  1.0) break;  /* r[j] is 1 halfwidth above x[i] */
	  kvalue = 1.0 - ABS(vij);
	  if(kvalue > 0)
	    f[j] += kvalue * wi;
	}
      }
    }
    /* no constant factor */
    return;
  } else if(*kerncode == EPANECHNIKOV) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      k = floor((xi - (double) 1.0)/dr); /* r[k] is 1 halfwidth below x[i] */
      if(k < 0) k = 0;
      if(k < Nr) {
	for(j = k; j < Nr; j++) {
	  vij = r[j] - xi;
	  if(vij >  1.0) break;  /* r[j] is 1 halfwidth above x[i] */
	  kvalue = (1.0 - vij * vij);
	  if(kvalue > 0) 
	    f[j] += kvalue * wi;
	}
      }
    }
    /* multiply by constant factor */
    for(j = 0; j < Nr; j++) 
      f[j] *= 3.0/4.0;
    return;
  } else if(*kerncode == BIWEIGHT) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      k = floor((xi - (double) 1.0)/dr); /* r[k] is 1 halfwidth below x[i] */
      if(k < 0) k = 0;
      if(k < Nr) {
	for(j = k; j < Nr; j++) {
	  vij = r[j] - xi;
	  if(vij >  1.0) break;  /* r[j] is 1 halfwidth above x[i] */
	  temp = (1 - vij * vij);
	  if(temp > 0.0) {
	    kvalue = temp * temp;
	    f[j] += kvalue * wi;
	  }
	}
      }
    }
    /* multiply by constant factor */
    for(j = 0; j < Nr; j++) 
      f[j] *= 15.0/16.0;
    return;
  } else if(*kerncode == COSINE) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      k = floor((xi - (double) 1.0)/dr); /* r[k] is 1 halfwidth below x[i] */
      if(k < 0) k = 0;
      if(k < Nr) {
	for(j = k; j < Nr; j++) {
	  vij = r[j] - xi;
	  if(vij > 1.0) break;  /* r[j] is 1 halfwidth above x[i] */
	  if(vij > -1.0) {
	    kvalue = 1.0 + cos(PI * vij);
	    f[j] += kvalue * wi;
	  }
	}
      }
    }
    /* divide by constant factor */
    for(j = 0; j < Nr; j++) 
      f[j] /= 2.0;
    return;
  } else if(*kerncode == OPTCOSINE) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      k = floor((xi - (double) 1.0)/dr); /* r[k] is 1 halfwidth below x[i] */
      if(k < 0) k = 0;
      if(k < Nr) {
	for(j = k; j < Nr; j++) {
	  vij = r[j] - xi;
	  if(vij > 1.0) break;  /* r[j] is 1 halfwidth above x[i] */
	  vij = ABS(vij);
	  if(vij < 1.0) {
	    kvalue = cos(vij * HALF_PI);
	    f[j] += kvalue * wi;
	  }
	}
      }
    }
    /* multiply by constant factor */
    for(j = 0; j < Nr; j++) 
      f[j] *= QUARTER_PI;
    return;
  } else {
    /* unrecognised kernel */
    *errcode = ERR_UNKNOWN_KERNEL;
    return;
  }
}

/* 
   ......................................................................
   bcolonel()  kernel density estimate with linear boundary correction at zero
   ......................................................................
*/

void bcolonel(
      int *kerncode,  /* integer code for kernel			  */
      int *nx,	      /* number of data values				  */
      double *x,      /* vector of data values				  */
      double *w,      /* vector of weights for data			  */
      int *nr,	      /* number of r values 		                  */
      double *r,      /* vector of r values (argument of density)         */
      double *nu0,    /* values of integral of kernel from -Inf to r	  */
      double *nu1,    /* values of integral of s k(s) from -Inf to r	  */
      double *nu2,    /* values of integral of s^2 k(s) from -Inf to r    */
      double *a,      /* scratch space 					  */
      double *b,      /* scratch space					  */
      double *f,      /* output values of density estimate f(r)           */
      int *errcode)   /* output integer code for errors                           */
{ 
  int i, j, Nx, Nr, bb, jbdry;
  double xi, wi, uij;
  double kvalue, denomj, temp, thresh, root2pi;

  /* extract arguments and validate */
  Nx = *nx;
  Nr = *nr;

  *errcode = OK;
  if(Nx < 0 || Nr <= 0) {
    *errcode = ERR_NEGATIVE_LENGTH;
    return;
  }

  /* initialise f(r) */
  for(j = 0; j < Nr; j++) 
    f[j] = 0.0;

  if(Nx == 0)
    return;

  /* compute coefficients a(r) and b(r) for boundary kernel */
  for(j = 0; j < Nr; j++) {
    denomj = nu0[j] * nu2[j] - nu1[j] * nu1[j];
    a[j] = nu2[j]/denomj;
    b[j] = nu1[j]/denomj;
  }

  /* determine jbdry (index of smallest r value outside support of kernel) */
  thresh = (*kerncode == GAUSSIAN) ? 3.0 : 1.0;
  for(jbdry = 0; jbdry < Nr; jbdry++) 
    if(r[jbdry] > thresh)
      break;

  /* go */
  if(*kerncode == GAUSSIAN) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      for(j = 0; j < jbdry; j++) {
	uij = r[j] - xi;
	kvalue = exp(- uij * uij/2);
	kvalue *= a[j] - b[j] * uij;
	f[j] += kvalue * wi;
      }
      for(j = jbdry; j < Nr; j++) {
	uij = r[j] - xi;
	kvalue = exp(- uij * uij/2);
	f[j] += kvalue * wi;
      }
    }
    /* divide by constant factor */
    for(j = 0; j < Nr; j++) 
      f[j] *= RECIPROCAL_SQRT_TWO_PI;
    return;
  } else if(*kerncode == RECTANGULAR) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      for(j = 0; j < jbdry; j++) {
	uij = r[j] - xi;
	if(INSIDE(uij)) {
	  kvalue = a[j] - b[j] * uij;
	  f[j] += kvalue * wi;
	}
      }
      for(j = jbdry; j < Nr; j++) {
	uij = r[j] - xi;
	if(INSIDE(uij))
	  f[j] += wi;
      }
    }
    /* divide by constant factor */
    for(j = 0; j < Nr; j++) 
      f[j] /= 2.0;
    return;
  } else if(*kerncode == TRIANGULAR) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      for(j = 0; j < jbdry; j++) {
	uij = r[j] - xi;
	if(INSIDE(uij)) {
	  kvalue = 1.0 - ABS(uij);
	  kvalue *= a[j] - b[j] * uij;
	  f[j] += kvalue * wi;
	}
      }
      for(j = jbdry; j < Nr; j++) {
	uij = r[j] - xi;
	if(INSIDE(uij)) {
	  kvalue = 1.0 - ABS(uij);
	  f[j] += kvalue * wi;
	}
      }
    }
    /* no constant factor */
    return;
  } else if(*kerncode == EPANECHNIKOV) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      for(j = 0; j < jbdry; j++) {
	uij = r[j] - xi;
	if(INSIDE(uij)) {
	  kvalue = (1 - uij * uij);
	  if(kvalue > 0) {
	    kvalue *= a[j] - b[j] * uij;
	    f[j] += kvalue * wi;
	  }
	}
      }
      for(j = jbdry; j < Nr; j++) {
	uij = r[j] - xi;
	if(INSIDE(uij)) {
	  kvalue = (1 - uij * uij);
	  if(kvalue > 0) 
	    f[j] += kvalue * wi;
	}
      }
    }
    /* multiply by constant factor */
    for(j = 0; j < Nr; j++) 
      f[j] *= 3.0/4.0;
    return;
  } else if(*kerncode == BIWEIGHT) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      for(j = 0; j < jbdry; j++) {
	uij = r[j] - xi;
	if(INSIDE(uij)) {
	  temp = (1 - uij * uij);
	  kvalue = temp * temp;
	  kvalue *= a[j] - b[j] * uij;
	  f[j] += kvalue * wi;
	}
      }
      for(j = jbdry; j < Nr; j++) {
	uij = r[j] - xi;
	if(INSIDE(uij)) {
	  temp = (1 - uij * uij);
	  kvalue = temp * temp;
	  f[j] += kvalue * wi;
	}
      }
    }
    /* multiply by constant factor */
    for(j = 0; j < Nr; j++) 
      f[j] *= 15.0/16.0;
    return;
  } else if(*kerncode == COSINE) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      for(j = 0; j < jbdry; j++) {
	uij = r[j] - xi;
	if(INSIDE(uij)) {
	  kvalue = 1.0 + cos(PI * uij);
	  kvalue *= a[j] - b[j] * uij;
	  f[j] += kvalue * wi;
	}
      }
      for(j = jbdry; j < Nr; j++) {
	uij = r[j] - xi;
	if(INSIDE(uij)) {
	  kvalue = 1.0 + cos(PI * uij);
	  f[j] += kvalue * wi;
	}
      }
    }
    /* divide by constant factor */
    for(j = 0; j < Nr; j++) 
      f[j] /= 2.0;
    return;
  } else if(*kerncode == OPTCOSINE) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      for(j = 0; j < jbdry; j++) {
	uij = r[j] - xi;
	if(INSIDE(uij)) {
	  kvalue = cos(uij * HALF_PI);
	  kvalue *= a[j] - b[j] * uij;
	  f[j] += kvalue * wi;
	}
      }
      for(j = jbdry; j < Nr; j++) {
	uij = r[j] - xi;
	if(INSIDE(uij)) {
	  kvalue = cos(uij * HALF_PI);
	  f[j] += kvalue * wi;
	}
      }
    }
    /* multiply by constant factor */
    for(j = 0; j < Nr; j++) 
      f[j] *= QUARTER_PI;
    return;
  } else {
    /* unrecognised kernel */
    *errcode = ERR_UNKNOWN_KERNEL;
    return;
  }
}


/* 
   ......................................................................
   fbcolonel()   Faster version of bcolonel() 
   ......................................................................
*/

void fbcolonel(
      int *kerncode,  /* integer code for kernel			  */
      int *nx,	      /* number of data values				  */
      double *x,      /* vector of data values				  */
      double *w,      /* vector of weights for data			  */
      int *nr,	      /* number of r values 		                  */
      double *r,      /* vector of r values (argument of density)         */
      double *nu0,    /* values of integral of kernel from -Inf to r	  */
      double *nu1,    /* values of integral of s k(s) from -Inf to r	  */
      double *nu2,    /* values of integral of s^2 k(s) from -Inf to r    */
      double *a,      /* scratch space 					  */
      double *b,      /* scratch space					  */
      double *f,      /* output values of density estimate f(r)           */
      int *errcode)   /* output integer code for errors                           */
{ 
  int i, j, Nx, Nr, bb, jbdry, jupperi;
  double xi, wi, uij;
  double kvalue, denomj, temp, thresh, root2pi;

  /* extract arguments and validate */
  Nx = *nx;
  Nr = *nr;

  *errcode = OK;
  if(Nx < 0 || Nr <= 0) {
    *errcode = ERR_NEGATIVE_LENGTH;
    return;
  }

  /* initialise f(r) */
  for(j = 0; j < Nr; j++) 
    f[j] = 0.0;

  if(Nx == 0)
    return;

  /* compute coefficients a(r) and b(r) for boundary kernel */
  for(j = 0; j < Nr; j++) {
    denomj = nu0[j] * nu2[j] - nu1[j] * nu1[j];
    a[j] = nu2[j]/denomj;
    b[j] = nu1[j]/denomj;
  }

  /* index of smallest r value outside support of kernel centred at zero */
  thresh = (*kerncode == GAUSSIAN) ? GAUSSTHRESH : 1.0;
  for(jbdry = 0; jbdry < Nr; jbdry++)
    if(r[jbdry] > thresh)
      break;
  
  /* go */
  if(*kerncode == GAUSSIAN) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      for(j = 0; j < jbdry; j++) {
	uij = r[j] - xi;
	kvalue = exp(- uij * uij/2);
	kvalue *= a[j] - b[j] * uij;
	f[j] += kvalue * wi;
      }
      for(j = jbdry; j < Nr; j++) {
	uij = r[j] - xi;
	if(uij > GAUSSTHRESH) break;  /* r[j] exceeds 8 sd above x[i] */
	kvalue = exp(- uij * uij/2);
	f[j] += kvalue * wi;
      }
    }
    /* divide by constant factor */
    for(j = 0; j < Nr; j++) 
      f[j] *= RECIPROCAL_SQRT_TWO_PI;
    return;
  } else if(*kerncode == RECTANGULAR) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      for(j = 0; j < jbdry; j++) {
	uij = r[j] - xi;
	if(INSIDE(uij)) {
	  kvalue = a[j] - b[j] * uij;
	  f[j] += kvalue * wi;
	}
      }
      for(j = jbdry; j < Nr; j++) {
	uij = r[j] - xi;
	if(uij > 1.0) break; /* r[j] exceeds 1 halfwidth above x[i] */
	if(uij > -1.0)
	  f[j] += wi;
      }
    }
    /* divide by constant factor */
    for(j = 0; j < Nr; j++) 
      f[j] /= 2.0;
    return;
  } else if(*kerncode == TRIANGULAR) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      for(j = 0; j < jbdry; j++) {
	uij = r[j] - xi;
	if(INSIDE(uij)) {
	  kvalue = 1.0 - ABS(uij);
	  kvalue *= a[j] - b[j] * uij;
	  f[j] += kvalue * wi;
	}
      }
      for(j = jbdry; j < Nr; j++) {
	uij = r[j] - xi;
	if(uij > 1.0) break; /* r[j] exceeds 1 halfwidth above x[i] */
	if(uij > -1.0) {
	  kvalue = 1.0 - ABS(uij);
	  f[j] += kvalue * wi;
	}
      }
    }
    /* no constant factor */
    return;
  } else if(*kerncode == EPANECHNIKOV) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      for(j = 0; j < jbdry; j++) {
	uij = r[j] - xi;
	if(INSIDE(uij)) {
	  kvalue = (1 - uij * uij);
	  if(kvalue > 0) {
	    kvalue *= a[j] - b[j] * uij;
	    f[j] += kvalue * wi;
	  }
	}
      }
      for(j = jbdry; j < Nr; j++) {
	uij = r[j] - xi;
	if(uij > 1.0) break; /* r[j] exceeds 1 halfwidth above x[i] */
	if(uij > -1.0) {
	  kvalue = (1 - uij * uij);
	  if(kvalue > 0) 
	    f[j] += kvalue * wi;
	}
      }
    }
    /* multiply by constant factor */
    for(j = 0; j < Nr; j++) 
      f[j] *= 3.0/4.0;
    return;
  } else if(*kerncode == BIWEIGHT) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      for(j = 0; j < jbdry; j++) {
	uij = r[j] - xi;
	if(INSIDE(uij)) {
	  temp = (1 - uij * uij);
	  kvalue = temp * temp;
	  kvalue *= a[j] - b[j] * uij;
	  f[j] += kvalue * wi;
	}
      }
      for(j = jbdry; j < Nr; j++) {
	uij = r[j] - xi;
	if(uij > 1.0) break; /* r[j] exceeds 1 halfwidth above x[i] */
	if(uij > -1.0) {
	  temp = (1 - uij * uij);
	  if(temp > 0.0) {
	    kvalue = temp * temp;
	    f[j] += kvalue * wi;
	  }
	}
      }
    }
    /* multiply by constant factor */
    for(j = 0; j < Nr; j++) 
      f[j] *= 15.0/16.0;
    return;
  } else if(*kerncode == COSINE) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      for(j = 0; j < jbdry; j++) {
	uij = r[j] - xi;
	if(INSIDE(uij)) {
	  kvalue = 1.0 + cos(PI * uij);
	  kvalue *= a[j] - b[j] * uij;
	  f[j] += kvalue * wi;
	}
      }
      for(j = jbdry; j < Nr; j++) {
	uij = r[j] - xi;
	if(uij > 1.0) break; /* r[j] exceeds 1 halfwidth above x[i] */
	if(uij > -1.0) {
	  kvalue = 1.0 + cos(PI * uij);
	  f[j] += kvalue * wi;
	}
      }
    }
    /* divide by constant factor */
    for(j = 0; j < Nr; j++) 
      f[j] /= 2.0;
    return;
  } else if(*kerncode == OPTCOSINE) {
    /* sum contributions */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      for(j = 0; j < jbdry; j++) {
	uij = r[j] - xi;
	if(INSIDE(uij)) {
	  kvalue = cos(uij * HALF_PI);
	  kvalue *= a[j] - b[j] * uij;
	  f[j] += kvalue * wi;
	}
      }
      for(j = jbdry; j < Nr; j++) {
	uij = r[j] - xi;
	if(uij > 1.0) break; /* r[j] exceeds 1 halfwidth above x[i] */
	if(uij > -1.0) {
	  kvalue = cos(uij * HALF_PI);
	  f[j] += kvalue * wi;
	}
      }
    }
    /* multiply by constant factor */
    for(j = 0; j < Nr; j++) 
      f[j] *= QUARTER_PI;
    return;
  } else {
    /* unrecognised kernel */
    *errcode = ERR_UNKNOWN_KERNEL;
    return;
  }
}



  
  

