/*
  access.c

  R interface functions
  just to check the validity of internal code

*/

#include "interfacecodes.h"
#include "kernels.h"

void kermom(int *nx,     
	    double *x,
	    double *mean,
	    double *sd,
	    int *m, 
	    int *kerncode,
	    double *y,
	    int *errcode) {
  int i, n;
  n = *nx;
  *errcode = OK;
  
  if(*m == 0) {
    /* cumulative distribution function */
    if(*kerncode == GAUSSIAN) {
      for(i = 0; i < n; i++) y[i] = pgaussian(x[i], mean[i], sd[i]);
    } else if(*kerncode == RECTANGULAR){ 
      for(i = 0; i < n; i++) y[i] = prectangular(x[i], mean[i], sd[i]);
    } else if(*kerncode == TRIANGULAR) {
      for(i = 0; i < n; i++) y[i] = ptriangular(x[i], mean[i], sd[i]);
    } else if(*kerncode == EPANECHNIKOV) {
      for(i = 0; i < n; i++) y[i] = pepanechnikov(x[i], mean[i], sd[i]);
    } else if(*kerncode == BIWEIGHT) {     
      for(i = 0; i < n; i++) y[i] = pbiweight(x[i], mean[i], sd[i]);
    } else if(*kerncode == COSINE) {       
      for(i = 0; i < n; i++) y[i] = pcosine(x[i], mean[i], sd[i]);
    } else if(*kerncode == OPTCOSINE) {
      for(i = 0; i < n; i++) y[i] = poptcosine(x[i], mean[i], sd[i]);
    }
  } else if(*m == 1) {
    /* partial first moment */
    if(*kerncode == GAUSSIAN) {
      for(i = 0; i < n; i++) y[i] = m1gaussian(x[i], mean[i], sd[i]);
    } else if(*kerncode == RECTANGULAR){ 
      for(i = 0; i < n; i++) y[i] = m1rectangular(x[i], mean[i], sd[i]);
    } else if(*kerncode == TRIANGULAR) {
      for(i = 0; i < n; i++) y[i] = m1triangular(x[i], mean[i], sd[i]);
    } else if(*kerncode == EPANECHNIKOV) {
      for(i = 0; i < n; i++) y[i] = m1epanechnikov(x[i], mean[i], sd[i]);
    } else if(*kerncode == BIWEIGHT) {     
      for(i = 0; i < n; i++) y[i] = m1biweight(x[i], mean[i], sd[i]);
    } else if(*kerncode == COSINE) {       
      for(i = 0; i < n; i++) y[i] = m1cosine(x[i], mean[i], sd[i]);
    } else if(*kerncode == OPTCOSINE) {
      for(i = 0; i < n; i++) y[i] = m1optcosine(x[i], mean[i], sd[i]);
    }
  } else if(*m == 2) {
    /* partial second moment */
    if(*kerncode == GAUSSIAN) {
      for(i = 0; i < n; i++) y[i] = m2gaussian(x[i], mean[i], sd[i]);
    } else if(*kerncode == RECTANGULAR){ 
      for(i = 0; i < n; i++) y[i] = m2rectangular(x[i], mean[i], sd[i]);
    } else if(*kerncode == TRIANGULAR) {
      for(i = 0; i < n; i++) y[i] = m2triangular(x[i], mean[i], sd[i]);
    } else if(*kerncode == EPANECHNIKOV) {
      for(i = 0; i < n; i++) y[i] = m2epanechnikov(x[i], mean[i], sd[i]);
    } else if(*kerncode == BIWEIGHT) {     
      for(i = 0; i < n; i++) y[i] = m2biweight(x[i], mean[i], sd[i]);
    } else if(*kerncode == COSINE) {       
      for(i = 0; i < n; i++) y[i] = m2cosine(x[i], mean[i], sd[i]);
    } else if(*kerncode == OPTCOSINE) {
      for(i = 0; i < n; i++) y[i] = m2optcosine(x[i], mean[i], sd[i]);
    }
  } else {
    *errcode = NOT_SUPPORTED;
    return;
  }
}
