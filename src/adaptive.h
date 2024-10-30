/* 

   adaptive.h

   Adaptive kernel density estimation, brute force, arbitrary weights,
   optional boundary correction for positive data

   This file is #included multiple times in adaptive.c 
   using different values of the macro 'ZEROCOR'

   values of ZEROCOR:

   NONE        no correction
   WEIGHTED    weighting boundary correction 
   CONVOLUTION convolution boundary correction
   REFLECTION  reflection boundary correction
   BDRYKERNEL  boundary kernel

   Copyright (c) 2008-2024 Adrian Baddeley, Tilman Davies and Martin Hazelton
   GNU Public Licence (>= 2.0)

*/


void FNAME(
	   int *kerncode, /* integer code for kernel                   */
	   int *nx,       /* number of data values	               */
	   double *x,	  /* vector of data values		       */
	   double *sd,	  /* vector of kernel bandwidths for data      */
	   double *w,	  /* vector of weights for data		       */
	   int *nr,	  /* number of r values			       */
	   double *r,	  /* vector of r values (argument of density)  */
	   double *f,	  /* resulting values of density estimate f(r) */
	   int *errcode)  /* integer code for errors                   */
{ 
  int i, j, Nx, Nr;
  double xi, si, wi, rj, kvalue;

#if (ZEROCOR == WEIGHTED)
  double kmass;
#endif
  
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
#undef KERNELNAME
#define KERNELNAME gaussian    
#include "adaptiveloop.h"
    return;
  } else if(*kerncode == RECTANGULAR) {
#undef KERNELNAME
#define KERNELNAME rectangular
#include "adaptiveloop.h"
    return;
  } else if(*kerncode == TRIANGULAR) {
#undef KERNELNAME
#define KERNELNAME triangular
#include "adaptiveloop.h"
    return;
  } else if(*kerncode == EPANECHNIKOV) {
#undef KERNELNAME
#define KERNELNAME epanechnikov
#include "adaptiveloop.h"
    return;
  } else if(*kerncode == BIWEIGHT) {
#undef KERNELNAME
#define KERNELNAME biweight
#include "adaptiveloop.h"
    return;
  } else if(*kerncode == COSINE) {
#undef KERNELNAME
#define KERNELNAME cosine
#include "adaptiveloop.h"
    return;
  } else if(*kerncode == OPTCOSINE) {
#undef KERNELNAME
#define KERNELNAME optcosine
#include "adaptiveloop.h"
    return;
  } else {
    /* unrecognised kernel */
    *errcode = ERR_UNKNOWN_KERNEL;
    return;
  }
}

#undef KERNELNAME

