#include <R.h>
#include <Rinternals.h>

/*
  Prototype declarations for all native routines in spatstat.univar package

  Automatically generated - do not edit! 

*/

/*
  
                  Functions invoked by .C

*/

void kermom(int *, double *, double *, double *, int *, int *, double *, int *);
void adaptiveKDE(int *, int *, double *, double *, double *, int *, double *, double *, int *); 
void adaptKDEweight(int *, int *, double *, double *, double *, int *, double *, double *, int *); 
void adaptKDEreflect(int *, int *, double *, double *, double *, int *, double *, double *, int *); 
void adaptKDEconvol(int *, int *, double *, double *, double *, int *, double *, double *, int *); 
void adaptKDEbdry(int *, int *, double *, double *, double *, int *, double *, double *, int *);
void taylorboot(double *, int *, double *, int *, double *);
void fcolonel(int *, int *, double *, double *, int *, double *, double *, int *); 
void colonel(int *, int *, double *, double *, int *, double *, double *, int *); 
void fbcolonel(int *, int *, double *, double *, int *, double *, double *, double *, double *, double *, double *, double *, int *); 
void bcolonel(int *, int *, double *, double *, int *, double *, double *, double *, double *, double *, double *, double *, int *);
void tabsumweight(int *, double *, double *, int *, double *, double *);
void hotrodInsul(int *, double *, double *, double *, double *, int *, double *); 
void hotrodAbsorb(int *, double *, double *, double *, double *, int *, double *);
/*

             Functions invoked by .Call

*/
SEXP Cwhist(SEXP, SEXP, SEXP);
