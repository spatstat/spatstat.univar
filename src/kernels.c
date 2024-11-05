/*
   kernels.c

   C support for the basic kernels recognised by default.density

   'gaussian'
   'rectangular'
   'triangular'
   'epanechnikov'
   'biweight'
   'cosine'
   'optcosine'

   Includes pdf, cdf, partial moments and boundary-corrected kernels.

     d: probability density
     p: cumulative distribution function
     m1: partial first moment
     m2: partial second moment
     b: linear boundary kernel at 0

   Prototypes for all these functions are declared in 'kernels.h'

   Constants are declared in 'kerconstants.h'

   Terminology:
       template kernel: kernel with mean 0, halfwidth 1
       standard kernel: kernel with mean 0, standard deviation 1

   See 'kerconstants.h' for further information.

   The programmer encodes the functions d, p, m1, m2 for the TEMPLATE kernel only
   and the general case is derived by standard transformations.

   Copyright (c) 2008-2024 Adrian Baddeley, Tilman Davies and Martin Hazelton
   GNU Public Licence (>= 2.0)

*/

#include <Rmath.h>
double sqrt(double), exp(double), cos(double), sin(double);

#include "kerconstants.h"

/* 
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
              T E M P L A T E    K E R N E L S
   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   Prefix 'TEM' indicates the template kernel (mean = 0, halfwidth = 1)

   The functions d, p, m1, m2 for the template kernel are encoded by hand.

*/

/* ------------ GAUSSIAN ---------------------------------  */

/* Gaussian density */
double dTEMgaussian(double x) {
  double fx;
  fx = M_1_SQRT_2PI * exp(- x * x/2);
  return(fx);
}

/* Gaussian cdf */
double pTEMgaussian(double x) {
  double Fx;
  Fx = pnorm(x, 0.0, 1.0, (int) 1, (int) 0);
  return(Fx);
}

/* Gaussian partial first moment */
double m1TEMgaussian(double x) {
  double z;
  z = -dnorm(x, 0.0, 1.0, (int) 0);
  return(z);
}

/* Gaussian partial second moment */
double m2TEMgaussian(double x) {
  double z;
  z = pnorm(x, 0.0, 1.0, (int) 1, (int) 0) - x * dnorm(x, 0.0, 1.0, (int) 0);
  return(z);
}

/* ------------ RECTANGULAR ---------------------------------  */

/* Rectangular density */
double dTEMrectangular(double x) {
  double fx;
  fx = (x < -1.0) ? 0.0 : (x > 1.0) ? 0.0 : 0.5;
  return(fx);
}

/* Rectangular, cdf */
double pTEMrectangular(double x) {
  double Fx;
  Fx = (x < -1.0) ? 0.0 : (x > 1.0) ? 1.0 :
    ((x + 1.0)/2);
  return(Fx);
}

/* Rectangular, partial first moment */
double m1TEMrectangular(double x) {
  double z;
  z = (x < -1.0) ? 0.0 : (x > 1.0) ? 0.0 :
    ((x * x - 1.0)/4);
  return(z);
}

/* Rectangular, partial second moment */
double m2TEMrectangular(double x) {
  double z;
  z = (x < -1.0) ? 0.0 : (x > 1.0) ? TEMPLATE_VARIANCE(rectangular) :
    ((x * x * x + 1.0)/6);
  return(z);
}


/* ------------ TRIANGULAR ---------------------------------  */

/* Triangular density */
double dTEMtriangular(double x) {
  double fx;
  if(x < 0.0) x = -x;
  fx = 1.0 - x;
  if(fx < 0.0) return(0.0);
  return(fx);
}

/* Triangular, cdf */
double pTEMtriangular(double x) {
  double Fx;
  Fx = (x < -1.0) ? 0.0 : (x > 1.0) ? 1.0 :
    (x < 0.0) ? (0.5 + x + x*x/2) : (0.5 + x - x*x/2);
  return(Fx);
}

/* Triangular, partial first moment */
double m1TEMtriangular(double x) {
  double z;
  z = (x < -1.0) ? 0.0 : (x > 1.0) ? 0.0 :
    (x < 0.0) ? (x*x/2 + x*x*x/3 - 1.0/6) : (x*x/2 - x*x*x/3 - 1.0/6);
  return(z);
}

/* Triangular, partial second moment */
double m2TEMtriangular(double x) {
  double z;
  z = (x < -1.0) ? 0.0 : (x > 1.0) ? TEMPLATE_VARIANCE(triangular) :
    (x < 0.0) ? (x*x*x/3 + x*x*x*x/4 + 1.0/12) : (x*x*x/3 - x*x*x*x/4 + 1.0/12);
  return(z);
}


/* ------------ EPANECHNIKOV ---------------------------------  */

/* Epanechnikov density */
double dTEMepanechnikov(double x) {
  double z, fx;
  z = 1.0 - x * x;
  if(z < 0.0) return(0.0);
  fx = (3.0/4) * z;
  return(fx);
}

/* Epanechnikov, cdf */
double pTEMepanechnikov(double x) {
  double Fx;
  Fx = (x < -1.0) ? 0.0 : (x > 1.0) ? 1.0 :
    ((2.0 + 3 * x - x*x*x)/4);
  return(Fx);
}

/* Epanechnikov, partial first moment */
double m1TEMepanechnikov(double x) {
  double z;
  z = (x < -1.0) ? 0.0 : (x > 1.0) ? 0.0 :
    ((-3.0 + 6 * x*x - 3*x*x*x*x)/16);
  return(z);
}

/* Epanechnikov, partial second moment */
double m2TEMepanechnikov(double x) {
  double z;
  z = (x < -1.0) ? 0.0 : (x > 1.0) ? TEMPLATE_VARIANCE(epanechnikov) :
    ((2.0 + 5 * x*x*x - 3*x*x*x*x*x)/20);
  return(z);
}

/* ------------ BIWEIGHT ---------------------------------  */

/* Biweight density */
double dTEMbiweight(double x) {
  double z, fx;
  z = 1.0 - x * x;
  if(z < 0.0)
    return(0.0);
  fx = (15.0/16) * z * z;
  return(fx);
}

/* biweight, cdf */
double pTEMbiweight(double x) {
  double Fx;
  Fx = (x < -1.0) ? 0.0 : (x > 1.0) ? 1.0 :
    ((15 * x - 10 * x*x*x + 3 * x*x*x*x*x + 8.0)/16);
  return(Fx);
}

/* biweight, partial first moment */
double m1TEMbiweight(double x) {
  double z;
  z = (x < -1.0) ? 0.0 : (x > 1.0) ? 0.0 :
    ((5 * R_pow(x, 6) - 15 * R_pow(x, 4) + 15 * x*x - 5.0)/32);
  return(z);
}

/* biweight, partial second moment */
double m2TEMbiweight(double x) {
  double z;
  z = (x < -1.0) ? 0.0 : (x > 1.0) ? TEMPLATE_VARIANCE(biweight) :
    ((15 * R_pow(x, 7) - 42 * R_pow(x, 5) + 35 * R_pow(x, 3) + 8.0)/112);
  return(z);
}

/* ------------ COSINE ---------------------------------  */

/* cosine density */
double dTEMcosine(double x) {
  double fx;
  if(x < -1.0 || x > 1.0) return(0.0);
  fx = (1.0 + cos(M_PI * x))/2;
  return(fx);
}

/* cosine, cdf */
double pTEMcosine(double x) {
  double Fx;
  Fx = (x < -1.0) ? 0.0 : (x > 1.0) ? 1.0 :
    ((x + sin(M_PI * x)/M_PI + 1.0)/2);
  return(Fx);
}

/* cosine, partial first moment */
double m1TEMcosine(double x) {
  double z, pix;
  pix = M_PI * x;
  z = (x < -1.0) ? 0.0 : (x > 1.0) ? 0.0 :
    ((x*x-1.0)/4 +
     (pix*sin(pix) + cos(pix) + 1.0)/(2*M_PI*M_PI)
     );
  return(z);
}

/* cosine, partial second moment */
double m2TEMcosine(double x) {
  double z, pix;
  pix = M_PI * x;
  z = (x < -1.0) ? 0.0 : (x > 1.0) ? TEMPLATE_VARIANCE(cosine):
    ((x * x * x + 1.0)/6 +
     (
      (pix * pix - 2.0) * sin(pix) +
      2 * pix * cos(pix) - 2 * M_PI
     )/(2 * M_PI * M_PI * M_PI)
    );
  return(z);
}

/* ------------ OPTIMAL COSINE ---------------------------------  */

/* optcosine density */
double dTEMoptcosine(double x) {
  double fx;
  if(x < -1.0 || x > 1.0) return(0.0);
  fx = M_PI_4 * cos(M_PI * x/2);
  return(fx);
}

/* optcosine, cdf */
double pTEMoptcosine(double x) {
  double Fx;
  Fx = (x < -1.0) ? 0.0 : (x > 1.0) ? 1.0 :
    ((sin(M_PI * x/2) + 1.0)/2);
  return(Fx);
}

/* optcosine, partial first moment */
double m1TEMoptcosine(double x) {
  double z, pi2x;
  pi2x = M_PI_2 * x;
  z = (x < -1.0) ? 0.0 : (x > 1.0) ? 0.0:
    ((pi2x * sin(pi2x) + cos(pi2x) - M_PI_2)/M_PI);
  return(z);
}

/* optcosine, partial second moment */
double m2TEMoptcosine(double x) {
  double z, pi2x;
  pi2x = M_PI_2 * x;
  z = (x < -1.0) ? 0.0 : (x > 1.0) ? TEMPLATE_VARIANCE(optcosine):
    ((2.0/(M_PI * M_PI)) * (
			    (pi2x * pi2x - 2.0) * sin(pi2x) +
			    2 * pi2x * cos(pi2x) +
			    M_PI_2 * M_PI_2 - 2.0
			    )
     );
  return(z);
}

/*
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          D E R I V E D    F U N C T I O N S
   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*/

#define TEMPLATE_PDF(KERNEL) dTEM ## KERNEL
#define TEMPLATE_CDF(KERNEL) pTEM ## KERNEL
#define TEMPLATE_FIRST(KERNEL) m1TEM ## KERNEL
#define TEMPLATE_SECOND(KERNEL) m2TEM ## KERNEL

#define GENERAL_PDF(KERNEL) d ## KERNEL
#define GENERAL_CDF(KERNEL) p ## KERNEL
#define GENERAL_FIRST(KERNEL) m1 ## KERNEL
#define GENERAL_SECOND(KERNEL) m2 ## KERNEL

#define GENERAL_BDRYKERN(KERNEL) b ## KERNEL

/* 
   Define probability densities, for general mean and standard deviation
*/

/* Gaussian density is a special case */
double dgaussian(double x, double mean, double sd) {
  double fx;
  fx = dnorm(x, mean, sd, (int) 0);
  return(fx);
}

/* define drectangular, dtriangular etc */

#define DEFINE_GENERAL_PDF(KERNEL) \
  double GENERAL_PDF(KERNEL)(double x, double mean, double sd) {	\
  double h, y, fx; \
  h = sd * HALFWIDTH_ON_SIGMA(KERNEL); \
  y = (x - mean)/h; \
  fx = TEMPLATE_PDF(KERNEL)(y)/h;		\
  return(fx); \
}

DEFINE_GENERAL_PDF(rectangular)

DEFINE_GENERAL_PDF(triangular)

DEFINE_GENERAL_PDF(epanechnikov)

DEFINE_GENERAL_PDF(biweight)

DEFINE_GENERAL_PDF(cosine)

DEFINE_GENERAL_PDF(optcosine)

/* 
    Cumulative distribution functions
*/

/* Gaussian CDF is special case */
  
double pgaussian(double x, double mean, double sd) {
  double Fx;
  Fx = pnorm(x, mean, sd, (int) 1, (int) 0);
  return(Fx);
}

/* define prectangular, ptriangular etc */

#define DEFINE_GENERAL_CDF(KERNEL)				 \
  double GENERAL_CDF(KERNEL)(double x, double mean, double sd) {	\
  double h, y, Fx; \
  h = sd * HALFWIDTH_ON_SIGMA(KERNEL) ; \
  y = (x - mean)/h; \
  Fx = TEMPLATE_CDF(KERNEL)(y);			\
  return(Fx); \
}


DEFINE_GENERAL_CDF(rectangular)

DEFINE_GENERAL_CDF(triangular)

DEFINE_GENERAL_CDF(epanechnikov)

DEFINE_GENERAL_CDF(biweight)

DEFINE_GENERAL_CDF(cosine)

DEFINE_GENERAL_CDF(optcosine)


/* ----------- PARTIAL MOMENTS  -------------------- */

/* 
  The partial moment of order m is the function
     a_m(x) = \int_{-\infty}^x t^m f(t) dt
   where f is the probability density.

*/


/* 
    Partial first moments for general case
*/

#define DEFINE_GENERAL_FIRST(KERNEL) \
  double GENERAL_FIRST(KERNEL)(double x, double mean, double sd) { \
  double h, y, z; \
  h = sd * HALFWIDTH_ON_SIGMA(KERNEL); \
  y = (x - mean)/h; \
  z = mean * TEMPLATE_CDF(KERNEL)(y) + h * TEMPLATE_FIRST(KERNEL)(y); \
  return(z); \
}

DEFINE_GENERAL_FIRST(gaussian)
  
DEFINE_GENERAL_FIRST(rectangular)

DEFINE_GENERAL_FIRST(triangular)

DEFINE_GENERAL_FIRST(epanechnikov)

DEFINE_GENERAL_FIRST(biweight)

DEFINE_GENERAL_FIRST(cosine)

DEFINE_GENERAL_FIRST(optcosine)

/* 
    Partial second moments for general case
*/

#define DEFINE_GENERAL_SECOND(KERNEL) \
  double GENERAL_SECOND(KERNEL)(double x, double mean, double sd) { \
  double h, y, z; \
  h = sd * HALFWIDTH_ON_SIGMA(KERNEL); \
  y = (x - mean)/h; \
  z = mean * mean * TEMPLATE_CDF(KERNEL)(y) + \
      2 * mean * h * TEMPLATE_FIRST(KERNEL)(y) + \
      h * h * TEMPLATE_SECOND(KERNEL)(y); \
  return(z); \
}

DEFINE_GENERAL_SECOND(gaussian)
  
DEFINE_GENERAL_SECOND(rectangular)

DEFINE_GENERAL_SECOND(triangular)

DEFINE_GENERAL_SECOND(epanechnikov)

DEFINE_GENERAL_SECOND(biweight)

DEFINE_GENERAL_SECOND(cosine)

DEFINE_GENERAL_SECOND(optcosine)

/* ------- LINEAR BOUNDARY KERNELS ---------------- */
/* query point = x, data point = mean */

#define DEFINE_GENERAL_BDRYKERN(KERNEL) \
  double GENERAL_BDRYKERN(KERNEL)(double x, double mean, double sd) { \
  double h, p, u, a0, a1, a2, fy; \
  fy = GENERAL_PDF(KERNEL)(x, mean, sd); \
  if(fy == 0.0) return(0.0); \
  h = HALFWIDTH_ON_SIGMA(KERNEL) * sd; \
  p = x/h; \
  u = (x-mean)/h; \
  a0 = TEMPLATE_CDF(KERNEL)(p); \
  a1 = TEMPLATE_FIRST(KERNEL)(p); \
  a2 = TEMPLATE_SECOND(KERNEL)(p); \
  fy *= (a2 - a1 * u)/(a0 * a2 - a1 * a1); \
  return(fy); \
}

DEFINE_GENERAL_BDRYKERN(gaussian)
  
DEFINE_GENERAL_BDRYKERN(rectangular)

DEFINE_GENERAL_BDRYKERN(triangular)

DEFINE_GENERAL_BDRYKERN(epanechnikov)

DEFINE_GENERAL_BDRYKERN(biweight)

DEFINE_GENERAL_BDRYKERN(cosine)

DEFINE_GENERAL_BDRYKERN(optcosine)

