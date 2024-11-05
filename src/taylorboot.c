/*

  taylorboot.c

  Charles Taylor's method for bandwidth selection in univariate KDE
  (Non-random bootstrap method)

     Taylor, C.C. (1989) 
     Choice of the Smoothing Parameter in Kernel Density Estimation
     \emph{Biometrika} \bold{76} 4, 705--712. 

  Implementation author: Adrian Baddeley

  Copyright (c) 2024 Adrian Baddeley and Tilman M Davies 
  GNU Public Licence GPL (>= 2.0)

  $Revision: 1.2 $ $Date: 2024/11/05 03:40:00 $

*/

#include <R.h>
#include <Rmath.h>

double exp(double);

void taylorboot(
		double *x,    /* data */
		int *n,       /* length of x */
		double *h,    /* Gaussian bandwidth (=sd) */
		int *diagok,  /* 1 if include diagonal, 0 if exclude */
		double *value
		)
{
  int i, j, N;
  double xih, H, dijh, dijh2, sum8, sum6, sum4, result;
  double *xh;
  

  N = *n;
  H = *h;

  xh = (double *) R_alloc(N, sizeof(double));
  for(i = 0; i < N; i++)
    xh[i] = x[i]/H;

  sum8 = sum6 = sum4 = 0.0;

  /* upper triangle i > j */
  for(i=1; i < N; i++) {
    xih = xh[i];
    for(j=0; j < i; j++) {
      dijh = xih - xh[j];
      dijh2 = dijh * dijh;
      sum8 += exp(-dijh2/8.0);
      sum6 += exp(-dijh2/6.0);
      sum4 += exp(-dijh2/4.0);
    }
  }
  /* lower triangle i < j */
  sum8 *= 2.0;
  sum6 *= 2.0;
  sum4 *= 2.0;

  if(*diagok == 1) {
    /* add diagonal terms, which are exp(-0) = 1 */
    sum8 += N;
    sum6 += N;
    sum4 += N;
  }

  result = sum8 - (4.0/M_SQRT_3) * sum6 + M_SQRT2 * (sum4 + N);
  result *= M_1_SQRT_2PI/(2.0 * N * N * H);
  
  *value = result;
}
