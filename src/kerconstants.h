/* 
   kerconstants.h

   Constants for smoothing kernels

   Terminology:
       template kernel: kernel with mean 0, halfwidth 1
       standard kernel: kernel with mean 0, standard deviation 1

   KFAC_XXX: ratio of halfwidth to standard deviation
   TVAR_XXX: variance of template kernel

   For densities with compact support (namely all except the Gaussian)
   the halfwidth 'h' is half the width of the support.

   The halfwidth is proportional to the standard deviation,
        h = c. sigma
   where 'c' is the constant 'KFAC' defined below.

   For the Gaussian density, h=sigma and c=1 by convention.

   The constants 'TVAR' give the variance of the template kernel.
   The variance of the kernel with halfwidth h is
      TVAR * h^2
   The variance of the standard kernel is TVAR * KFACTOR^2.
   The variance of the kernel with standard deviation sigma is 
      TVAR * sigma^2 * KFACTOR^2
   
   Copyright (c) 2008-2024 Adrian Baddeley, Tilman Davies and Martin Hazelton
   GNU Public Licence (>= 2.0)

*/

#define KFAC_gaussian     1.0
#define KFAC_rectangular  sqrt(3.0)
#define KFAC_triangular   sqrt(6.0)
#define KFAC_epanechnikov sqrt(5.0)
#define KFAC_biweight     sqrt(7.0)
#define KFAC_cosine       (1.0/sqrt(1.0/3 - 2.0/(M_PI * M_PI)))
#define KFAC_optcosine    (1.0/sqrt(1.0 - 8.0/(M_PI * M_PI)))

#define TVAR_gaussian     1.0
#define TVAR_rectangular  (1.0/3)
#define TVAR_triangular   (1.0/6)
#define TVAR_epanechnikov (1.0/5)
#define TVAR_biweight     (1.0/7)
#define TVAR_cosine       (1.0/3 - 2.0/(M_PI * M_PI))
#define TVAR_optcosine    (1.0 - 8.0/(M_PI * M_PI))

#define HALFWIDTH_ON_SIGMA(KERNELNAME) KFAC_ ## KERNELNAME
#define TEMPLATE_VARIANCE(KERNELNAME) TVAR_ ## KERNELNAME
