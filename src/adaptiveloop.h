/* 
   adaptiveloop.h

   Code for the actual computation of adaptive kernel estimates

   This file is #included multiple times in 'adaptive.h',
   once for each kernel.

   Uses macros ZEROCOR and KERNELNAME

   Copyright (c) 2008-2024 Adrian Baddeley, Tilman Davies and Martin Hazelton
   GNU Public Licence (>= 2.0)

*/

#define PREFIX(A, KNAME) A ## KNAME
#define CDF(KNAME)         PREFIX(p, KNAME)
#define DENSITY(KNAME)     PREFIX(d, KNAME)
#define BDRYDENSITY(KNAME) PREFIX(b, KNAME)

  /* sum contributions from each data point x[i] */
    for(i = 0; i < Nx; i++) {
      xi = x[i];
      wi = w[i];
      si = sd[i];
#if (ZEROCOR == WEIGHTED)
      /* 
	 divide by mass of kernel on positive half-line 
	 so that total mass is conserved 
      */
      kmass = 1.0 - CDF(KERNELNAME) (0.0, xi, si);
      wi /= kmass;
#endif
      /* visit each query point r[j] */
      for(j = 0; j < Nr; j++) {
	rj = r[j];
#if (ZEROCOR == BDRYKERNEL)	
	kvalue = BDRYDENSITY(KERNELNAME) (rj, xi, si);
#else	
	kvalue = DENSITY(KERNELNAME) (rj, xi, si);
#endif	
#if (ZEROCOR == REFLECTION)
	/* add contribution from reflected data point -xi */
	kvalue += DENSITY(KERNELNAME) (rj, -xi, si);
#elif (ZEROCOR == CONVOLUTION)
	/* renormalise to give unbiased estimator for uniform density */
	kvalue /= 1.0 - CDF(KERNELNAME) (0.0, rj, si);
#endif	
	f[j] += kvalue * wi;
      }
    }

#undef CDF
#undef DENSITY
#undef BDRYDENSITY
#undef PREFIX

