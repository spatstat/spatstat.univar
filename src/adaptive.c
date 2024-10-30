/* 

   adaptive.c

   Adaptive kernel density estimation, brute force, 
   arbitrary weights, boundary correction

   Copyright (c) 2008-2024 Adrian Baddeley, Tilman Davies and Martin Hazelton
   GNU Public Licence (>= 2.0)

*/

#include "kernels.h"
#include "interfacecodes.h"

#define FNAME adaptiveKDE
#define ZEROCOR NONE
#include "adaptive.h"
#undef FNAME
#undef ZEROCOR

#define FNAME adaptKDEweight
#define ZEROCOR WEIGHTED
#include "adaptive.h"
#undef FNAME
#undef ZEROCOR

#define FNAME adaptKDEreflect
#define ZEROCOR REFLECTION
#include "adaptive.h"
#undef FNAME
#undef ZEROCOR

#define FNAME adaptKDEconvol
#define ZEROCOR CONVOLUTION
#include "adaptive.h"
#undef FNAME
#undef ZEROCOR

#define FNAME adaptKDEbdry
#define ZEROCOR BDRYKERNEL
#include "adaptive.h"
#undef FNAME
#undef ZEROCOR



