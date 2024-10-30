/*
  interfacecodes.h

  Numerical codes for the interface between R and C

   Copyright (c) 2008-2024 Adrian Baddeley, Tilman Davies and Martin Hazelton
   GNU Public Licence (>= 2.0)

*/

/* kernel codes */
#define GAUSSIAN     1
#define RECTANGULAR  2
#define TRIANGULAR   3
#define EPANECHNIKOV 4
#define BIWEIGHT     5
#define COSINE       6
#define OPTCOSINE    7

/* error codes */
#define OK                         0
#define ERR_NEGATIVE_LENGTH        1
#define ERR_UNKNOWN_KERNEL         2
#define NOT_SUPPORTED              3

/* codes for zero boundary correction (density on positive half-line) */
#define NONE        0
#define WEIGHTED    1
#define CONVOLUTION 2
#define REFLECTION  3
#define BDRYKERNEL  4

