#ifndef __MATH_SOBOL__
#define __MATH_SOBOL__

#include "IntArray.h"
#include "MathVector.h"

#define POLY_COUNT   36
#define SOBOL_BITS   30
#define SOBOL_FACTOR (1.0 / (1L << SOBOL_BITS))

class SobolSequence
   {
   public:
      IntArray * bits;
      IntArray   x;
      int        dim;
      long       counter;

      SobolSequence();
      ~SobolSequence();

      void     Init(int dimensions);
      Vector & Next(Vector & point);

   private:
      static int poly_integers[POLY_COUNT];
      static int poly_degrees[POLY_COUNT];
   };

#endif

