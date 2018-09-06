#ifndef __MATH_MISER__
#define __MATH_MISER__

#include "Random.h"
#include "MathMatrix.h"
#include "MathVector.h"
#include "MathSobol.h"

// Monte Carlo Samples a user-supplied function in a rectangular volume
// specified by region[2][dim]. The total of ncalls are made to the function
// The integral of the function is returned in trgral and the standard
// deviation of this estimate is in stdev.
//

struct MiserStack
   {
   int    points;
   double weight;
   Matrix region;
   };

class MathMiser
   {
   public:
      SobolSequence sobol;

      long   ncall;

      double tgral;
      double stdev;

      VectorFunc * vfunc;

      MathMiser() : sobol()
         {
         ncall = 1000;
         }

      double Integrate(Matrix & region);

   protected:
      // local variables for integration
      // are here... to save on new / delete
      // calls

      void RandomPoint(Matrix & region, Vector & point);

      double func (Vector & v)
         { return vfunc->Evaluate(v); }

   private:
      MiserStack stack[32];      // should be good for at least 2^31 points

      Vector midpoint, point, minl, minr, maxl, maxr;

   };

#endif

