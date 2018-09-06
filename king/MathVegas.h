#ifndef __MATHVEGAS_H__
#define __MATHVEGAS_H__

#include "MathVector.h"
#include "MathMatrix.h"
#include "Random.h"

#define ALPH   1.5
#define NDMX   50       // Maximum number of increments for each axis
#define MXDIM  10       // Maximum number of dimensions

// Monte-carlo integration of user supplied ndimensional
// function, in a rectangular volume specified by matrix
// Volume[2][ndim], consisting of lower and upper bounds
// itmx iterations each with about ncall function calls
// The sampling grid is refined iteratively. Produces
// the integral tgral, with standard deviation sd, and
// an indicator of integrity chi2a (should be less than 1).
class Vegas
   {
   public:
      int    itmx, ncall;
      double tgral, sd, chi2a;
      static Random rand;
      VectorFunc * vfunc;

   Vegas();
   ~Vegas();

   double func(Vector & point)
      { return vfunc->Evaluate(point); }

   // Three levels of initialization possible
   // 0 - Total reset
   // 1 - Keep Grid, clear Estimates
   // 2 - Keep Grid and Estimates
   // 3 - Do additional iterations, no changes
   void Init(Matrix & Volume, int level = 0);

   // Integrate the function
   double Integrate(Matrix & Volume);

   private:
      void   Rebin(double rc, Vector & xi);

      int    mds, nd, ndo, ng, npg, * ia, * kg;
      double calls, dv2g, dxg, rc;
      double wgt, xjac, xn, xnd, xo, schi, si, swgt;
      Vector dt, dx, r, x, xin;
      Matrix d, di, xi;
   };

#endif
