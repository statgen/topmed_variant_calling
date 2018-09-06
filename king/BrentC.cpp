#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "analysis.h"

/// <summary>Brent's method for minimizing a 1d function</summary>
// Machine eps
double Engine::MACHEPS = (double)2.2204460492503131e-016;
double Engine::MACHEPS_SQRT = sqrt(MACHEPS);
double Engine::cbrent = ((double)3.0 - sqrt((double)5.0)) / (double)2.0;
   
/// <summary>Minimize the function over the interval [a, b]</summary>
/// <param name="f">Function to minimize</param>
/// <param name="a">Left side of the bracket</param>
/// <param name="b">Right side of the bracket</param>
/// <param name="eps">Stopping tolerance</param>
/// <param name="funcx">Function evaluated at the minimum</param>
/// <param name="numiter">number of function evaluations</param>
/// <param name="maxIter">maximum number of function evaluations allowed</param>
/// <param name="quiet">print out function evaluations?</param>
/// <returns>Point that minimizes the function</returns>
/// <remarks>This implements the algorithm from Brent's book, "Algorithms for Minimization without Derivatives"</remarks>
//double BrentC::minimize(BrentFunctor &f, double a, double b, double eps, double &funcx, size_t &numiter, size_t maxIter, bool quiet)
double Engine::minimize(double a, double b, double eps, double &funcx, int &numiter, int maxIter, bool quiet)
{
   if (a >= b){
      printf("Exception: a must be < b");
      throw(1);
   }

   double x = a + cbrent * (b - a);
   double v = x;
   double w = x;
   double e = 0;

   double fx = fLL(x);
   double fv = fx;
   double fw = fx;

   numiter = 0;

   while (true)
   {
      double m = (double)0.5 * (a + b);
      double tol = MACHEPS_SQRT * abs(x) + eps;
      double tol2 = (double)2 * tol;

      // Check the stopping criterion
      if (abs(x - m) <= tol2 - 0.5 * (b - a)){ break; }

      // Stop if we've exceeded the maximum number of iterations
      numiter++;
      if (numiter > maxIter){
         printf("Exception: Exceeded maximum number of iterations.");
         throw(2);
      }
      double p = 0.0, q = 0.0, r = 0.0;
      double d = 0.0;
      double u = 0.0;

      if (abs(e) > tol)
      {
         // Fit parabola
         r = (x - w) * (fx - fv);
         q = (x - v) * (fx - fw);
         p = (x - v) * q - (x - w)*r;
         q = (double)2.0 * (q - r);
         if (q > (double)0.0)
            p = -p;
         else
            q = -q;
         r = e;
         e = d;
      }

      if ((abs(p) < abs((double)0.5*q*r)) && (p < q*(a-x)) && (p < q*(b-x)))
      {
         // Parabolic interpolation step
         d = p / q;
         u = x + d;
         // f must not be evaluated too close to a or b
         if (u - a < tol2 || b - u < tol2)
            d = (x < m) ? tol : -tol;
      }
      else
      {
         // Golden section step
         e = (x < m) ? b - x : a - x;
         d = cbrent * e;
      }

      // f must not be evaluated too close to x
      if (abs(d) >= tol)
         u = x + d;
      else if (d > 0.0)
         u = x + tol;
      else
         u = x - tol;
      double fu = fLL(u);

      // Update
      if (fu <= fx)
      {
         if (u < x)
            b = x;
         else
            a = x;
         v = w; fv = fw;
         w = x; fw = fx;
         x = u; fx = fu;
      }
      else
      {
         if (u < x)
            a = u;
         else
            b = u;

         if (fu <= fw || w == x)
         {
            v = w; fv = fw;
            w = u; fw = fu;
         }
         else if (fu <= fv || v == x || v == w)
         {
            v = u; fv = fu;
         }
      }

      if ( !quiet ){
         printf("Iteration %d, min_x = %.4lf, f(min_x) = %.4lf\n", numiter, x, fx);

      /*
#if defined( _MSC_VER )    // Windows/VC uses a %Iu specifier for size_t
            const char *szFmt1 = "Iteration %Iu, min_x = %f, f(min_x) = %f, ";
#else                      // Linux/g++ uses a %zu specifier for size_t
            const char *szFmt1 = "Iteration %zu, min_x = %f, f(min_x) = %f, ";
#endif
         printf( szFmt1, numiter, x, fx );
         */
      }
   }
   funcx = fLL(x);
   return x;
}

