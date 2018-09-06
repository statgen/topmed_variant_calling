#ifndef __MATH_CHOLESKY__
#define __MATH_CHOLESKY__

#include "MathMatrix.h"
#include "MathVector.h"

class Cholesky
   {
   public:
      Matrix      L, inv;
      Vector      x;

   Cholesky() : L("cholesky.L"), inv("cholesky.inverse"), x("cholesky.x")
      { }

   ~Cholesky()
      { }

   // Given a symmetric positive definite matrix A finds
   // a lower triangular matrix L such that L * transpose(L) = A
   // Only the upper triangle of A need be given
   void Decompose(Matrix & A);

   // If you call fast decompose the upper triangle of U is
   // undefined (as opposed to zero). This is often okay and
   // allows for a little more speed...
   void FastDecompose(Matrix & A);

   // Tries to decompose matrix A, returning true on success
   // or zero on failure ... you should also check that
   // determinant is not zero before using results if this
   // is a concern
   bool TryDecompose(Matrix & A);

   // solve Y = X b
   void BackSubst(Vector & b);
   void BackSubst0(Vector & b);   
   void Invert();

   // determinant functions
   double lnDeterminantL();
   double DeterminantL();

   double lnDeterminant()
      {
      return 2 * lnDeterminantL();
      }
   double Determinant()
      {
      double temp = DeterminantL();
      return temp * temp;
      }
   };

#endif
