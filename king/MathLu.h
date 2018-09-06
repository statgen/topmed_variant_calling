#ifndef __MATH_LU__
#define __MATH_LU__

#include "MathMatrix.h"
#include "MathVector.h"
#include "IntArray.h"

class LU
   {
   public:
      Matrix   lu, inv;
      Vector   x;
      IntArray permutation;
      double   d;

   LU() : lu("LU.LU"), x("LU.x"), inv("LU.inv") { }
   ~LU();

   // Given a square matrix a, decomposes a permutation of a
   // into an LU product, stored in LU as follows:
   //    Lij = LUij when i > j; 1.0 when i == j; 0.0 otherwise
   //    Uij = LUij when i <= j; 0.0 otherwise
   // permutation[1..n] records the permutation effected by
   // partial pivoting
   // d is output as +1 or -1 depending on whether the number
   // of row interchanges was even or odd
   // (for calculating determinants)
   void Decompose(Matrix & a);

   // Solves LU*X = B, taking b as the right hand side vector
   // and storing the solution in x.
   void BackSubst(Vector & b);

   // Calculate matrix inverse by backsubstituting basis vectors
   void Invert();

   // Calculate determinant
   double Determinant();

   // Calculate log of determinant
   double lnDeterminant();

   };

#endif
