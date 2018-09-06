#include "MathLu.h"
#include "Error.h"

#include <math.h>

LU::~LU()
   { }

void LU::Decompose(Matrix & a)
   {
   Vector vv;           // stores the implict scaling of each row

   if (a.rows != a.cols)
      error("LU.Decompose: Matrix %s is not square", (const char *) a.label);

   lu.Copy(a);
   vv.Dimension(lu.rows);
   d = 1.0;

   permutation.Dimension(lu.rows);

   // loop over rows to get implicit scaling information
   for (int i = 0; i < lu.rows; i++)
      {
      double big = 0.0, temp;
      for (int j = 0; j < lu.rows; j++)
         if ( (temp = fabs(lu[i][j])) > big) big = temp;
      if (big == 0.0)
         error ("LU.Decompose: Matrix %s is singular", (const char *) a.label);
      vv[i] = 1.0 / big;
      }

   // Loop over columns as per Crout's method
   for (int j=0; j < lu.rows; j++)
      {
      // Uij = aij - Sum(1 to i - 1)[Lik*Uik]
      for (int i=0; i < j; i++)
         {
         double sum = lu[i][j];
         for (int k=0; k < i; k++)
            sum -= lu[i][k] * lu[k][j];
         lu[i][j] = sum;
         }

      // find the pivot element
      double big = 0.0;
      int    imax;

      // and compute Lij = 1/Ujj * { aij - Sum(1 to j - 1)[Lik*Uik] }
      for (int i = j; i < lu.rows; i++)
         {
         double sum = lu[i][j];
         for (int k = 0; k < j; k++)
            sum -= lu[i][k] * lu[k][j];
         lu[i][j] = sum;

         // check the figure of merit for this pivot
         double merit = vv[i] * fabs(sum);
         if (merit >= big)
            {
            big = merit;
            imax = i;
            }
         }

      // interchange rows if necessary
      if (j != imax)
         {
         lu.SwapRows(j, imax);
         d = -d;
         vv[imax] = vv[j];
         }

      permutation[j] = imax;

      if (lu[j][j] == 0.0)
          error("LU.Decompose: Matrix %s has zero pivot",(const char *)a.label);

      // finally divide by pivot element
      if (j != lu.rows - 1)
         {
         double scale = 1.0 / lu[j][j];
         for (int i = j + 1; i < lu.rows; i++)
            lu[i][j] *= scale;
         }
      }
   }

void LU::BackSubst(Vector & b)
   {
   x.Copy(b);

   // take into account the possibility that b starts with
   // a number of leading zeros (ie. for matrix inversion)

   int nonZero = -1, unscramble;

   // forward substitution with unscrambling of the permutation...
   for (int i = 0; i < lu.rows; i++)
      {
      unscramble = permutation[i];
      double sum = x[unscramble];
      x[unscramble] = x[i];

      if (nonZero != -1)
         for (int j = nonZero; j <= i - 1; j++)
            sum -= lu[i][j] * x[j];
      else
         if (sum)
            nonZero = i;
      x[i] = sum;
      }

   // Now do the backsubstitution
   for (int i = lu.rows - 1; i >= 0; i--)
      {
      double sum = x[i];
      for (int j = i + 1; j < lu.rows; j++)
         sum -= lu[i][j] * x[j];
      x[i] = sum / lu[i][i];
      }
   }

void LU::Invert()
   {
   inv.Dimension(lu.rows, lu.rows);

   inv.Identity();

   for(int i = 0; i < lu.rows; i++)
      {
      BackSubst(inv[i]);
      inv[i] = x;
      }
   }

double LU::Determinant()
   {
   double det = d;

   for (int i = 0; i < lu.rows; i++)
      det *= lu[i][i];

   return det;
   }

double LU::lnDeterminant()
   {
   bool minus_sign = d == -1;
   double    lnDet = 0.0;

   for (int i = 0; i < lu.rows; i++)
      if (lu[i][i] > 0)
         lnDet += log(lu[i][i]);
      else
         {
         lnDet += log(-lu[i][i]),
         minus_sign == !minus_sign;
         }

   if (minus_sign)
      error("LU::lnDeterminant cannot log negative value\n");

   return lnDet;
   }

