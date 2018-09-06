#include "MathCholesky.h"
#include "Error.h"

#include <math.h>

void Cholesky::Decompose(Matrix & A)
   {
   L.Dimension(A.rows, A.rows);
   L.Zero();
   FastDecompose(A);
   }

void Cholesky::FastDecompose(Matrix & A)
   {
   if (A.rows != A.cols)
      error("Cholesky.Decompose: Matrix %s is not square",
            (const char *) A.label);

   L.Dimension(A.rows, A.rows);

   for (int i=0; i<L.rows; i++)
      for (int j=i; j<L.rows; j++)
         {
         double sum = A.data[i]->data[j];
         for (int k = i - 1; k >= 0; k--)
            sum -= L.data[i]->data[k] * L.data[j]->data[k];
         if (i == j)
            if (sum <= 0.0)
               error("Cholesky - matrix %s is not positive definite",
                     (const char *) A.label);
            else
               L.data[i]->data[i] = sqrt(sum);
         else
            L.data[j]->data[i] = sum / L.data[i]->data[i];
         }
   }

bool Cholesky::TryDecompose(Matrix & A)
   {
   L.Dimension(A.rows, A.rows);
   L.Zero();

   if (A.rows != A.cols)
      return false;

   L.Dimension(A.rows, A.rows);

   for (int i=0; i<L.rows; i++)
      for (int j=i; j<L.rows; j++)
         {
         double sum = A[i][j];
         for (int k = i - 1; k >= 0; k--)
            sum -= L.data[i]->data[k] * L.data[j]->data[k];
         if (i == j)
            if (sum <= 0.0)
               return false;
            else
               L.data[i]->data[i] = sqrt(sum);
         else
            L.data[j]->data[i] = sum / L.data[i]->data[i];
         }

   return true;
   }
   
void Cholesky::BackSubst0(Vector & b)
   {
   x.Dimension(L.rows);

   // Solve L*v = b (store v in x)
   for (int i = 0; i < L.rows; i++)
      {
      double sum = b.data[i];
      for (int k = i-1; k>=0; k--)
         sum -= L.data[i]->data[k] * x.data[k];
      x.data[i] = sum / L.data[i]->data[i];
      }
   }

void Cholesky::BackSubst(Vector & b)
   {
   x.Dimension(L.rows);

   // Solve L*v = b (store v in x)
   for (int i = 0; i < L.rows; i++)
      {
      double sum = b[i];
      for (int k = i-1; k>=0; k--)
         sum -= L.data[i]->data[k] * x.data[k];
      x.data[i] = sum / L.data[i]->data[i];
      }

   // Solve transpose(L)*x = v
   // End result is ... A*x = L*t(L)*x = L*v = b
   for (int i=L.rows-1; i>=0; i--)
      {
      double sum = x[i];
      for (int k = i+1; k < L.rows; k++)
         sum -= L.data[k]->data[i] * x.data[k];
      x.data[i] = sum / L.data[i]->data[i];
      }

   // Done!
   }

void Cholesky::Invert()
   {
   inv.Dimension(L.rows, L.rows);

   inv.Identity();

   for(int i = 0; i < L.rows; i++)
      {
      BackSubst(inv[i]);
      inv[i] = x;
      }
   }

double Cholesky::lnDeterminantL()
   {
   double sum = 0;
   for (int i = 0; i < L.rows; i++)
      sum += log(L[i][i]);
   return sum;
   }

double Cholesky::DeterminantL()
   {
   double product = 1;
   for (int i=0; i<L.rows; i++)
      product *= L[i][i];
   return product;
   }

