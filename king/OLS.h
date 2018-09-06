#ifndef __OLS_h__
#define __OLS_h__

#include "Pedigree.h"
#include "IntArray.h"
#include "MathMatrix.h"
#include "MathVector.h"
#include "MathCholesky.h"

class OLS_REGRESSION{
      Matrix MatrixOne;
      Matrix L, Linverse;
      double Q;
      Matrix tMatrix;
      Vector tVector;
   public:
   // Input
      Vector Y;
      Matrix X;

   // Output
      int N;         // sample size
      int P;         // # covariates
      int testCount;
      StringArray covariateNames;
      int nuisanceCount;
      Vector beta;   // regression coefficient
      Vector SE;
      Matrix Cov;
      double loglik; // log likelihood
      Vector t_statistic;
      Vector pvalue;
      Vector R2;     // r-square: a Cov(X, Y) / Var(Y)
//      Vector R2_alt; // r-square: a^2 Var(X) / Var(Y)
      bool failure;

      OLS_REGRESSION();
      void run();
      void run(Vector y, Matrix X);
      void run(Vector y, Vector X);
      void Print();
      void Print(const char* title);
};

#endif


