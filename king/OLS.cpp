// Wei-Min Chen, 10/18/2006

#include <math.h>
#include "OLS.h"
#include "MathStats.h"
#include "MathSVD.h"

OLS_REGRESSION::OLS_REGRESSION()
{
   nuisanceCount = 0;
}

void OLS_REGRESSION::Print()
{
   printf("%10s%7s%10s%10s%7s", "covariate", "t", "pvalue", "effect", "h2");
   if(nuisanceCount) printf("%7s", "h2_adj");
   printf("\n");

   printf("%10s%7.2f%10.2G%10.2f\n", "Mu", t_statistic[0], pvalue[0], beta[0]);
   double sum = 0.0;
   for(int c = 0; c < covariateNames.Length(); c++){
      printf("%10s%7.2f%10.2G%10.2f%6.1f%%",
         (const char*)covariateNames[c], t_statistic[c+1],
         pvalue[c+1], beta[c+1], R2[c]*100);
      if(nuisanceCount){
         if(c < nuisanceCount)
            sum += R2[c];
         else
            printf("%6.1f%%", R2[c] / (1-sum) * 100);
      }
      printf("\n");
   }
}

void OLS_REGRESSION::Print(const char* title)
{
   printf("%s (N=%d, loglik=%.2f)\n", title, Y.Length(), loglik);
   int length = 0;
   for(; title[length]; length++);
   for(int i = 0; i < length + 25; i++) printf("=");
   printf("\n");
   Print();
}

void OLS_REGRESSION::run()
{
   failure = true;
   Cholesky chol;
   SVD svd;
   N = Y.Length();
   P = X.cols - 1;

   t_statistic.Dimension(P+1);
   pvalue.Dimension(P+1);
   R2.Dimension(P);
   R2.Zero();
   t_statistic.Zero();
   pvalue.Set(1.0);
   Vector x(P);
   tMatrix.Transpose(X);
   for(int i = 0; i < P; i++)
      x[i] = tMatrix[i+1].Sum() / N;

   svd.Decompose(X);
   svd.Edit();
   svd.BackSubst(Y);
   beta = svd.x;

   IntArray pheno(0);
   for(int i = 0; i < P; i++)
      if(svd.w[i+1] != 0.0)
         pheno.Push(i);
      else beta[i+1] = 0.0;

   L.Dimension(pheno.Length(), pheno.Length());
   L.Zero();
   for(int i = 0; i < pheno.Length(); i++)
      for(int j = 0; j < pheno.Length(); j++){
         for(int k = 0; k < N; k++)
            L[i][j] += X[k][pheno[i]+1] * X[k][pheno[j]+1];
         L[i][j] -= N * x[pheno[i]] * x[pheno[j]];
      }
   if(pheno.Length()==0){/*printf("Regression failed.\n");*/ return;}
   if(!chol.TryDecompose(L)) return;
   chol.Decompose(L);
   chol.Invert();
   Linverse = chol.inv;
   tVector.Dimension(N);
   Y.Add(-Y.Sum()/N);
   for(int i = 0; i < N; i++){
      tVector[i] = Y[i];
      for(int j = 0; j < P+1; j++)
         tVector[i] -= X[i][j] * beta[j];
   }
   Q = tVector.SumSquares();
   Cov.Dimension(pheno.Length()+1, pheno.Length()+1);
   Cov.Zero();
   Cov[0][0] = 1.0 / N;
   for(int i = 0; i < pheno.Length(); i++)
      for(int j = 0; j < pheno.Length(); j++)
         Cov[0][0] += x[pheno[i]]*Linverse[i][j]*x[pheno[j]];
   for(int i = 0; i < pheno.Length(); i++){
      for(int j = 0; j < pheno.Length(); j++)
         Cov[0][1+pheno[i]] -= x[pheno[j]]*Linverse[j][i];
      Cov[1+pheno[i]][0] = Cov[0][1+pheno[i]];
   }
   for(int i = 0; i < pheno.Length(); i++)
      for(int j = 0; j < pheno.Length(); j++)
         Cov[pheno[i]+1][pheno[j]+1] = Linverse[i][j];
   for(int i = 0; i < P+1; i++)
      for(int j = 0; j < P+1; j++)
         Cov[i][j] *= Q / (N-P-1.0);

   SE.Dimension(P+1);
   for(int i = 0; i < P+1; i++) SE[i] = sqrt(Cov[i][i]);
   t_statistic[0] = beta[0] / SE[0];
   pvalue[0] = tdist(fabs(t_statistic[0]), N-P-1);
   for(int i = 0; i < pheno.Length(); i++){
      t_statistic[1+pheno[i]] = beta[pheno[i]+1] / SE[pheno[i]+1];
      pvalue[1+pheno[i]] = tdist(fabs(t_statistic[1+pheno[i]]), N-P-1);
   }
   loglik = -0.5 * N * log(Q/N);

   // R2 = beta * Sxy /Syy;
   tMatrix.Transpose(X);
   for(int i = 0; i < P; i++)
      tMatrix[i+1] -= x[i];
   for(int i = 0; i < P; i++){
     R2[i] = tMatrix[i+1].InnerProduct(Y) * beta[i+1] / Y.SumSquares();
   }
   failure = false;
}

void OLS_REGRESSION::run(Vector y, Vector x)
{
   N = y.Length();
   Y = y;
   X.Dimension(N, 2);
   for(int i = 0; i < N; i++){
      X[i][0] = 1;
      X[i][1] = x[i];
   }
   run();
}


void OLS_REGRESSION::run(Vector y, Matrix x)
{
   Y = y;
   MatrixOne.Dimension(y.Length(),1);
   MatrixOne.Set(1.0);
   X = MatrixOne;
   X.StackLeft(x);
   run();
}






//   R2_alt.Dimension(P);
//   R2_alt.Zero();
//     R2_alt[i] = beta[i+1] * beta[i+1] * tMatrix[i+1].SumSquares() / Y.SumSquares();


/*
   tVector.Dimension(P);
   tMatrix.Transpose(X);
   for(int i = 0; i < P; i++)
      tVector[i] = tMatrix[i+1].Sum() / N;
   for(int i = 0; i < N; i++)
      for(int j = 0; j < P; j++)
         X[i][j+1] -= tVector[j];
*/


   /*
   svd.Covariances();
   for(int i = 0; i < P; i++){
      if( (svd.w[i+1]!=0.0) && svd.cov[i+1][i+1]>0){
         t_statistic[i] = beta[i+1] / sqrt(svd.cov[i+1][i+1]);
         pvalue[i] = tdist(fabs(t_statistic[i]), N-P-1);
      }
   }
   */
/*   tMatrix.Transpose(X);
   tVector.Dimension(pheno.Length());
   for(int i = 0; i < pheno.Length(); i++)
      tVector[i] = tMatrix[pheno[i]+1].Sum() / N;
*/
/*   for(int i = 0; i < P; i++){
      tVector = tMatrix[i+1];
      tVector.Add(-tVector.Sum()/N);
      R2[i] = tVector.InnerProduct(Y) * beta[i+1] / Y.SumSquares();
   }
   */


/*
   chol.Decompose(XX);
   chol.BackSubst(XY);
   beta = chol.x;
   tVector.Dimension(P);
   for(int i = 0; i < P; i++)
      tVector[i] = tMatrix[i+1].Sum() / N;
   L.Dimension(P, P);
   L.Zero();
   for(int i = 0; i < P; i++)
      for(int j = 0; j < P; j++)
         for(int k = 0; k < N; k++)
            L[i][j] += (X[k][i+1]-tVector[i]) * (X[k][j+1]-tVector[j]);
   chol.Decompose(L);
   chol.Invert();
   Linverse = chol.inv;

   tVector.Dimension(N);
   for(int i = 0; i < N; i++){
      tVector[i] = Y[i];
      for(int j = 0; j < P+1; j++)
         tVector[i] -= X[i][j] * beta[j];
   }
   Q = tVector.SumSquares();

   t_statistic.Dimension(P);
   pvalue.Dimension(P);
   R2.Dimension(P);
   for(int i = 0; i < P; i++){
      t_statistic[i] = beta[i+1] / sqrt(Linverse[i][i] * Q / (N-P-1.0));
      pvalue[i] = tdist(fabs(t_statistic[i]), N-P-1);
   }

   loglik = -0.5 * N * log(Q/N);    

   // R2 = beta * Sxy /Syy;
   Y.Add(-Y.Sum()/N);
   tMatrix.Transpose(X);
   for(int i = 0; i < P; i++){
      tVector = tMatrix[i+1];
      tVector.Add(-tVector.Sum()/N);
      R2[i] = tVector.InnerProduct(Y) * beta[i+1] / Y.SumSquares();
   }
  */
/*   tMatrix.Transpose(X);
   XX.Product(tMatrix, X);
   XY.Dimension(N);
   XY.Zero();
   for(int i = 0; i < P+1; i++)
      for(int j = 0; j < N; j++)
         XY[i] += X[j][i] * Y[j];
*/

