// analysis.cpp
// 9/26/2006 Wei-Min Chen

#include "analysis.h"
#include "MathGenMin.h"
#include "MathStats.h"
#include "MapFunction.h"
#include "QuickIndex.h"
#include "Kinship.h"
#include "diseaseGEE.h"
#include "OLS.h"
#include "VCLinear.h"

void Engine::polygenic()
{
   means.Dimension(traits.Length(), covariates.Length()+1);
   variances.Dimension(traits.Length());
   heritabilities.Dimension(traits.Length());
   SampleSize.Dimension(traits.Length());

   if(diseases.Length()>0){
      GEE_DIS gee(ped);
      gee.mCovariate = covariates;
      for(int i = 0; i < diseases.Length(); i++){
         gee.disease = diseases[i];
         gee.solve();
         gee.print();
      }
      return;
   }
   if(unrelatedFlag){
      heritabilities.Zero();
      Vector Y;
      Matrix X(ped.count, 1+covariates.Length());
      Vector XY(1+covariates.Length());
      Matrix XX(1+covariates.Length(),1+covariates.Length());
      Cholesky chol;
      means.Zero();

      for(int i = 0; i < ped.count; i++)
         X[i][0] = 1.0;
      for(int t = 0; t < traits.Length(); t++){
         Y.Dimension(0);
         for(int i = 0; i < ped.count; i++)
            if(ped[i].isPhenotyped(traits[t]) && CheckCovariates(ped[i])){
               for(int j = 0; j < covariates.Length(); j++)
                  X[Y.Length()][j+1] = ped[i].covariates[covariates[j]];
               Y.Push(ped[i].traits[traits[t]]);
             }
         if(covariates.Length()==0){
            variances[t] = Y.Var();
            means[t][0] = Y.Sum() / Y.Length();
         }else{
            XX.Zero();
            XY.Zero();
            for(int i = 0; i < 1+covariates.Length(); i++)
               for(int j = 0; j < Y.Length(); j++)
                  XY[i] += X[j][i] * Y[j];
            for(int i = 0; i < 1+covariates.Length(); i++)
               for(int j = 0; j < 1+covariates.Length(); j++)
                  for(int k = 0; k < Y.Length(); k++)
                     XX[i][j] += X[k][i] * X[k][j];
            chol.Decompose(XX);
            chol.Invert();
            for(int i = 0; i < covariates.Length()+1; i++)
               for(int j = 0; j < covariates.Length()+1; j++)
                  means[t][i] += chol.inv[i][j] * XY[j];
            for(int i = 0; i < Y.Length(); i++)
               for(int j = 0; j < covariates.Length()+1; j++)
                  Y[i] -= means[t][j] * X[i][j];
            variances[t] = Y.Var();
         }
      }
      return;
   }
   
   POLY *GEEpoly = new POLY(ped);
   GEEpoly->mCovariate = covariates;

   for(int i = 0; i < traits.Length(); i++){
      GEEpoly->trait = traits[i];
      GEEpoly->solve();
      if(GEEpoly->AtBorder){
         double oldLoglik = GEEpoly->loglik;
         double oldScale = GEEpoly->deltaScale;
         GEEpoly->deltaScale *= 0.5;
         GEEpoly->solve();
         if(GEEpoly->loglik < oldLoglik - 0.001){
            GEEpoly->deltaScale = oldScale;
            GEEpoly->AtBorder = 0;
            GEEpoly->solve();
         }
         GEEpoly->deltaScale = oldScale;
         GEEpoly->AtBorder = 0;
      }
      if(GEEpoly->coef.Min() < -1e100 || GEEpoly->coef.Max() > 1E100){
         means[i] = GEEpoly->coef;
         means[i].Zero();
      }else
         means[i] = GEEpoly->coef;
      if(GEEpoly->totalVariance > 1e100 || GEEpoly->totalVariance < -1e100)
         variances[i] = 0;
      else
         variances[i] = GEEpoly->totalVariance;
      if(GEEpoly->H2 > 1)
         heritabilities[i] = 1;
      else if(GEEpoly->H2 < 0)
         heritabilities[i] = 0;
      else
         heritabilities[i] = GEEpoly->H2;
      SampleSize[i] = GEEpoly->ValidPersons;
   }
   PrintPolygenic();
}

bool Engine::CheckCovariates(Person & p)
{
   for (int i = 0; i < covariates.Length(); i++)
      if (p.covariates[covariates[i]] == _NAN_)
         return false;
   return true;
}

void Engine::PrintPolygenic()
{
   printf("\nFitted Polygenic Models\n");
   for(int i = 0; i < 48 + 9*covariates.Length(); i++) printf("=");
   printf("\n");
   printf("%15s%9s%9s%6s", "Trait", "Herit", "Variance", "N");
   if(covariates.Length()) printf("%9s", "Mu");
   for(int i = 0; i < covariates.Length(); i++)
      printf("%9s", (const char*)ped.covariateNames[covariates[i]]);
   printf("\n");
   for(int trait = 0; trait < traits.Length(); trait++){
      printf("%15s%8.1f%%%9.3f", (const char*)ped.traitNames[traits[trait]],
         heritabilities[trait]*100, variances[trait]);
      printf("%6d", SampleSize[trait]);
      for(int i = 0; i < means[trait].Length(); i++)
         printf("%9.3f", means[trait][i]);
      printf("\n");
   }
}

//   POLYloglik.Dimension(traits.Length());
//   if(moreFlag) GEEpoly->moreFlag = 1;


