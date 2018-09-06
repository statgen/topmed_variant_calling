//////////////////////////////////////////////////////////////////////
// VCPOLY.cpp
// Author: Wei-Min Chen
// March 16, 2005

#include "diseaseGEE.h"
#include "Pedigree.h"
#include "Kinship.h"
#include "MathStats.h"
#include "QuickIndex.h"
#include <math.h>

GEE_DIS::GEE_DIS(Pedigree & pedigree):GEE(pedigree)
{
   LoopCount = 20; Epsilon = 0.0001;
   diseases = NULL;
   resid = NULL;
}

GEE_DIS::~GEE_DIS()
{
   if(diseases) delete[]diseases;
   if(resid) delete[]resid;
}

void GEE_DIS::RefreshD(int f)
{
   D.Dimension(size, parCount);
   D.Zero();
   for(int i = 0; i < size; i++){
      double t = resid[f][i] * (1-resid[f][i]);
      D[i][0] = t;
      for(int j = 1; j < parCount; j++)
         D[i][j] = covariates[f][j][i] * t;
   }
}

void GEE_DIS::InitCoef()
{
   if(resid==NULL) resid = new Vector[ped.familyCount];
   if(pheno==NULL) pheno = new IntArray[ped.familyCount];
   if(diseases==NULL) diseases = new IntArray[ped.familyCount];
   if(covariates==NULL) covariates = new Matrix[ped.familyCount];
   for(int f = 0; f < ped.familyCount; f++) {
      pheno[f].Dimension(0);
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
         int missing = 0;
         for(int k = 0; k < mCovariate.Length(); k++){
            if(!ped[i].isControlled(mCovariate[k])) {
               missing=1;
               break;
            }
         }
         if(!missing && ped[i].isDiagnosed(disease))
            pheno[f].Push(i);
      }
   }

   coef.Dimension(mCovariate.Length()+1);
   coef.Zero();
   parCount = coef.Length();
   ValidFamilies = ValidPersons = 0;
   for(int f = 0; f < ped.familyCount; f++){
      int size = pheno[f].Length();
      diseases[f].Dimension(size);
      covariates[f].Dimension(parCount, size);
      if(size == 0) continue;
      ValidFamilies++;
      ValidPersons += size;
      for(int i = 0; i < size; i++){
         diseases[f][i] = ped[pheno[f][i]].affections[disease]-1;
         covariates[f][0][i] = 1;
         for(int j = 1; j < parCount; j++)
            covariates[f][j][i] = ped[pheno[f][i]].covariates[mCovariate[j-1]];
      }
   }
   int sampleSize = 0;
   int affectedSize = 0;
   for (int f = 0; f < ped.familyCount; f++)
      for(int i = 0; i < pheno[f].Length(); i++)
         if(ped[pheno[f][i]].affections[disease]){
            sampleSize ++;
            affectedSize += (ped[pheno[f][i]].affections[disease]==2);
         }
   double K = (double)affectedSize / sampleSize;
   if(K <= 0.0) error("log(0)");
   coef[0] = log(K / (1 - K));
   for (int f = 0; f < ped.familyCount; f++){
      resid[f].Dimension(pheno[f].Length());
      for(int i = 0; i < pheno[f].Length(); i++)
         resid[f][i] = K;
   }
}

void GEE_DIS::solve()
{
   InitCoef();
   Matrix DVD(parCount, parCount);
   Vector DVS(parCount);
   Matrix DVS_FAM(parCount, ped.familyCount);
   for(int loop = 0; loop < LoopCount; loop++){
//      printf("Iteration %d...\n", loop+1);
      DVD.Zero();
      DVS_FAM.Zero();
      for(int f = 0; f < ped.familyCount; f++){
         size = pheno[f].Length();
         if(size==0) continue;
         RefreshD(f);
         for(int i = 0; i < parCount; i++)
            for(int j = 0; j < parCount; j++)
               for(int u = 0; u < size; u++)
                  DVD[i][j] += D[u][i] * D[u][j];
         for(int i = 0; i < parCount; i++)
            for(int u = 0; u < size; u++)
               DVS_FAM[i][f] += D[u][i] * (diseases[f][u] - resid[f][u]);
      }
      for(int i = 0; i < parCount; i++)
         DVS[i] = DVS_FAM[i].Sum();
      covariance.CholeskyInvert(DVD);
      Vector delta(parCount);
      delta.Zero();
      for(int i = 0; i < parCount; i++)
         for(int j = 0; j < parCount; j++)
            delta[i] += covariance[i][j] * DVS[j];
      coef.Add(delta);
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = 0; i < pheno[f].Length(); i++){
            double t = coef[0];
            for(int j = 1; j < parCount; j++)
               t += coef[j] * covariates[f][j][i];
            if(t>20) resid[f][i] = 1;
            else resid[f][i] = exp(t)/(1+exp(t));
         }
      if(delta.SumSquares() < Epsilon) break;
   }
   covariance_R.Dimension(parCount, parCount);
   covariance_R.Zero();
   Matrix DVD2(parCount, parCount);
   DVD2.Zero();
   for(int i = 0; i < parCount; i++)
      for(int j = 0; j < parCount; j++)
         for(int f = 0; f < ped.familyCount; f++)
            DVD2[i][j] += DVS_FAM[i][f] * DVS_FAM[j][f];
   for(int i = 0; i < parCount; i++)
      for(int j = i; j < parCount; j++){
         for(int u = 0; u < parCount; u++)
            for(int v = 0; v < parCount; v++)
               covariance_R[i][j] += covariance[i][u] * DVD2[u][v] * covariance[v][j];
         covariance_R[j][i] = covariance_R[i][j];
      }
   SEcoef_R.Dimension(parCount);
   for(int i = 0; i < parCount; i++)
      if(covariance_R[i][i]>0) SEcoef_R[i] = sqrt(covariance_R[i][i]);
      else continue;
}


void GEE_DIS::print()
{
   printf("\nFitted Polygenic Models (GEE)\n");
   for(int i = 0; i < 48 + 9*mCovariate.Length(); i++) printf("=");
   printf("\n");
   printf("%15s%6s", "Trait", "N");
   printf("%11s", "Mean");
   for(int i = 0; i < mCovariate.Length(); i++)
      printf("%11s", (const char*)ped.covariateNames[mCovariate[i]]);
   printf("\n");
   printf("%15s", (const char*)ped.affectionNames[disease]);
   printf("%6d", ValidPersons);
   for(int i = 0; i < parCount; i++)
      printf("%11.3f", coef[i]);
   printf("\n");
   printf("%15s%6s", "(se)","");
   for(int i = 0; i < parCount; i++)
      printf("%11.3f", SEcoef_R[i]);
   printf("\n\n");
//   coef.Print();
//   SEcoef_R.Print();
//   covariance_R.Print(stdout);
}






