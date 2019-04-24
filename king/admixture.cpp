//////////////////////////////////////////////////////////////////////
// admixture.cpp
// (c) 2010-2019 Wei-Min Chen
//
// This file is distributed as part of the KING source code package
// and may not be redistributed in any form, without prior written
// permission from the author. Permission is granted for you to
// modify this file for your own personal use, but modified versions
// must retain this copyright notice and must not be distributed.
//
// Permission is granted for you to use this file to compile KING.
//
// All computer programs have bugs. Use this file at your own risk.
//
// Feb 22, 2019

#include "Random.h"
#include "analysis.h"
#include "MathStats.h"
#include "Kinship.h"
#include "QuickIndex.h"
#include "Davies.h"
#include "MathSVD.h"
#include "QuickIndex.h"
#include "MathCholesky.h"
#include <math.h>
#include "diseaseGEE.h"
#include "VCLinear.h"

#ifdef WITH_LAPACK
extern "C" void dgesdd_(char*, int*, int*, double*, int*, double*, double *, int*, double*, int*, double*, int*, int*, int*);
extern "C" void dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);
#endif

void Engine::AdmixtureMapping(double windowSize, double moveSize)
{
   if(geno.Length()==0) BuildBinary();
   bool isDisease = false;
   if(traitList.Length()==0 && ped.affectionCount==0) {printf("No phenotype specified\n"); return;}
   printf("Admixture mapping starts at %s", currentTime());

   unrelatedFlag = false;
   int relatedCount = 0;
   Kinship kin;
   for(int f = 0; f < ped.familyCount; f++){
      kin.Setup(*ped.families[f]);
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         for(int j = i+1; j <= ped.families[f]->last; j++)
            if(CheckCovariates(ped[i]) && geno[i]!=-1 &&
               CheckCovariates(ped[j]) && geno[j]!=-1 &&
               kin(ped[i], ped[j])>0)
               relatedCount ++;
   }
   if(relatedCount < 5) {
      unrelatedFlag = true;
      if(relatedCount)
         printf("Samples are treated as unrelated. It's recommanded to remove the %d relative pairs.\n", relatedCount);
      else
         printf("All samples are unrelated\n");
   }else
      printf("There are %d relative pairs in the known pedigrees.\n", relatedCount);

   int trait = 0;
   int nid = 0;
   if(traitList.Length()==0) {   // Dichotomous trait
      isDisease = true;
      for(int i = 0; i < ped.count; i++)
         if(ped[i].affections[0]!=0 && CheckCovariates(ped[i]) && geno[i]!=-1)
            nid++;
      printf("The dichotomous trait is %s (N=%d)\n", (const char*)ped.affectionNames[0], nid);
   }else{   // quantitative trait
      for(int i = 0; i < ped.count; i++)
         if(ped[i].traits[traits[trait]]!=_NAN_ && CheckCovariates(ped[i]) && geno[i]!=-1)
            nid++;
      printf("The quantitative trait is %s (N=%d)\n", (const char*)ped.traitNames[traits[trait]], nid);
   }

   int b;
   int *vid = new int[nid];
   Vector adjResid(nid);
   adjResid.Zero();
   Matrix Omega;
   IntArray ids;

   if(isDisease)
      ResidualDisease(vid, adjResid);
   else
      ResidualQuantitative(0, vid, adjResid);

   double pvalue, Q, pvalue2, Q2, mean0;
   Matrix kernel;
   kernel.Dimension(nid, nid);
   Vector mean(nid);
   int HetHetCount, IBS0Count, het1Count, het2Count, notMissingCount, HomHomCount;
   int id1, id2;
   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++){
      oneoneCount[i] = 0;
      for(int j = 0; j < 16; j++)
         if(i & shortbase[j]) oneoneCount[i]++;
   }

   if(windowSize > 9999){ // test of global ancestry
      printf("Calculating pair-wise distances among %d individuals...", nid);
      kernel.Zero();
      for(int i = 0; i < nid; i++){
      id1 = geno[vid[i]];
      for(int j = i+1; j < nid; j++){
         id2 = geno[vid[j]];
         HetHetCount = IBS0Count = het1Count = het2Count = HomHomCount = 0;
         for(int m = 0; m < shortCount; m++){
            HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
            IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
            HomHomCount += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & ~(GG[1][id1][m]^GG[1][id2][m])];
            het1Count += oneoneCount[(GG[0][id2][m] | GG[1][id2][m]) & (~GG[0][id1][m]) & GG[1][id1][m]];
            het2Count += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
         }
         notMissingCount = HomHomCount + IBS0Count + het1Count + het2Count - HetHetCount;
         kernel[i][j] = kernel[j][i] = -(het1Count + het2Count - 2*HetHetCount + 4*IBS0Count) * 1.0 / notMissingCount;
         }
      }
      printf("Done\n");
      mean.Zero();
      for(int i = 0; i < nid; i++){
         for(int j = 0; j < nid; j++)
            mean[i] += kernel[i][j];
         mean[i] /= nid;
      }
      mean0 = mean.Sum() / nid;
      for(int i = 0; i < nid; i++)
         for(int j = 0; j < nid; j++)
            kernel[i][j] -= (mean[i] + mean[j] - mean0);
      Q = 0.0;
      for(int i = 0; i < nid; i++)
         for(int j = 0; j < nid; j++)
            Q += adjResid[i] * kernel[i][j] * adjResid[j];
      printf("SVD starts at %s", currentTime());
      pvalue = ComputeDaviesP(Q, kernel);
      if(pvalue<0) pvalue=0;
      printf("Association between the global ancestry and the trait is Q = %.1lf with P value = %.1G\n",
         Q, pvalue);
      printf("\nGlobal ancestry association analysis ends at %s", currentTime());
      return;
   }

   printf("Local ancestry inference and test starts at %s", currentTime());
   double firstPos, lastPos, length;
   IntArray vbyte, vbit;
   unsigned short int G1[2], G2[2];
   int SNPCount;
   double phi;
   double D[3];
   int nD[3];

   if(unrelatedFlag){
      if(isDisease)
         printf("\n%5s%8s%8s%8s%7s%7s%7s%8s%8s\n",
            "Chr", "Start", "Stop", "N_SNP", "D_UU", "D_AU", "D_AA", "Q", "P");
      else
         printf("\n%5s%8s%8s%8s%8s%8s\n",
            "Chr", "Start", "Stop", "N_SNP", "Q", "P");
   }else{
      if(isDisease)
         printf("\n%5s%8s%8s%8s%7s%7s%7s%8s%8s%8s%8s\n",
            "Chr", "Start", "Stop", "N_SNP", "D_UU", "D_AU", "D_AA", "Q_unadj", "P_unadj", "Q_adj", "P_adj");
      else
         printf("\n%5s%8s%8s%8s%8s%8s%8s%8s\n",
            "Chr", "Start", "Stop", "N_SNP", "Q_unadj", "P_unadj", "Q_adj", "P_adj");
   }
   for(int chr=0; chr<SEXCHR; chr++){   // by chromosome
      firstPos=1000000000.0; lastPos=0;
      for(int m = 0; m < markerCount; m++){
         if(chromosomes[m]!=chr) continue;
         if(bp[m] < firstPos) firstPos = bp[m];
         if(bp[m] > lastPos) lastPos = bp[m];
      }
      length = lastPos - firstPos;
      int windowCount =
         length < windowSize? 1: int((length - windowSize) / moveSize) + 2;

      for(int i = 0; i < windowCount; i++){  // by window
         double start = firstPos + i*moveSize;
         double stop = firstPos + i*moveSize + windowSize;
         vbyte.Dimension(0);
         vbit.Dimension(0);
         SNPCount = 0;
         for(int m = 0; m < markerCount; m++){
            if(chromosomes[m]!= chr) continue;
            if(bp[m] < start || bp[m] >= stop) continue;
            SNPCount ++;
            int proposedByte = m/16;
            int proposedBit = m%16;
            if(vbyte.Find(proposedByte)>-1){ // in previous byte
               vbit[vbyte.Find(proposedByte)] |= (1<<proposedBit);
            }else{   // a new byte
               vbyte.Push(proposedByte);
               vbit.Push(1<<proposedBit);
            }
         }

         if(SNPCount < 1000) continue; // local ancestry inference must have at least 500 SNPs
         kernel.Zero();
         int id1, id2;
         for(int i = 0; i < nid; i++){
            id1 = geno[vid[i]];
            for(int j = i+1; j < nid; j++){
               id2 = geno[vid[j]];
               HetHetCount = IBS0Count = het1Count = het2Count = HomHomCount = 0;
               for(int v = 0; v < vbyte.Length(); v++){
                  b = vbyte[v];
                  G1[0] = GG[0][id1][b] & vbit[v];
                  G1[1] = GG[1][id1][b] & vbit[v];
                  G2[0] = GG[0][id2][b] & vbit[v];
                  G2[1] = GG[1][id2][b] & vbit[v];
                  HetHetCount += oneoneCount[(~G1[0]) & G1[1] & (~G2[0]) & G2[1]];
                  IBS0Count += oneoneCount[G1[0] & G2[0] & (G1[1] ^ G2[1])];
                  HomHomCount += oneoneCount[G1[0] & G2[0] & ~(G1[1]^G2[1])];
                  het1Count += oneoneCount[(G2[0] | G2[1]) & (~G1[0]) & G1[1]];
                  het2Count += oneoneCount[(G1[0] | G1[1]) & (~G2[0]) & G2[1]];
               }
               notMissingCount = HomHomCount + IBS0Count + het1Count + het2Count - HetHetCount;
               kernel.data[i]->data[j] = kernel.data[j]->data[i] =
                  -(het1Count + het2Count - 2*HetHetCount + 4*IBS0Count) * 1.0 / notMissingCount;
            }
         }

         D[0]=D[1]=D[2]=0.0;
         nD[0]=nD[1]=nD[2]=0;
         if(isDisease){
            for(int i = 0; i < nid; i++)
               for(int j = i+1; j < nid; j++){
                  if(adjResid[i] < 0 && adjResid[j] < 0){ // UU
                     nD[0]++;
                     D[0] -= kernel[i][j];
                  }else if(adjResid[i] > 0 && adjResid[j] > 0){ // DD
                     nD[2]++;
                     D[2] -= kernel[i][j];
                  }else{ // DU
                     nD[1]++;
                     D[1] -= kernel[i][j];
                  }
               }
            for(int i = 0; i < 3; i++) D[i] /= nD[i];
         }

         mean.Zero();
         for(int i = 0; i < nid; i++){
            for(int j = 0; j < nid; j++)
               mean[i] += kernel[i][j];
            mean[i] /= nid;
         }
         mean0 = mean.Sum() / nid;
         for(int i = 0; i < nid; i++)
            for(int j = 0; j < nid; j++)
               kernel[i][j] -= (mean[i] + mean[j] - mean0);
         Q = 0.0;
         for(int i = 0; i < nid; i++)
            for(int j = 0; j < nid; j++)
               Q += adjResid[i] * kernel[i][j] * adjResid[j];

         pvalue = Q>0? ComputeDaviesP(Q, kernel): 1;
         if(pvalue<0) pvalue=0;

         for(int i = 0; i < nid; i++)
            for(int j = 0; j < nid; j++)
               kernel[i][j] += (mean[i] + mean[j] - mean0);

         if(!unrelatedFlag){
            int baseN=0;
            for(int f = 0; f < ped.familyCount; f++){
               ids.Dimension(0);
               for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
               if(!isDisease){
                  if(ped[i].traits[traits[trait]]!=_NAN_ && CheckCovariates(ped[i]) && geno[i]!=-1)
                     ids.Push(i);
               }else{
                  if(ped[i].affections[0]!=0 && CheckCovariates(ped[i]) && geno[i]!=-1)
                     ids.Push(i);
               }
               int N = ids.Length();
               if(N==0) continue;
               else if(N==1){
                  baseN++;
                  continue;
               }
               kin.Setup(*ped.families[f]);
               for(int i = 0; i < N; i++)
               for(int j = i+1; j < N; j++){
                  phi = kin(ped[ids[i]],ped[ids[j]]);
                  if(phi < 0.499)
                     kernel[baseN+i][baseN+j] = kernel[baseN+j][baseN+i]
                        = kernel[baseN+i][baseN+j] / (1-2*phi);
               }
               baseN += N;
            }
            mean.Zero();
            for(int i = 0; i < nid; i++){
               for(int j = 0; j < nid; j++)
                  mean[i] += kernel[i][j];
               mean[i] /= nid;
            }
            mean0 = mean.Sum() / nid;
            for(int i = 0; i < nid; i++)
               for(int j = 0; j < nid; j++)
               kernel[i][j] -= (mean[i] + mean[j] - mean0);
            Q2 = 0.0;
            for(int i = 0; i < nid; i++)
               for(int j = 0; j < nid; j++)
               Q2 += adjResid[i] * kernel[i][j] * adjResid[j];

            pvalue2 = Q2>0? ComputeDaviesP(Q2, kernel): 1;
            if(pvalue2<0) pvalue2=0;
            if(isDisease)
               printf("%5d%8.3lf%8.3lf%8d%7.4lf%7.4lf%7.4lf%8.1lf%8.1G%8.1f%8.1G\n",
                  chr, start, stop, SNPCount, D[0], D[1], D[2], Q, pvalue, Q2, pvalue2);
            else
               printf("%5d%8.3lf%8.3lf%8d%8.1lf%8.1G%8.1f%8.1G\n",
                  chr, start, stop, SNPCount, Q, pvalue, Q2, pvalue2);
         }else{ // unrelated
            if(isDisease)
               printf("%5d%8.3lf%8.3lf%8d%7.4lf%7.4lf%7.4lf%8.1lf%8.1G\n",
                  chr, start, stop, SNPCount, D[0], D[1], D[2], Q, pvalue);
            else
               printf("%5d%8.3lf%8.3lf%8d%8.1lf%8.1G\n",
                  chr, start, stop, SNPCount, Q, pvalue);
         }
      }  // end of window loop
  }   // end of chr loop
   delete []vid;
   printf("\nAdmixture mapping ends at %s", currentTime());
}





double Engine::ComputeDaviesP(double Q, Matrix & kernel)
{
   Vector lambda;
   lambda.Dimension(0);
   int nid = kernel[0].Length();
#ifdef WITH_LAPACK
   double *A, *S, *U, *VT, *WORK;
   int *IWORK;
   int LWORK = 8*nid + 4*nid*nid;
   char JOBZ = 'A';
   int dimLA = nid;
   A = new double[nid*nid];
   S = new double[nid];
   U = new double[nid*nid];
   VT = new double[nid* nid];
   IWORK = new int[nid*8];
   WORK = new double[LWORK];
   int info;
   for(int i = 0; i < nid; i++)
      for(int j = 0; j < nid; j++)
         A[i*nid+j] = kernel[j][i];
   dgesdd_(&JOBZ, &dimLA, &dimLA, A, &dimLA, S, U, &dimLA, VT, &dimLA, WORK, &LWORK, IWORK, &info);
   for(int i = 0; i < nid; i++)
      if(S[i] > 1E-100) lambda.Push(S[i]);

   delete U;
   delete IWORK;
   delete WORK;
   delete A;
   delete S;
   delete VT;
#else
   SVD svd;
   svd.Decompose(kernel);
   if(svd.n == 0) return _NAN_;
   for(int i = 0; i < nid; i++)
      if(svd.w[i] > 1E-200) lambda.Push(svd.w[i]);
#endif

   double pvalue;
   pvalue = lambda.Length()? 1-Davies(Q, lambda.data, lambda.Length()): 1;
   if(pvalue > 1) pvalue = 1;
//      lambda.Print(20);
   return pvalue;
}

void Engine::ResidualQuantitative(int trait, int * vid, Vector & adjResid)
{
   int nid = adjResid.Length();
   Vector qtrait(nid);
   Matrix qcov(nid, covariates.Length()+1);
   qtrait.Zero();
   qcov.Zero();
   int baseN = 0;
   Kinship kin;
   Matrix Omega;
   IntArray ids;
   Cholesky chol;

   int count = 0;
   Vector resid(nid);
   diseases.Dimension(0);

   polygenic();
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(ped[i].traits[traits[trait]]!=_NAN_ && CheckCovariates(ped[i]) && geno[i]!=-1){
         vid[count] = i;
         qtrait[count] = ped[i].traits[traits[trait]];
         qcov[count][0] = 1;
         for(int j = 0; j < covariates.Length(); j++)
            qcov[count][j+1] = ped[i].covariates[covariates[j]];
         count++;
         }
      for(int i = 0; i < nid; i++){
         resid[i] = qtrait[i] - means[trait][0];
         for(int j = 0; j < covariates.Length(); j++)
            resid[i] -= means[trait][j+1] * qcov[i][j+1];
      }
      if(!unrelatedFlag)
      for(int f = 0; f < ped.familyCount; f++){
         ids.Dimension(0);
         for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(ped[i].traits[traits[trait]]!=_NAN_ && CheckCovariates(ped[i]) && geno[i]!=-1)
            ids.Push(i);
         int N = ids.Length();
         if(N==0) continue;
         kin.Setup(*ped.families[f]);
         Omega.Dimension(N, N);
         for(int i = 0; i < N; i++){
            Omega[i][i] = variances[trait];
            for(int j = i+1; j < N; j++)
            Omega[i][j] = Omega[j][i] = 2*kin(ped[ids[i]], ped[ids[j]])*heritabilities[trait]*variances[trait];
         }
         chol.Decompose(Omega);
         Vector temp(N);
         for(int i = 0; i < N; i++)
            temp[i] = resid[baseN+i];
         chol.BackSubst0(temp);
         for(int i = 0; i < N; i++)
            adjResid[baseN+i] = chol.x[i];
         baseN += N;
      }else{   // unrelated
         for(int i = 0; i < nid; i++)
            adjResid[i] = resid[i] / sqrt(variances[trait]);
      }
}

void Engine::ResidualDisease(int * vid, Vector & adjResid)
{
   int nid = adjResid.Length();
   Vector qtrait(nid);
   Matrix qcov(nid, covariates.Length()+1);
   qtrait.Zero();
   qcov.Zero();
   Kinship kin;
   Matrix Omega;
   IntArray ids;
   Cholesky chol;
   int baseN = 0;
   int count = 0;
   
      Vector resid(nid);
      GEE_DIS gee(ped);
      gee.mCovariate = covariates;
      gee.disease = 0;
      gee.solve();
      gee.print();
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
            if(ped[i].affections[0]!=0 && CheckCovariates(ped[i]) && geno[i]!=-1){
               vid[count] = i;
               qtrait[count] = ped[i].affections[0]-1;
               qcov[count][0] = 1;
               for(int j = 0; j < covariates.Length(); j++)
                  qcov[count][j+1] = ped[i].covariates[covariates[j]];
               count++;
            }
      for(int i = 0; i < nid; i++){
         double coef = gee.coef[0];
         for(int c = 0; c < covariates.Length(); c++)
            coef += qcov[i][c+1] * gee.coef[c+1];
         coef = exp(coef);
         resid[i] = qtrait[i] - coef/(1+coef);
      }

      if(!unrelatedFlag)  // family data: normalize within each family
         for(int f = 0; f < ped.familyCount; f++){
            ids.Dimension(0);
            for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
               if(ped[i].affections[0]!=0  && CheckCovariates(ped[i]) && geno[i]!=-1)
                  ids.Push(i);
            int N = ids.Length();
            if(N==0) continue;
            double var = 0;
            for(int i = 0; i < N; i++)
               var += resid[baseN+i] * resid[baseN+i];
            var /= N;
            kin.Setup(*ped.families[f]);
            Omega.Dimension(N, N);
            for(int i = 0; i < N; i++){
               Omega[i][i] = 2*kin(ped[ids[i]], ped[ids[i]]);
               for(int j = i+1; j < N; j++)
               Omega[i][j] = Omega[j][i] = 2*kin(ped[ids[i]], ped[ids[j]]);
            }
            chol.Decompose(Omega);
            Vector temp(N);
            for(int i = 0; i < N; i++)
               temp[i] = resid[baseN+i]/sqrt(var);
            chol.BackSubst0(temp);
            for(int i = 0; i < N; i++)
               adjResid[baseN+i] = chol.x[i];
            baseN += N;
         }
      else{ // unrelated
         double sd = sqrt(resid.SumSquares()/nid);
         if(sd < 1E-200) error("There is no variation in the disease phenotype.");
         for(int i = 0; i < nid; i++)
            adjResid[i] = resid[i] / sd;
      }
}


//adjResid[i] = resid[i] * (1+coef) / sqrt(coef);    // var = coef / (1+coef)^2
  /*
      if(!unrelatedFlag)
      for(int f = 0; f < ped.familyCount; f++){
         ids.Dimension(0);
         for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(ped[i].affections[0]!=0  && CheckCovariates(ped[i]) && geno[i]!=-1)
            ids.Push(i);
         int N = ids.Length();
         if(N==0) continue;
         kin.Setup(*ped.families[f]);
         Omega.Dimension(N, N);
         for(int i = 0; i < N; i++){
            Omega[i][i] = 2*kin(ped[ids[i]], ped[ids[i]]);
            for(int j = i+1; j < N; j++)
            Omega[i][j] = Omega[j][i] = 2*kin(ped[ids[i]], ped[ids[j]]);
         }
         chol.Decompose(Omega);
         Vector temp(N);
         for(int i = 0; i < N; i++)
            temp[i] = resid[baseN+i];
         chol.BackSubst0(temp);
         for(int i = 0; i < N; i++)
            adjResid[baseN+i] = chol.x[i] / sd;
         baseN += N;
*/
