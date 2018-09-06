//////////////////////////////////////////////////////////////////////
// VC.cpp
// (c) 2010-2018 Wei-Min Chen
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
// May 9, 2018

#include "analysis.h"
#include <math.h>
#include "Kinship.h"
#include "KinshipX.h"
#include "MathStats.h"
#include "MathSVD.h"
#include "QuickIndex.h"
#include "MathCholesky.h"

#ifdef WITH_LAPACK
extern "C" void dgesdd_(char*, int*, int*, double*, int*, double*, double *, int*, double*, int*, double*, int*, int*, int*);
extern "C" void dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);
#endif

const double PI = 4.0*atan(1.0);

void Engine::PreVC()
{
   int N_Trait = traits.Length();
   if(N_Trait == 0) error("Quantitative traits are not found for variance component analysis.");
   int N_Covariate = covariates.Length();
   IntArray validFlag(ped.count);
   validFlag.Set(1);
   for(int i = 0; i < ped.count; i++){
      bool somemissing = false;
      for(int t = 0; t < N_Trait; t++)
         if(!ped[i].isPhenotyped(traits[t])) {
            somemissing = true;
            break;
         }
      if(somemissing) validFlag[i] = 0;
      somemissing = false;
      for(int j = 0; j < N_Covariate; j++)
         if(ped[i].covariates[covariates[j]] == _NAN_){
            somemissing = true;
            break;
         }
      if(somemissing) validFlag[i] = 0;
      if(geno[i]<0 || ped[i].ngeno == 0) validFlag[i] = 0;
   }
   ID.Dimension(0);
   for(int i = 0; i < ped.count; i++)
      if(validFlag[i]) ID.Push(i);
   int N_ID = ID.Length();
   if(N_Trait>1)
      printf("Individuals with less than %d phenotypes are removed.\n", N_Trait);
   printf("Genotypes stored in %d integers for each of %d individuals used in analysis.\n",
      shortCount, N_ID);
   EV.Dimension(N_ID);
   UT = new double * [N_ID];
   for(int i = 0; i < N_ID; i++)
      UT[i] = new double [N_ID];
   if(svdinfile!="")
      ReadSVD();
   else
      ComputeSVD();
   if(normalization)
      printf("Inverse normal transformation is applied to phenotypes.\n");

   chisqs = new Vector[N_Trait];
   for(int t = 0; t < N_Trait; t++)
      chisqs[t].Dimension(0);
}

void Engine::VC()
{
   if(detailFlag) countGenotype();
   if(geno.Length()==0) {
      individualInfo = true;
      if(shortFlag) BuildShortBinary();
      else BuildBinary();
   }
   const int max_iter = 20;
   PreVC();

   int N_Trait = traits.Length();
   int N_Covariate = covariates.Length();
   int N_ID = ID.Length();
   Matrix Y(N_ID, N_Trait);
   Matrix X(N_ID, 1+N_Covariate);
   for(int i = 0; i < N_ID; i++)
      X[i][0] = 1.0;
   for(int t = 0; t < N_Trait; t++)
      for(int i = 0; i < N_ID; i++)
         Y[i][t] = ped[ID[i]].traits[traits[t]];
   for(int i = 0; i < N_ID; i++)
      for(int j = 0; j < N_Covariate; j++)
         X[i][j+1] = ped[ID[i]].covariates[covariates[j]];

   Vector meanY(N_Trait);
   meanY.Set(0.0);
   int n;
   for(int t = 0; t < N_Trait; t++){
      n = 0;
      for(int i = 0; i < N_ID; i++)
         if(Y[i][t] != _NAN_) {
            meanY[t] += Y[i][t];
            n ++;
         }
      meanY[t] /= n;
      for(int i = 0; i < N_ID; i++)
         if(Y[i][t] == _NAN_)
            Y[i][t] = meanY[t];
   }

   Matrix UYs(N_ID, N_Trait); // U'Y += UT[i][k] * Y[k]
   UYs.Zero();
   UX.Dimension(N_ID, N_Covariate+1);
   UX.Zero();
   Vector tV(N_ID);

   Cholesky chol;
   QuickIndex idx;
   for(int t = 0; t < N_Trait; t++){
      if(normalization){
         for(int i = 0; i < N_ID; i++)
            tV[i] = Y[i][t];
         idx.Index(tV);
         for(int i = 0; i < N_ID;){
            int start = i, end = i + 1;
            while(end < N_ID && tV[idx[end]] == tV[idx[start]]) end++;
            end --;
            double q = ninv((start + (end-start)/2.0 + 0.5)/N_ID);
            for(int j = start; j <= end; j++)
               Y[idx[j]][t] = q;
            i = end + 1;
         }
      }
      for(int i = 0; i < N_ID; i++)
         for(int k = 0; k < N_ID; k++)
            UYs[i][t] += UT[i][k] * Y[k][t];
   }
   for(int i = 0; i < N_ID; i++)
      for(int j = 0; j < N_Covariate+1; j++)
         for(int k = 0; k < N_ID; k++)
            UX[i][j] += UT[i][k] * X[k][j];

   UY.Dimension(N_ID);
   double a, b, LB, T, x, temp, funcx;
   int numiter;
   int maxIter=10000;
   bool quiet=true;
   T = 0.001;
   double R1R;

   lambda0.Dimension(N_Trait);
   lambda0.Zero();
   Matrix score0(N_Trait, N_Covariate+1);
   score0.Zero();
   Matrix *Info0 = new Matrix[N_Trait];
   printf("\nPolygenic parameter estimates\n");
   printf("%-15s %7s %7s %7s %7s %9s",
      "TraitName", "Herit", "LogLik", "Lambda", "Tau", "Mu");
   for(int i = 0; i < N_Covariate; i++)
      printf(" %9s", (const char*)ped.covariateNames[covariates[i]]);
   printf("\n");

   Vector beta(N_Covariate+2);
   Vector beta0(N_Covariate+1);
   double Resid0;
   Vector Resid(N_ID);
   tau0.Dimension(N_Trait);
   double lambda;
   double loglik;
   Vector loglik0(N_Trait);

   for(int t = 0; t < N_Trait; t++){
      currentT = t;
      for(int i = 0; i < N_ID; i++)
         UY[i] = UYs[i][currentT];
      loglik0[currentT] = 1E10;
      for(int i = 0; i < 100; i++){
         a = exp(i*0.1 - 5);
         b = exp(i*0.1 - 4.9);
         x = minimize(a, b, T, funcx, numiter, maxIter, quiet);
         if(funcx < loglik0[currentT]) {
            loglik0[currentT] = funcx;
            lambda0[currentT] = 1/x;
         }
      }
      if(lambda0[currentT]<0.0067){  // border
         lambda0[currentT] = 0;
      }
      Info0[currentT].Dimension(N_Covariate+1, N_Covariate+1);
      Info0[currentT].Zero();
      for(int i = 0; i < N_Covariate+1; i++)
         for(int j = 0; j < N_Covariate+1; j++)
            for(int k = 0; k < N_ID; k++)
               Info0[currentT][i][j] += UX[k][i] * UX[k][j] / (lambda0[currentT] * EV[k] + 1);
      for(int i = 0; i < N_Covariate+1; i++)
         for(int k = 0; k < N_ID; k++)
            score0[currentT][i] += UX[k][i] * UY[k] / (lambda0[currentT] * EV[k] + 1);
      if(N_Covariate==0)
         beta0[0] = score0[currentT][0]/Info0[currentT][0][0];
      else if(chol.TryDecompose(Info0[currentT])==0)
         beta0.Zero();
      else{
         chol.Decompose(Info0[currentT]);
         chol.BackSubst(score0[currentT]);
         beta0 = chol.x;
      }
      R1R = 0;
      for(int i = 0; i < N_ID; i++){
         Resid0 = UY[i];
         for(int j = 0; j < N_Covariate+1; j++)
            Resid0 -= UX[i][j] * beta0[j];
         R1R += Resid0 * Resid0 / (lambda0[currentT] * EV[i]+1);
      }
      tau0[currentT] = N_ID / R1R;
      printf("%-15s %7.4lf %7.2lf %7.4lf %7.3lf",
         (const char*)ped.traitNames[traits[currentT]], lambda0[currentT]/(1+lambda0[currentT]),
         -loglik0[currentT], lambda0[currentT], tau0[currentT]);
      for(int i = 0; i < N_Covariate+1; i++)
         printf(" %9.3lf", beta0[i]);
      printf("\n");
   }
   if(HeritFlag) {
      if(Info0) delete []Info0;
      for(int i = 0; i < N_ID; i++)
         delete UT[i];
      delete []UT;
      return;
   }
   /*
   printf("Variance component genome scan starts at %s", currentTime());

   int id, AA, Aa, missing, lowb, highb;
   double Info1;
   Vector Info01(N_Covariate+1);
   double se, zstat, pvalue, h2;
   Matrix Info(N_Covariate+2, N_Covariate+2);
   Vector score(N_Covariate+2);
   double **UG;
   UG = new double * [N_ID];
   for(int i = 0; i < N_ID; i++)
      UG[i] = new double [16];

   int MaskArray[256][9];
   for(int i = 0; i < 256; i++){
      MaskArray[i][0] = 0;
      for(int k = 0; k < 8; k++)
         if(i & (1<<k)){
            MaskArray[i][0]++;
            MaskArray[i][MaskArray[i][0]] = k;
         }
   }
   double **UT2 = new double * [N_ID];
   for(int i = 0; i < N_ID; i++)
      UT2[i] = new double [N_ID];
   for(int i = 0; i < N_ID; i++)
      for(int j = 0; j < N_ID; j++)
         UT2[i][j] = UT[i][j] * 2;
   double freq[16];
   int AACount[16], AaCount[16], missingCount[16];
   int pos;
   String lmmfile = prefix;
   lmmfile.Add("vc.txt");
   FILE *fp = fopen(lmmfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)lmmfile);
   fprintf(fp, "SNP");
   if(N_Trait > 1)
      fprintf(fp, "\tTraitName");
   if(chromosomes.Length()) fprintf(fp, "\tChr");
   if(positions.Length()) fprintf(fp, "\tPos");
   fprintf(fp, "\tLabelA\tLabela\tFreqA\tN\tNIter\tLogLik\tP_MLE\tLambda\tTao\tBeta\tSE\tZ\tH2\tPvalue\n");

   int iter=0;
   int k, mm;
   double deriv1, deriv2, R2R, R3R;
   double pMLE;
   bool failureFlag;
   for(int b = 0; b < shortCount; b++){ // loop over short integer
      for(int i = 0; i < 16; i++){
         AACount[i] = AaCount[i] = missingCount[i] = 0;
         freq[i] = 0.0;
      }
      for(int i = 0; i < N_ID; i++){
         id = geno[ID[i]];
         AA = GG[0][id][b] & GG[1][id][b];
         Aa = (~GG[0][id][b]) & GG[1][id][b];
         missing = (~GG[0][id][b]) & (~GG[1][id][b]) & 65535;
         for(int j = 0; j < 16; j++){
            if(AA & shortbase[j])
               AACount[j] ++;
            else if(Aa & shortbase[j])
               AaCount[j] ++;
            else if(missing & shortbase[j])
               missingCount[j] ++;
         }
      }
      for(int j = 0; j < 16; j++)
         if(missingCount[j] < N_ID)  // frequency of allele A
            freq[j] = (AACount[j] + AaCount[j]*0.5) / (N_ID - missingCount[j]);

      for(int i = 0; i < N_ID; i++)
         for(int m = 0; m < 16; m++)
            UG[i][m] = 0;
      for(int i = 0; i < N_ID; i++){ // O(N^2): most time-consuming loop in genome scan
         id = geno[ID[i]];
         AA = GG[0][id][b] & GG[1][id][b];
         Aa = (~GG[0][id][b]) & GG[1][id][b];
         missing = (~GG[0][id][b]) & (~GG[1][id][b]) & 65535;
         if(Aa & 255){
            lowb = Aa & 255;
            for(int m = 0; m < MaskArray[lowb][0]; m++){
               mm = MaskArray[lowb][m+1];
               for(k = 0; k < N_ID; k++)
                  UG[k][mm] += UT[k][i];
            }
         }
         if(Aa & 65280){
            highb = Aa >> 8;
            for(int m = 0; m < MaskArray[highb][0]; m++){
               mm = MaskArray[highb][m+1]+8;
               for(k = 0; k < N_ID; k++)
                  UG[k][mm] += UT[k][i];
            }
         }
         if(AA & 255){
            lowb = AA & 255;
            for(int m = 0; m < MaskArray[lowb][0]; m++){
               mm = MaskArray[lowb][m+1];
               for(k = 0; k < N_ID; k++)
                  UG[k][mm] += UT2[k][i];
            }
         }
         if(AA & 65280){
            highb = AA >> 8;
            for(int m = 0; m < MaskArray[highb][0]; m++){
               mm = MaskArray[highb][m+1]+8;
               for(k = 0; k < N_ID; k++)
                  UG[k][mm] += UT2[k][i];
            }
         }

         if(missing & 255){
            lowb = missing & 255;
            for(int m = 0; m < MaskArray[lowb][0]; m++){
               mm = MaskArray[lowb][m+1];
               for(k = 0; k < N_ID; k++)
                  UG[k][mm] += UT2[k][i] * freq[mm];
            }
         }
         if(missing & 65280){
            highb = missing >> 8;
            for(int m = 0; m < MaskArray[highb][0]; m++){
               mm = MaskArray[highb][m+1]+8;
               for(k = 0; k < N_ID; k++)
                  UG[k][mm] += UT2[k][i] * freq[mm];
            }
         }
      }
      for(int m = 0; m < 16; m++){
         pos = b*16 + m;
         if(pos >= markerCount) break;
         for(int t = 0; t < N_Trait; t++){
            failureFlag = false;
            lambda = lambda0[t];
            Info1 = 0.0;
            for(int k = 0; k < N_ID; k++)
               Info1 += UG[k][m] * UG[k][m] / (EV[k]*lambda + 1);
            Info01.Zero();
            for(int j = 0; j < N_Covariate+1; j++)
               for(int k = 0; k < N_ID; k++)
                  Info01[j] += UX[k][j] * UG[k][m] / (EV[k]*lambda + 1);
            for(int i = 0; i < N_Covariate+1; i++)
               for(int j = 0; j < N_Covariate+1; j++)
                  Info[i][j] = Info0[t][i][j];
            Info[N_Covariate+1][N_Covariate+1] = Info1;
            for(int i = 0; i < N_Covariate+1; i++)
               Info[i][N_Covariate+1] = Info[N_Covariate+1][i] = Info01[i];
            score.Zero();
            for(int i = 0; i < N_Covariate+1; i++)
               score[i] = score0[t][i];
            for(int k = 0; k < N_ID; k++)
               score[N_Covariate+1] += UG[k][m] * UYs[k][t] / (EV[k]*lambda + 1);
            if(chol.TryDecompose(Info)==0){
               beta[N_Covariate+1] = se = zstat = pvalue = _NAN_;
               failureFlag = true;
            }else{
               chol.Decompose(Info);
               chol.BackSubst(score);
               beta = chol.x;
               for(int i = 0; i < N_ID; i++){
                  Resid[i] = UYs[i][t] - UG[i][m] * beta[N_Covariate+1];
                  for(int j = 0; j < N_Covariate+1; j++)
                     Resid[i] -= UX[i][j] * beta[j];
               }
               R1R = 0;
               for(int i = 0; i < N_ID; i++)
                  R1R += Resid[i] * Resid[i] / (lambda*EV[i]+1);
               loglik = N_ID*log(2*PI*R1R/N_ID) + N_ID;
               for(int i = 0; i < N_ID; i++)
                  loglik += log(lambda*EV[i]+1);
               loglik *= 0.5;

               if(!noiterFlag){
                  double loglik_init = loglik;
                  double delta_local = 100;
                  for(iter = 0; (iter < max_iter) && ((delta_local < -0.0001) || (delta_local > 0.0001)); iter++){
                     R2R = R3R = deriv1 = deriv2 = 0.0;
                     for(int i = 0; i < N_ID; i++){
                        double temp = lambda*EV[i]+1;
                        R2R += Resid.data[i] * Resid.data[i] * EV[i] / temp / temp;
                        R3R += Resid.data[i] * Resid.data[i] * EV[i] * EV[i] / temp / temp / temp;
                        deriv1 += EV[i] / temp;
                        deriv2 += EV[i] * EV[i] / temp / temp;
                     }
                     deriv1 = N_ID * 0.5 * R2R / R1R - 0.5 * deriv1;
                     deriv2 = N_ID * 0.5 * (R2R*R2R-2*R3R*R1R) / R1R / R1R + 0.5 * deriv2;
                     delta_local = -deriv1 / deriv2;
                     lambda += delta_local;
                     Info.Zero();
                     score.Zero();
                     for(int k = 0; k < N_ID; k++){
                        for(int i = 0; i < N_Covariate+1; i++)
                           for(int j = 0; j < N_Covariate+1; j++)
                              Info.data[i]->data[j] += UX.data[k]->data[i] * UX.data[k]->data[j] / (lambda*EV[k]+1);
                        Info.data[N_Covariate+1]->data[N_Covariate+1] += UG[k][m] * UG[k][m] / (lambda*EV[k]+1);
                        for(int i = 0; i < N_Covariate+1; i++)
                           Info.data[i]->data[N_Covariate+1] += UG[k][m] * UX.data[k]->data[i] / (lambda*EV[k]+1);
                     }
                     for(int i = 0; i < N_Covariate+1; i++)
                        Info.data[N_Covariate+1]->data[i] = Info.data[i]->data[N_Covariate+1];
                     if(chol.TryDecompose(Info)==0){
                        failureFlag = true;
                        //beta[N_Covariate+1] = se = zstat = pvalue = _NAN_;
                        break;
                     }
                     chol.Decompose(Info);
                     for(int k = 0; k < N_ID; k++){
                        for(int i = 0; i < N_Covariate+1; i++)
                           score[i] += UX[k][i] * UYs[k][t] / (lambda*EV[k]+1);
                        score[N_Covariate+1] += UG[k][m] * UYs[k][t] / (lambda*EV[k]+1);
                     }
                     chol.BackSubst(score);
                     beta = chol.x;
                     for(int i = 0; i < N_ID; i++){
                        Resid.data[i] = UYs.data[i]->data[t] - UG[i][m] * beta.data[N_Covariate+1];
                        for(int j = 0; j < N_Covariate+1; j++)
                           Resid.data[i] -= UX.data[i]->data[j] * beta.data[j];
                     }
                     R1R=0;
                     for(int i = 0; i < N_ID; i++)
                        R1R += Resid.data[i] * Resid.data[i] / (lambda*EV[i]+1);
                  }
                  loglik = N_ID*log(2*PI*R1R/N_ID) + N_ID;
                  for(int i = 0; i < N_ID; i++)
                     loglik += log(lambda*EV[i]+1);
                  loglik *= 0.5;
                  if(loglik > loglik_init || failureFlag){ // MLE failed
                     loglik = loglik_init;
                     lambda = lambda0[t];
                     iter = -1;
                     Info1 = 0.0;
                     for(int k = 0; k < N_ID; k++)
                        Info1 += UG[k][m] * UG[k][m] / (EV[k]*lambda + 1);
                     Info01.Zero();
                     for(int j = 0; j < N_Covariate+1; j++)
                        for(int k = 0; k < N_ID; k++)
                           Info01.data[j] += UX.data[k]->data[j] * UG[k][m] / (EV[k]*lambda + 1);
                     for(int i = 0; i < N_Covariate+1; i++)
                        for(int j = 0; j < N_Covariate+1; j++)
                           Info[i][j] = Info0[t][i][j];
                     Info[N_Covariate+1][N_Covariate+1] = Info1;
                     for(int i = 0; i < N_Covariate+1; i++)
                        Info[i][N_Covariate+1] = Info[N_Covariate+1][i] = Info01[i];
                     score.Zero();
                     for(int i = 0; i < N_Covariate+1; i++)
                        score[i] = score0[t][i];
                     for(int k = 0; k < N_ID; k++)
                        score.data[N_Covariate+1] += UG[k][m] * UYs.data[k]->data[t] / (EV[k]*lambda + 1);
                     chol.Decompose(Info);
                     chol.BackSubst(score);
                     beta = chol.x;
                     for(int i = 0; i < N_ID; i++){
                        Resid[i] = UYs[i][t] - UG[i][m] * beta[N_Covariate+1];
                        for(int j = 0; j < N_Covariate+1; j++)
                           Resid[i] -= UX[i][j] * beta[j];
                     }
                     R1R = 0;
                     for(int i = 0; i < N_ID; i++)
                        R1R += Resid[i] * Resid[i] / (lambda*EV[i]+1);
                  }
               }
               h2 = beta[N_Covariate+1] * beta[N_Covariate+1] * 2 * freq[m] * (1-freq[m]);
               h2 = h2 / ((1+lambda)*R1R/N_ID+h2);
               chol.Invert();
               se = sqrt(chol.inv[N_Covariate+1][N_Covariate+1] * R1R / N_ID);
               if(se < 1E-200 || failureFlag){
                  beta[N_Covariate+1] = se = zstat = pvalue = pMLE = _NAN_;
               }else{
                  zstat = beta[N_Covariate+1] / se;
                  pvalue = ndist(fabs(zstat))*2;
                  chisqs[t].Push(zstat*zstat);
                  pMLE = ndist(sqrt((loglik0[t] - loglik)*2))*2;
               }
            }
            if(snpName.Length())
               fprintf(fp, "%s", (const char*)snpName[pos]);
            else
               fprintf(fp, "SNP%d", pos+1);
            if(N_Trait>1)
               fprintf(fp, "\t%s", (const char*)ped.traitNames[traits[t]]);
            if(chromosomes.Length())
               fprintf(fp, "\t%d", chromosomes[pos]);
            if(positions.Length())
               fprintf(fp, "\t%.6lf", positions[pos]);
            fprintf(fp, "\t%s\t%s\t%.3lf\t%d",
               (const char*)alleleLabel[0][pos], (const char*)alleleLabel[1][pos],
               freq[m], N_ID-missingCount[m]);
            if(zstat != _NAN_)
               fprintf(fp, "\t%d\t%.2lf\t%.3G\t%.3lf\t%.3lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.3G\n",
                  iter, -loglik, pMLE,lambda, N_ID / R1R, beta[N_Covariate+1], se, zstat, h2, pvalue);
            else
               fprintf(fp, "\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s\n",
                  "NA", "NA", "NA", "NA", "NA", "NA", "NA");
         }
      }
   }
   fclose(fp);
   printf("Variance component genome scan ends at %s", currentTime());
   printf("Variance component genome scan results saved in file %s\n\n", (const char*)lmmfile);

   PostVC();
   for(int i = 0; i < N_ID; i++)
      delete UT[i];
   delete []UT;
   if(Info0) delete []Info0;
   for(int i = 0; i < N_ID; i++)
      delete UG[i];
   delete []UG;
   for(int i = 0; i < N_ID; i++)
      delete UT2[i];
   delete []UT2;
   */
}



void Engine::ReadSVD()
{
   StringArray tokens;
   String line;
   int N_ID = ID.Length();
   IFILE input = ifopen(svdinfile, "rt");
   if(input == NULL)
      error("Cannot open %s to read", (const char*)svdinfile);
   tokens.Clear();
   line.ReadLine(input);
   tokens.AddTokens(line);
   if(tokens.Length() != N_ID)
      error("%d samples in SVD file %s; however there are %d samples in the data.\n",
            tokens.Length(), (const char*)svdinfile, N_ID);
   String tempID;
   IntArray key(N_ID);
   for(int i = 0; i < N_ID; i++){
         tempID = ped[ID[i]].famid;
         tempID.Add("_");
         tempID.Add(ped[ID[i]].pid);
         if(tempID == tokens[i])
            key[i] = i;
         else{
            int k = tokens.Find(tempID);
            if(k < 0) error("Family %s ID %s cannot be found in SVD file",
               (const char*)ped[ID[i]].famid, (const char*)ped[ID[i]].pid);
            key[k] = i;
         }
   }
   tokens.Clear();
   line.ReadLine(input);
   tokens.AddTokens(line);
   if(tokens.Length() < N_ID)
      error("There are only %d Eigenvalues (not %d) in the 2nd row of the SVD file",
         tokens.Length(), N_ID);
   for(int i = 0; i < N_ID; i++)
      EV[key[i]] = tokens[i].AsDouble();
   int row = 0;
   while(!ifeof(input)){
      if(row == N_ID) break;
      tokens.Clear();
      line.ReadLine(input);
      tokens.AddTokens(line);
      if(tokens.Length() < N_ID) continue;
      for(int i = 0; i < N_ID; i++)
         UT[key[row]][key[i]] = tokens[i].AsDouble();
      row++;
   }
   if(row < N_ID)
      error("There are only %d rows (not %d) in the Eigenvector data",
            row, N_ID);
   printf("Eigenvalues and Eigenvectors are successfully loaded from file %s\n",
         (const char*)svdinfile);
   ifclose(input);
}

void Engine::ComputeSVD()
{
   int IBS0Count, het1Count, notMissingCount;
   int id1, id2, count;
   int N_ID = ID.Length();

   printf("Calculating similarity matrix starts at %s", currentTime());
   Matrix D(N_ID, N_ID);
   D.Zero();
   if(Bit64==64){
      bool IBDvalidFlag = PreSegment();
      if(!IBDvalidFlag){
         printf("%s\n", (const char*)segmessage);
         printf("  Note chromosomal positions can be sorted conveniently using other tools such as PLINK.\n");
         return;
      }                 
      IntArray pairIndex[2];
      for(int i = 0; i < 2; i++)
         pairIndex[i].Dimension(0);
      IntArray allpairs(0);
      for(int i = 0; i < N_ID; i++)
         for(int j = i+1; j < N_ID; j++){
            if(geno[ID[i]] < 0 || geno[ID[j]] < 0) continue;
            allpairs.Push(geno[ID[i]]);
            allpairs.Push(geno[ID[j]]);
            pairIndex[0].Push(i);
            pairIndex[1].Push(j);
         }
      IntArray ibdsegStorage[2], ibdsegIndex[2];
      int segCount = (chrSeg.Length()>>2);
      for(int seg = 0; seg < segCount; seg++){
         IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], true);
         int pairCount = pairIndex[0].Length();
         for(int p = 0; p < pairCount; p++){
            int i = pairIndex[0][p];
            int j = pairIndex[1][p];
            D[i][j] += (ibdsegStorage[0][p]*0.5 + ibdsegStorage[1][p])*0.000001;
         }  // end of pair p loop
      }  // end of seg loop
      for(int i = 0; i < N_ID; i++)
         for(int j = i+1; j < N_ID; j++){
            D[i][j] /= totalLength;
            D[j][i] = D[i][j];
         }
//      D.Print(stdout, 20,20);
   }else{
      char oneoneCount[65536];
      for(int i = 0; i < 65536; i++){
         oneoneCount[i] = 0;
         for(int j = 0; j < 16; j++)
            if(i & shortbase[j]) oneoneCount[i]++;
      }
      for(int i = 0; i < N_ID; i++)
         for(int j = i+1; j < N_ID; j++){
            id1 = geno[ID[i]]; id2 = geno[ID[j]];
            if(id1 < 0 || id2 < 0) continue;
            IBS0Count = het1Count = notMissingCount = 0;
            for(int m = 0; m < shortCount; m++){
               IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
               notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
               het1Count += oneoneCount[(GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m]) |
                     (GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]) ];
            }
            if(notMissingCount)
               D[i][j] = D[j][i] = (het1Count + 4.0*IBS0Count) / notMissingCount;
         }
         Vector tempV(N_ID);
         tempV.Zero();
         for(int j = 0; j < N_ID; j++)
            for(int i = 0; i < N_ID; i++)
               tempV[j] += D[i][j];
         tempV.Multiply(1.0/N_ID);
         // (I-11'/N) * D
         for(int i = 0; i < N_ID; i++)
            for(int j = 0; j < N_ID; j++)
               D[i][j] -= tempV[j];
         tempV.Zero();
         for(int i = 0; i < N_ID; i++)
            for(int j = 0; j < N_ID; j++)
               tempV[i] += D[i][j];
         tempV.Multiply(1.0/N_ID);
         // D * (I-11'/N)
         for(int i = 0; i < N_ID; i++)
            for(int j = 0; j < N_ID; j++)
               D[i][j] -= tempV[i];
         D.Multiply(-0.5);
         for(int i = 0; i < N_ID; i++)
            D[i][i] = sqrt(D[i][i]); // save computation for below
         for(int i = 0; i < N_ID; i++)
            for(int j = i+1; j < N_ID; j++)
               D[i][j] = D[j][i] = D[i][j] / (D[i][i]*D[j][j]);
   }
   for(int i = 0; i < N_ID; i++)
      D[i][i] = 1.0;

   printf("Similarity matrix consists of sample covariance estimates, standardized by diagonal elements\n");
   int degree[3];
   for(int d = 0; d < 3; d++)
      degree[d] = 0;
   for(int i = 0; i < N_ID; i++)
      for(int j = i+1; j < N_ID; j++)
         if(D[i][j] > 0.8) // MZ twins
            degree[0] ++;
         else if(D[i][j] > 0.354) // 1st-degree relatives
            degree[1] ++;
         else if(D[i][j] > 0.177) // 2nd-degree relatives
            degree[2] ++;
   if(degree[0] || degree[1] || degree[2])
      printf("Close relatives identified: %d MZ twins/duplicates, %d 1st-degree and %d 2nd-degree relative pairs\n",
         degree[0], degree[1], degree[2]);
   else
      printf("No second-degree or closer relationships were found\n");
   printf("SVD starts at %s...", currentTime());
#ifdef WITH_LAPACK
   printf("  LAPACK is used.\n");
   char JOBZ = 'A';
   int info;
   double *A = new double[N_ID*N_ID];
   int dimLA = N_ID;
   for(int i = 0; i < N_ID; i++)
     for(int j = 0; j < N_ID; j++)
       A[i*N_ID+j] = D[j][i];
   double *S = new double[N_ID];
   double *U = new double[N_ID*N_ID];
   double *VT = new double[N_ID*N_ID];
   int *IWORK = new int[N_ID*8];
   int LWORK = 8*N_ID + 4*N_ID*N_ID;
   double *WORK = new double[LWORK];
   dgesdd_(&JOBZ, &dimLA, &dimLA, A, &dimLA, S, U, &dimLA, VT, &dimLA, WORK, &LWORK, IWORK, &info);

   delete []U;
   delete []IWORK;
   delete []WORK;
   delete []A;
   for(int i = 0; i < N_ID; i++) EV[i] = S[i];
   delete S;
   // VT[k*N_ID+j] stores jth eigenvector. k: id; j: marker
   for(int j = 0; j < N_ID; j++)
      for(int k = 0; k < N_ID; k++)
         UT[j][k] = VT[k*N_ID+j];
   delete []VT;
#else
   SVD svd;
   svd.Decompose(D);
   printf("done\n");
   if(svd.n == 0) return;
   for(int i = 0; i < N_ID; i++) EV[i] = svd.w[i];
   // svd.v[k][idx[N_ID-1-j]] stores jth eigenvector
   for(int j = 0; j < N_ID; j++)
      for(int k = 0; k < N_ID; k++)
         UT[j][k] = svd.v[k][j];
#endif
   printf("SVD ends at %s", currentTime());
   printf("Eigenvalues range from %G to %G.\n", EV.Min(), EV.Max());

   if(svdoutFlag){
      String svdoutfile = prefix;
      svdoutfile.Add(".svd");
      FILE *fp = fopen(svdoutfile, "wt");
      if(fp == NULL) error("Cannot open %s to write.", (const char*)svdoutfile);
      fprintf(fp, "%s_%s",
         (const char*)ped[ID[0]].famid, (const char*)ped[ID[0]].pid);
      for(int i = 1; i < N_ID; i++)
         fprintf(fp, "\t%s_%s", (const char*)ped[ID[i]].famid, (const char*)ped[ID[i]].pid);
      fprintf(fp, "\n");
      fprintf(fp, "%lf", EV[0]);
      for(int i = 1; i < N_ID; i++)
         if(EV[i] > 1E-6 || EV[i] < -1E-6)
            fprintf(fp, "\t%lf", EV[i]);
         else
            fprintf(fp, "\t%G", EV[i]);
      fprintf(fp, "\n");
      for(int i = 0; i < N_ID; i++){
         fprintf(fp, "%G", UT[i][0]);
         for(int j = 1; j < N_ID; j++)
            fprintf(fp, "\t%G", UT[i][j]);
         fprintf(fp, "\n");
      }
      for(int i = 0; i < N_ID; i++){
         fprintf(fp, "%G", D[i][0]);
         for(int j = 1; j < N_ID; j++)
            fprintf(fp, "\t%G", D[i][j]);
         fprintf(fp, "\n");
      }
      fclose(fp);
      printf("Eigenvalues, eigenvectors and genetic similarity are saved in file %s\n",
         (const char*)svdoutfile);
      exit(0);
   }
}

   /*
      for(int i = 0; i < N_ID; i++)
         for(int j = i+1; j < N_ID; j++){
            id1 = geno[ID[i]]; id2 = geno[ID[j]];
            if(id1 < 0 || id2 < 0) continue;
            IBS0Count = het1Count = notMissingCount = 0;
            for(int m = 0; m < longCount; m++){
               IBS0Count += popcount(LG[0][id1][m] & LG[0][id2][m] & (LG[1][id1][m] ^ LG[1][id2][m]));
               notMissingCount += popcount((LG[0][id1][m] | LG[1][id1][m]) & (LG[0][id2][m] | LG[1][id2][m]));
               het1Count += popcount((LG[0][id2][m] & (~LG[0][id1][m]) & LG[1][id1][m]) |
                     (LG[0][id1][m] & (~LG[0][id2][m]) & LG[1][id2][m]) );
            }
            if(notMissingCount)
               D[i][j] = D[j][i] = (het1Count + 4.0*IBS0Count) / notMissingCount;
         }
*/

void Engine::PostVC()
{
   int N_Trait = traits.Length();
   QuickIndex idx;
   if(N_Trait>1){
      printf("\nTraitName\tN_SNP\tSmallestP\tGC_lambda\n");
      for(int t = 0; t < N_Trait; t++){
         idx.Index(chisqs[t]);
         printf("%s\t%d\t%.3G\t%.3lf\n",
         (const char*)ped.traitNames[traits[t]],
         chisqs[t].Length(),
         chidist(chisqs[t].Max(),1),
         chisqs[t][idx[chisqs[t].Length()/2]]/.456);
      }
   }else{
      idx.Index(chisqs[0]);
      printf("\nGC lambda for the genome scan of trait %s in %d autosome SNPS is %.3lf\n",
         (const char*)ped.traitNames[traits[0]],
         chisqs[0].Length(), chisqs[0][idx[chisqs[0].Length()/2]]/.456);
   }
   printf("\n");
   if(chisqs) delete []chisqs;
}

