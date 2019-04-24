//////////////////////////////////////////////////////////////////////
// LMM.cpp
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

double Engine::fLL(double x)// -log(likelihood)
{
   double N = UY.Length();
   double J = UX.cols; // # covariates
   Matrix M(J, J);
   for(int i = 0; i < J; i++)
      for(int j = 0; j < J; j++){
         double temp = 0.0;
         for(int k = 0; k < N; k++)
            temp += UX[k][i] * UX[k][j] / (EV[k] + x);
         M[i][j] = temp;
      }
   Cholesky chol;
   chol.Decompose(M);
   Vector tV(J);
   for(int i = 0; i < J; i++){
      double temp = 0.0;
      for(int k = 0; k < N; k++)
         temp += UX[k][i] * UY[k] / (EV[k] + x);
      tV[i] = temp;
   }
   chol.BackSubst(tV);
   double result = 0.0;
   double sigmaG = 0.0;
   for(int i = 0; i < N; i++){
      double temp = 0.0;
      for(int j = 0; j < J; j++)
         temp += UX[i][j] * chol.x[j];
      temp = UY[i] - temp;
      double temp2 = EV[i] + x;
      sigmaG += temp * temp / temp2;
      result += log(temp2);
   }
   sigmaG /= N;
   result += N * (log(2*PI*sigmaG) + 1);
   return result*0.5;
}

void Engine::LMM()
{
   printf("Calculating similarity matrix starts at %s", currentTime());
   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++){
      oneoneCount[i] = 0;
      for(int j = 0; j < 16; j++)
         if(i & shortbase[j]) oneoneCount[i]++;
   }
   if(detailFlag) countGenotype();
   if(geno.Length()==0) {
      individualInfo = true;
      if(shortFlag) BuildShortBinary();
      else BuildBinary();
   }

   int N_Trait = traits.Length();
   if(N_Trait == 0) error("Quantitative traits are not found for LMM analysis.");
   int N_Covariate = covariates.Length();
   IntArray validFlag(ped.count);
   validFlag.Set(1);

   for(int i = 0; i < ped.count; i++){
      bool allmissing = true;
      for(int t = 0; t < N_Trait; t++)
         if(ped[i].isPhenotyped(traits[t])) {
            allmissing = false;
            break;
         }
      if(allmissing) validFlag[i] = 0;
      bool somemissing = false;
      for(int j = 0; j < N_Covariate; j++)
         if(ped[i].covariates[covariates[j]] == _NAN_){
            somemissing = true;
            break;
         }
      if(somemissing) validFlag[i] = 0;
      if(geno[i]<0 || ped[i].ngeno < MINSNPCOUNT) validFlag[i] = 0;
   }
   IntArray ID(0);
   for(int i = 0; i < ped.count; i++)
      if(validFlag[i]) ID.Push(i);
   int N_ID = ID.Length();

   printf("Genotypes stored in %d integers for each of %d individuals used in analysis.\n",
      shortCount, N_ID);
   int id1, id2, count;

   EV.Dimension(N_ID);
   int IBS0Count, het1Count, notMissingCount;
   Matrix D(N_ID, N_ID);
   D.Zero();
   double **UT = new double * [N_ID];
   for(int i = 0; i < N_ID; i++)
      UT[i] = new double [N_ID];
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

   if(normalization)
      printf("Inverse normal transformation is applied to phenotypes.\n");

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

   Vector delta(N_Trait);
   delta.Zero();
   Matrix score0(N_Trait, N_Covariate+1);
   score0.Zero();
   Matrix *Info0 = new Matrix[N_Trait];
   printf("\nPolygenic parameter estimates\n");
   printf("%-15s %7s %7s %7s %7s %9s",
      "TraitName", "Herit", "LogLik", "Delta", "SigmaG", "Mu");
   for(int i = 0; i < N_Covariate; i++)
      printf(" %9s", (const char*)ped.covariateNames[covariates[i]]);
   printf("\n");

   Vector beta(N_Covariate+2);
   Vector beta0(N_Covariate+1);
   double Resid0;
   double Resid;
   Vector sigmaG0(N_Trait);
   double sigmaG;
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
            delta[currentT] = x;
         }
      }

      Info0[currentT].Dimension(N_Covariate+1, N_Covariate+1);
      Info0[currentT].Zero();
      for(int i = 0; i < N_Covariate+1; i++)
         for(int j = 0; j < N_Covariate+1; j++)
            for(int k = 0; k < N_ID; k++)
               Info0[currentT][i][j] += UX[k][i] * UX[k][j] / (EV[k] + delta[currentT]);
      chol.Decompose(Info0[currentT]);
      for(int i = 0; i < N_Covariate+1; i++)
         for(int k = 0; k < N_ID; k++)
            score0[currentT][i] += UX[k][i] * UY[k] / (EV[k] + delta[currentT]);
      chol.BackSubst(score0[currentT]);
      beta0 = chol.x;

      sigmaG0[currentT] = 0.0;
      for(int i = 0; i < N_ID; i++){
         Resid0 = UY[i];
         for(int j = 0; j < N_Covariate+1; j++)
            Resid0 -= UX[i][j] * beta0[j];
         sigmaG0[currentT] += Resid0 * Resid0 / (EV[i] + delta[currentT]);
      }
      sigmaG0[currentT] /= N_ID;
      printf("%-15s %7.4lf %7.2lf %7.4lf %7.3lf",
         (const char*)ped.traitNames[traits[currentT]], 1/(1+delta[currentT]),
         -loglik0[currentT], delta[currentT], sigmaG0[currentT]);
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
   printf("LMM genome scan start...");

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
   Vector *chisq = new Vector[N_Trait];
   for(int t = 0; t < N_Trait; t++)
      chisq[t].Dimension(0);
   String lmmfile = prefix;
   lmmfile.Add("lmm.txt");
   FILE *fp = fopen(lmmfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)lmmfile);
   fprintf(fp, "SNP");
   if(N_Trait > 1)
      fprintf(fp, "\tTraitName");
   if(chromosomes.Length()) fprintf(fp, "\tChr");
   if(bp.Length()) fprintf(fp, "\tPos");
   fprintf(fp, "\tLabelA\tLabela\tFreqA\tLogLik\tBeta\tSE\tZ\tSigmaG\tH2\tPvalue\n");

   int k, mm;
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
      for(int i = 0; i < N_ID; i++){ // most time-consuming loop
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
            Info1 = 0.0;
            for(int k = 0; k < N_ID; k++)
               Info1 += UG[k][m] * UG[k][m] / (EV[k] + delta[t]);
            Info01.Zero();
            for(int j = 0; j < N_Covariate+1; j++)
               for(int k = 0; k < N_ID; k++)
                  Info01[j] += UX[k][j] * UG[k][m] / (EV[k] + delta[t]);
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
               score[N_Covariate+1] += UG[k][m] * UYs[k][t] / (EV[k] + delta[t]);
            if(chol.TryDecompose(Info)==0){
               beta[N_Covariate+1] = se = zstat = pvalue = _NAN_;
            }else{
               chol.Decompose(Info);
               chol.BackSubst(score);
               beta = chol.x;
               sigmaG = 0.0;
               for(int i = 0; i < N_ID; i++){
                  Resid = UYs[i][t] - UG[i][m] * beta[N_Covariate+1];
                  for(int j = 0; j < N_Covariate+1; j++)
                     Resid -= UX[i][j] * beta[j];
                  sigmaG += Resid * Resid / (EV[i] + delta[t]);
               }
               sigmaG /= N_ID;
               loglik = loglik0[t] + 0.5 * N_ID * log(sigmaG/sigmaG0[t]);
               h2 = beta[N_Covariate+1] * beta[N_Covariate+1] * 2 * freq[m] * (1-freq[m]);
               h2 = h2 / (sigmaG*(1+delta[t])+h2);
               chol.Invert();
               se = sqrt(chol.inv[N_Covariate+1][N_Covariate+1] * sigmaG);
               if(se > 1E-200){
                  zstat = beta[N_Covariate+1] / se;
                  pvalue = ndist(fabs(zstat))*2;
                  chisq[t].Push(zstat*zstat);
               }else{
                  zstat = _NAN_;
                  pvalue = 1;
               }
            }
            //output results here
            if(snpName.Length())
               fprintf(fp, "%s", (const char*)snpName[pos]);
            else
               fprintf(fp, "SNP%d", pos+1);
            if(N_Trait>1)
               fprintf(fp, "\t%s", (const char*)ped.traitNames[traits[t]]);
            if(chromosomes.Length())
               fprintf(fp, "\t%d", chromosomes[pos]);
            if(bp.Length())
               fprintf(fp, "\t%.6lf", bp[pos]*0.000001);
            fprintf(fp, "\t%s\t%s\t%.3lf",
               (const char*)alleleLabel[0][pos], (const char*)alleleLabel[1][pos], freq[m]);
            if(zstat != _NAN_)
               fprintf(fp, "\t%.2lf\t%.4lf\t%.4lf\t%.4lf\t%.3lf\t%.4lf\t%.2G\n",
                  loglik, beta[N_Covariate+1], se, zstat, sigmaG, h2, pvalue);
            else
               fprintf(fp, "\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s\n",
                  "NA", "NA", "NA", "NA", "NA", "NA", "NA");
         }
      }
   }
   fclose(fp);
   printf("done\n");
   printf("LMM genome scan ends at %s", currentTime());
   printf("LMM genome scan results saved in file %s\n\n", (const char*)lmmfile);

   if(N_Trait>1){
      printf("\nTraitName\tN_SNP\tSmallestP\tGC_lambda\n");
      for(int t = 0; t < N_Trait; t++){
         idx.Index(chisq[t]);
         printf("%s\t%d\t%.3G\t%.3lf\n",
         (const char*)ped.traitNames[traits[t]],
         chisq[t].Length(),
         chidist(chisq[t].Max(),1),
         chisq[t][idx[chisq[t].Length()/2]]/.456);
      }
   }else{
      idx.Index(chisq[0]);
      printf("\nGC lambda for the genome scan of trait %s in %d autosome SNPS is %.3lf\n",
         (const char*)ped.traitNames[traits[0]],
         chisq[0].Length(), chisq[0][idx[chisq[0].Length()/2]]/.456);
   }
   printf("\n");
   if(Info0) delete []Info0;
   if(chisq) delete []chisq;
   for(int i = 0; i < N_ID; i++)
      delete UG[i];
   delete []UG;
   for(int i = 0; i < N_ID; i++)
      delete UT[i];
   delete []UT;
   for(int i = 0; i < N_ID; i++)
      delete UT2[i];
   delete []UT2;
}


               /*
               loglik = 0.0;
               for(int i = 0; i < N_ID; i++)
                  loglik += log(EV[i] + delta[t]);
               loglik += N_ID * (log(2*PI*sigmaG) + 1);
               loglik = -loglik*0.5;
                 */

         /*
         if(Aa) // genotype: 1
            for(int m = 0; m < MaskArray[Aa][0]; m++){
               mm = MaskArray[Aa][m+1];
               for(int k = 0; k < N_ID; k++)
                  UG.data[k]->data[mm] += UT.data[k]->data[i];
            }
         if(AA) // genotype: 2
            for(int m = 0; m < MaskArray[AA][0]; m++){
               mm = MaskArray[AA][m+1];
               for(int k = 0; k < N_ID; k++)
                  UG.data[k]->data[mm] += UT2.data[k]->data[i];
            }
         if(missing) // genotype: 2*MAF
            for(int m = 0; m < MaskArray[missing][0]; m++){
               mm = MaskArray[missing][m+1];
               for(int k = 0; k < N_ID; k++)
                  UG.data[k]->data[mm] += UT2.data[k]->data[i] * freq[m];
            }
            */

               /*
   IntArray MaskArray[65536];
   for(int i = 0; i < 65536; i++){
      MaskArray[i].Dimension(1);
      MaskArray[i][0] = 0;
      for(int k = 0; k < 16; k++)
         if(i & (1<<k)){
            MaskArray[i].Push(k);
            MaskArray[i][0]++;
         }
   } */

