//////////////////////////////////////////////////////////////////////
// ancestry.cpp
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
// August 7, 2018

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

void Engine::SlidingWindows()
{
   double sizeMainWindow = 20;
   double sizeOffset = 2;

   double mv_start, mv_stop;
   Vector ancestry0, tancestry;
   Matrix allancestry;
   IntArray ancestryOneWindow;
   for(int chr = 0; chr < SEXCHR; chr++){
      int nW = 10; // (maxPos - sizeOffset) / sizeMainWindow + 1;

      // for each window
      // mv_start, mv_stop
      tancestry = ancestry0;
      OneWindowProjection(chr, mv_start, mv_stop, tancestry);
      // derive ancestryOneWindow from tancestry
      // majority vote here

      // allancestry 

   }


}

void Engine::OneWindowProjection(int mv_chr, double mv_start, double mv_stop, Vector & ancestry)
{
   IntArray ID(0), ID_UN(0);
   int N = ancestry.Length();
   for(int i = 0; i < N; i++){
      if(ancestry[i] == 0 || ancestry[i] == 1) ID_UN.Push(i);
      ID.Push(i);
   }
   int dimN = ID_UN.Length();
   if(dimN < 2) return;
   int dimM = markerCount;
   Matrix X(dimM, dimN);
   Matrix Z(dimM, ID.Length());
   X.Zero(); Z.Zero();
   int minorAllele;
   int validSNP=0;
   int byte, offset;
   double p, g, mean, se;
   int k, freqCount;
   bool flipFlag;

   for(int m = 0; m < markerCount; m ++){
      if(chromosomes[m]!=mv_chr || positions[m] < mv_start || positions[m] >= mv_stop)
         continue;
      byte = m/16;
      offset = m%16;

      p = 0.0;
      freqCount = 0;
      for(int i = 0; i < dimN; i++){
         k = ID_UN[i];
         if(GG[1][k][byte] & shortbase[offset]){ // AA or Aa
            p ++;
            if(GG[0][k][byte] & shortbase[offset]) p ++; // AA
            freqCount += 2;
         }else if(GG[0][k][byte] & shortbase[offset]) // aa
            freqCount += 2;
      }
      if(freqCount == 0) continue;
      p /= freqCount;
      flipFlag = false;
      if(p > 0.5) {
         flipFlag = true;
         p = 1-p;
      }
      if(p < 0.001) continue;
      mean = p*2;
      se = sqrt(2*p*(1-p));
      if(se < 1E-100) continue;

      for(int i = 0; i < dimN; i++){
         k = ID_UN[i];
         g = -1.0;
         if(GG[0][k][byte] & shortbase[offset])  // homozygote
            g = (GG[1][k][byte] & shortbase[offset])? 2.0: 0.0;
         else if(GG[1][k][byte] & shortbase[offset])   // Aa
            g = 1.0;
         if(g > -0.9){
            if(flipFlag)
               X[validSNP][i] = (2 - g - mean) / se;
            else
               X[validSNP][i] = (g - mean) / se;
         }
      }
      for(int i = 0; i < ID.Length(); i++){
         k = ID[i];
         g = -1.0;
         if(GG[0][k][byte] & shortbase[offset])  // homozygote
            g = (GG[1][k][byte] & shortbase[offset])? 2.0: 0.0;
         else if(GG[1][k][byte] & shortbase[offset])   // Aa
            g = 1.0;
         if(g > -0.9) {
            if(flipFlag)
               Z[validSNP][i] = (2 - g - mean) / se;
            else
               Z[validSNP][i] = (g - mean) / se;
         }
      }
      validSNP++;
   }
   dimM = validSNP;
   if(dimM == 0) return; // no SNPs involved
   X.Dimension(dimM, dimN);
   Z.Dimension(dimM, ID.Length());

#ifdef WITH_LAPACK
   char JOBU = 'S';
   char JOBVT = 'O';
   int LDU = dimM;
   int LDVT = 1;
   int dimS = dimM>dimN? dimN: dimM;
   int info;
   double *A = (double *)malloc(dimM*dimN*sizeof(double));
   for(int i = 0; i < dimM; i++)
     for(int j = 0; j < dimN; j++)
       A[j*dimM+i] = X[i][j];
   double *S = (double *)malloc(dimS*sizeof(double));
   double *U = (double *)malloc(LDU*dimS*sizeof(double));
   double *VT = (double *)malloc(LDVT*sizeof(double));
   int LWORK = dimM>dimN? dimM*2+8*dimN: dimN+8*dimM;
   double *WORK = (double *)malloc(LWORK*sizeof(double));
   dgesvd_(&JOBU, &JOBVT, &dimM, &dimN, A, &dimM, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &info);
   if(info!=0) error("SVD failed.");
   free(WORK);
   free(A);
   free(VT);

   int dimPC = 2;

   Matrix EV(ID.Length(), dimPC);
   EV.Zero();
   for(int i = 0; i < ID.Length(); i++)
      for(int j = 0; j < dimPC; j++){
         for(int m = 0; m < dimM; m++)
            EV.data[i]->data[j] += Z.data[m]->data[i] * U[j*dimM+m];
         if(S[j] > 1E-100)
            EV[i][j] /= S[j];
      }
   free(U);
   free(S);
#else
   SVD svd;
   svd.Decompose(X);
   if(svd.n == 0) return;
   QuickIndex idx;
   idx.Index(svd.w);

   int dimPC = 2;
   Matrix EV(ID.Length(), dimPC);
   EV.Zero();
   for(int i = 0; i < ID.Length(); i++)
      for(int j = 0; j < dimPC; j++){
         for(int m = 0; m < dimM; m++)
            EV[i][j] += Z[m][i] * svd.u[m][idx[dimN-1-j]];
         if(svd.w[idx[dimN-1-j]] > 1E-10)
            EV[i][j] /= svd.w[idx[dimN-1-j]];
      }
#endif
   Vector ref0(0), ref1(0);
   for(int i = 0; i < N; i++){
      if(ancestry[i]==0) ref0.Push(EV[i][0]);
      else if(ancestry[i]==1) ref1.Push(EV[i][0]);
   }
   double x0, x1;
   if(ref0.Max() > ref1.Max()) {
      x0 = ref0.Max();
      x1 = ref1.Min();
   }else{
      x0 = ref0.Min();
      x1 = ref1.Max();
   }
   for(int i = 0; i < N; i++){
      if(ancestry[i] == 0 || ancestry[i] == 1) continue;
      ancestry[i] = (EV[i][0] - x0) / (x1 - x0);
   }                                     
}


   /*
   IntArray segIndex(0), segPartialIndex(0), segMask(0);
   for(int m = 0; m < markerCount; m+=16){
      int b = m / 16;
      int mask = 0;
      for(int j = 0; j < 16 && m+j < markerCount; j++)
         if(chromosomes[m+j]==mv_chr && positions[m+j] >= mv_start && positions[m+j] < mv_stop)
            mask |= (1<<j);
      if(mask == 65535) segIndex.Push(b);
      else if(mask !=0){
         segPartialIndex.Push(b);
         segMask.Push(mask);
      }
   } */

/*
   for(int i = 0; i < ped.count; i++){
      if(typed[i] == -1) continue;
      fprintf(fp, "%s %s %s %s %d",
         (const char*)ped[i].famid, (const char*)ped[i].pid,
         (const char*)ped[i].fatid, (const char*)ped[i].motid,
         ped[i].sex);

      if(inSVD[i]!=-1) fprintf(fp, " 1"); // unrelated in SVD
      else fprintf(fp, " 2"); // related in PCA
      for(int j = 0; j < dimPC; j++)
         fprintf(fp, " %.4lf", EV[typed[i]][j]);
   }
*/


/*
   IntArray typed(ped.count);
   typed.Set(-1);
   for(int i = 0; i < ID.Length(); i++)
      typed[ID[i]] = i;
   IntArray inSVD(ped.count);
   inSVD.Set(-1);
   for(int i = 0; i < ID_UN.Length(); i++)
      inSVD[ID_UN[i]] = i;
      */
   /*
   for(int i = 0; i < ped.count; i++)
      if(ped[i].ngeno >= MINSNPCOUNT){
         if(ped[i].affections[0]!=2)
            ID_UN.Push(i);
         ID.Push(i);
      }*/

