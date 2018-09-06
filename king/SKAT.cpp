//////////////////////////////////////////////////////////////////////
// SKAT.cpp
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

#define MAXSNP 50000

void Engine::SKAT_Batch_WeightedLinear(const char* skatfile)
{
/*
   IFILE input = ifopen(skatfile, "rt");
   if(input == NULL)
      error("SKAT parameter file %s cannot be openned.", skatfile);

   StringArray skatGene(0);
   IntArray skatChr(0);
   IntArray skatStart(0);
   IntArray skatStop(0);

   String line;
   StringArray tokens;
   tokens.Clear();
   line.ReadLine(input);
   tokens.AddTokens(line);
   // gene chr start stop
   if(tokens.Length() < 4)
      error("Format for SKAT parameter file: GeneName Chr Start Stop");
   if(tokens[0] == "GENE" || tokens[0]=="gene" || tokens[1]=="CHR" || tokens[1]=="chr"); // skip
   else ifrewind(input);
   while(!ifeof(input)){
      tokens.Clear();
      line.ReadLine(input);
      tokens.AddTokens(line);
      if(tokens.Length() < 4) continue;
      skatGene.Push(tokens[0]);
      if(tokens[1]=="X") skatChr.Push(SEXCHR);
      else if(tokens[1]=="Y") skatChr.Push(SEXCHR+1);
      else if(tokens[1]=="XY") skatChr.Push(SEXCHR+2);
      else if(tokens[1]=="MT") skatChr.Push(SEXCHR+3);
      else skatChr.Push(tokens[1].AsInteger());
      skatStart.Push(tokens[2].AsInteger());
      skatStop.Push(tokens[3].AsInteger());
   }
   ifclose(input);

   int testCount = skatChr.Length();
   if(testCount==0) {
      error("Format for SKAT parameter file: GeneName Chr Start Stop");
      return;
   }
   printf("SKAT with weighted linear kernel.\n");
   if(geno.Length()==0) BuildBinary();

   SVD svd;
   unsigned char AA, Aa, missing;
   int b, id;
   int *isAff = new int[idCount];
   for(int i = 0; i < idCount; i++)
      isAff[i] = -1;
   for(int i = 0; i < ped.count; i++){
      id = geno[i];
      if(id == -1) continue;
      if(ped[i].affections[0]==0) continue;
      isAff[id] = ped[i].affections[0]-1;
   }
   int nid = 0;
   for(int i = 0; i < idCount; i++)
      if(isAff[i]!=-1) nid++;
   if(nid==0) return;

   int *vid = new int[nid];
   int count=0;
   for(int i = 0; i < idCount; i++)
      if(isAff[i]!=-1) {
         vid[count] = i;
         count++;
      }

   int Na = 0;
   for(int i = 0; i < nid; i++)
      Na += isAff[vid[i]];
   int Nu = nid - Na;
   int *affvid = new int[nid];
   for(int i = 0; i < nid; i++)
      affvid[i] = isAff[vid[i]];
   double diseaseCoef = sqrt(Na*Nu*1.0)/nid;

   double weight[MAXSNP];
   double frequencies[MAXSNP];
   Matrix Geno(nid, 1000);
   int pos;
   char *GAA[256], *GAa[256], *GMissing[256];
   for(int i = 0; i < 256; i++){
      GAA[i] = new char[8];
      GAa[i] = new char[8];
      GMissing[i] = new char[8];
      for(int j = 0; j < 8; j++)
         GAA[i][j] = GAa[i][j] = GMissing[i][j] = 0;
   }
   for(int i = 0; i < 256; i++)
      for(int j = 0; j < 8; j++)
         if(i & base[j]){
            GAA[i][j] = 2;
            GAa[i][j] = 1;
            GMissing[i][j] = -1;
         }

   String skatout = prefix;
   skatout.Add("skat.tbl");
   FILE *fp = fopen((const char*)skatout, "wt");
   if(fp==NULL) error("Cannot open %s to write.", (const char*)skatout);
   fprintf(fp, "Gene\tChr\tStart\tStop\tSNPCount\tRareCount\tChisq\tDistribution\tP\n");

#ifdef WITH_LAPACK
   printf("  LAPACK is used.\n");
#endif
   printf("Scanning regions/genes");
   double numerator, denominator;
   IntArray vbyte, vbit;
   int validtestCount = 0;
//   Matrix ZPZ;
   Matrix PK;
   double pvalue;
   Vector lambda;
   Vector YG;

   for(int t = 0; t < testCount; t++){ // loop over genes/regions
      Geno.Zero();
      vbyte.Dimension(0);
      vbit.Dimension(0);
      int indexCount = 0;
      for(int m = 0; m < positions.Length(); m++)
         if(chromosomes[m] == skatChr[t] && positions[m] >= skatStart[t]*0.000001
            && positions[m] <= skatStop[t]*0.000001){
            int proposedByte = m/8;
            int proposedBit = m%8;
            if(vbyte.Find(proposedByte)>-1){ // in previous byte
               vbit[vbyte.Find(proposedByte)] |= (1<<proposedBit);
            }else{   // a new byte
               vbyte.Push(proposedByte);
               vbit.Push(1<<proposedBit);
            }
            indexCount ++;
         }
      if(indexCount==0 || indexCount > MAXSNP) continue;
      if(indexCount>1000) Geno.Dimension(nid, indexCount);
      int cIndex = 0;
      int offset;
      int rareCount = 0;
      for(int v = 0; v < vbyte.Length(); v++){ // loop over vbyte
         b = vbyte[v];
         for(int i = 0; i < nid; i++){
            id = vid[i];
            AA = G[0][id][b] & G[1][id][b];
            Aa = (~G[0][id][b]) & G[1][id][b];
            missing = (~G[0][id][b]) & (~G[1][id][b]);
            offset = 0;
            for(int j = 0; j < 8; j++){
               if(vbit[v] & base[j]){
                  pos = cIndex+offset;
                  Geno[i][pos] += GAA[AA][j];
                  Geno[i][pos] += GAa[Aa][j];
                  Geno[i][pos] += GMissing[missing][j];
                  offset ++;
               }
            }
         }
         offset = 0;
         for(int j = 0; j < 8; j++)
            if(vbit[v] & base[j]){
               pos = cIndex+offset;
               numerator = denominator = 0;
               for(int i = 0; i < nid; i++){
                  if(Geno[i][pos]==2) {
                     numerator+=2;
                     denominator+=2;
                  }else if(Geno[i][pos]==1) {
                     numerator+=1;
                     denominator+=2;
                  }else if(Geno[i][pos]==0)
                     denominator+=2;
               }
               if(denominator > 0){
                  frequencies[pos] = numerator / denominator;
                  if(frequencies[pos] < 0.5)
                     weight[pos] = pow(1-frequencies[pos], 24)*25; // sqrt of weight
                  else
                     weight[pos] = pow(frequencies[pos], 24)*25;
               }else{
                  frequencies[pos] = 0;
                  weight[pos] = 0;
               }
               for(int i = 0; i < nid; i++)
                  if(Geno[i][pos] == -1)
                     Geno[i][pos] = frequencies[pos]*2;
               offset++;
               if(frequencies[pos]<0.05 || frequencies[pos]>0.95) rareCount++;
            }
         cIndex += offset;
      } // loop over vbyte


      double stat = 0;
      for(int i = 0; i < nid; i++)
         for(int j = 0; j < indexCount; j++)
            Geno.data[i]->data[j] -= frequencies[j]*2;

      YG.Dimension(indexCount);
      YG.Zero();
      for(int j = 0; j < indexCount; j++)
         for(int i = 0; i < nid; i++)
            if(affvid[i]) YG[j] += Geno[i][j];
      for(int j = 0; j < indexCount; j++)
         YG[j] *= weight[j];
      for(int j = 0; j < indexCount; j++)
         stat += YG[j] * YG[j];

      PK.Dimension(nid, indexCount);
      for(int i = 0; i < nid; i++)
         for(int j = 0; j < indexCount; j++)
            PK.data[i]->data[j] = Geno.data[i]->data[j] * weight[j] * diseaseCoef;
   lambda.Dimension(0);
#ifdef WITH_LAPACK
   if(indexCount == 1){
      double temp = 0;
      for(int i = 0; i < nid; i++)
         temp += PK[i][0] * PK[i][0];
      if(temp > 0.0001) lambda.Push(temp);
   }else if(indexCount > 1){
      int dimM = nid;
      int dimN = indexCount;
      char JOBU = 'S';
      char JOBVT = 'O';
      int LDU = dimM;
      int LDVT = 1;
      int dimS = dimM>dimN? dimN: dimM;
      int info;
      double *A = (double *)malloc(dimM*dimN*sizeof(double));
      for(int i = 0; i < dimM; i++)
         for(int j = 0; j < dimN; j++)
            A[j*dimM+i] = PK.data[i]->data[j];
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
      free(U);
      for(int i = 0; i < dimS; i++)
         if(S[i] > 0.01) lambda.Push(S[i]*S[i]);
      free(S);
   }else
      return;
#else
   if(indexCount == 1){
      double temp = 0;
      for(int i = 0; i < nid; i++)
         temp += PK[i][0] * PK[i][0];
      if(temp > 0.0001) lambda.Push(temp);
   }else if(indexCount > 1){
      svd.Decompose(PK);
      if(svd.n == 0) return;
      for(int i = 0; i < (nid>indexCount? indexCount: nid); i++)
         if(svd.w[i] > 0.01) lambda.Push(svd.w[i]*svd.w[i]);
   }else
      continue;
#endif
      if(stat>0 && lambda.Length()>0){
         if(lambda.Length()==1)
           pvalue = chidist(stat/lambda[0], 1);
         else
            pvalue = 1-Davies(stat, lambda.data, lambda.Length());
      }else
         pvalue = 1;
      if(pvalue<0) pvalue=0;
      if(lambda.Length()>0)
      fprintf(fp, "%s\t%d\t%d\t%d\t%d\t%d\t%.3lf\t%.1lf",
         (const char*)skatGene[t], skatChr[t], skatStart[t], skatStop[t],
            indexCount, rareCount, stat, lambda[0]);
      else
      fprintf(fp, "%s\t%d\t%d\t%d\t%d\t%d\t%.3lf\t%.1lf",
         (const char*)skatGene[t], skatChr[t], skatStart[t], skatStop[t],
            indexCount, rareCount, stat, 0);
      for(int k = 1; k < lambda.Length(); k++)
         fprintf(fp, ",%.1lf", lambda[k]);
      fprintf(fp, "\t%.3G\n", pvalue);
      validtestCount ++;
      printf(".");
   } // loop over genes/regions

   delete []affvid;
   delete []isAff;
   delete []vid;
   for(int i = 0; i < 256; i++){
      delete []GAA[i];
      delete []GAa[i];
      delete []GMissing[i];
   }
   printf("\nSKAT results for %d genes/regions saved in file %s\n",
      validtestCount, (const char*)skatout);
*/
}

/*
   if(indexCount == 1){
      if(ZPZ[0][0] > 0.0001)
         lambda.Push(ZPZ[0][0]);
   }else if(indexCount > 1){
      int dimN = indexCount;
      char JOBZ = 'A';
      int info;
      double *A = new double[dimN*dimN];
      int dimLA = dimN;
      for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimN; j++)
          A[i*dimN+j] = ZPZ[j][i];
      double *S = new double[dimN];
      double *U = new double[dimN*dimN];
      double *VT = new double[dimN*dimN];
      int *IWORK = new int[dimN*8];
      int LWORK = 8*dimN + 4*dimN*dimN;
      double *WORK = new double[LWORK];
      dgesdd_(&JOBZ, &dimLA, &dimLA, A, &dimLA, S, U, &dimLA, VT, &dimLA, WORK, &LWORK, IWORK, &info);
      delete []U;
      delete []IWORK;
      delete []WORK;
      delete []A;
      delete []VT;
      for(int i = 0; i < dimN; i++)
         if(S[i] > 0.0001) lambda.Push(S[i]);
      delete []S;
   }else
      continue;
*/
/*
   if(indexCount == 1){
      if(ZPZ[0][0] > 0.0001)
         lambda.Push(ZPZ[0][0]);
   }else if(indexCount > 1){
      svd.Decompose(ZPZ);
      if(svd.n == 0) return;
      for(int i = 0; i < indexCount; i++)
         if(svd.w[i] > 0.001) lambda.Push(svd.w[i]);
   }else
      continue;
*/

void Engine::SKAT_WeightedLinear()
{
/*
   printf("SKAT with weighted linear kernel.\n");
   if(geno.Length()==0) BuildBinary();
   if(byteCount > 10000) error("Too many SNPs included.\n");

   SVD svd;
   unsigned char AA, Aa, missing;
   int b, id;
   int *isAff = new int[idCount];
   for(int i = 0; i < idCount; i++)
      isAff[i] = -1;
   for(int i = 0; i < ped.count; i++){
      id = geno[i];
      if(id == -1) continue;
      if(ped[i].affections[0]==0) continue;
      isAff[id] = ped[i].affections[0]-1;
   }
   int nid = 0;
   for(int i = 0; i < idCount; i++)
      if(isAff[i]!=-1) nid++;
   int *vid = new int[nid];
   int count=0;
   for(int i = 0; i < idCount; i++)
      if(isAff[i]!=-1) {
         vid[count] = i;
         count++;
      }
   Vector weight(byteCount*8);
   Vector frequencies(byteCount*8);
   Matrix Geno(nid, byteCount*8);
   Geno.Zero();
   int pos;
   char *GAA[256], *GAa[256], *GMissing[256];
   for(int i = 0; i < 256; i++){
      GAA[i] = new char[8];
      GAa[i] = new char[8];
      GMissing[i] = new char[8];
      for(int j = 0; j < 8; j++)
         GAA[i][j] = GAa[i][j] = GMissing[i][j] = 0;
   }
   for(int i = 0; i < 256; i++)
      for(int j = 0; j < 8; j++)
         if(i & base[j]){
            GAA[i][j] = 2;
            GAa[i][j] = 1;
            GMissing[i][j] = -1;
         }

   double numerator, denominator;
   for(int b = 0; b < byteCount; b++){
      for(int i = 0; i < nid; i++){
         id = vid[i];
         AA = G[0][id][b] & G[1][id][b];
         Aa = (~G[0][id][b]) & G[1][id][b];
         missing = (~G[0][id][b]) & (~G[1][id][b]);
         for(int j = 0; j < 8; j++){
            Geno[i][b*8+j] += GAA[AA][j];
            Geno[i][b*8+j] += GAa[Aa][j];
            Geno[i][b*8+j] += GMissing[missing][j];
         }
      }
      for(int j = 0; j < 8; j++){
         pos = b*8+j;
         numerator = denominator = 0;
         for(int i = 0; i < nid; i++){
            if(Geno[i][pos]==2) {
               numerator+=2;
               denominator+=2;
            }else if(Geno[i][pos]==1) {
               numerator+=1;
               denominator+=2;
            }else if(Geno[i][pos]==0)
               denominator+=2;
         }
         if(denominator > 0){
            frequencies[pos] = numerator / denominator;
            if(frequencies[pos] < 0.5)
               weight[pos] = pow(1-frequencies[pos], 24)*25; // sqrt of weight
            else
               weight[pos] = pow(frequencies[pos], 24)*25;
          }else{
            frequencies[pos] = 0;
            weight[pos] = 0;
         }
         for(int i = 0; i < nid; i++)
            if(Geno[i][pos] == -1)
               Geno[i][pos] = frequencies[pos]*2;
      }
   }
   int Na = 0;
   for(int i = 0; i < nid; i++)
      Na += isAff[vid[i]];
   int Nu = nid - Na;
   int *affvid = new int[nid];
   for(int i = 0; i < nid; i++)
      affvid[i] = isAff[vid[i]];
   double stat = 0;
   for(int i = 0; i < nid; i++)
      for(int j = 0; j < markerCount; j++)
         Geno[i][j] -= frequencies[j]*2;

   Vector YG(markerCount);
   YG.Zero();
   for(int j = 0; j < markerCount; j++)
      for(int i = 0; i < nid; i++)
         if(affvid[i]) YG[j] += Geno[i][j];
   for(int j = 0; j < markerCount; j++)
      YG[j] *= weight[j];
   for(int j = 0; j < markerCount; j++)
      stat += YG[j] * YG[j];

   Matrix ZPZ(markerCount, markerCount);
   ZPZ.Zero();
   double coef = Na*Nu*1.0/(nid*nid);
   for(int i = 0; i < markerCount; i++)
      for(int j = 0; j < markerCount; j++)
         for(int k = 0; k < nid; k++)
            ZPZ.data[i]->data[j] += Geno.data[k]->data[i] * weight.data[i] *
            Geno.data[k]->data[j] * weight.data[j] * coef;

   printf("%d affected and %d unaffected\n", Na, Nu);
   printf("%20s%10s%10s\n", "SNP", "MAF", "Weight");
   for(int j = 0; j < markerCount; j++)
      printf("%20s%10.3lf%10.3lf\n",
         (const char*)snpName[j], frequencies[j], weight[j]);
   printf("Observed statistic Q = %.3lf\n", stat);

   double pvalue;
   Vector lambda(0);
#ifdef WITH_LAPACK
   if(markerCount == 1){
      if(ZPZ[0][0] > 0.0001)
         lambda.Push(ZPZ[0][0]);
   }else if(markerCount > 1){
      int dimN = markerCount;
      char JOBZ = 'A';
      int info;
      double *A = new double[dimN*dimN];
      int dimLA = dimN;
      for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimN; j++)
          A[i*dimN+j] = ZPZ[j][i];
      double *S = new double[dimN];
      double *U = new double[dimN*dimN];
      double *VT = new double[dimN*dimN];
      int *IWORK = new int[dimN*8];
      int LWORK = 8*dimN + 4*dimN*dimN;
      double *WORK = new double[LWORK];
      dgesdd_(&JOBZ, &dimLA, &dimLA, A, &dimLA, S, U, &dimLA, VT, &dimLA, WORK, &LWORK, IWORK, &info);
      delete U;
      delete IWORK;
      delete WORK;
      delete A;
      delete VT;
      for(int i = 0; i < dimN; i++)
         if(S[i] > 0.0001) lambda.Push(S[i]);
      delete S;
   }else
      return;
#else
   if(markerCount == 1){
      lambda.Push(ZPZ[0][0]);
   }else if(markerCount > 1){
      svd.Decompose(ZPZ);
      if(svd.n == 0) return;
      for(int i = 0; i < markerCount; i++)
         if(svd.w[i] > 0.0001) lambda.Push(svd.w[i]);
   }else // no marker
      return;
#endif
   printf("The null distribution of Q is a mixture of %d independent chisq(1) with weights\n",
      lambda.Length());
   for(int i = 0; i < lambda.Length(); i++)
      printf(" %.2lf", lambda[i]);
   printf("\n");

   if(stat>0 && lambda.Length()>0){
      if(lambda.Length()==1)
         pvalue = chidist(stat/lambda[0], 1);
      else
         pvalue = 1-Davies(stat, lambda.data, lambda.Length());
   }else
      pvalue = 1.0;
   printf("Analytical (Davies) P-value = %G\n", pvalue);

   delete []isAff;
   delete []vid;
   delete []affvid;
   for(int i = 0; i < 256; i++){
      delete []GAA[i];
      delete []GAa[i];
      delete []GMissing[i];
   }
   */
}
/*
   int dimM = nid;
   int dimN = markerCount;
   printf("  LAPACK is used.\n");
   char JOBU = 'S';
   char JOBVT = 'O';
   int LDU = dimM;
   int LDVT = 1;
   int dimS = dimM>dimN? dimN: dimM;
   int info;
   double *A = (double *)malloc(dimM*dimN*sizeof(double));
   for(int i = 0; i < dimM; i++)
     for(int j = 0; j < dimN; j++)
       A[j*dimM+i] = PK[i][j];
   double *S = (double *)malloc(dimS*sizeof(double));
   double *U = (double *)malloc(LDU*dimS*sizeof(double));
   double *VT = (double *)malloc(LDVT*sizeof(double));
   int LWORK = dimM>dimN? dimM*2+8*dimN: dimN+8*dimM;
   double *WORK = (double *)malloc(LWORK*sizeof(double));
   dgesvd_(&JOBU, &JOBVT, &dimM, &dimN, A, &dimM, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &info);
   if(info!=0) error("SVD failed.");
   delete WORK;
   delete A;
   delete VT;
   delete U;
   for(int i = 0; i < dimS; i++)
      if(S[i] > 0.01) lambda.Push(S[i]*S[i]);
   delete S;
      */
   /*
   svd.Decompose(PK);
   if(svd.n == 0) return;
   for(int i = 0; i < (nid>markerCount? markerCount: nid); i++)
      if(svd.w[i] > 0.01) lambda.Push(svd.w[i]*svd.w[i]);
      */

void Engine::SKAT()
{
/*
   int nperm = permuteCount;
   if(ibsFlag)
      printf("SKAT with IBS kernel.\n");
   else
      printf("SKAT with Euclidean distance kernel.\n");
   int missingCount[8], AACount[8], AaCount[8];
   unsigned char AA, Aa, missing;
   int b, id1, id2;
   double D[3];
   if(geno.Length()==0) BuildBinary();

   IntArray bmask(0);
   for(int b = 0; b < byteCount; b++)
      if( (b < byteCount-1) || (markerCount%8==0) )
         bmask.Push(255);
      else
         bmask.Push((1<<(markerCount%8))-1);
   double frequencies[8];
   int *isAff = new int[idCount];
   for(int i = 0; i < idCount; i++)
      isAff[i] = -1;
   for(int i = 0; i < ped.count; i++){
      id1 = geno[i];
      if(id1 == -1) continue;
      if(ped[i].affections[0]==0) continue;
      isAff[id1] = ped[i].affections[0]-1;
   }

   int nid = 0;
   for(int i = 0; i < idCount; i++)
      if(isAff[i]!=-1) nid++;
   int *vid = new int[nid];
   int count=0;
   for(int i = 0; i < idCount; i++)
      if(isAff[i]!=-1) {
         vid[count] = i;
         count++;
      }
   double *WS[256], *WM[256], weight[8];
   int tmask[8];
   for(int i = 0; i < 256; i++){
      WS[i] = new double[byteCount];
      WM[i] = new double[byteCount];
   }
   for(int b = 0; b < byteCount; b++){
      for(int i = 0; i < 8; i++){
         AACount[i] = AaCount[i] = missingCount[i] = 0;
         frequencies[i] = 0.0;
      }
      for(int i = 0; i < idCount; i++){
         AA = G[0][i][b] & G[1][i][b];
         Aa = (~G[0][i][b]) & G[1][i][b];
         missing = (~G[0][i][b]) & (~G[1][i][b]);
         for(int j = 0; j < 8; j++){
            if(AA & base[j])
               AACount[j] ++;
            else if(Aa & base[j])
               AaCount[j] ++;
            else if(missing & base[j])
               missingCount[j] ++;
         }
      }
      for(int i = 0; i < 256; i++)
         WS[i][b] = WM[i][b] = 0;
      for(int j = 0; j < 8; j++){
         if(missingCount[j] < idCount){
            frequencies[j] = (AACount[j] + AaCount[j]*0.5) / (idCount - missingCount[j]);
            if(frequencies[j] > 0.5) frequencies[j] = 1-frequencies[j];
            weight[j] = pow(1-frequencies[j], 48)*625;
         }else{
            frequencies[j] = 0;
            weight[j] = 0;
         }
         tmask[j] = base[j] & bmask[b];
         for(int i = 0; i < 256; i++)
            if(i & tmask[j]){
               WS[i][b] += weight[j];
               WM[i][b] += weight[j]*4*frequencies[j]*(1-frequencies[j]);
            }
      }
   }
   int Na = 0;
   for(int i = 0; i < nid; i++)
      Na += isAff[vid[i]];
   int Nu = nid - Na;
   int notMissingCount[3];
   notMissingCount[0] = Nu*(Nu-1)/2;
   notMissingCount[1] = Na*Nu;
   notMissingCount[2] = Na*(Na-1)/2;
   double DistCount[3];
   for(int i = 0; i < 3; i++)
      DistCount[i] = 0;

   double *pDistCount[3];
   for(int i = 0; i < 3; i++){
      pDistCount[i] = new double[nperm];
      for(int j = 0; j < nperm; j++)
         pDistCount[i][j] = 0;
   }

   Matrix Dist(nid, nid);
   Dist.Zero();
   int *affvid = new int[nid];
   for(int i = 0; i < nid; i++)
      affvid[i] = isAff[vid[i]];

   int T = nid / 10;
   if(T==0) T=1;
   double tD;
   printf("\nScan");

   if(ibsFlag) // IBS kernel
   for(int i = 0; i < nid; i++){
      id1 = vid[i];
      if(i % T == 0) printf(".");
      for(int j = i+1; j < nid; j++){
         id2 = vid[j];
         tD = 0;
         for(int b = 0; b < byteCount; b++){
            missing = ~((G[0][id1][b] | G[1][id1][b]) & (G[0][id2][b] | G[1][id2][b]));
            tD +=
            2*WS[G[0][id1][b] & G[0][id2][b] & (G[1][id1][b] ^ G[1][id2][b])][b]
            + WS[(G[0][id2][b] & (~G[0][id1][b]) & G[1][id1][b]) |
                  (G[0][id1][b] & (~G[0][id2][b]) & G[1][id2][b])][b]
            + WM[missing][b];         // Dist = IBS0*2 + IBS1
         }
         DistCount[affvid[i]+affvid[j]] += tD;
         Dist.data[i]->data[j] = tD;
      }
   }
   else // Euclidean Distance kernel
   for(int i = 0; i < nid; i++){
      id1 = vid[i];
      if(i % T == 0) printf(".");
      for(int j = i+1; j < nid; j++){
         id2 = vid[j];
         tD = 0;
         for(int b = 0; b < byteCount; b++){
            missing = ~((G[0][id1][b] | G[1][id1][b]) & (G[0][id2][b] | G[1][id2][b]));
            tD +=
            4*WS[G[0][id1][b] & G[0][id2][b] & (G[1][id1][b] ^ G[1][id2][b])][b]
            + WS[(G[0][id2][b] & (~G[0][id1][b]) & G[1][id1][b]) |
                  (G[0][id1][b] & (~G[0][id2][b]) & G[1][id2][b])][b]
            + WM[missing][b];         // Dist = IBS0*4 + IBS1
         }
         DistCount[affvid[i]+affvid[j]] += tD;
         Dist.data[i]->data[j] = tD;
      }
   }
   printf("done\n");

   StringArray sCat(3);
   sCat[0]="Unaffected"; sCat[1]="Discordant"; sCat[2] = "Affected";
   printf("\n%15s%15s%15s%10s\n", "Category", "Dist", "Count", "AveDist");
   for(int c = 0; c < 3; c++){
      D[c] = DistCount[c] / notMissingCount[c];
      printf("%15s%15.1lf%15d%10.4lf\n",
         (const char*)sCat[c], DistCount[c], notMissingCount[c], D[c]);
   }
   printf("\n");

   Matrix PKP = Dist;
   for(int i = 0; i < nid; i++)
      for(int j = i+1; j < nid; j++)
         PKP[j][i] = Dist[i][j];
   PKP *= (-Na*Nu*1.0/nid/nid);
   Vector mean(nid);
   mean.Zero();
   for(int i = 0; i < nid; i++){
      for(int j = 0; j < nid; j++)
         mean[i] += PKP[i][j];
      mean[i] /= nid;
   }
   double mean0 = mean.Sum() / nid;
   for(int i = 0; i < nid; i++)
      for(int j = 0; j < nid; j++)
         PKP[i][j] = PKP[i][j] - mean[i] - mean[j] + mean0;

   double stat = (D[1] - ((Na-1.0)/Na*D[0]+(Nu-1.0)/Nu*D[2])*0.5)*Na*Na*Nu*Nu*2.0/nid/nid;
   printf("Observed statistic Q = %.4lf\n", stat);

   double pvalue;
   Vector lambda(0);

#ifdef WITH_LAPACK
   printf("  LAPACK is used.\n");
   char JOBZ = 'A';
   int info;
   double *A = new double[nid*nid];
   int dimLA = nid;
   for(int i = 0; i < nid; i++)
     for(int j = 0; j < nid; j++)
       A[i*nid+j] = PKP[j][i];
   double *S = new double[nid];
   double *U = new double[nid*nid];
   double *VT = new double[nid*nid];
   int *IWORK = new int[nid*8];
   int LWORK = 8*nid + 4*nid*nid;
   double *WORK = new double[LWORK];
   dgesdd_(&JOBZ, &dimLA, &dimLA, A, &dimLA, S, U, &dimLA, VT, &dimLA, WORK, &LWORK, IWORK, &info);

   delete U;
   delete IWORK;
   delete WORK;
   delete A;

   for(int i = 0; i < nid; i++)
      if(S[i] > 0.01) lambda.Push(S[i]);

#else
   SVD svd;
   svd.Decompose(PKP);
   if(svd.n == 0) return;
   for(int i = 0; i < nid; i++)
      if(svd.w[i] > 0.01) lambda.Push(svd.w[i]);
#endif
   printf("The null distribution of Q is a mixture of %d independent chisq(1) with weights\n",
      lambda.Length());
   for(int i = 0; i < lambda.Length(); i++){
      printf(" %.1lf", lambda[i]);
   }
   printf("\n");
   pvalue = 1-Davies(stat, lambda.data, lambda.Length());
   printf("Analytical (Davies) P-value = %G\n", pvalue);

   if(nperm==0){
      for(int i = 0; i < 256; i++){
         delete []WS[i];
         delete []WM[i];
      }
      delete []isAff;
      delete []vid;
      delete []affvid;
      return;
   }

   printf("%d permutations start", nperm);
   Random rand;
   int *perm = new int[nid];
   int pos, swap;
   for(int k = 0; k < nperm; k++){
      if(k%100==0) printf(".");
      for(int j = 0; j < nid; j++) perm[j] = j;
      for(int j = nid; j > 1; j--){
         pos = rand.NextInt()%j;
         swap = perm[j-1];
         perm[j-1] = perm[pos];
         perm[pos] = swap;
      }
      for(int i = 0; i < nid; i++){
         count = affvid[perm[i]];
         for(int j = i+1; j < nid; j++) // computationally intensive
            pDistCount[count + affvid[perm[j]]][k] += Dist.data[i]->data[j];
      }
   }
   printf("done\n");
   delete []perm;
   for(int i = 0; i < 256; i++){
      delete []WS[i];
      delete []WM[i];
   }
   delete []isAff;
   delete []vid;
   delete []affvid;

   Vector pstat(nperm);
   String permfile = prefix;
   permfile.Add("perm.txt");
   FILE *fp = fopen(permfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)permfile);
   fprintf(fp, "N\tD_UU\tD_Au\tD_AA\tStat\n");
   fprintf(fp, "%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\n",
      0, D[0], D[1], D[2], stat);
   for(int k = 0; k < nperm; k++){
      for(int c = 0; c < 3; c++)
         D[c] = pDistCount[c][k] / notMissingCount[c];
      pstat[k] = (D[1] - ((Na-1.0)/Na*D[0]+(Nu-1.0)/Nu*D[2])*0.5)*Na*Na*Nu*Nu*2.0/nid/nid;
      fprintf(fp, "%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\n",
         k+1, D[0], D[1], D[2], pstat[k]);
   }
   fclose(fp);
   printf("\nPermutation results saved in file %s\n", (const char*)permfile);
   for(int i = 0; i < 3; i++)
      delete []pDistCount[i];

   QuickIndex idx;
   idx.Index(pstat);
   for(pos = 0; pos < nperm && pstat[idx[pos]] < stat; pos++);
   pvalue = 1-pos * 1.0 / nperm;
   printf("Permuted quantiles: %.4lf %.4lf %.4lf %.4lf %.4lf\n",
      pstat[idx[0]], pstat[idx[int(nperm*0.25-0.5)]], pstat[idx[int(nperm*0.5-0.5)]],
      pstat[idx[int(nperm*0.75-0.5)]], pstat[idx[nperm-1]]);
   printf("P value is %G\n", pvalue);
*/
}

/*
void Engine::SKAT_family_quantitative(int trait)
{
   printf("SKAT with weighted linear kernel for quantitative traits in families.\n");
   if(geno.Length()==0) BuildBinary();

   int nid = 0;
   for(int i = 0; i < ped.count; i++)
      if(ped[i].traits[trait]!=_NAN_ && CheckCovariates(ped[i]) && geno[i]!=-1)
         nid++;
   int *vid = new int[nid];
   double *resid = new double[nid];
   Vector qtrait(nid);
   Matrix qcov(nid, covariates.Length()+1);
   qtrait.Zero();
   qcov.Zero();
   int count = 0;

   for(int f = 0; f < ped.familyCount; f++)
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
      if(ped[i].traits[trait]!=_NAN_ && CheckCovariates(ped[i]) && geno[i]!=-1){
         vid[count] = i;
         qtrait[count] = ped[i].traits[trait];
         qcov[count][1] = 1;
         for(int j = 0; j < covariates.Length(); j++)
            qcov[count][j+1] = ped[i].covariates[covariates[j]];
         count++;
       }

   for(int i = 0; i < nid; i++){
      resid[i] = means[trait][0];
      for(int j = 0; j < covariates.Length(); j++)
         resid[i] -= means[trait][j+1] * qtrait[i];
   }
   unsigned char AA, Aa, missing;
   int b, id;
   Vector weight(byteCount*8);
   Vector frequencies(byteCount*8);
   Matrix Geno(nid, byteCount*8);
   Geno.Zero();
   int pos;
   char *GAA[256], *GAa[256], *GMissing[256];
   for(int i = 0; i < 256; i++){
      GAA[i] = new char[8];
      GAa[i] = new char[8];
      GMissing[i] = new char[8];
      for(int j = 0; j < 8; j++)
         GAA[i][j] = GAa[i][j] = GMissing[i][j] = 0;
   }
   for(int i = 0; i < 256; i++)
      for(int j = 0; j < 8; j++)
         if(i & base[j]){
            GAA[i][j] = 2;
            GAa[i][j] = 1;
            GMissing[i][j] = -1;
         }

   double numerator, denominator;
   for(int b = 0; b < byteCount; b++){
      for(int i = 0; i < nid; i++){
         id = vid[i];
         AA = G[0][id][b] & G[1][id][b];
         Aa = (~G[0][id][b]) & G[1][id][b];
         missing = (~G[0][id][b]) & (~G[1][id][b]);
         for(int j = 0; j < 8; j++){
            Geno[i][b*8+j] += GAA[AA][j];
            Geno[i][b*8+j] += GAa[Aa][j];
            Geno[i][b*8+j] += GMissing[missing][j];
         }
      }
      for(int j = 0; j < 8; j++){
         pos = b*8+j;
         numerator = denominator = 0;
         for(int i = 0; i < nid; i++){
            if(Geno[i][pos]==2) {
               numerator+=2;
               denominator+=2;
            }else if(Geno[i][pos]==1) {
               numerator+=1;
               denominator+=2;
            }else if(Geno[i][pos]==0)
               denominator+=2;
         }
         if(denominator > 0){
            frequencies[pos] = numerator / denominator;
            if(frequencies[pos] < 0.5)
               weight[pos] = pow(1-frequencies[pos], 24)*25; // sqrt of weight
            else
               weight[pos] = pow(frequencies[pos], 24)*25;
          }else{
            frequencies[pos] = 0;
            weight[pos] = 0;
         }
         for(int i = 0; i < nid; i++)
            if(Geno[i][pos] == -1)
               Geno[i][pos] = frequencies[pos]*2;
      }
   }
   for(int i = 0; i < nid; i++)
      for(int j = 0; j < markerCount; j++)
         Geno[i][j] -= frequencies[j]*2;

   Kinship kin;
   Matrix Omega;
   IntArray ids;
   Cholesky chol;
   int N;

   Matrix kernel(nid, markerCount);
   int baseN = 0;
   for(int f = 0; f < ped.familyCount; f++){
      ids.Dimension(0);
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(ped[i].traits[trait]!=_NAN_ && CheckCovariates(ped[i]) && geno[i]!=-1)
            ids.Push(i);
      N = ids.Length();
      if(N==0) continue;
      kin.Setup(*ped.families[f]);
      for(int i = 0; i < N; i++)
         for(int j = i; j < N; j++) // Upper triangle matrix
            Omega[i][j] = variances[trait] -
            (1-2*kin(ped[ids[i]], ped[ids[j]]))*heritabilities[trait]*variances[trait];
      chol.Decompose(Omega);

      Vector temp;
      for(int j = 0; j < markerCount; j++){
         for(int i = 0; i < N; i++){
            temp.Dimension(N);
            temp[i] = Geno[ids[i]][j] * weight[j];
         }
         chol.BackSubst0(temp);
         for(int i = 0; i < N; i++)
            kernel[baseN+i][j] = chol.x[i];
      }
   }

   Vector lambda(0);
   Matrix svdU(nid, markerCount);
#ifdef WITH_LAPACK
   int dimM = nid;
   int dimN = markerCount;
   printf("  LAPACK is used.\n");
   char JOBU = 'S';
   char JOBVT = 'O';
   int LDU = dimM;
   int LDVT = 1;
   int dimS = dimM>dimN? dimN: dimM;
   int info;
   double *A = (double *)malloc(dimM*dimN*sizeof(double));
   for(int i = 0; i < dimM; i++)
     for(int j = 0; j < dimN; j++)
       A[j*dimM+i] = kernel[i][j];
   double *S = (double *)malloc(dimS*sizeof(double));
   double *U = (double *)malloc(LDU*dimS*sizeof(double));
   double *VT = (double *)malloc(LDVT*sizeof(double));
   int LWORK = dimM>dimN? dimM*2+8*dimN: dimN+8*dimM;
   double *WORK = (double *)malloc(LWORK*sizeof(double));
   dgesvd_(&JOBU, &JOBVT, &dimM, &dimN, A, &dimM, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &info);
   if(info!=0) error("SVD failed.");
   delete WORK;
   delete A;
   delete VT;
   count = 0;
   for(int i = 0; i < dimS; i++){
      if(S[i] > 0.01) lambda.Push(S[i]*S[i]);
      for(int k = 0; k < nid; k++)
         svdU[k][count] = U[k*nid+i];
      count++;
   }
#else
   SVD svd;
   svd.Decompose(kernel);
   if(svd.n == 0) return;
   count = 0;
   for(int i = 0; i < (nid>markerCount? markerCount: nid); i++){
      if(svd.w[i] > 0.01) lambda.Push(svd.w[i]*svd.w[i]);
      for(int k = 0; k < nid; k++)
         svdU[k][count] = svd.u[k][i];
      count++;
   }
#endif

   count = 0;
   double temp1, temp2;
   double Q = 0;
   for(int f = 0; f < ped.familyCount; f++){
      ids.Dimension(0);
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(ped[i].traits[trait]!=_NAN_ && CheckCovariates(ped[i]) && geno[i]!=-1)
            ids.Push(i);
      N = ids.Length();
      if(N==0) continue;
      if(N==1){
         for(int m = 0; m < lambda.Length(); m++){
            temp1 = resid[count] * svdU[count][m] / variances[m];
            Q += temp1 * temp1 * lambda[m];
         }
         continue;
      }
      kin.Setup(*ped.families[f]);
      for(int i = 0; i < N; i++)
         for(int j = i; j < N; j++) // Upper triangle matrix
            Omega[i][j] = variances[trait] -
            (1-2*kin(ped[ids[i]], ped[ids[j]]))*heritabilities[trait]*variances[trait];
      chol.Decompose(Omega);
      chol.Invert();
      for(int m = 0; m < lambda.Length(); m++){
         temp1 = 0;
         for(int i = 0; i < N; i++)
            for(int j = 0; j < N; j++)
               temp1 += resid[count+i] * chol.inv[i][j] * Geno[count+j][m];
         Q += temp1*temp1*lambda[m];
      }
      count += ids.Length();
   }

   printf("The null distribution of Q is a mixture of %d independent chisq(1) with weights\n",
      lambda.Length());
   for(int i = 0; i < lambda.Length(); i++)
      printf(" %.1lf", lambda[i]);
   printf("\n");

   double pvalue = 1-Davies(Q, lambda.data, lambda.Length());
   printf("Analytical (Davies) P-value = %G\n", pvalue);

   delete []vid;
   delete []resid;
   for(int i = 0; i < 256; i++){
      delete []GAA[i];
      delete []GAa[i];
      delete []GMissing[i];
   }
}

*/

