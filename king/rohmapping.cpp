//////////////////////////////////////////////////////////////////////
// rohmapping.cpp
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

#include <math.h>
#include <time.h>
#include "analysis.h"
#include "Random.h"
#include "MathCholesky.h"
#include "MathStats.h"
#include "QuickIndex.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

void Engine::HomozygosityMappingForQTMH(const char *popName)
{
   printf("\nOptions in effect:\n");
   printf("\t--mthomo\n");
   printf("\t--strat %s\n", popName);
   if(Bit64Flag)
      printf("\t--sysbit 64\n");
   if(CoreCount)
      printf("\t--cpus %d\n", CoreCount);
   if(lessmemFlag)
      printf("\t--lessmem\n");
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");

   const double MINLOD=1.0;
   const int MINROHCOUNT=5;
   const int MAXROHLENGTH=10000000;
   bool IBDvalidFlag = PreSegment();
   if(!IBDvalidFlag){
      printf("%s\n", (const char*)segmessage);
      printf("  Note chromosomal positions can be sorted conveniently using other tools such as PLINK.\n");
      return;
   }

   int traitCount = traits.Length();
   if(traitCount==0){
      printf("No quantitative trait data.\n");
      return;
   }
   printf("Stratified homozygosity mapping for quantitative traits starts at %s", currentTime());
   int covariateCount = covariates.Length();
   IntArray validFlag(ped.count);
   validFlag.Set(1);
   for(int i = 0; i < ped.count; i++){
      bool allmissing = true;
      for(int t = 0; t < traitCount; t++)
         if(ped[i].isPhenotyped(traits[t])) {
            allmissing = false;
            break;
         }
      if(allmissing) validFlag[i] = 0;
      bool somemissing = false;
      for(int j = 0; j < covariateCount; j++)
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
   int NID = ID.Length();
   IntArray *validpheno = new IntArray[traitCount];
   for(int t = 0; t < traitCount; t++)
      validpheno[t].Dimension(0);
   Matrix Y(NID, traitCount);
   Matrix X(NID, 1+covariateCount);
   for(int i = 0; i < NID; i++)
      X[i][0] = 1.0;
   for(int t = 0; t < traitCount; t++)
      for(int i = 0; i < NID; i++)
         Y[i][t] = ped[ID[i]].traits[traits[t]];
   for(int i = 0; i < NID; i++)
      for(int j = 0; j < covariateCount; j++)
         X[i][j+1] = ped[ID[i]].covariates[covariates[j]];
   Vector meanY(traitCount);
   meanY.Set(0.0);
   for(int t = 0; t < traitCount; t++){
      int n = 0;
      for(int i = 0; i < NID; i++)
         if(Y[i][t] != _NAN_) {
            meanY[t] += Y[i][t];
            n ++;
         }
      meanY[t] /= n;
      for(int i = 0; i < NID; i++)
         if(Y[i][t] != _NAN_)
            validpheno[t].Push(i);
   }
   IntArray validphenoCount(traitCount);
   for(int t = 0; t < traitCount; t++)
      validphenoCount[t] = validpheno[t].Length();
   if(covariateCount==0) // chol is not needed
      for(int t = 0; t < traitCount; t++){
         double temp = 0.0;
         for(int k = 0; k < validphenoCount[t]; k++)
            temp += Y[validpheno[t][k]][t];
         temp /= validphenoCount[t];
         for(int k = 0; k < validphenoCount[t]; k++)
            Y[validpheno[t][k]][t] -= temp;
      }
   else{
      Vector tV(NID);
      Cholesky chol;
      Matrix Info0;
      Vector beta0(covariateCount+1);
      for(int t = 0; t < traitCount; t++){
         Info0.Dimension(covariateCount+1, covariateCount+1);
         Info0.Zero();
         for(int i = 0; i < covariateCount+1; i++)
            for(int j = 0; j < covariateCount+1; j++){
               double sum = 0.0;
               for(int k = 0; k < validphenoCount[t]; k++)
                  sum += X[validpheno[t][k]][i] * X[validpheno[t][k]][j];
               Info0[i][j] = sum;
            }
         if(chol.TryDecompose(Info0)==0)
            beta0.Zero();
         else{
            chol.Decompose(Info0);
            tV.Dimension(covariateCount+1);
            tV.Zero();
            for(int i = 0; i < covariateCount+1; i++)
               for(int k = 0; k < validphenoCount[t]; k++)
                  tV[i] += X[validpheno[t][k]][i] * Y[validpheno[t][k]][t];
            chol.BackSubst(tV);
            beta0 = chol.x;
         }
         for(int k = 0; k < validphenoCount[t]; k++){
            int index = validpheno[t][k];
            double sum = Y[index][t];
            for(int j = 0; j < covariateCount+1; j++)
               sum -= beta0[j] * X[index][j];
            Y[index][t] = sum;
         }
      }
   }
   StringArray traitNames(traitCount);
   for(int t = 0; t < traitCount; t++)
      traitNames[t] = ped.traitNames[traits[t]];
   IntArray membership(NID);
   membership.Set(-1);
   int covRef = -1;
   int popCount = 0;
   if(ped.covariateCount){
      covRef = ped.covariateNames.Find(popName);
      if(covRef > -1){
         for(int i = 0; i < NID; i++){
            int pop = int(ped[ID[i]].covariates[covRef]+0.5)-1;
            if(pop > -1){
               membership[i] = pop;
               if(pop >= popCount) popCount = pop+1;
            }
         }
         printf("There are %d populations in the stratified analysis\n", popCount);
      }
   }
   if(covRef==-1){
      popCount = 1;
      for(int i = 0; i < NID; i++)
         membership[i] = 0;
      printf("Covariate %s cannot be found. Non-stratified analysis is running.\n", popName);
   }
   IntArray *traitCounts = new IntArray [popCount];
   for(int pop = 0; pop < popCount; pop++){
      traitCounts[pop].Dimension(traitCount);
      traitCounts[pop].Zero();
   }
   IntArray SS(popCount);
   SS.Zero();
   for(int i = 0; i < NID; i++){
      int pop = membership[i];
      if(pop < 0) continue;
      SS[pop]++;
      for(int t = 0; t < traitCount; t++)
         if(Y[i][t] != _NAN_)
            traitCounts[pop][t] ++;
   }
   Matrix *popY = new Matrix[popCount];
   for(int pop = 0; pop < popCount; pop++)
      popY[pop].Dimension(SS[pop], traitCount);
   IntArray tempIndex(popCount);
   tempIndex.Zero();
   for(int i = 0; i < NID; i++){
      int pop = membership[i];
      if(pop < 0) continue;
      for(int t = 0; t < traitCount; t++)
         popY[pop][tempIndex[pop]][t] = Y[i][t];
      tempIndex[pop]++;
   }
   if(normalization){
      Vector tV;
      QuickIndex idx;
      IntArray validQ;
      for(int pop = 0; pop < popCount; pop++)
         for(int t = 0; t < traitCount; t++){
            validQ.Dimension(0);
            for(int i = 0; i < SS[pop]; i++)
               if(popY[pop][i][t] != _NAN_)
                  validQ.Push(i);
            tV.Dimension(traitCounts[pop][t]);
            for(int i = 0; i < traitCounts[pop][t]; i++)
               tV[i] = popY[pop][validQ[i]][t];
            idx.Index(tV);
            for(int i = 0; i < traitCounts[pop][t];){
               int start = i, end = i + 1;
               while(end < traitCounts[pop][t] && tV[idx[end]] == tV[idx[start]])
                  end++;
               end --;
               double q = ninv((start + (end-start)/2.0 + 0.5)/traitCounts[pop][t]);
               for(int j = start; j <= end; j++)
                  popY[pop][validQ[idx[j]]][t] = q;
               i = end + 1;
            }  // end of i
         }  // end of t for trait
   }  // end of if InvNorm

   Matrix traitMean(popCount, traitCount);
   traitMean.Zero();
   Matrix traitVar(popCount, traitCount);
   traitVar.Zero();
   for(int pop = 0; pop < popCount; pop++)
      for(int t = 0; t < traitCount; t++){
        for(int i = 0; i < SS[pop]; i++)
            if(popY[pop][i][t] != _NAN_){
               traitMean[pop][t] += popY[pop][i][t];
               traitVar[pop][t] += popY[pop][i][t]*popY[pop][i][t];
            }
         if(traitCounts[pop][t] > 1){
            traitMean[pop][t] /= traitCounts[pop][t];
            traitVar[pop][t] -= traitCounts[pop][t]*traitMean[pop][t]*traitMean[pop][t];
            traitVar[pop][t] /= (traitCounts[pop][t]-1);
         }
      }
   printf("The following traits are used for homozygosity mapping:\n");
   printf("Trait");
   for(int pop = 0; pop < popCount; pop++)
      printf("\tN_%d\tMean_%d\tSD_%d", pop+1, pop+1, pop+1);
   printf("\n");
   for(int t = 0; t < traitCount; t++){
      printf("%s", (const char*)traitNames[t]);
      for(int pop = 0; pop < popCount; pop++)
         printf("\t%d\t%.3lf\t%.3lf",
            traitCounts[pop][t], traitMean[pop][t], sqrt(traitVar[pop][t]));
      printf("\n");
   }
   printf("Scanning genome...\n");

   IntArray idList, rohStorage, rohIndex;
   const double LODfactor = 0.5 / log(10.0);
   Matrix *R, *R2, *extraR, *extraR2;
   R = new Matrix [popCount];
   R2 = new Matrix [popCount];
   extraR = new Matrix [popCount];
   extraR2 = new Matrix [popCount];
   IntArray **ROH = new IntArray *[popCount];
   IntArray **extraROH = new IntArray *[popCount];
   for(int pop = 0; pop < popCount; pop++){
      ROH[pop] = new IntArray [traitCount];
      extraROH[pop] = new IntArray [traitCount];
   }
   char buffer[0x20000];
   int pbuffer = 0;
   IntArray *ibuffer = new IntArray[defaultMaxCoreCount];
   Vector *rbuffer = new Vector[defaultMaxCoreCount];
   IntArray **starts = new IntArray *[popCount];
   IntArray **stops = new IntArray *[popCount];
   String outfile(prefix);
   outfile.Add(".mthomo");
   FILE *fp = fopen(outfile, "wb");
   pbuffer += sprintf(&buffer[pbuffer], "Chr\tPos\tFlankSNP1\tFlankSNP2\tTrait");
   for(int pop = 0; pop < popCount; pop ++)
      pbuffer += sprintf(&buffer[pbuffer], "\tLOD_%d", pop+1);
   pbuffer += sprintf(&buffer[pbuffer], "\tN_ROH\tN_Trait\tDiff\tSE\tLOD\n");
   int segCount = (chrSeg.Length()>>2);
   for(int seg = 0; seg < segCount; seg++){
      int chrsegMin = chrSeg[seg<<2];
      int chrsegMax = chrSeg[(seg<<2)|1];
      int ndim = chrsegMax - chrsegMin + 1;
      int chr = chromosomes[chrsegMin<<6];
      for(int pop = 0; pop < popCount; pop++){
         R[pop].Dimension(traitCount, ndim);
         R[pop].Zero();
         R2[pop].Dimension(traitCount, ndim);
         R2[pop].Zero();
         extraR[pop].Dimension(traitCount, ndim);
         extraR[pop].Zero();
         extraR2[pop].Dimension(traitCount, ndim);
         extraR2[pop].Zero();
         for(int t = 0; t < traitCount; t++){
            ROH[pop][t].Dimension(ndim);
            ROH[pop][t].Zero();
            extraROH[pop][t].Dimension(ndim);
            extraROH[pop][t].Zero();
         }
         starts[pop] = new IntArray [ndim];
         stops[pop] = new IntArray [ndim];
         for(int w = 0; w < ndim; w++){
            starts[pop][w].Dimension(0);
            stops[pop][w].Dimension(0);
         }
         idList.Dimension(0);
         for(int i = 0; i < NID; i++)
            if(membership[i]==pop)
               idList.Push(geno[ID[i]]);
         ROHOnly(idList, seg, rohStorage, rohIndex);
         int rohsegCount = (rohStorage.Length()>>1);
         for(int i = 0; i < rohsegCount; i++){
            int startPos = rohStorage[i*2];
            int stopPos = rohStorage[i*2+1];
            if(bp[stopPos] - bp[startPos] > MAXROHLENGTH) continue;
            int w = (startPos>>6)-chrsegMin;
            int offset = startPos&0x3F;
            if(w >= 0 && offset)
               starts[pop][w].Push((rohIndex[i]<<6) | offset);
            w = (stopPos>>6)-chrsegMin;
            offset = (stopPos+1)&0x3F;
            if(w < ndim && offset)
               stops[pop][w].Push((rohIndex[i]<<6) | offset);
         }
         for(int i = 0; i < rohsegCount; i++){
            int startPos = rohStorage[i*2];
            int stopPos = rohStorage[i*2+1];
            int startword = ((startPos+63)>>6);
            int stopword = ((stopPos+1)>>6)-1;
            int id = rohIndex[i];
            for(int t = 0; t < traitCount; t++){
               double localtrait = popY[pop][id][t];
               if(localtrait == _NAN_) continue;
               double localtrait2 = localtrait*localtrait;
               int w = startword-chrsegMin-1;
               if((w >= 0) && (startPos <= (startword<<6)-33)){
                  extraR[pop][t][w] += localtrait;
                  extraR2[pop][t][w] += localtrait2;
                  extraROH[pop][t][w] ++;
               }
               int localcount = stopword-chrsegMin+1;
               for(w++; w < localcount; w++){
                  R[pop][t][w] += localtrait;
                  R2[pop][t][w] += localtrait2;
                  ROH[pop][t][w] ++;
               }
               if((w < ndim) && (stopPos&0x3F) >= 32){
                  extraR[pop][t][w] += localtrait;
                  extraR2[pop][t][w] += localtrait2;
                  extraROH[pop][t][w] ++;
               }
            }  // end of trait t loop
         }  // end of ROH i loop
      }  // end of pop loop
      for(int c = 0; c < defaultMaxCoreCount; c++){
         ibuffer[c].Dimension(0);
         rbuffer[c].Dimension(0);
      }
      int thread = 0;
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) private(thread)
{
#endif
      Vector bitR[64], bitR2[64], bitDiff[64], bitVar[64];
      IntArray bitROH[64], bitROHsum[64], bitSSsum[64];
      Matrix bitLOD[64];
      Vector localLOD(popCount);
      for(int b = 0; b < 64; b++){
         bitR[b].Dimension(traitCount);
         bitR2[b].Dimension(traitCount);
         bitROH[b].Dimension(traitCount);
         bitDiff[b].Dimension(traitCount);
         bitVar[b].Dimension(traitCount);
         bitROHsum[b].Dimension(traitCount);
         bitSSsum[b].Dimension(traitCount);
         bitLOD[b].Dimension(popCount, traitCount);
      }
      IntArray detailedtraits, nodetails;
#ifdef _OPENMP
   thread = omp_get_thread_num();
   #pragma omp for
#endif
      for(int w = 0; w < ndim; w++){
         int detailCount = 0;
         detailedtraits.Dimension(0);
         nodetails.Dimension(0);
         int basePos = ((w+chrsegMin)<<6);
         for(int t = 0; t < traitCount; t++){
            double diff = 0.0;
            double var = 0.0;
            for(int pop = 0; pop < popCount; pop++){
               int localROH = ROH[pop][t][w] + extraROH[pop][t][w];
               if(localROH < MINROHCOUNT) continue;
               double localR = R[pop][t][w] + extraR[pop][t][w];
               double localR2 = R2[pop][t][w] + extraR2[pop][t][w];
               localR2 -= localR*localR/localROH;
               localR2 /= (localROH-1);
               diff += localR - traitMean[pop][t]*localROH;
               var += (localR2 * (traitCounts[pop][t]*1.0 / localROH-2) + traitVar[pop][t]) * localROH * localROH / traitCounts[pop][t];
            }
            double LOD = diff * diff / var * LODfactor;
            if(LOD >= MINLOD){
               detailCount++;
               detailedtraits.Push(t);
            }else
               nodetails.Push(t);
         }
         if(detailCount){
            for(int b = 0; b < 64; b++){
               bitDiff[b].Zero();
               bitVar[b].Zero();
               bitROHsum[b].Zero();
               bitSSsum[b].Zero();
            }
            for(int pop = 0; pop < popCount; pop++){
               for(int b = 0; b < 64; b++){
                  bitR[b].Zero();
                  bitR2[b].Zero();
                  bitROH[b].Zero();
               }
               int startsCount = starts[pop][w].Length();
               for(int i = 0; i < startsCount; i++){
                  int offset = (starts[pop][w][i]&0x3F);
                  int id = (starts[pop][w][i]>>6);
                  for(int t = 0; t < detailCount; t++){
                     int trait = detailedtraits[t];
                     double localtrait = popY[pop][id][trait];
                     if(localtrait == _NAN_) continue;
                     double localtrait2 = localtrait*localtrait;
                     for(int b = offset; b < 64; b++){
                        bitR[b][trait] += localtrait;
                        bitR2[b][trait] += localtrait2;
                        bitROH[b][trait] ++;
                     }
                  }
               }
               int stopsCount = stops[pop][w].Length();
               for(int i = 0; i < stopsCount; i++){
                  int offset = (stops[pop][w][i]&0x3F);
                  int id = (stops[pop][w][i]>>6);
                  for(int t = 0; t < detailCount; t++){
                     int trait = detailedtraits[t];
                     double localtrait = popY[pop][id][trait];
                     if(localtrait == _NAN_) continue;
                     double localtrait2 = localtrait*localtrait;
                     for(int b = 0; b < offset; b++){
                        bitR[b][trait] += localtrait;
                        bitR2[b][trait] += localtrait2;
                        bitROH[b][trait] ++;
                     }  // end of b loop
                  }  // end of t loop for detailedtraits
               }  // end of i loop for stops
               for(int b = 0; b < 64; b++)
                  for(int t = 0; t < detailCount; t++){
                     int trait = detailedtraits[t];
                     int localROH = ROH[pop][trait][w] + bitROH[b][trait];
                     if(localROH < MINROHCOUNT) {
                        bitLOD[b][pop][trait] = _NAN_;
                        continue;
                     }
                     double localR = R[pop][trait][w] + bitR[b][trait];
                     double localR2 = R2[pop][trait][w] + bitR2[b][trait];
                     localR2 -= localR*localR/localROH;
                     localR2 /= (localROH-1);
                     double tempDiff = localR - traitMean[pop][trait]*localROH;
                     double tempVar = (localR2 * (traitCounts[pop][trait]*1.0 / localROH-2)
                     + traitVar[pop][trait]) * localROH * localROH / traitCounts[pop][trait];
                     bitDiff[b][trait] += tempDiff;
                     bitVar[b][trait] += tempVar;
                     bitLOD[b][pop][trait] = tempDiff * tempDiff / tempVar * LODfactor;
                     if(tempDiff < 0) bitLOD[b][pop][trait] = -bitLOD[b][pop][trait];
                     bitROHsum[b][trait] += localROH;
                     bitSSsum[b][trait] += traitCounts[pop][trait];
                  }
            }  // end of pop loop
            for(int b = 0; b < 64; b++)
               for(int t = 0; t < detailCount; t++){
                  int trait = detailedtraits[t];
                  int localROH = bitROHsum[b][trait];
                  if(localROH < 7) continue;
                  ibuffer[thread].Push(basePos+b);
                  ibuffer[thread].Push(basePos+b);
                  ibuffer[thread].Push(trait);
                  ibuffer[thread].Push(localROH);
                  ibuffer[thread].Push(bitSSsum[b][trait]);
                  for(int pop = 0; pop < popCount; pop++)
                     rbuffer[thread].Push(bitLOD[b][pop][trait]);
                  rbuffer[thread].Push(bitDiff[b][trait]);
                  rbuffer[thread].Push(bitVar[b][trait]);
               }
         }  // end of if detailCount
         if(detailCount < traitCount){
            int nodetailsCount = nodetails.Length();
            for(int t = 0; t < nodetailsCount; t++){
               int trait = nodetails[t];
               double diff = 0.0;
               double var = 0.0;
               int ROHsum = 0;
               int SSsum = 0;
               for(int pop = 0; pop < popCount; pop++){
                  int localROH = ROH[pop][trait][w] + extraROH[pop][trait][w];
                  if(localROH < MINROHCOUNT) {
                     localLOD[pop] = _NAN_;
                     continue;
                  }
                  double localR = R[pop][trait][w] + extraR[pop][trait][w];
                  double localR2 = R2[pop][trait][w] + extraR2[pop][trait][w];
                  localR2 -= localR*localR/localROH;
                  localR2 /= (localROH-1);
                  double tempDiff = localR - traitMean[pop][trait]*localROH;
                  double tempVar = (localR2 * (traitCounts[pop][trait]*1.0 / localROH-2)
                  + traitVar[pop][trait]) * localROH * localROH / traitCounts[pop][trait];
                  localLOD[pop] = tempDiff * tempDiff / tempVar * LODfactor;
                  if(tempDiff < 0) localLOD[pop] = -localLOD[pop];
                  diff += tempDiff;
                  var += tempVar;
                  ROHsum += localROH;
                  SSsum += traitCounts[pop][trait];
               }  // end of pop loop
               if(ROHsum < 7) continue;
               ibuffer[thread].Push(basePos+31);
               ibuffer[thread].Push(basePos+32);
               ibuffer[thread].Push(trait);
               ibuffer[thread].Push(ROHsum);
               ibuffer[thread].Push(SSsum);
               for(int pop = 0; pop < popCount; pop++)
                  rbuffer[thread].Push(localLOD[pop]);
               rbuffer[thread].Push(diff);
               rbuffer[thread].Push(var);
            }
         }
      }  // end of word w loop
#ifdef _OPENMP
}  // extra bracket for omp
#endif
      for(int c = 0; c < defaultMaxCoreCount; c++){
         int count = ibuffer[c].Length()/5;
         for(int i = 0; i < count; i++){
            int base = i*5;
            int pos1 = ibuffer[c][base];
            int pos2 = ibuffer[c][base+1];
            int trait = ibuffer[c][base+2];
            int localROH = ibuffer[c][base+3];
            int localSS = ibuffer[c][base+4];
            base = i*(2+popCount);
            double diff = rbuffer[c][base+popCount] / localROH;
            double var = rbuffer[c][base+popCount+1] / localROH / localROH;
            double LOD = diff * diff / var * LODfactor;
            if(diff < 0) LOD = -LOD;
            pbuffer += sprintf(&buffer[pbuffer], "%d\t%.3lf\t%s\t%s\t%s",
               chr, (pos1==pos2? bp[pos1]:(bp[pos1]+bp[pos2])/2)*0.000001,
               (const char*)snpName[pos1], (const char*)snpName[pos2],
               (const char*)traitNames[trait]);
            for(int pop = 0; pop < popCount; pop++)
               if(rbuffer[c][base+pop] == _NAN_)
                  pbuffer += sprintf(&buffer[pbuffer], "\t%s", "NA");
               else
                  pbuffer += sprintf(&buffer[pbuffer], "\t%.2lf", rbuffer[c][base+pop]);
            pbuffer += sprintf(&buffer[pbuffer], "\t%d\t%d\t%.3lf\t%.3lf\t%.2lf\n",
               localROH, localSS, diff, sqrt(var), LOD);
            if(pbuffer > 0xFFFF){  // buffer big enough for writing
               fwrite(buffer, 1, pbuffer, fp);
               pbuffer = 0;
            }
         }
      }
      for(int pop = 0; pop < popCount; pop++){
         delete []starts[pop];
         delete []stops[pop];
      }
   }  // end of seg loop
   for(int pop = 0;pop < popCount; pop++){
      delete []ROH[pop];
      delete []extraROH[pop];
   }
   delete []ROH;
   delete []extraROH;
   delete []starts;
   delete []stops;
   delete []ibuffer;
   delete []rbuffer;
   delete []traitCounts;
   delete []R;
   delete []R2;
   delete []extraR;
   delete []extraR2;
   delete []popY;
   if(pbuffer>0)
      fwrite(buffer, 1, pbuffer, fp);
   fclose(fp);
   printf("Stratified homozygosity mapping for quantitative traits ends at %s", currentTime());
   printf("Stratified homozygosity mapping for quantitative traits scan results saved in file %s\n", (const char*)outfile);
}

   /*
   for(int i = 0; i < NID; i++){
      int pop = membership[i];
      if(pop < 0) continue;
      for(int t = 0; t < traitCount; t++){
         double temp = Y[i][t];
         if(temp != _NAN_){
            traitMean[pop][t] += temp;
            traitVar[pop][t] += temp*temp;
         }
      }
   }
   for(int pop = 0; pop < popCount; pop++)
      for(int t = 0; t < traitCount; t++)
         if(traitCounts[pop][t] > 1){
            traitMean[pop][t] /= traitCounts[pop][t];
            traitVar[pop][t] -= traitCounts[pop][t]*traitMean[pop][t]*traitMean[pop][t];
            traitVar[pop][t] /= (traitCounts[pop][t]-1);
         }
*/

void Engine::HomozygosityMappingMH(const char *popName)
{
   printf("\nOptions in effect:\n");
   printf("\t--homomap\n");
   printf("\t--strat %s\n", popName);
   if(Bit64Flag)
      printf("\t--sysbit 64\n");
   if(CoreCount)
      printf("\t--cpus %d\n", CoreCount);
   if(lessmemFlag)
      printf("\t--lessmem\n");
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");

   const double MINLOD=1.0;
   const int MINMARGINCOUNT=5;
   const double MAXROHLENGTH=10.0;
   bool IBDvalidFlag = PreSegment();
   if(!IBDvalidFlag){
      printf("%s\n", (const char*)segmessage);
      printf("  Note chromosomal positions can be sorted conveniently using other tools such as PLINK.\n");
      return;
   }

   if(ped.affectionCount==0){
      printf("No disease data.\n");
      return;
   }
   printf("Stratified homozygosity mapping starts at %s", currentTime());

   int covRef = -1;
   int popCount = 0;
   if(ped.covariateCount){
      covRef = ped.covariateNames.Find(popName);
      if(covRef > -1){
         for(int i = 0; i < idCount; i++)
            if(ped[phenoid[i]].covariates[covRef] > popCount) popCount = int(ped[phenoid[i]].covariates[covRef]+0.5);
      }else{
         popCount = 1;
         printf("Covariate %s cannot be found. Non-stratified analysis is running.\n", popName);
      }
   }else{
      popCount = 1;
      printf("Covariate %s cannot be found. Non-stratified analysis is running.\n", popName);
   }
   IntArray *affection[2];
   for(int aff = 0; aff < 2; aff++){
      affection[aff] = new IntArray[popCount];
      for(int pop = 0; pop < popCount; pop++)
         affection[aff][pop].Dimension(0);
   }
   if(covRef > -1)
      for(int i = 0; i < idCount; i++){
         int pop = int(ped[phenoid[i]].covariates[covRef]+0.5)-1;
         int aff = ped[phenoid[i]].affections[0]-1;
         if(pop > -1 && aff > -1) affection[aff][pop].Push(i);
      }
   else
      for(int i = 0; i < idCount; i++){
         int aff = ped[phenoid[i]].affections[0]-1;
         if(aff > -1) affection[aff][0].Push(i);
      }
   IntArray affCount[2];
   for(int a = 0; a < 2; a++){
      affCount[a].Dimension(popCount);
      for(int pop = 0; pop < popCount; pop++)
         affCount[a][pop] = affection[a][pop].Length();
   }
   IntArray validpop(0);
   for(int pop = 0; pop < popCount; pop++)
      if(affCount[0][pop] >= MINMARGINCOUNT && affCount[1][pop] >= MINMARGINCOUNT)
         validpop.Push(pop);
   int validpopCount = validpop.Length();
   if(!validpopCount) {
      printf("No populations have enough cases/controls.\n");
      return;
   }
   if(validpopCount < popCount) {
      printf("The following populations are excluded from homozygosity mapping for too few cases or controls.\n");
      printf("Pop\tN_Contr\tN_Case\n");
      for(int pop = 0; pop < popCount; pop++)
         if(affCount[0][pop] < MINMARGINCOUNT || affCount[1][pop] < MINMARGINCOUNT)
            printf("%d\t%d\t%d\n", pop+1, affCount[0][pop], affCount[1][pop]);
   }
   IntArray idList, rohStorage, rohIndex;
   const double LODfactor = 0.5 / log(10.0);
   int bitPi[2][64];
   int validLoci = 0;
   Vector OR, ORdeno, Var, LOD;
   IntArray SS[3];
   Matrix eachOR, eachSE;
   double localLOD[64], localVar[64], localOR[64], localORdeno[64];
   Matrix localeachOR(popCount, 64);
   Matrix localeachSE(popCount, 64);
   IntArray *pi[2], *extraPi[2], **starts[2], **stops[2];
   for(int aff = 0; aff < 2; aff++){
      starts[aff] = new IntArray *[popCount];
      stops[aff] = new IntArray *[popCount];
      pi[aff] = new IntArray [popCount];
      extraPi[aff] = new IntArray [popCount];
   }
   printf("Scanning genome...\n");
   String outfile(prefix);
   outfile.Add(".homomapMH");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "Chr\tPos\tFlankSNP1\tFlankSNP2");
   for(int p = 0; p < validpopCount; p++)
      fprintf(fp, "\tOR_%d\tSE_%d", validpop[p]+1, validpop[p]+1);
   fprintf(fp, "\tN_Contr\tN_Case\tN_ROH\tOR_MH\tLOD\n");
   int segCount = (chrSeg.Length()>>2);
   for(int seg = 0; seg < segCount; seg++){
      int chrsegMin = chrSeg[seg<<2];
      int chrsegMax = chrSeg[(seg<<2)|1];
      int ndim = chrsegMax - chrsegMin + 1;
      OR.Dimension(ndim);
      ORdeno.Dimension(ndim);
      Var.Dimension(ndim);
      LOD.Dimension(ndim);
      eachOR.Dimension(popCount, ndim);
      eachSE.Dimension(popCount, ndim);
      for(int aff = 0; aff < 3; aff++){
         SS[aff].Dimension(ndim);
         SS[aff].Zero();
      }
      int chr = chromosomes[chrsegMin<<6];
      for(int p = 0; p < validpopCount; p++){
         int pop = validpop[p];
         for(int aff = 0; aff < 2; aff++){
            pi[aff][pop].Dimension(ndim);
            pi[aff][pop].Zero();
            extraPi[aff][pop].Dimension(ndim);
            extraPi[aff][pop].Zero();
            idList.Dimension(0);
            for(int s = 0; s < affCount[aff][pop]; s++)
               idList.Push(affection[aff][pop][s]);
            ROHOnly(idList, seg, rohStorage, rohIndex);
            starts[aff][pop] = new IntArray [ndim];
            stops[aff][pop] = new IntArray [ndim];
            for(int w = 0; w < ndim; w++){
               starts[aff][pop][w].Dimension(0);
               stops[aff][pop][w].Dimension(0);
            }
            int rohsegCount = (rohStorage.Length()>>1);
            for(int i = 0; i < rohsegCount; i++){
               int startPos = rohStorage[i*2];
               int stopPos = rohStorage[i*2+1];
               if(bp[stopPos] - bp[startPos] > MAXROHLENGTH) continue;
               int w = (startPos>>6)-chrsegMin;
               int offset = startPos&0x3F;
               if(w >= 0 && offset)
                  starts[aff][pop][w].Push(offset);
               w = (stopPos>>6)-chrsegMin;
               offset = (stopPos+1)&0x3F;
               if(w < ndim && offset)
                  stops[aff][pop][w].Push(offset);
               int startword = ((startPos-1)>>6)+1;
               int stopword = ((stopPos+1)>>6)-1;
               w = startword-chrsegMin-1;
               if(w >= 0 && startPos <= (startword<<6)-33)
                  extraPi[aff][pop][w] ++;
               int localcount = stopword-chrsegMin+1;
               for(w++; w < localcount; w++)
                  pi[aff][pop][w] ++;
               if(w < ndim && (stopPos&0x3F) >= 32)
                  extraPi[aff][pop][w] ++;
            }  // end of ith ROH loop
         }  // end of aff loop
      }  // end of pop loop
      for(int w = 0; w < ndim; w++){
         LOD[w] = Var[w] = OR[w] = ORdeno[w] = 0.0;
         for(int p = 0; p < validpopCount; p++){
            int pop = validpop[p];
            int a = pi[1][pop][w] + extraPi[1][pop][w];
            int b = pi[0][pop][w] + extraPi[0][pop][w];
            int c = affCount[1][pop] - a;
            int d = affCount[0][pop] - b;
            int n1 = a + c;
            int n2 = b + d;
            int m1 = a + b;
            int m2 = c + d;
            int n = n1 + n2;
            if(n1 < MINMARGINCOUNT || n2 < MINMARGINCOUNT || m1 < MINMARGINCOUNT || m2 < MINMARGINCOUNT) {
               eachOR[pop][w] = _NAN_;
               continue;
            }
            SS[1][w] += n1;
            SS[0][w] += n2;
            SS[2][w] += m1;
            LOD[w] += (a - 1.0 * n1 * m1/ n);
            Var[w] += (1.0 * n1 * n2 * m1 * m2/ n / n / (n-1));
            OR[w] += (a+0.5) * (d+0.5) / n;
            ORdeno[w] += (b+0.5) * (c+0.5) / n;
            eachOR[pop][w] = (a+0.5)*(d+0.5)/(b+0.5)/(c+0.5);
            eachSE[pop][w] = sqrt(1/(a+0.5)+1/(d+0.5)+1/(b+0.5)+1/(c+0.5));
         }
         if(Var[w] == 0.0){
            OR[w] = _NAN_;
            LOD[w] = 0.0;
         }else{
            OR[w] /= ORdeno[w];
            LOD[w] = LOD[w] * LOD[w] * LODfactor / Var[w];
            validLoci++;
         }
      }  // end of word w loop
      for(int w = 0; w < ndim; w++)
         if(LOD[w] >= MINLOD){
            if(w > 0 && LOD[w-1]!=_NAN_) LOD[w-1] = _NAN_;
            LOD[w] = _NAN_;
            if(w < ndim-1 && LOD[w+1] < MINLOD) LOD[w+1] = _NAN_;
         }
      for(int w = 0; w < ndim; w++){
         if(OR[w] == _NAN_) continue;
         int basePos = ((w+chrsegMin)<<6);
         if(LOD[w] == _NAN_){
            for(int bit = 0; bit < 64; bit++)
               localLOD[bit] = localVar[bit] = localOR[bit] = localORdeno[bit] = 0.0;
            for(int p = 0; p < validpopCount; p++){
               int pop = validpop[p];
               for(int aff = 0; aff < 2; aff++){
                  for(int b = 0; b < 64; b++)
                     bitPi[aff][b] = 0;
                  int startsCount = starts[aff][pop][w].Length();
                  for(int i = 0; i < startsCount; i++)
                     for(int b = starts[aff][pop][w][i]; b < 64; b++)
                     bitPi[aff][b] ++;
                  int stopsCount = stops[aff][pop][w].Length();
                  for(int i = 0; i < stopsCount; i++)
                     for(int b = 0; b < stops[aff][pop][w][i]; b++)
                        bitPi[aff][b] ++;
               }
               for(int bit = 0; bit < 64; bit++){
                  int a = pi[1][pop][w] + bitPi[1][bit];
                  int b = pi[0][pop][w] + bitPi[0][bit];
                  int c = affCount[1][pop] - a;
                  int d = affCount[0][pop] - b;
                  int n1 = a + c;
                  int n2 = b + d;
                  int m1 = a + b;
                  int m2 = c + d;
                  int n = n1 + n2;
                  if(n1 < MINMARGINCOUNT || n2 < MINMARGINCOUNT || m1 < MINMARGINCOUNT || m2 < MINMARGINCOUNT){
                     localeachOR[pop][bit] = _NAN_;
                     continue;
                  }
                  localLOD[bit] += a - 1.0 * n1 * m1 / n;
                  localVar[bit] += 1.0 * n1 * n2 * m1 * m2 / n / n / (n-1);
                  localOR[bit] += (a+0.5) * (d+0.5) / n;
                  localORdeno[bit] += (b+0.5) * (c+0.5) / n;
                  localeachOR[pop][bit] = (a+0.5)*(d+0.5)/(b+0.5)/(c+0.5);
                  localeachSE[pop][bit] = sqrt(1/(a+0.5)+1/(d+0.5)+1/(b+0.5)+1/(c+0.5));
               }
            }
            for(int bit = 0; bit < 64; bit++){
               if(localVar[bit] == 0.0) continue;
               localLOD[bit] = localLOD[bit] * localLOD[bit] * LODfactor / localVar[bit];
               localOR[bit] /= localORdeno[bit];
               if(localOR[bit] < 1) localLOD[bit] = -localLOD[bit];
               fprintf(fp, "%d\t%.3lf\t%s\t%s",
                  chr, bp[basePos+bit]*0.000001,
                  (const char*)snpName[basePos+bit],
                  (const char*)snpName[basePos+bit]);
               for(int p = 0; p < validpopCount; p++)
                  if(localeachOR[validpop[p]][bit] == _NAN_)
                     fprintf(fp, "\tNA\tNA");
                  else
                     fprintf(fp, "\t%.3lf\t%.3lf",
                        localeachOR[validpop[p]][bit], localeachSE[validpop[p]][bit]);
               fprintf(fp, "\t%d\t%d\t%d\t%.3lf\t%.2lf\n",
                  SS[0][w], SS[1][w], SS[2][w], localOR[bit], localLOD[bit]);
            }
         }else{
            if(OR[w] < 1) LOD[w] = -LOD[w];
            fprintf(fp, "%d\t%.3lf\t%s\t%s",
               chr, (bp[basePos+31]+bp[basePos+32])*0.0000005,
               (const char*)snpName[basePos+31], (const char*)snpName[basePos+32]);
            for(int p = 0; p < validpopCount; p++)
                  if(eachOR[validpop[p]][w] == _NAN_)
                     fprintf(fp, "\tNA\tNA");
                  else
                     fprintf(fp, "\t%.3lf\t%.3lf",
                        eachOR[validpop[p]][w], eachSE[validpop[p]][w]);
            fprintf(fp, "\t%d\t%d\t%d\t%.3lf\t%.2lf\n",
               SS[0][w], SS[1][w], SS[2][w], OR[w], LOD[w]);
         }
      }  // end of w
      for(int aff = 0; aff < 2; aff++)
         for(int p = 0; p < validpopCount; p++){
            int pop = validpop[p];
            delete []starts[aff][pop];
            delete []stops[aff][pop];
         }
   }  // end of seg loop
   fclose(fp);
   for(int aff = 0; aff < 2; aff++){
      delete []starts[aff];
      delete []stops[aff];
      delete []pi[aff];
      delete []extraPi[aff];
   }
   printf("\nStratified homozygosity mapping ends at %s", currentTime());
   printf("%d loci are tested\n", validLoci);
   printf("Stratified homozygosity mapping results saved in file %s\n",
      (const char*)outfile);
}

void Engine::HomozygosityMapping()
{
   printf("\nOptions in effect:\n");
   printf("\t--homomap\n");
   if(Bit64Flag)
      printf("\t--sysbit 64\n");
   if(CoreCount)
      printf("\t--cpus %d\n", CoreCount);
   if(lessmemFlag)
      printf("\t--lessmem\n");
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");

   const int MINMARGINCOUNT = 10;
   const int MAXROHLENGTH = 10000000;
   bool IBDvalidFlag = PreSegment();
   if(!IBDvalidFlag){
      printf("%s\n", (const char*)segmessage);
      printf("  Note chromosomal positions can be sorted conveniently using other tools such as PLINK.\n");
      return;
   }

   if(ped.affectionCount==0){
      printf("No disease data.\n");
      return;
   }
   printf("Homozygosity mapping starts at %s", currentTime());
   IntArray aff[2];
   for(int p = 0; p < 2; p++)
      aff[p].Dimension(0);
   for(int i = 0; i < idCount; i++){
      int pop = ped[phenoid[i]].affections[0]-1;
      if(pop > -1) aff[pop].Push(i);
   }
   int affCount[2];
   for(int a = 0; a < 2; a++)
      affCount[a] = aff[a].Length();
   if(!affCount[0] || !affCount[1]) {
      printf("Both affected and unaffected are required.\n");
      return;
   }
   printf("%d affected and %d unaffected individuals are used for homozygosity mapping.\n",
      affCount[1], affCount[0]);
   IntArray idList, rohStorage, rohIndex;
   const double LODfactor = 0.5 / log(10.0);
   Vector FROH[2];
   for(int pop = 0; pop < 2; pop++){
      FROH[pop].Dimension(affCount[pop]);
      FROH[pop].Zero();
   }
   IntArray pi[2], extraPi[2];
   int bitPi[2][64];
   Vector OR, se, LOD;
   IntArray *starts[2], *stops[2];
   printf("Scanning genome...\n");
   String outfile(prefix);
   outfile.Add(".homomap");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "Chr\tPos\tFlankSNP1\tFlankSNP2\tN_Con\tHomoCon\tN_Cas\tHomoCas\tOR\tSE\tLOD\n");
   int validLoci = 0;
   int segCount = (chrSeg.Length()>>2);
   for(int seg = 0; seg < segCount; seg++){
      int chrsegMin = chrSeg[seg<<2];
      int chrsegMax = chrSeg[(seg<<2)|1];
      int ndim = chrsegMax - chrsegMin + 1;
      OR.Dimension(ndim);
      se.Dimension(ndim);
      LOD.Dimension(ndim);
      int chr = chromosomes[chrsegMin<<6];
      for(int pop = 0; pop < 2; pop++){
         pi[pop].Dimension(ndim);
         pi[pop].Zero();
         extraPi[pop].Dimension(ndim);
         extraPi[pop].Zero();
         idList.Dimension(0);
         for(int s = 0; s < affCount[pop]; s ++)
            idList.Push(aff[pop][s]);
         ROHOnly(idList, seg, rohStorage, rohIndex);
         starts[pop] = new IntArray [ndim];
         stops[pop] = new IntArray [ndim];
         for(int w = 0; w < ndim; w++){
            starts[pop][w].Dimension(0);
            stops[pop][w].Dimension(0);
         }
         int rohsegCount = (rohStorage.Length()>>1);
         for(int i = 0; i < rohsegCount; i++){
            int startPos = rohStorage[i*2];
            int stopPos = rohStorage[i*2+1];
            int length = bp[stopPos] - bp[startPos];
            if(length > MAXROHLENGTH) continue;
            int w = (startPos>>6)-chrsegMin;
            int offset = startPos&0x3F;
            if(w >= 0 && offset)
               starts[pop][w].Push(offset);
            w = (stopPos>>6)-chrsegMin;
            offset = (stopPos+1)&0x3F;
            if(w < ndim && offset)
               stops[pop][w].Push(offset);
            int id = rohIndex[i];
            FROH[pop][id] += length;
            int startword = ((startPos-1)>>6)+1;
            int stopword = ((stopPos+1)>>6)-1;
            w = startword-chrsegMin-1;
            if(w >= 0 && startPos <= (startword<<6)-33)
               extraPi[pop][w] ++;
            int localcount = stopword-chrsegMin+1;
            for(w++; w < localcount; w++)
               pi[pop][w] ++;
            if(w < ndim && (stopPos&0x3F) >= 32)
               extraPi[pop][w] ++;
         }  // end of ith ROH loop
      }  // end of pop loop
      for(int w = 0; w < ndim; w++){
         int a = pi[1][w] + extraPi[1][w];
         int b = pi[0][w] + extraPi[0][w];
         int c = affCount[1] - a;
         int d = affCount[0] - b;
         int n1 = a + c;
         int n2 = b + d;
         int m1 = a + b;
         int m2 = c + d;
         int n = n1 + n2;
         if(n1 < MINMARGINCOUNT || n2 < MINMARGINCOUNT || m1 < MINMARGINCOUNT || m2 < MINMARGINCOUNT){
            OR[w] = _NAN_;
            LOD[w] = 0.0;
         }else{
            OR[w] = (a+0.5)*(d+0.5)/(b+0.5)/(c+0.5);
            se[w] = sqrt(1/(a+0.5) + 1/(b+0.5) + 1/(c+0.5) + 1/(d+0.5));
            double temp = a*d-b*c;
            LOD[w] = temp*temp*n*1.0/n1/n2/m1/m2*LODfactor;
            validLoci++;
         }
      }  // end of word w loop
      for(int w = 0; w < ndim; w++)
         if(LOD[w] >= 2.0){
            if(w > 0 && LOD[w-1]!=_NAN_) LOD[w-1] = _NAN_;
            LOD[w] = _NAN_;
            if(w < ndim-1 && LOD[w+1] < 2.0) LOD[w+1] = _NAN_;
         }
      for(int w = 0; w < ndim; w++){
         if(OR[w] == _NAN_) continue;
         int basePos = ((w+chrsegMin)<<6);
         if(LOD[w] == _NAN_){
            for(int pop = 0; pop < 2; pop++){
               for(int b = 0; b < 64; b++)
                  bitPi[pop][b] = 0;
               int startsCount = starts[pop][w].Length();
               for(int i = 0; i < startsCount; i++)
                  for(int b = starts[pop][w][i]; b < 64; b++)
                     bitPi[pop][b] ++;
               int stopsCount = stops[pop][w].Length();
               for(int i = 0; i < stopsCount; i++)
                  for(int b = 0; b < stops[pop][w][i]; b++)
                     bitPi[pop][b] ++;
            }
            for(int bit = 0; bit < 64; bit++){
               int a = pi[1][w] + bitPi[1][bit];
               int b = pi[0][w] + bitPi[0][bit];
               int c = affCount[1] - a;
               int d = affCount[0] - b;
               int n1 = a + c;
               int n2 = b + d;
               int m1 = a + b;
               int m2 = c + d;
               int n = n1 + n2;
               if(n1 < MINMARGINCOUNT || n2 < MINMARGINCOUNT || m1 < MINMARGINCOUNT || m2 < MINMARGINCOUNT) continue;
               double localOR = (a+0.5)*(d+0.5)/(b+0.5)/(c+0.5);
               double localse = sqrt(1/(a+0.5) + 1/(b+0.5) + 1/(c+0.5) + 1/(d+0.5));
               double temp = a*d-b*c;
               double localLOD = temp*temp*n*1.0/n1/n2/m1/m2*LODfactor;
               if(localOR < 1) localLOD = -localLOD;
               fprintf(fp, "%d\t%.3lf\t%s\t%s\t%d\t%d\t%d\t%d\t%.3lf\t%.3lf\t%.2lf\n",
                  chr, bp[basePos+bit]*0.000001,
                  (const char*)snpName[basePos+bit], (const char*)snpName[basePos+bit],
                  affCount[0], b, affCount[1], a, localOR, localse, localLOD);
            }
         }else{
            if(OR[w] < 1) LOD[w] = -LOD[w];
            fprintf(fp, "%d\t%.3lf\t%s\t%s\t%d\t%d\t%d\t%d\t%.3lf\t%.3lf\t%.2lf\n",
               chr, (bp[basePos+31]+bp[basePos+32])*0.0000005,
               (const char*)snpName[basePos+31], (const char*)snpName[basePos+32],
               affCount[0], pi[0][w]+extraPi[0][w], affCount[1],
               pi[1][w]+extraPi[1][w], OR[w], se[w], LOD[w]);
         }
      }
      for(int pop = 0; pop < 2; pop++){
         delete []starts[pop];
         delete []stops[pop];
      }
   }  // end of seg loop
   fclose(fp);
   printf("\nAffect\tF_ROH\tSD_FROH\n");
   for(int pop = 0; pop < 2; pop++){
      double factor = 1.0 / totalLength;
      FROH[pop].Multiply(factor);
      printf("%d\t%.4lf\t%.4lf\n",
         pop, FROH[pop].Average(), sqrt(FROH[pop].Var()));
   }
   printf("\nHomozygosity mapping ends at %s", currentTime());
   printf("%d loci are tested\n", validLoci);
   printf("Homozygosity mapping scan results saved in file %s\n", (const char*)outfile);
}

void Engine::HomozygosityMappingForQT()
{
   printf("\nOptions in effect:\n");
   printf("\t--mthomo\n");
   if(Bit64Flag)
      printf("\t--sysbit 64\n");
   if(CoreCount)
      printf("\t--cpus %d\n", CoreCount);
   if(lessmemFlag)
      printf("\t--lessmem\n");
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");

   const double MAXROHLENGTH=10.0;
   bool IBDvalidFlag = PreSegment();
   if(!IBDvalidFlag){
      printf("%s\n", (const char*)segmessage);
      printf("  Note chromosomal positions can be sorted conveniently using other tools such as PLINK.\n");
      return;
   }

   int traitCount = traits.Length();
   if(traitCount==0){
      printf("No quantitative trait data.\n");
      return;
   }
   printf("Homozygosity mapping for quantitative traits starts at %s", currentTime());
   int covariateCount = covariates.Length();
   IntArray validFlag(ped.count);
   validFlag.Set(1);
   for(int i = 0; i < ped.count; i++){
      bool allmissing = true;
      for(int t = 0; t < traitCount; t++)
         if(ped[i].isPhenotyped(traits[t])) {
            allmissing = false;
            break;
         }
      if(allmissing) validFlag[i] = 0;
      bool somemissing = false;
      for(int j = 0; j < covariateCount; j++)
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
   int NID = ID.Length();
   IntArray *validpheno = new IntArray[traitCount];
   for(int t = 0; t < traitCount; t++)
      validpheno[t].Dimension(0);
   Matrix Y(NID, traitCount);
   Matrix X(NID, 1+covariateCount);
   for(int i = 0; i < NID; i++)
      X[i][0] = 1.0;
   for(int t = 0; t < traitCount; t++)
      for(int i = 0; i < NID; i++)
         Y[i][t] = ped[ID[i]].traits[traits[t]];
   for(int i = 0; i < NID; i++)
      for(int j = 0; j < covariateCount; j++)
         X[i][j+1] = ped[ID[i]].covariates[covariates[j]];
   Vector meanY(traitCount);
   meanY.Set(0.0);
   for(int t = 0; t < traitCount; t++){
      int n = 0;
      for(int i = 0; i < NID; i++)
         if(Y[i][t] != _NAN_) {
            meanY[t] += Y[i][t];
            n ++;
         }
      meanY[t] /= n;
      for(int i = 0; i < NID; i++)
         if(Y[i][t] != _NAN_)
            validpheno[t].Push(i);
   }
   IntArray validphenoCount(traitCount);
   for(int t = 0; t < traitCount; t++)
      validphenoCount[t] = validpheno[t].Length();
   if(normalization){
      Vector tV(NID);
      QuickIndex idx;
      for(int t = 0; t < traitCount; t++){
         tV.Dimension(validpheno[t].Length());
         for(int i = 0; i < validphenoCount[t]; i++)
            tV[i] = Y[validpheno[t][i]][t];
         idx.Index(tV);
         for(int i = 0; i < validphenoCount[t];){
            int start = i, end = i + 1;
            while(end < validphenoCount[t] && tV[idx[end]] == tV[idx[start]])
               end++;
            end --;
            double q = ninv((start + (end-start)/2.0 + 0.5)/validpheno[t].Length());
            for(int j = start; j <= end; j++)
               Y[validpheno[t][idx[j]]][t] = q;
            i = end + 1;
         }
      }
   }
   if(covariateCount==0) // chol is not needed
      for(int t = 0; t < traitCount; t++){
         double temp = 0.0;
         for(int k = 0; k < validphenoCount[t]; k++)
            temp += Y[validpheno[t][k]][t];
         temp /= validphenoCount[t];
         for(int k = 0; k < validphenoCount[t]; k++)
            Y[validpheno[t][k]][t] -= temp;
      }
   else{
      Vector tV(NID);
      Cholesky chol;
      Matrix Info0;
      Vector beta0(covariateCount+1);
      for(int t = 0; t < traitCount; t++){
         Info0.Dimension(covariateCount+1, covariateCount+1);
         Info0.Zero();
         for(int i = 0; i < covariateCount+1; i++)
            for(int j = 0; j < covariateCount+1; j++){
               double sum = 0.0;
               for(int k = 0; k < validphenoCount[t]; k++)
                  sum += X[validpheno[t][k]][i] * X[validpheno[t][k]][j];
               Info0[i][j] = sum;
            }
         if(chol.TryDecompose(Info0)==0)
            beta0.Zero();
         else{
            chol.Decompose(Info0);
            tV.Dimension(covariateCount+1);
            tV.Zero();
            for(int i = 0; i < covariateCount+1; i++)
               for(int k = 0; k < validphenoCount[t]; k++)
                  tV[i] += X[validpheno[t][k]][i] * Y[validpheno[t][k]][t];
            chol.BackSubst(tV);
            beta0 = chol.x;
         }
         for(int k = 0; k < validphenoCount[t]; k++){
            int index = validpheno[t][k];
            double sum = Y[index][t];
            for(int j = 0; j < covariateCount+1; j++)
               sum -= beta0[j] * X[index][j];
            Y[index][t] = sum;
         }
      }
   }
   StringArray traitNames(traitCount);
   for(int t = 0; t < traitCount; t++)
      traitNames[t] = ped.traitNames[traits[t]];
   IntArray traitCounts(traitCount);
   traitCounts.Zero();
   Vector traitMean(traitCount);
   traitMean.Zero();
   Vector traitVar(traitCount);
   traitVar.Zero();
   for(int i = 0; i < NID; i++)
      for(int t = 0; t < traitCount; t++){
         double temp = Y[i][t];
         if(temp != _NAN_){
            traitCounts[t] ++;
            traitMean[t] += temp;
            traitVar[t] += temp*temp;
         }
      }
   for(int t = 0; t < traitCount; t++)
      if(traitCounts[t] > 1){
         traitMean[t] /= traitCounts[t];
         traitVar[t] -= traitCounts[t]*traitMean[t]*traitMean[t];
         traitVar[t] /= (traitCounts[t]-1);
      }
   printf("The following traits are used for homozygosity mapping:\n");
   printf("Trait\tN\tMean\tSD\n");
   for(int t = 0; t < traitCount; t++)
      printf("%s\t%d\t%.3lf\t%.3lf\n",
         (const char*)ped.traitNames[traits[t]], traitCounts[t], traitMean[t], sqrt(traitVar[t]));
   printf("Scanning genome...\n");

   IntArray idList(0), rohStorage, rohIndex;
   for(int i = 0; i < NID; i++)
      idList.Push(geno[ID[i]]);
   const double LODfactor = 0.5 / log(10.0);
   Matrix R, R2, extraR, extraR2;
   IntArray *ROH = new IntArray[traitCount];
   IntArray *extraROH = new IntArray[traitCount];
   char buffer[0x20000];
   int pbuffer = 0;

   String outfile(prefix);
   outfile.Add(".mthomo");
   FILE *fp = fopen(outfile, "wb");
   pbuffer += sprintf(&buffer[pbuffer], "Chr\tPos\tFlankSNP1\tFlankSNP2\tTrait\tN_ROH\tN_Trait\tAve_ROH\tAve_All\tDiff\tSE\tLOD\n");
   int segCount = (chrSeg.Length()>>2);
   for(int seg = 0; seg < segCount; seg++){
      int chrsegMin = chrSeg[seg<<2];
      int chrsegMax = chrSeg[(seg<<2)|1];
      int ndim = chrsegMax - chrsegMin + 1;
      int chr = chromosomes[chrsegMin<<6];
      R.Dimension(traitCount, ndim);
      R.Zero();
      R2.Dimension(traitCount, ndim);
      R2.Zero();
      extraR.Dimension(traitCount, ndim);
      extraR.Zero();
      extraR2.Dimension(traitCount, ndim);
      extraR2.Zero();
      for(int t = 0; t < traitCount; t++){
         ROH[t].Dimension(ndim);
         ROH[t].Zero();
         extraROH[t].Dimension(ndim);
         extraROH[t].Zero();
      }
      IntArray *starts = new IntArray [ndim];
      IntArray *stops = new IntArray [ndim];
      for(int w = 0; w < ndim; w++){
         starts[w].Dimension(0);
         stops[w].Dimension(0);
      }
      ROHOnly(idList, seg, rohStorage, rohIndex);
      int rohsegCount = (rohStorage.Length()>>1);
      for(int i = 0; i < rohsegCount; i++){
         int startPos = rohStorage[i*2];
         int stopPos = rohStorage[i*2+1];
         if(bp[stopPos] - bp[startPos] > MAXROHLENGTH) continue;
         int w = (startPos>>6)-chrsegMin;
         int offset = startPos&0x3F;
         if(w >= 0 && offset)
            starts[w].Push((rohIndex[i]<<6) | offset);
         w = (stopPos>>6)-chrsegMin;
         offset = (stopPos+1)&0x3F;
         if(w < ndim && offset)
            stops[w].Push((rohIndex[i]<<6) | offset);
      }
      for(int i = 0; i < rohsegCount; i++){
         int startPos = rohStorage[i*2];
         int stopPos = rohStorage[i*2+1];
         int startword = ((startPos+63)>>6);
         int stopword = ((stopPos+1)>>6)-1;
         int id = rohIndex[i];
         for(int t = 0; t < traitCount; t++){
            double localtrait = Y[id][t];
            if(localtrait == _NAN_) continue;
            double localtrait2 = localtrait*localtrait;
            int w = startword-chrsegMin-1;
            if((w >= 0) && (startPos <= (startword<<6)-33)){
               extraR[t][w] += localtrait;
               extraR2[t][w] += localtrait2;
               extraROH[t][w] ++;
            }
            int localcount = stopword-chrsegMin+1;
            for(w++; w < localcount; w++){
               R[t][w] += localtrait;
               R2[t][w] += localtrait2;
               ROH[t][w] ++;
            }
            if((w < ndim) && (stopPos&0x3F) >= 32){
               extraR[t][w] += localtrait;
               extraR2[t][w] += localtrait2;
               extraROH[t][w] ++;
            }
         }  // end of trait t loop
      }  // end of ROH i loop
      IntArray *ibuffer = new IntArray[defaultMaxCoreCount];
      Vector *rbuffer = new Vector[defaultMaxCoreCount];
      for(int c = 0; c < defaultMaxCoreCount; c++){
         ibuffer[c].Dimension(0);
         rbuffer[c].Dimension(0);
      }
      int thread = 0;
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(thread)
{
#endif
      Vector bitR[64], bitR2[64];
      IntArray bitROH[64];
      for(int b = 0; b < 64; b++){
         bitR[b].Dimension(traitCount);
         bitR2[b].Dimension(traitCount);
         bitROH[b].Dimension(traitCount);
      }
      IntArray detailedtraits, nodetails;
#ifdef _OPENMP
   thread = omp_get_thread_num();
   #pragma omp for
#endif
      for(int w = 0; w < ndim; w++){
         int detailCount = 0;
         detailedtraits.Dimension(0);
         nodetails.Dimension(0);
         int basePos = ((w+chrsegMin)<<6);
         for(int t = 0; t < traitCount; t++){
            int localROH = ROH[t][w] + extraROH[t][w];
            if(localROH < 7) continue;
            double localR = R[t][w] + extraR[t][w];
            double localR2 = R2[t][w] + extraR2[t][w];
            localR /= localROH;
            localR2 -= localR*localR*localROH;
            localR2 /= (localROH-1);
            double diff = localR - traitMean[t];
            double var = (localR2 * (traitCounts[t]*1.0 / localROH-2) + traitVar[t]) / traitCounts[t];
            double LOD = diff * diff / var * LODfactor;
            if(LOD >= 2.0){
               detailCount++;
               detailedtraits.Push(t);
            }else
               nodetails.Push(t);
         }
         if(detailCount){
            for(int b = 0; b < 64; b++){
               bitR[b].Zero();
               bitR2[b].Zero();
               bitROH[b].Zero();
            }
            int startsCount = starts[w].Length();
            for(int i = 0; i < startsCount; i++){
               int offset = (starts[w][i]&0x3F);
               int id = (starts[w][i]>>6);
               for(int t = 0; t < detailCount; t++){
                  int trait = detailedtraits[t];
                  double localtrait = Y[id][trait];
                  if(localtrait == _NAN_) continue;
                  double localtrait2 = localtrait*localtrait;
                  for(int b = offset; b < 64; b++){
                     bitR[b][trait] += localtrait;
                     bitR2[b][trait] += localtrait2;
                     bitROH[b][trait] ++;
                  }
               }
            }
            int stopsCount = stops[w].Length();
            for(int i = 0; i < stopsCount; i++){
               int offset = (stops[w][i]&0x3F);
               int id = (stops[w][i]>>6);
               for(int t = 0; t < detailCount; t++){
                  int trait = detailedtraits[t];
                  double localtrait = Y[id][trait];
                  if(localtrait == _NAN_) continue;
                  double localtrait2 = localtrait*localtrait;
                  for(int b = 0; b < offset; b++){
                     bitR[b][trait] += localtrait;
                     bitR2[b][trait] += localtrait2;
                     bitROH[b][trait] ++;
                  }
               }
            }
            for(int b = 0; b < 64; b++)
               for(int t = 0; t < detailCount; t++){
                  int trait = detailedtraits[t];
                  int localROH = ROH[trait][w] + bitROH[b][trait];
                  if(localROH < 7) continue;
                  double localR = R[trait][w] + bitR[b][trait];
                  double localR2 = R2[trait][w] + bitR2[b][trait];
                  localR /= localROH;
                  localR2 -= localR*localR*localROH;
                  localR2 /= (localROH-1);
                  ibuffer[thread].Push(basePos+b);
                  ibuffer[thread].Push(basePos+b);
                  ibuffer[thread].Push(trait);
                  ibuffer[thread].Push(localROH);
                  rbuffer[thread].Push(localR);
                  rbuffer[thread].Push(localR2);
               }
         }
         if(detailCount < traitCount){
            int nodetailsCount = nodetails.Length();
            for(int t = 0; t < nodetailsCount; t++){
               int trait = nodetails[t];
               int localROH = ROH[trait][w] + extraROH[trait][w];
               if(localROH < 7) continue;
               double localR = R[trait][w] + extraR[trait][w];
               double localR2 = R2[trait][w] + extraR2[trait][w];
               localR /= localROH;
               localR2 -= localR*localR*localROH;
               localR2 /= (localROH-1);
               ibuffer[thread].Push(basePos+31);
               ibuffer[thread].Push(basePos+32);
               ibuffer[thread].Push(trait);
               ibuffer[thread].Push(localROH);
               rbuffer[thread].Push(localR);
               rbuffer[thread].Push(localR2);
            }
         }
      }  // end of word w loop
#ifdef _OPENMP
}  // extra bracket for omp
#endif
      for(int c = 0; c < defaultMaxCoreCount; c++){
         int count = ibuffer[c].Length()/4;
         for(int i = 0; i < count; i++){
            int pos1 = ibuffer[c][i*4];
            int pos2 = ibuffer[c][i*4+1];
            int trait = ibuffer[c][i*4+2];
            int localROH = ibuffer[c][i*4+3];
            double localR = rbuffer[c][i*2];
            double localR2 = rbuffer[c][i*2+1];
            double diff = localR - traitMean[trait];
            double var = (localR2 * (traitCounts[trait]*1.0 / localROH-2) + traitVar[trait]) / traitCounts[trait];
            double LOD = diff * diff / var * LODfactor;
            if(diff < 0) LOD = -LOD;
            pbuffer += sprintf(&buffer[pbuffer], "%d\t%.3lf\t%s\t%s\t%s\t%d\t%d\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.2lf\n",
               chr, (pos1==pos2? bp[pos1]:(bp[pos1]+bp[pos2])/2)*0.000001,
               (const char*)snpName[pos1], (const char*)snpName[pos2],
               (const char*)traitNames[trait], localROH, traitCounts[trait],
               localR, traitMean[trait], diff, sqrt(var), LOD);
            if(pbuffer > 0xFFFF){  // buffer big enough for writing
               fwrite(buffer, 1, pbuffer, fp);
               pbuffer = 0;
            }
         }
      }
      delete []starts;
      delete []stops;
      delete []ibuffer;
      delete []rbuffer;
   }  // end of seg loop
   if(pbuffer>0)
      fwrite(buffer, 1, pbuffer, fp);
   fclose(fp);

   printf("Homozygosity mapping for quantitative traits ends at %s", currentTime());
   printf("Homozygosity mapping for quantitative traits scan results saved in file %s\n", (const char*)outfile);
}

void Engine::PopulationROH()
{
   printf("\nOptions in effect:\n");
   printf("\t--poproh\n");
   if(Bit64Flag)
      printf("\t--sysbit 64\n");
   if(CoreCount)
      printf("\t--cpus %d\n", CoreCount);
   if(lessmemFlag)
      printf("\t--lessmem\n");
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");

   bool IBDvalidFlag = PreSegment();
   if(!IBDvalidFlag){
      printf("%s\n", (const char*)segmessage);
      printf("  Note chromosomal positions can be sorted conveniently using other tools such as PLINK.\n");
      return;
   }

   if(ped.affectionCount==0 && ped.covariateCount==0){
      printf("Reference populations need to be specified as covariate Reference.\n");
      return;
   }
   printf("Population IBD analysis starts at %s", currentTime());
   int covRef = -1;
   int popCount = 0;
   if(ped.covariateCount){
      covRef = ped.covariateNames.Find("Reference");
      if(covRef > -1)
         for(int i = 0; i < idCount; i++)
            if(ped[phenoid[i]].covariates[covRef] > popCount) popCount = int(ped[phenoid[i]].covariates[covRef]+0.5);
   }
   if(covRef == -1){
      for(int i = 0; i < idCount; i++)
         if(ped[phenoid[i]].affections[0] > popCount) popCount = ped[phenoid[i]].affections[0];
      if(popCount > 0) printf("Populations are specified through affection status.\n");
      else{
         printf("Reference populations need to be specified as covariate Reference or affection status.\n");
         return;
      }
   }else
      printf("Populations are specified through covariate Reference.\n");
   IntArray *aff = new IntArray[popCount];
   for(int p = 0; p < popCount; p++)
      aff[p].Dimension(0);
   int NAcount = 0;
   if(covRef > -1)
      for(int i = 0; i < idCount; i++){
         int pop = int(ped[phenoid[i]].covariates[covRef]+0.5);
         if(pop < 1) NAcount++;
         else aff[pop-1].Push(i);
      }
   else
      for(int i = 0; i < idCount; i++){
         int pop = ped[phenoid[i]].affections[0]-1;
         if(pop > -1) aff[pop].Push(i);
         else NAcount++;
      }
   if(NAcount)
      printf("%d samples are skipped for unspecified population memberships.\n", NAcount);
   IntArray affCount(popCount);
   for(int a = 0; a < popCount; a++)
      affCount[a] = aff[a].Length();

   IntArray idList, rohStorage, rohIndex;
   Vector *FROH = new Vector [popCount];
   Vector *FROH10 = new Vector [popCount];
   Vector *FROH25 = new Vector [popCount];
   for(int pop = 0; pop < popCount; pop++){
      FROH[pop].Dimension(affCount[pop]);
      FROH[pop].Zero();
      FROH10[pop].Dimension(affCount[pop]);
      FROH10[pop].Zero();
      FROH25[pop].Dimension(affCount[pop]);
      FROH25[pop].Zero();
   }
   const double LODfactor = 0.5 / log(10.0);
   IntArray *pi = new IntArray [popCount];
   printf("Scanning genome...\n");
   String outfile(prefix);
   outfile.Add(".poproh");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "Chr\tPos\tFlankSNP1\tFlankSNP2");
   for(int pop = 0; pop < popCount; pop++)
      fprintf(fp, "\tROH_%d", pop+1);
   fprintf(fp, "\n");
   String outfile2(prefix);
   outfile2.Add(".rohdiff");
   FILE *fp2 = fopen(outfile2, "wt");
   fprintf(fp2, "Chr\tPos\tFlankSNP1\tFlankSNP2\tPop_Neg\tPop_Pos\tN_Neg\tROH_Neg\tN_Pos\tROH_Pos\tOR\tSE\tLOD\n");
   int segCount = (chrSeg.Length()>>2);
   for(int seg = 0; seg < segCount; seg++){
      int chrsegMin = chrSeg[seg<<2];
      int chrsegMax = chrSeg[(seg<<2)|1];
      int ndim = chrsegMax - chrsegMin + 1;
      int chr = chromosomes[chrsegMin<<6];
      for(int pop = 0; pop < popCount; pop++){
         pi[pop].Dimension(ndim);
         pi[pop].Zero();
         idList.Dimension(0);
         for(int s = 0; s < affCount[pop]; s ++)
            idList.Push(aff[pop][s]);
         ROHOnly(idList, seg, rohStorage, rohIndex);
         int tempcount = (rohStorage.Length()>>1);
         for(int i = 0; i < tempcount; i++){
            int startPos = rohStorage[i*2];
            int stopPos = rohStorage[i*2+1];
            int length = bp[stopPos] - bp[startPos];
            int id = rohIndex[i];
            FROH[pop][id] += length;
            if(length > 10000000){
               FROH10[pop][id] += length;
               FROH25[pop][id] += length;
            }else if(length > 2500000)
               FROH25[pop][id] += length;
            int startword = ((startPos-1)>>6)+1;
            int stopword = ((stopPos+1)>>6)-1;
            int w = startword-chrsegMin-1;
            if(w >= 0 && startPos <= (startword<<6)-33)
               pi[pop][w] ++;
            int localcount = stopword-chrsegMin+1;
            for(w++; w < localcount; w++)
               pi[pop][w] ++;
            if(w < ndim && (stopPos&0x3F) >= 32)
               pi[pop][w] ++;
         }  // end of ith ROH loop
      }  // end of pop loop
      char buffer[256];
      for(int w = 0; w < ndim; w++){
         int base = ((w+chrsegMin)<<6);
         sprintf(buffer, "%d\t%.3lf\t%s\t%s",
            chr, (bp[base+31]+bp[base+32])*0.0000005,
            (const char*)snpName[base+31], (const char*)snpName[base+32]);
         fprintf(fp, "%s", buffer);
         for(int pop = 0; pop < popCount; pop++)
            if(affCount[pop])
               fprintf(fp, "\t%.3lf", pi[pop][w] * 1.0 / affCount[pop]);
            else
               fprintf(fp, "\t0");
         fprintf(fp, "\n");
         for(int p1 = 0; p1 < popCount; p1++)
            for(int p2 = p1+1;  p2 < popCount; p2++){
               if(!pi[p1][w] || !pi[p2][w]) continue;
               int a = affCount[p1] - pi[p1][w];
               int b = affCount[p2] - pi[p2][w];
               if(!a || !b) continue;
               double OR = pi[p2][w] * a * 1.0 / (pi[p1][w] * b);
               double se = sqrt(1.0/pi[p1][w] + 1.0/pi[p2][w] + 1.0/a + 1.0/b);
               double LOD = log(OR) / se;
               LOD = LOD*LOD*LODfactor;
               if(OR < 1) LOD = -LOD;
               fprintf(fp2, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.3lf\t%.3lf\t%.2lf\n",
                  buffer, p1+1, p2+1, affCount[p1], pi[p1][w],
                  affCount[p2], pi[p2][w], OR, se, LOD);
            }
      }  // end of word w loop
   }  // end of seg loop
   fclose(fp);
   fclose(fp2);
   printf("\nPop\tF_ROH10\tF_ROH2\tF_ROH\tSD_FROH\n");
   for(int pop = 0; pop < popCount; pop++){
      double factor = 1.0 / totalLength;
      FROH[pop].Multiply(factor);
      FROH10[pop].Multiply(factor);
      FROH25[pop].Multiply(factor);
      printf("%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\n",
         pop+1, FROH10[pop].Average(), FROH25[pop].Average(),
         FROH[pop].Average(), sqrt(FROH[pop].Var()));
   }
   printf("\nPopulation ROH analysis ends at %s", currentTime());
   printf("ROH by population saved in file %s\n", (const char*)outfile);
   printf("ROH difference between populations saved in file %s\n", (const char*)outfile2);
}

void Engine::ROHOnly(IntArray & idList, int segment, IntArray & rohStorage, IntArray & rohIndex, bool LengthOnly)
{
   const int MINCCOUNT = 10; // N_CS >= 100
   int cCount, icCount, segstart, segstop;
   unsigned long long int word;
   int listCount = idList.Length();
   IntArray *allsegments = new IntArray[listCount];
   for(int id = 0; id < listCount; id++)
      allsegments[id].Dimension(0);
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
   private(word, segstart, segstop, cCount, icCount)
{
#endif
   IntArray tempStart, tempStop, mergedStart, mergedStop;
   IntArray startPos, stopPos, startExtraBit, stopExtraBit, cCounts;
#ifdef _OPENMP
   #pragma omp for
#endif
   for(int i = 0; i < listCount; i++){
      int id = idList[i];
      startPos.Dimension(0);
      stopPos.Dimension(0);
      startExtraBit.Dimension(0);
      stopExtraBit.Dimension(0);
      int chrsegMin = chrSeg[segment<<2];
      int chrsegMax = chrSeg[(segment<<2)|1];
      int minExtraBit = chrSeg[(segment<<2)|2];
      int maxExtraBit = chrSeg[(segment<<2)|3];
      tempStart.Dimension(0);
      tempStop.Dimension(0);
      cCounts.Dimension(0);
      for(int m = chrsegMin; m <= chrsegMax; m++){
         for(; (m <= chrsegMax) &&  // keep passing if m word includes IC (Aa) or does not include any C (AA)
         (((~LG[0][id][m]) & LG[1][id][m]) || (LG[0][id][m] & LG[1][id][m])==0);
            m++); //includes IC (Aa) OR does not include any C (AA)
         for(cCount = 0, segstart = m; m <= chrsegMax &&
            ((~LG[0][id][m]) & LG[1][id][m])==0; m++)
            cCount += popcount(LG[0][id][m] & LG[1][id][m]);
         if( (cCount < 5) || // for sparse array
            (cCount < 10 && bp[((m-1)<<6)|0x3F] - bp[segstart<<6] < 20000) ) // for dense array or WGS
            continue;// continue only when cCount>=10 AND > 20Kb
         for(segstart--; segstart >= chrsegMin &&  // icCount==0
            ((~LG[0][id][segstart]) & LG[1][id][segstart])==0; segstart--);
         tempStart.Push(segstart+1);
         tempStop.Push(m-1);
         cCounts.Push(cCount);
      }  // end of scan indexed by m
      int tempcount = tempStart.Length();
      if(tempcount == 0) continue;
      mergedStart.Dimension(0);
      mergedStop.Dimension(0);
      for(int t = 0; t < tempcount-1; t++){
         int gap = bp[tempStart[t+1]<<6]-bp[(tempStop[t]<<6)|0x3F];
         if(tempStart[t+1] - tempStop[t] < 3){
            tempStart[t+1] = tempStart[t];// merge if 1 word in-between
         }else if((gap < 5000000 && tempStart[t+1] - tempStop[t] < 100) || (gap < 1.0)){
            cCount = 0; // consistency C (AA) count
            icCount = -2;   // inconsistency IC (Aa) count
            for(int m = tempStop[t]+1; m < tempStart[t+1]; m++){
               icCount += popcount((~LG[0][id][m]) & LG[1][id][m]);
               cCount += popcount(LG[0][id][m] & LG[1][id][m]);
            }
            if(cCount > icCount*3){ // merge if C_roh > 75%
               tempStart[t+1] = tempStart[t];
               cCounts[t+1] += cCounts[t] + cCount;
            }else if((bp[(tempStop[t]<<6)|0x3F] - bp[tempStart[t]<<6] > 2500000) ||
               (cCounts[t] >= MINCCOUNT) ){
               mergedStart.Push(tempStart[t]);  // ROH segments need to be > 2.5MB
               mergedStop.Push(tempStop[t]);
            } // else discard the left interval
         }else if((bp[(tempStop[t]<<6)|0x3F] - bp[tempStart[t]<<6] > 2500000) ||
            (cCounts[t] >= MINCCOUNT) ){
            mergedStart.Push(tempStart[t]);  // ROH segments need to be > 2.5MB
            mergedStop.Push(tempStop[t]);    // No gap to consider
         } // else discard the left interval
      }
      if((bp[(tempStop[tempcount-1]<<6)|0x3F] - bp[tempStart[tempcount-1]<<6] > 2500000) ||
         (cCounts[tempcount-1] >= MINCCOUNT) ){
         mergedStart.Push(tempStart[tempcount-1]);  // ROH segments need to be > 2.5MB
         mergedStop.Push(tempStop[tempcount-1]);
      }
      tempcount = mergedStart.Length();
      for(int t = 0; t < tempcount; t++){
         startPos.Push(mergedStart[t]);
         stopPos.Push(mergedStop[t]);
         int m = mergedStart[t];
         if(m > chrsegMin){
            m --; // Aa
            word = LG[1][id][m] & (~LG[0][id][m]);
            int bit = 0;
            for(; bit < 63 && (word & (1<<(63-bit))) == 0; bit++);
            startExtraBit.Push(bit);
         }else
            startExtraBit.Push(minExtraBit);
         m = mergedStop[t];
         if(m < chrsegMax){
            m ++; // Aa
            word = LG[1][id][m] & (~LG[0][id][m]);
            int bit = 0;
            for(; bit < 63 && (word & (1<<bit)) == 0; bit++);
               stopExtraBit.Push(bit);
         }else
            stopExtraBit.Push(maxExtraBit);
      }// end of a ROH segment indexed by t
      int newsegCount = startPos.Length();
      for(int seg = 0; seg < newsegCount; seg++){
         allsegments[i].Push((startPos[seg]<<6)-startExtraBit[seg]);
         allsegments[i].Push((stopPos[seg] == longCount-1)? (markerCount-1):((stopPos[seg]<<6)|0x3F)+stopExtraBit[seg]);
      }
   }  // end of i loop
#ifdef _OPENMP
}
#endif
   rohStorage.Dimension(0);
   rohIndex.Dimension(0);
   for(int i = 0; i < listCount; i++){
      int segmentLength = allsegments[i].Length()/2;
      for(int k = 0; k < segmentLength; k++){
         rohStorage.Push(allsegments[i][k*2]);
         rohStorage.Push(allsegments[i][k*2+1]);
         rohIndex.Push(i);
      }
   }
   delete []allsegments;
}


