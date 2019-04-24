//////////////////////////////////////////////////////////////////////
// ibdmapping.cpp
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
#include "MathSVD.h"
#include "QuickIndex.h"
#ifdef _OPENMP
  #include <omp.h>
#endif
#ifdef __ZLIB_AVAILABLE__
  #include <zlib.h>
#endif

#ifdef WITH_LAPACK
extern "C" void dgesdd_(char*, int*, int*, double*, int*, double*, double *, int*, double*, int*, double*, int*, int*, int*);
extern "C" void dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);
#endif

void Engine::IBDMDS_Projection()
{
   printf("\nOptions in effect:\n");
   printf("\t--ibdmds\n");
   printf("\t--projection\n");
   if(mincons)
      printf("\t--mincons %d\n", mincons);
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
   printf("IBD-segment-based MDS Projection analysis starts at %s", currentTime());
   IntArray subset[2];
   for(int i = 0; i < 2; i++) subset[i].Dimension(0);
   for(int i = 0; i < ped.count; i++)
      if(ped[i].ngeno >= MINSNPCOUNT){
         if(ped[i].affections[0]!=2)
            subset[0].Push(geno[i]);
         else subset[1].Push(geno[i]);
      }
   int dimN = subset[0].Length();
   if(dimN < 2) {
      printf("The number of unaffected individuals is < 2.\n");
      return;
   }
   int dimN2 = subset[1].Length();
   printf("MDS on a subset of %d reference samples...\n", dimN);
   Vector EigenValue;
   Matrix EigenVector;
   IBDMDS_Internal(subset[0], EigenValue, EigenVector);
   printf("Each of %d samples is now projected to reference PCs...\n", dimN2);
   int dimPC = EigenValue.Length();
   IntArray allpairs(0);
   for(int i = 0; i < dimN2; i++)
      for(int j = 0; j < dimN; j++){
         allpairs.Push(subset[1][i]);
         allpairs.Push(subset[0][j]);
      }
   int allpairCount = allpairs.Length()/2;
   Vector pi(allpairCount);
   pi.Zero();
   IntArray ibdsegStorage[2], ibdsegIndex[2];
   int segCount = (chrSeg.Length()>>2);
   for(int seg = 0; seg < segCount; seg++){
      if(mincons)
         IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], true, 2500000, mincons);
      else
         IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], true);
      #pragma omp parallel for num_threads(defaultMaxCoreCount)
      for(int pair = 0; pair < allpairCount; pair++)
         pi[pair] += (ibdsegStorage[0][pair]*0.5 + ibdsegStorage[1][pair])*0.000001;
   }  // end of seg loop
   for(int pair = 0; pair < allpairCount; pair++)
      pi[pair] /= totalLength;
   String pedfile = prefix;
   pedfile.Add("pc.txt");
   FILE *fp = fopen(pedfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
   fprintf(fp, "FID IID FA MO SEX AFF");
   for(int j = 0; j < dimPC; j++)
      fprintf(fp, " PC%d", j+1);
   fprintf(fp, "\n");
   for(int i = 0; i < dimN; i++){
      int id = phenoid[subset[0][i]];
      fprintf(fp, "%s %s %s %s %d 1",
         (const char*)ped[id].famid, (const char*)ped[id].pid,
         (const char*)ped[id].fatid, (const char*)ped[id].motid,
         ped[id].sex);
      for(int j = 0; j < dimPC; j++)
         fprintf(fp, " %.4lf", EigenVector[i][j]);
      fprintf(fp, "\n");
   }
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimPC; j++)
         EigenVector[i][j] /= EigenValue[j];
   double PC[20];
   for(int i = 0; i < dimN2; i++){
      for(int k = 0; k < dimPC; k++)
         PC[k] = 0;
      for(int j = 0; j < dimN; j++){
         double p = pi[i*dimN+j];
         for(int k = 0; k < dimPC; k++)
            PC[k] += p * EigenVector[j][k];
      }
      int id = phenoid[subset[1][i]];
      fprintf(fp, "%s %s %s %s %d 2",
         (const char*)ped[id].famid, (const char*)ped[id].pid,
         (const char*)ped[id].fatid, (const char*)ped[id].motid,
         ped[id].sex);
      for(int k = 0; k < dimPC; k++)
         fprintf(fp, " %.4lf", PC[k]);
      fprintf(fp, "\n");
   }
   fclose(fp);
   printf("%d principal components saved in file %s\n", dimPC, (const char*)pedfile);
   printf("IBD-based MDS ends at %s", currentTime());
}

void Engine::IBDMDS_Internal(IntArray & subset, Vector & EigenValue, Matrix & EigenVector)
{
   int subsetCount = subset.Length();
   if(subsetCount > 0x7FFF) error("Internal ibdmds cannot handle %d samples at the moment", subsetCount);
   IntArray allpairs(0), index(0);
   const int BLOCKSIZE = 4;
   for(int s = 0; s < subsetCount; s += BLOCKSIZE)
      for(int s2 = s; s2 < s+BLOCKSIZE && s2 < subsetCount; s2++)
         for(int t = s; t < subsetCount; t += BLOCKSIZE)
            for(int t2 = t; t2 < t+BLOCKSIZE && t2 < subsetCount; t2++){
               if(t==s && t2 <= s2) continue;
               allpairs.Push(subset[s2]);
               allpairs.Push(subset[t2]);
               index.Push(s2);
               index.Push(t2);
            }
   int allpairCount = allpairs.Length()/2;
   Matrix D(subsetCount, subsetCount);
   D.Zero();
   IntArray ibdsegStorage[2], ibdsegIndex[2];
   int segCount = (chrSeg.Length()>>2);
   for(int seg = 0; seg < segCount; seg++){
      if(mincons)
         IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], true, 2500000, mincons);
      else
         IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], true);
      #pragma omp parallel for num_threads(defaultMaxCoreCount)
      for(int pair = 0; pair < allpairCount; pair++)
         D[index[pair*2]][index[pair*2+1]] += (ibdsegStorage[0][pair]*0.5 + ibdsegStorage[1][pair])*0.000001;
   }  // end of seg loop
   for(int i = 0; i < subsetCount; i++)
      for(int j = i+1; j < subsetCount; j++)
         D[j][i] = D[i][j] = D[i][j] / totalLength;
   for(int i = 0; i < subsetCount; i++)
      D[i][i] = 1.0;
   int dimN = subsetCount;
   int dimPC = dimN>20?20:dimN;
   EigenValue.Dimension(dimPC);
   EigenVector.Dimension(dimN, dimPC);
#ifdef WITH_LAPACK
   char JOBZ = 'A';
   int info;
   double *A = new double[dimN*dimN];
   int dimLA = dimN;
   for(int i = 0; i < dimN; i++)
     for(int j = 0; j < dimN; j++)
       A[i*dimN+j] = D[j][i];
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
   for(int i = 0; i < dimPC; i++)
      EigenValue[i] = S[i];
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimPC; j++)
         EigenVector[i][j] = VT[i*dimN+j];
   delete []S;
   delete []VT;
#else
   SVD svd;
   svd.Decompose(D);
   if(svd.n == 0) return;
   QuickIndex idx;
   idx.Index(svd.w);
   for(int i = 0; i < dimPC; i++)
      EigenValue[i] = svd.w[idx[dimN-1-i]];
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimPC; j++)
         EigenVector[i][j] = svd.v[i][idx[dimN-1-j]];
#endif
}



   /*
   Vector tempV(dimN);
   tempV.Zero();
   for(int j = 0; j < dimN; j++)
      for(int i = 0; i < dimN; i++)
         tempV[j] += D[i][j];
   tempV.Multiply(1.0/dimN);
   // (I-11'/N) * D
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimN; j++)
         D[i][j] -= tempV[j];
   tempV.Zero();
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimN; j++)
         tempV[i] += D[i][j];
   tempV.Multiply(1.0/dimN);
   // D * (I-11'/N)
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimN; j++)
         D[i][j] -= tempV[i];
     */

void Engine::HEreg()
{
   printf("\nOptions in effect:\n");
   printf("\t--HEreg\n");
   if(Bit64Flag)
      printf("\t--sysbit 64\n");
   if(CoreCount)
      printf("\t--cpus %d\n", CoreCount);
   if(lessmemFlag)
      printf("\t--lessmem\n");
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");
   printf("Haseman-Elston regression scan starts at %s", currentTime());
   bool IBDvalidFlag = PreSegment();
   if(!IBDvalidFlag){
      printf("%s\n", (const char*)segmessage);
      printf("  Note chromosomal positions can be sorted conveniently using other tools such as PLINK.\n");
      return;
   }
   IntArray allpairs(0);
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(ped[i].sibCount > 1 && ped[i].sibs[0]->serial == i)
            for(int s = 0; s < ped[i].sibCount; s++)
               for(int t = s+1; t < ped[i].sibCount; t++){
                  allpairs.Push(geno[ped[i].sibs[s]->serial]);
                  allpairs.Push(geno[ped[i].sibs[t]->serial]);
               }
   int pairCount = allpairs.Length()>>1;
   if(pairCount == 0) {printf("No sib pairs are found.\n"); return;}
   printf("%d sib pairs are used for Haseman-Elston regression.\n", pairCount);
   int traitCount = traits.Length();
   Matrix allY(pairCount, traitCount);
   allY.Set(_NAN_);
   Vector meanYs(traitCount);
   meanYs.Zero();
   Vector meanY2s(traitCount);
   meanY2s.Zero();
   IntArray yCount(traitCount);
   yCount.Zero();
   IntArray *sibs = new IntArray[traitCount];
   for(int t = 0; t < traitCount; t++)
      sibs[t].Dimension(0);
   for(int p = 0; p < pairCount; p++){
      int p1 = phenoid[allpairs[p*2]];
      int p2 = phenoid[allpairs[p*2+1]];
      for(int t = 0; t < traitCount; t++){
         double y1 = ped[p1].traits[traits[t]];
         double y2 = ped[p2].traits[traits[t]];
         if(y1 != _NAN_ && y2 != _NAN_){
            sibs[t].Push(p1);
            sibs[t].Push(p2);
            double temp = y1 - y2;
            temp = temp*temp;
            allY[p][t] = temp;
            meanYs[t] += temp;
            meanY2s[t] += temp*temp;
            yCount[t]++;
         }
      }
   }
   for(int t = 0; t < traitCount; t++)
      if(yCount[t] >= 5){
         meanYs[t] /= yCount[t];
         meanY2s[t] /= yCount[t];
      }
   Vector vars(traitCount);
   IntArray Ns(traitCount);
   Vector H2(traitCount);
   Vector Vg(traitCount);
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount)
#endif
   for(int t = 0; t < traitCount; t++){
      sibs[t].Sort();
      double sum = ped[sibs[t][0]].traits[traits[t]];
      double sqsum = sum * sum;
      Ns[t] = 1;
      int sibsCount = sibs[t].Length();
      for(int i = 1; i < sibsCount; i++){
         if(sibs[t][i] != sibs[t][i-1]){
            double temp = ped[sibs[t][i]].traits[traits[t]];
            sum += temp;
            sqsum += temp*temp;
            Ns[t]++;
         }
      }
      vars[t] = (sqsum - sum*sum/Ns[t])/(Ns[t]-1);
      Vg[t] = vars[t]*2 - meanYs[t];
      H2[t] = Vg[t] > 0? Vg[t] / vars[t]: 0;
      if(H2[t] > 1) H2[t] = 1;
   }
   printf("Heritability estimates using phenotypes of sibpairs only:\n");
   printf("%15s %7s %7s %10s %10s %10s\n", "Trait", "N_pairs", "N", "Var", "Var_g", "H2");
   for(int t = 0; t < traitCount; t++)
      printf("%15s %7d %7d %10.2lf %10.2lf %10.3lf\n",
         (const char*)ped.traitNames[traits[t]],
         yCount[t], Ns[t], vars[t], Vg[t], H2[t]);
   printf("\n");

   IntArray invalid(traitCount);
   invalid.Zero();
   for(int t = 0; t < traitCount; t++){
      if(yCount[t] < 5){
         printf("Trait %s is skipped for having less than 5 pairs of sibpairs.\n",
            (const char*)ped.traitNames[traits[t]]);
         invalid[t] = 1;
      }else if(meanYs[t] < 1E-10){
         printf("Trait %s is skipped for no variation within any sibpair.\n",
            (const char*)ped.traitNames[traits[t]]);
         invalid[t] = 1;
      }else if(vars[t] < 1E-10){
         printf("Trait %s is skipped for no variation across individuals.\n",
            (const char*)ped.traitNames[traits[t]]);
         invalid[t] = 1;
      }
   }
   IntArray *NbyIBD[2];
   Vector *meanYbyIBD[2];
   for(int k = 0; k < 2; k++){
      NbyIBD[k] = new IntArray[traitCount];
      meanYbyIBD[k] = new Vector[traitCount];
   }
   char header[256];
   char **buffer = new char *[traitCount];
   for(int t = 0; t < traitCount; t++){
      buffer[t] = new char[256];
      buffer[t][0]='\0';
   }
   String buffers;
   sprintf(header, "Chr\tPos\tFlankSNP1\tFlankSNP2\tN_IBD0\tN_IBD1\tN_IBD2\tSqDiff0\tSqDiff1\tSqdiff2\th2\tAlpha\tNegBeta\tSE\tT\tPvalue\tLOD\n");
   String *outfiles = new String[traitCount];
   FILE **fps = new FILE *[traitCount];
   for(int t = 0; t < traitCount; t++){
      if(invalid[t]) continue;
      outfiles[t] = prefix;
      outfiles[t].Add('_');
      outfiles[t].Add(ped.traitNames[traits[t]]);
      outfiles[t].Add(".her");
      fps[t] = fopen(outfiles[t], "wt");
      fprintf(fps[t], "%s", header);
   }
   bool headerPrinted=false;
   IntArray ibdsegStorage[2], ibdsegIndex[2];
   const double LODfactor = 0.5 / log(10.0);
   int segCount = (chrSeg.Length()>>2);
   for(int seg = 0; seg < segCount; seg++){
      int chrsegMin = chrSeg[seg<<2];
      int chrsegMax = chrSeg[(seg<<2)|1];
      int ndim = chrsegMax - chrsegMin + 1;
      IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0],
         ibdsegStorage[1], ibdsegIndex[1], false, 2500000, 10000);
      for(int k = 0; k < 2; k++){
         for(int t = 0; t < traitCount; t++){
            NbyIBD[k][t].Dimension(ndim);
            NbyIBD[k][t].Zero();
            meanYbyIBD[k][t].Dimension(ndim);
            meanYbyIBD[k][t].Zero();
         }
         int tempcount = ibdsegIndex[k].Length();
         for(int i = 0; i < tempcount; i++){
            int pair = ibdsegIndex[k][i];
            int startword = ((ibdsegStorage[k][i*2]-1)>>6)+1;
            int stopword = ((ibdsegStorage[k][i*2+1]+1)>>6)-1;
            int localcount = stopword-chrsegMin;
            int startw = startword-chrsegMin;
            if(startword-chrsegMin-1 >= 0 && ibdsegStorage[k][i*2] <= (startword<<6)-33)
               startw = startword-chrsegMin-1;
            int stopw = localcount;
            if(localcount+1 < ndim && (ibdsegStorage[k][i*2+1]&0x3F) >= 32)
               stopw = localcount + 1;
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount)
#endif
            for(int t = 0; t < traitCount; t++){
               double y = allY[pair][t];
               if(y != _NAN_)
                  for(int w = startw; w <= stopw; w++){
                     NbyIBD[k][t][w] ++;
                     meanYbyIBD[k][t][w] += y;
                  }
            }  // end of loop t for traits
         }  // end of loop i for pairs
      }  // end of loop k for IBD1 or IBD2
      int chr = chromosomes[chrsegMin<<6];
      buffers.Clear();
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount)
#endif
      for(int t = 0; t < traitCount; t++){
         if(invalid[t]) continue;
         int N = yCount[t];
         double meanY = meanYs[t];
         double meanY2 = meanY2s[t];
         double var = vars[t];
         double maxLOD = 3.6;
         buffer[t][0]='\0';
         String traitName = ped.traitNames[t];
         for(int w = 0; w <= chrsegMax - chrsegMin; w++){
            double meanX = (NbyIBD[0][t][w] * 0.5 + NbyIBD[1][t][w]) / N;
            double lxx = NbyIBD[0][t][w] * 0.25 + NbyIBD[1][t][w] - meanX * meanX * N;
            double lxy = meanYbyIBD[0][t][w] * 0.5 + meanYbyIBD[1][t][w] - meanX * meanY * N;
            meanX -= 0.5;
            for(int k = 0; k < 2; k++) meanYbyIBD[k][t][w] /= NbyIBD[k][t][w];
            double beta = lxy / lxx;
            double alpha = meanY - meanX * beta;
            double h2 = -beta*0.5 / var;
            if(h2>1) h2 = 1;
            double Q = (meanY2 + alpha*alpha - 2*alpha*meanY + 2*alpha*beta*meanX) * N
               + beta * beta * (lxx+meanX*meanX*N) - 2*beta*(lxy+meanX*meanY*N);// Q = Sum(Y - a - b X)^2
            double se = sqrt(Q / (N-1) / lxx);
            double Tstat = beta / se;
            double pvalue = beta < 0? tdist(fabs(Tstat), N-1)*0.5: 1;
            double LOD = Tstat * Tstat * LODfactor;
            if(beta > 0) LOD = -LOD;
            int N0 = N - NbyIBD[0][t][w] - NbyIBD[1][t][w];
            int base = ((w+chrsegMin)<<6);
            fprintf(fps[t], "%d\t%.3lf\t%s\t%s\t%d\t%d\t%d\t%.2lf\t%.2lf\t%.2lf\t%.3lf\t%.2lf\t%.2lf\t%.3lf\t%.3lf\t%.2G\t%.2lf\n",
               chr, (bp[base+31]+bp[base+32])*0.0000005,
               (const char*)snpName[base+31], (const char*)snpName[base+32],
               N0, NbyIBD[0][t][w], NbyIBD[1][t][w],
               (meanY*N - meanYbyIBD[0][t][w] * NbyIBD[0][t][w] - meanYbyIBD[1][t][w] * NbyIBD[1][t][w]) / N0,
               meanYbyIBD[0][t][w],
               meanYbyIBD[1][t][w],
               h2, alpha, -beta, se, Tstat, pvalue, LOD);
            if(LOD > maxLOD){
               maxLOD = LOD;
               sprintf(buffer[t], "%s\t%d\t%.3lf\t%s\t%s\t%d\t%d\t%d\t%.2lf\t%.2lf\t%.2lf\t%.3lf\t%.2lf\t%.2lf\t%.3lf\t%.3lf\t%.2G\t%.2lf\n",
               (const char*)traitName, chr, (bp[base+31]+bp[base+32])*0.0000005,
               (const char*)snpName[base+31], (const char*)snpName[base+32],
               N0, NbyIBD[0][t][w], NbyIBD[1][t][w],
               (meanY*N - meanYbyIBD[0][t][w] * NbyIBD[0][t][w] - meanYbyIBD[1][t][w] * NbyIBD[1][t][w]) / N0,
               meanYbyIBD[0][t][w],
               meanYbyIBD[1][t][w],
               h2, alpha, -beta, se, Tstat, pvalue, LOD);
            }
         }  // end of position
      }  // end of trait
      for(int t = 0; t < traitCount; t++)
         if(buffer[t][0]){
            if(!headerPrinted){
               printf("The following linkage regions reach genome-wide significance:\n");
               printf("Trait\t%s", header);
               headerPrinted = true;
            }
            printf("%s", buffer[t]);
         }
   }  // end of seg
   for(int k = 0; k < 2; k++){
      delete []NbyIBD[k];
      delete []meanYbyIBD[k];
   }
   for(int t = 0; t < traitCount; t++)
      delete []buffer[t];
   delete []buffer;
   printf("\nHaseman-Elston regression results saved in files:\n");
   for(int t = 0; t < traitCount; t++){
      if(invalid[t]) continue;
      fclose(fps[t]);
      printf("%s\t", (const char*)outfiles[t]);
   }
   String outfile(prefix);
   outfile.Add(".her");
   printf("\nHaseman-Elston regression results are also available in a single file %s.\n", (const char*)outfile);
   FILE *fp = fopen((const char*)outfile, "wt");
   fprintf(fp, "Trait\t%s", header);
   String line;
   for(int t = 0; t < traitCount; t++){
      if(invalid[t]) continue;
      FILE *input=fopen((const char*)outfiles[t], "rt");
      if(input==NULL) continue;
      line.ReadLine(input);
      while(!feof(input)){
         line.ReadLine(input);
         if(line.Length()>1)
            fprintf(fp, "%s\t%s\n", (const char*)ped.traitNames[traits[t]], (const char*)line);
      }
      fclose(input);
   }
   fclose(fp);
   delete []outfiles;
   delete []fps;
   printf("Haseman-Elston regression scan ends at %s\n", currentTime());
}

void Engine::IBDVC()
{
   printf("\nOptions in effect:\n");
   printf("\t--ibdvc\n");
   if(normalization)
      printf("\t--invnorm\n");
   if(mincons)
      printf("\t--mincons %d\n", mincons);
   if(Bit64Flag)
      printf("\t--sysbit 64\n");
   if(CoreCount)
      printf("\t--cpus %d\n", CoreCount);
   if(lessmemFlag)
      printf("\t--lessmem\n");
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");

   Cholesky chol;
   QuickIndex idx;
   bool IBDvalidFlag = PreSegment();
   if(!IBDvalidFlag){
      printf("%s\n", (const char*)segmessage);
      printf("  Note chromosomal positions can be sorted conveniently using other tools such as PLINK.\n");
      return;
   }
   IntArray allpairs(0);
   for(int i = 0; i < idCount; i++)
      for(int j = i+1; j < idCount; j++){
         allpairs.Push(i);
         allpairs.Push(j);
      }
   int allpairCount = allpairs.Length()/2;
   Matrix pi(idCount, idCount);
   pi.Zero();
   IntArray ibdsegStorage[2], ibdsegIndex[2];
   int segCount = (chrSeg.Length()>>2);
   for(int seg = 0; seg < segCount; seg++){
      if(mincons)
         IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], true, 2500000, mincons);
      else
         IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], true);
      for(int p = 0; p < allpairCount; p++)
         pi[allpairs[p*2]][allpairs[p*2+1]] += (ibdsegStorage[0][p]*0.5 + ibdsegStorage[1][p])*0.000001;
   }  // end of seg loop
   for(int i = 0; i < idCount; i++)
      for(int j = i+1; j < idCount; j++)
         pi[i][j] /= totalLength;
   for(int i = 0; i < idCount; i++) pi[i][i] = 1.0;
   printf("Similarity matrix consists of sample covariance estimates, standardized by diagonal elements\n");
   int degree[3];
   for(int d = 0; d < 3; d++)
      degree[d] = 0;
   for(int i = 0; i < idCount; i++)
      for(int j = i+1; j < idCount; j++)
         if(pi[i][j] > 0.8) // MZ twins
            degree[0] ++;
         else if(pi[i][j] > 0.354) // 1st-degree relatives
            degree[1] ++;
         else if(pi[i][j] > 0.177) // 2nd-degree relatives
            degree[2] ++;
   if(degree[0] || degree[1] || degree[2])
      printf("Close relatives identified: %d MZ twins/duplicates, %d 1st-degree and %d 2nd-degree relative pairs\n",
         degree[0], degree[1], degree[2]);
   else
      printf("No second-degree or closer relationships were found\n");

   int N_Trait = traits.Length();
   if(N_Trait == 0) error("Quantitative traits are not found for variance component analysis.");
   int N_Covariate = covariates.Length();
   lambda0.Dimension(N_Trait);
   Vector loglik0(N_Trait);
   loglik0.Set(1E10);
   Vector score0(N_Covariate+1);
   Matrix Info0(N_Covariate+1, N_Covariate+1);
   tau0.Dimension(N_Trait);
   Matrix InvX(N_Covariate+1, idCount);

   printf("\nPolygenic parameter estimates\n");
   printf("%-15s %7s %7s %7s %7s %7s %9s",
         "TraitName", "N", "Herit", "LogLik", "Lambda", "Tau", "Mu");
   for(int i= 0; i < N_Covariate; i++)
      printf(" %9s", (const char*)ped.covariateNames[covariates[i]]);
   printf("\n");
   IntArray validFlag(idCount);
   for(int t = 0; t < N_Trait; t++){
      validFlag.Set(1);
      for(int i = 0; i < idCount; i++){
         int id = phenoid[i];
         if(!ped[id].isPhenotyped(traits[t])) {
            validFlag[i] = 0;
            continue;
         }
         bool missing = false;
         for(int j = 0; j < N_Covariate; j++)
            if(ped[id].covariates[covariates[j]] == _NAN_)
               missing = true;
         if(missing)
            validFlag[i] = 0;
      }
      ID.Dimension(0);
      for(int i = 0; i < idCount; i++)
         if(validFlag[i]) ID.Push(i);
      int N_ID = ID.Length();
      for(int i = 0; i < N_ID; i++)
         for(int j = 0; j < i; j++)
            pi[i][j] = pi[ID[j]][ID[i]];
      Matrix X(N_ID, 1+N_Covariate);
      Vector Y(N_ID);
      for(int i = 0; i < N_ID; i++){
         int id = phenoid[ID[i]];
         Y[i] = ped[id].traits[traits[t]];
         X[i][0] = 1.0;
         for(int j = 0; j < N_Covariate; j++)
            X[i][j+1] = ped[id].covariates[covariates[j]];
      }
      double meanY = 0.0;
      int n=0;
      for(int i = 0; i < N_ID; i++)
         if(Y[i] != _NAN_) {
            meanY += Y[i];
            n ++;
         }
      meanY /= n;
      for(int i = 0; i < N_ID; i++)
         if(Y[i] == _NAN_)
            Y[i] = meanY;
      Vector tV(N_ID);
      if(normalization){
         for(int i = 0; i < N_ID; i++)
            tV[i] = Y[i];
         idx.Index(tV);
         for(int i = 0; i < N_ID;){
            int start = i, end = i + 1;
            while(end < N_ID && tV[idx[end]] == tV[idx[start]]) end++;
            end --;
            double q = ninv((start + (end-start)/2.0 + 0.5)/N_ID);
            for(int j = start; j <= end; j++)
               Y[idx[j]] = q;
            i = end + 1;
         }
      }
#ifdef WITH_LAPACK
      char JOBZ = 'A';
      int info;
      double *A = new double[N_ID*N_ID];
      int dimLA = N_ID;
      for(int i = 0; i < N_ID; i++)
         for(int j = 0; j < i; j++)
            A[i*N_ID+j] = pi[i][j];
      for(int i = 0; i < N_ID; i++)
         for(int j = i; j < N_ID; j++)
            A[i*N_ID+j] = pi[j][i];
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
      EV.Dimension(N_ID);
      for(int i = 0; i < N_ID; i++) EV[i] = S[i];
      delete []S;
      // VT[k*N_ID+j] stores jth eigenvector. k: id; j: marker
      UT = new double * [N_ID];
      for(int i = 0; i < N_ID; i++)
         UT[i] = new double [N_ID];
      for(int j = 0; j < N_ID; j++)
         for(int k = 0; k < N_ID; k++)
            UT[j][k] = VT[k*N_ID+j];
      delete []VT;
#else
      Matrix D(N_ID, N_ID);
      for(int i = 0; i < N_ID; i++)
         for(int j = i+1; j < N_ID; j++)
            D[i][j] = D[j][i];
      for(int i = 0; i < N_ID; i++)
         D[i][i] = 1.0;
      SVD svd;
      svd.Decompose(D);
      if(svd.n == 0) return;
      EV.Dimension(N_ID);
      for(int i = 0; i < N_ID; i++) EV[i] = svd.w[i];
      // svd.v[k][idx[N_ID-1-j]] stores jth eigenvector
      UT = new double * [N_ID];
      for(int i = 0; i < N_ID; i++)
         UT[i] = new double [N_ID];
      for(int j = 0; j < N_ID; j++)
         for(int k = 0; k < N_ID; k++)
            UT[j][k] = svd.v[k][j];
#endif
      UY.Dimension(N_ID);     // U'Y += UT[i][k] * Y[k]
      for(int i = 0; i < N_ID; i++){
         double temp = 0.0;
         for(int k = 0; k < N_ID; k++)
            temp += UT[i][k] * Y[k];
         UY[i] = temp;
      }
      for(int i = 0; i < N_ID; i++)
         for(int j = 0; j < N_Covariate+1; j++)
            InvX[j][i] = X[i][j];
      UX.Dimension(N_ID, N_Covariate+1);
      for(int i = 0; i < N_ID; i++)
         for(int j = 0; j < N_Covariate+1; j++){
            double temp = 0.0;
            for(int k = 0; k < N_ID; k++)
               temp += UT[i][k] * InvX[j][k];
            UX[i][j] = temp;
         }

      double LB, T, temp;
      int numiter;
      int maxIter=10000;
      bool quiet=true;
      T = 0.001;
      double R1R;
      Vector beta(N_Covariate+2);
      Vector beta0(N_Covariate+1);
      double Resid0;
      Vector Resid(N_ID);
      double lambda;
      double loglik;
      const int TOTALCUT=200;
      double funcx[TOTALCUT], x[TOTALCUT];
      currentT = t;
#pragma omp parallel for num_threads(defaultMaxCoreCount) private(numiter)
      for(int i = 0; i < TOTALCUT; i++){
         double a = exp(i*0.1 - 5);
         double b = exp(i*0.1 - 4.9);
         x[i] = minimize(a, b, T, funcx[i], numiter, maxIter, quiet);
      }
      lambda0.Zero();
      for(int i = 0; i < TOTALCUT; i++){
         if(funcx[i] < loglik0[currentT]) {
            loglik0[currentT] = funcx[i];
            lambda0[currentT] = 1/x[i];
         }
      }
      if(lambda0[currentT]<0.0067){  // border
         lambda0[currentT] = 0;
      }
      for(int i = 0; i < N_Covariate+1; i++)
         for(int j = 0; j < N_Covariate+1; j++){
            double temp = 0.0;
            for(int k = 0; k < N_ID; k++)
               temp += UX[k][i] * UX[k][j] / (lambda0[currentT] * EV[k] + 1);
            Info0[i][j] = temp;
         }
      for(int i = 0; i < N_Covariate+1; i++){
         double temp = 0.0;
         for(int k = 0; k < N_ID; k++)
            temp += UX[k][i] * UY[k] / (lambda0[currentT] * EV[k] + 1);
         score0[i] = temp;
      }
      if(N_Covariate==0)
         beta0[0] = score0[0]/Info0[0][0];
      else if(chol.TryDecompose(Info0)==0)
         beta0.Zero();
      else{
         chol.Decompose(Info0);
         chol.BackSubst(score0);
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
      printf("%-15s %7d %7.4lf %7.1lf %7.2lf %7.2lf",
         (const char*)ped.traitNames[traits[t]], N_ID, lambda0[t]/(1+lambda0[t]),
         -loglik0[t], lambda0[t], tau0[t]);
      for(int i = 0; i < N_Covariate+1; i++)
         printf(" %9.4lf", beta0[i]);
      printf("\n");

      for(int i = 0; i < N_ID; i++)
         delete UT[i];
      delete []UT;
   }  // end of for loop t

   // GWAS scan here to be implemented


}

void Engine::AllIBDSegments()
{
#ifndef __ZLIB_AVAILABLE__
   error("--ibdall cannot run without ZLIB");
#else
   if(idCount > 0x7FFF) error("--ibdall cannot handle %d samples at the moment", idCount);
   printf("\nOptions in effect:\n");
   printf("\t--ibdall\n");
   if(mincons)
      printf("\t--mincons %d\n", mincons);
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
   if(!IBDvalidFlag && segmessage!="Segments too short."){
      printf("%s\n", (const char*)segmessage);
      printf("  Note chromosomal positions can be sorted conveniently using other tools such as PLINK.\n");
      return;
   }
   printf("IBD segment inference starts at %s", currentTime());
   IntArray allpairs(0);
   const int BLOCKSIZE = 4;
   for(int s = 0; s < idCount; s += BLOCKSIZE)
      for(int s2 = s; s2 < s+BLOCKSIZE && s2 < idCount; s2++)
         for(int t = s; t < idCount; t += BLOCKSIZE)
            for(int t2 = t; t2 < t+BLOCKSIZE && t2 < idCount; t2++){
               if(t==s && t2 <= s2) continue;
               allpairs.Push(s2);
               allpairs.Push(t2);
            }
   int pairCount = allpairs.Length()>>1;
   IntArray ibdsegStorage[2], ibdsegIndex[2];
   int pbuffer0 = 0;
   const int coreCount = defaultMaxCoreCount;
   gzFile fps[coreCount];
   char buffer[coreCount][0x10000];
   int pbuffer[coreCount];
   String outfile(prefix);
   outfile.Add(".segments.gz");
   fps[0] = gzopen((const char*)outfile, "wb");
   pbuffer0 += sprintf(&buffer[0][pbuffer0], "FID1\tID1\tFID2\tID2\tIBDType\tChr\tStartMB\tStopMB\tStartSNP\tStopSNP\t\tN_SNP\tLength\n");
   gzwrite(fps[0], buffer[0], pbuffer0);
   gzclose(fps[0]);
   StringArray outfiles(defaultMaxCoreCount);
   for(int c = 1; c < defaultMaxCoreCount; c++){
      outfiles[c].Copy(prefix);
      outfiles[c] += (c+1);
      outfiles[c].Add("$$$.segments.gz");
   }
   int segCount = (chrSeg.Length()>>2);
   for(int seg = 0; seg < segCount; seg++){
      if(mincons)
         IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], false, 2500000, mincons);
      else
         IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1]);

      fps[0] = gzopen((const char*)outfile, "ab");
      int k = 1;   // IBD2
      pbuffer0=0;
      int tempcount = (ibdsegStorage[k].Length()>>1);
      for(int i = 0; i < tempcount; i++){
         int index = ibdsegIndex[k][i];
         int id1 = allpairs[index*2];
         int id2 = allpairs[index*2+1];
         int segstart = ibdsegStorage[k][i*2];
         int segstop = ibdsegStorage[k][i*2+1];
         pbuffer0 += sprintf(&buffer[0][pbuffer0], "%s\t%s\t%s\t%s\t%s\t%d\t%.3lf\t%.3lf\t%s\t%s\t%d\t%.1lf\n",
            (const char*)ped[phenoid[id1]].famid, (const char*)ped[phenoid[id1]].pid,
            (const char*)ped[phenoid[id2]].famid, (const char*)ped[phenoid[id2]].pid,
            "IBD2",
            chromosomes[segstart], bp[segstart]*0.000001, bp[segstop]*0.000001,
            (const char*)snpName[segstart], (const char*)snpName[segstop],
            segstop-segstart+1, (bp[segstop]-bp[segstart])*0.000001);
         if(pbuffer0 > 0xFEFF){  // buffer big enough for writing
            gzwrite(fps[0], buffer[0], pbuffer0);
            pbuffer0 = 0;
         }
      }  // end of i
      if(pbuffer0)
         gzwrite(fps[0], buffer[0], pbuffer0);

      k = 0;  // IBD1
      tempcount = (ibdsegStorage[k].Length()>>1);
   int thread = 0;
   pbuffer[0] = 0;
#ifdef _OPENMP
   for(int c = 1; c < defaultMaxCoreCount; c++)
      fps[c] = NULL;
   #pragma omp parallel num_threads(defaultMaxCoreCount) private(thread)
{
#endif
   char tbuffer[1024];
#ifdef _OPENMP
   thread = omp_get_thread_num();
   pbuffer[thread] = 0;
   #pragma omp for
#endif
      for(int i = 0; i < tempcount; i++){
         int index = ibdsegIndex[k][i];
         int id1 = allpairs[index*2];
         int id2 = allpairs[index*2+1];
         int segstart = ibdsegStorage[k][i*2];
         int segstop = ibdsegStorage[k][i*2+1];
         sprintf(tbuffer, "%s\t%s\t%s\t%s\t%s\t%d\t%.3lf\t%.3lf\t%s\t%s\t%d\t%.1lf\n",
            (const char*)ped[phenoid[id1]].famid, (const char*)ped[phenoid[id1]].pid,
            (const char*)ped[phenoid[id2]].famid, (const char*)ped[phenoid[id2]].pid,
            "IBD1",
            chromosomes[segstart], bp[segstart]*0.000001, bp[segstop]*0.000001,
            (const char*)snpName[segstart], (const char*)snpName[segstop],
            segstop-segstart+1, (bp[segstop]-bp[segstart])*0.000001);
         pbuffer[thread] += sprintf(&buffer[thread][pbuffer[thread]], "%s", tbuffer);
         if(pbuffer[thread] > 0xFEFF){  // buffer big enough for writing
            if(fps[thread] == NULL)
               fps[thread] = gzopen((const char*)outfiles[thread], "wb");
            gzwrite(fps[thread], buffer[thread], pbuffer[thread]);
            pbuffer[thread] = 0;
         }
      }  // end of i
#ifdef _OPENMP
      if(thread && fps[thread]!=NULL) gzclose(fps[thread]); // fps[0] still open
}  // end of OpenMP
#endif
      for(int c = 0; c < defaultMaxCoreCount; c++)
         if(pbuffer[c])
            gzwrite(fps[0], buffer[c], pbuffer[c]);
      gzclose(fps[0]);
#ifdef _OPENMP
      FILE *fp0 = fopen(outfile, "ab");
      for(int c = 1; c < defaultMaxCoreCount; c++){
         FILE *fp = fopen(outfiles[c], "rb");
         if(fp!=NULL){
            int count = fread(buffer[0], 1, 0x10000, fp);
            for(; count == 0x10000; count = fread(buffer, 1, 0x10000, fp))
               fwrite(buffer[0], 1, 0x10000, fp0);
            if(count)
               fwrite(buffer[0], 1, count, fp0);
            fclose(fp);
            remove(outfiles[c]);
         }
      }
      fclose(fp0);
#endif
   }  // end of seg
   printf("                        ends at %s", currentTime());
   printf("All IBD segments saved in a gzipped file %s\n", (const char*)outfile);
#endif
}

void Engine::IBDMDS()
{
   if(idCount > 0x7FFF) error("--ibdmds cannot handle %d samples at the moment", idCount);
   printf("\nOptions in effect:\n");
   printf("\t--ibdmds\n");
   if(mincons)
      printf("\t--mincons %d\n", mincons);
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
   printf("IBD-segment-based MDS analysis starts at %s", currentTime());
   IntArray allpairs(0);
   const int BLOCKSIZE = 4;
   for(int s = 0; s < idCount; s += BLOCKSIZE)
      for(int s2 = s; s2 < s+BLOCKSIZE && s2 < idCount; s2++)
         for(int t = s; t < idCount; t += BLOCKSIZE)
            for(int t2 = t; t2 < t+BLOCKSIZE && t2 < idCount; t2++){
               if(t==s && t2 <= s2) continue;
               allpairs.Push(s2);
               allpairs.Push(t2);
            }
   int allpairCount = allpairs.Length()/2;
   Matrix D(idCount, idCount);
   D.Zero();
   IntArray ibdsegStorage[2], ibdsegIndex[2];
   int segCount = (chrSeg.Length()>>2);
   for(int seg = 0; seg < segCount; seg++){
      if(mincons)
         IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], true, 2500000, mincons);
      else
         IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], true);
      #pragma omp parallel for num_threads(defaultMaxCoreCount)
      for(int pair = 0; pair < allpairCount; pair++)
         D[allpairs[pair*2]][allpairs[pair*2+1]] += (ibdsegStorage[0][pair]*0.5 + ibdsegStorage[1][pair])*0.000001;
   }  // end of seg loop
   for(int i = 0; i < idCount; i++)
      for(int j = i+1; j < idCount; j++)
         D[i][j] /= totalLength;
   IntArray toberemoved(idCount);
   toberemoved.Zero();
   const double T=0.3535534;
   for(int i = 0; i < idCount; i++)
      for(int j = i+1; j < idCount; j++)
         if(D[i][j] > T) { // close relative, i to be removed
            toberemoved[i] = 1;
            break;
         }
   IntArray ID(0);
   for(int i = 0; i < toberemoved.Length(); i++)
      if(toberemoved[i]==0) ID.Push(i);
   int dimN = ID.Length();
   for(int i = 0; i < dimN; i++)
      for(int j = i+1; j < dimN; j++)
         D[j][i] = D[ID[i]][ID[j]];
   for(int i = 0; i < dimN; i++)
      for(int j = i+1; j < dimN; j++)
         D[i][j] = D[j][i];
   for(int i = 0; i < dimN; i++)
      D[i][i] = 1.0;

   if(dimN < idCount){
      String removefile = prefix;
      removefile.Add("_relative_removed.txt");
      FILE *fp = fopen(removefile, "wt");
      if(fp == NULL) error("Cannot open %s to write.", (const char*)removefile);
      printf("  %d samples as in file %s are removed due to close relatedness.\n",
         idCount - dimN, (const char*)removefile);
      for(int i = 0; i < idCount; i++)
         if(toberemoved[i])
            fprintf(fp, "%s %s\n",
            (const char*)ped[phenoid[i]].famid, (const char*)ped[phenoid[i]].pid);
      fclose(fp);
   }
   /*
   Vector tempV(dimN);
   tempV.Zero();
   for(int j = 0; j < dimN; j++)
      for(int i = 0; i < dimN; i++)
         tempV[j] += D[i][j];
   tempV.Multiply(1.0/dimN);
   // (I-11'/N) * D
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimN; j++)
         D[i][j] -= tempV[j];
   tempV.Zero();
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimN; j++)
         tempV[i] += D[i][j];
   tempV.Multiply(1.0/dimN);
   // D * (I-11'/N)
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimN; j++)
         D[i][j] -= tempV[i];
   */
   printf("SVD starts at %s", currentTime());

#ifdef WITH_LAPACK
   printf("  LAPACK is used.\n");
   char JOBZ = 'A';
   int info;
   double *A = new double[dimN*dimN];
   int dimLA = dimN;
   for(int i = 0; i < dimN; i++)
     for(int j = 0; j < dimN; j++)
       A[i*dimN+j] = D[j][i];
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

   int dimPC = dimN>20?20:dimN;
   printf("Largest %d eigenvalues:", dimPC);
   for(int i = 0; i < dimPC; i++)
      printf(" %.2lf", S[i]);
   printf("\n");
   double totalVariance = 0;
   double variance20 = 0;
   for(int i = 0; i < dimPC; i++)
      variance20 += S[i];
   for(int i = 0; (S[i] > 1E-10) && (i < dimN-1); i++)
      totalVariance += S[i];
   printf("The first %d PCs are able to explain %.2lf / %.2lf = %.1lf%% of total variance.\n",
      dimPC, variance20, totalVariance, variance20/totalVariance*100);
   printf("The proportion of total variance explained (%%) by each PC is:\n  ");
   for(int i = 0; i < dimPC; i++)
      printf(" %.1lf", S[i]/totalVariance*100);
   printf("\n");
   delete []S;
#else
   printf("  Please re-compile KING with LAPACK library.\n");
   SVD svd;
   svd.Decompose(D);
   printf("done\n");
   if(svd.n == 0) return;
   QuickIndex idx;
   idx.Index(svd.w);

   int dimPC = dimN>20?20:dimN;
   printf("Largest %d eigenvalues:", dimPC);
   for(int i = 0; i < dimPC; i++)
      printf(" %.2lf", svd.w[idx[dimN-1-i]]);
   printf("\n");
#endif
   String pedfile = prefix;
   pedfile.Add("pc.txt");
   FILE *fp = fopen(pedfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
   fprintf(fp, "FID IID FA MO SEX AFF");
   for(int j = 0; j < dimPC; j++)
      fprintf(fp, " PC%d", j+1);
   fprintf(fp, "\n");
   for(int i = 0; i < dimN; i++){
      int id = ID[i];
      fprintf(fp, "%s %s %s %s %d",
         (const char*)ped[phenoid[id]].famid, (const char*)ped[phenoid[id]].pid,
         (const char*)ped[phenoid[id]].fatid, (const char*)ped[phenoid[id]].motid,
         ped[phenoid[i]].sex);
         fprintf(fp, " 1");
         for(int j = 0; j < dimPC; j++)
#ifdef WITH_LAPACK
            fprintf(fp, " %.4lf", VT[i*dimN+j]);
#else
            fprintf(fp, " %.4lf", svd.v[i][idx[dimN-1-j]]);
#endif
      fprintf(fp, "\n");
   }
   fclose(fp);
   printf("MDS ends at %s", currentTime());
   printf("%d principal components saved in file %s\n",
      dimPC, (const char*)pedfile);
#ifdef WITH_LAPACK
   delete []VT;
#endif
}

void Engine::IBDGRM()
{
#ifndef __ZLIB_AVAILABLE__
   error("--ibdGRM cannot run without ZLIB");
#else
   if(idCount > 0x7FFF) error("--ibdGRM cannot handle %d samples at the moment", idCount);
   printf("\nOptions in effect:\n");
   printf("\t--ibdGRM\n");
   if(mincons)
      printf("\t--mincons %d\n", mincons);
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
   printf("IBD-segment-based GRM inference starts at %s", currentTime());
   IntArray allpairs(0);
   const int BLOCKSIZE = 4;
   for(int s = 0; s < idCount; s += BLOCKSIZE)
      for(int s2 = s; s2 < s+BLOCKSIZE && s2 < idCount; s2++)
         for(int t = s; t < idCount; t += BLOCKSIZE)
            for(int t2 = t; t2 < t+BLOCKSIZE && t2 < idCount; t2++){
               if(t==s && t2 <= s2) continue;
               allpairs.Push(s2);
               allpairs.Push(t2);
            }
   int allpairCount = allpairs.Length()/2;
   Matrix D(idCount, idCount);
   D.Zero();
   IntArray ibdsegStorage[2], ibdsegIndex[2];
   int segCount = (chrSeg.Length()>>2);
   for(int seg = 0; seg < segCount; seg++){
      if(mincons)
         IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], true, 2500000, mincons);
      else
         IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], true);
      #pragma omp parallel for num_threads(defaultMaxCoreCount)
      for(int pair = 0; pair < allpairCount; pair++)
         D[allpairs[pair*2]][allpairs[pair*2+1]] += (ibdsegStorage[0][pair]*0.5 + ibdsegStorage[1][pair])*0.000001;
   }  // end of seg loop
   for(int i = 0; i < idCount; i++)
      D[i][i] = totalLength;
   printf("Saving the IBD-segment-based GRM at %s", currentTime());

   int validmarkerCount = 0;
   for(int seg = 0; seg < segCount; seg++)
      validmarkerCount += ((chrSeg[(seg<<2)|1] - chrSeg[seg<<2] + 1)<<6)+chrSeg[(seg<<2)|2]+chrSeg[(seg<<2)|3];

   allpairs.Dimension(0);
   for(int i = 0; i < idCount; i++)
      for(int j = 0; j <= i; j++){
         allpairs.Push(i);
         allpairs.Push(j);
      }
   allpairCount = allpairs.Length()/2;
   const int coreCount = defaultMaxCoreCount;
   gzFile fps[coreCount];
   char buffer[coreCount][0x10000];
   int pbuffer[coreCount];
   String outfile(prefix);
   outfile.Add(".grm.gz");
   StringArray outfiles(defaultMaxCoreCount);
   for(int c = 1; c < defaultMaxCoreCount; c++){
      outfiles[c].Copy(prefix);
      outfiles[c] += (c+1);
      outfiles[c].Add("$$$.grm.gz");
   }
   outfiles[0] = outfile;
   for(int c = 0; c < defaultMaxCoreCount; c++)
      fps[c] = NULL;
   int thread = 0;
   pbuffer[0] = 0;
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) private(thread)
{
#endif
   char tbuffer[1024];
#ifdef _OPENMP
   thread = omp_get_thread_num();
   pbuffer[thread] = 0;
   #pragma omp for
#endif
   for(int pair = 0; pair < allpairCount; pair++){
      int i = allpairs[pair*2];
      int j = allpairs[pair*2+1];
      D[j][i] /= totalLength; // i >= j
      sprintf(tbuffer, "%d\t%d\t%d\t%.5lf\n", i+1, j+1, validmarkerCount, D[j][i]);
      pbuffer[thread] += sprintf(&buffer[thread][pbuffer[thread]], "%s", tbuffer);
      if(pbuffer[thread] > 0xFEFF){  // buffer big enough for writing
         if(fps[thread] == NULL)
            fps[thread] = gzopen((const char*)outfiles[thread], "wb");
         gzwrite(fps[thread], buffer[thread], pbuffer[thread]);
         pbuffer[thread] = 0;
      }
   }
   if(thread==0 && pbuffer[thread]){  // buffer big enough for writing
      if(fps[thread] == NULL)
         fps[thread] = gzopen((const char*)outfiles[thread], "wb");
      gzwrite(fps[thread], buffer[thread], pbuffer[thread]);
      gzclose(fps[thread]);
   }
#ifdef _OPENMP
   if(thread && fps[thread])
      gzclose(fps[thread]);
}  // end of OpenMP
   for(int c = 1; c < defaultMaxCoreCount; c++){
      FILE *fp = fopen(outfiles[c], "rb");
      if(fp){
         FILE *fp0 = fopen(outfile, "ab");
         int count = fread(buffer[0], 1, 0x10000, fp);
         for(; count == 0x10000; count = fread(buffer, 1, 0x10000, fp))
            fwrite(buffer[0], 1, 0x10000, fp0);
         if(count)
            fwrite(buffer[0], 1, count, fp0);
         fclose(fp);
         remove(outfiles[c]);
         fclose(fp0);
      }
      if(pbuffer[c]){
         fps[0] = gzopen((const char*)outfiles[0], "ab");
         gzwrite(fps[0], buffer[c], pbuffer[c]);
         gzclose(fps[0]);
      }
   }
#endif
   String outfile2(prefix);
   outfile2.Add(".grm.id");
   FILE *fp=fopen((const char*)outfile2, "wt");
   for(int i = 0; i < idCount; i++)
      fprintf(fp, "%s\t%s\n",
         (const char*)ped[phenoid[i]].famid, (const char*)ped[phenoid[i]].pid);
   fclose(fp);
   printf("IBD-segment-based GRM inference ends at %s", currentTime());
   printf("GRM files are saved in %s and %s\n", (const char*)outfile, (const char*)outfile2);
#endif
}


void Engine::IBDmapping(int nperm)
{
   printf("\nOptions in effect:\n");
   printf("\t--ibdmap\n");
   printf("\t--nperm %d\n", nperm);
   if(mincons)
      printf("\t--mincons %d\n", mincons);
   if(Bit64Flag)
      printf("\t--sysbit 64\n");
   if(CoreCount)
      printf("\t--cpus %d\n", CoreCount);
   if(lessmemFlag)
      printf("\t--lessmem\n");
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");

   const int FREEPERMCOUNT=10;
   const int FIXEDPERMCOUNT=200;
   int permCount = nperm / FREEPERMCOUNT;
   if(nperm==1000) printf("Note number of permutations can be specified through --nperm.\n");
   if(nperm%FREEPERMCOUNT)
      printf("Permutation count is rounded to %d\n", permCount * FREEPERMCOUNT);

   bool IBDvalidFlag = PreSegment();
   if(!IBDvalidFlag){
      printf("%s\n", (const char*)segmessage);
      printf("  Note chromosomal positions can be sorted conveniently using other tools such as PLINK.\n");
      return;
   }

   printf("IBD mapping starts at %s", currentTime());

   IntArray allpairs[3], aff[2];
   for(int k = 0; k < 3; k++)
      allpairs[k].Dimension(0);
   for(int k = 0; k < 2; k++)
      aff[k].Dimension(0);
   IntArray validgeno(0);
   IntArray inverseid(idCount);
   inverseid.Set(-1);
   for(int i = 0; i < idCount; i++)
      if(ped[phenoid[i]].affections[0]==1 || ped[phenoid[i]].affections[0] == 2){
         inverseid[i] = validgeno.Length();
         validgeno.Push(i);
      }
   int validCount = validgeno.Length();
   IntArray afftype(validCount);
   for(int i = 0; i < validCount; i++){
      int tempaff = ped[phenoid[validgeno[i]]].affections[0];
      if(tempaff == 2){
         aff[1].Push(i);
         afftype[i] = 1;
      }else if(tempaff == 1){
         aff[0].Push(i);
         afftype[i] = 0;
      }
   }
   int pairCount[3], affcount[2];
   const int BLOCKSIZE = 8;
   for(int k = 0; k < 2; k++){
      affcount[k] = aff[k].Length();
      for(int s = 0; s < affcount[k]; s += BLOCKSIZE)
         for(int s2 = s; s2 < s+BLOCKSIZE && s2 < affcount[k]; s2++)
            for(int t = s; t < affcount[k]; t += BLOCKSIZE)
               for(int t2 = t; t2 < t+BLOCKSIZE && t2 < affcount[k]; t2++){
                  if(t==s && t2 <= s2) continue;
                  allpairs[k*2].Push(validgeno[aff[k][s2]]);
                  allpairs[k*2].Push(validgeno[aff[k][t2]]);
               }
   }
   for(int s = 0; s < affcount[0]; s += BLOCKSIZE)
      for(int s2 = s; s2 < s+BLOCKSIZE && s2 < affcount[0]; s2++)
         for(int t = 0; t < affcount[1]; t += BLOCKSIZE)
            for(int t2 = t; t2 < t+BLOCKSIZE && t2 < affcount[1]; t2++){
               allpairs[1].Push(validgeno[aff[0][s2]]);
               allpairs[1].Push(validgeno[aff[1][t2]]);
            }
   for(int k = 0; k < 3; k++)
      pairCount[k] = allpairs[k].Length()>>1;
   if(pairCount[2] == 0) {printf("No ARPs are found.\n"); return;}
   else printf("#%d URPs, #%d DRPs, #%d ARPs are used for IBD mapping.\n",
      pairCount[0], pairCount[1], pairCount[2]);
   char header[256];
   sprintf(header, "Chr\tPosMb\tFlank1\tFlank2\tPI_URP\tPI_DRP\tPI_ARP\tPI_Diff\tP\n");
   String outfile(prefix);
   outfile.Add(".ibdmap");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "%s", header);

   Random rand/*((long)time(&t))*/;  // set a random seed. To be improved later
   long *seeds = new long[permCount];
   for(int p = 0; p < permCount; p++)
      seeds[p] = rand.NextInt();
   IntArray ibdsegStorage[3][2], ibdsegIndex[3][2], pi[3];
   Vector delta, allnewdelta[FIXEDPERMCOUNT];
   IntArray *extremeCount = new IntArray[defaultMaxCoreCount];
   IntArray *allextremeCount[FIXEDPERMCOUNT];
   for(int p = 0; p < FIXEDPERMCOUNT; p++)
      allextremeCount[p] = new IntArray[defaultMaxCoreCount];
   Vector smallestPV(FIXEDPERMCOUNT);
   smallestPV.Set(1);
   printf("Scanning genome...\n");
   int segCount = (chrSeg.Length()>>2);
   for(int seg = 0; seg < segCount; seg++){
      int chrsegMin = chrSeg[seg<<2];
      int chrsegMax = chrSeg[(seg<<2)|1];
      if(seg == 0 || chromosomes[chrsegMin<<6]!=chromosomes[chrSeg[(seg-1)<<2]<<6])
         printf("  Chr %d...\n", chromosomes[chrsegMin<<6]);
      int ndim = chrsegMax - chrsegMin + 1;
      for(int type = 0; type < 3; type++){
         if(mincons)
            IBDSegOnly(allpairs[type], seg, ibdsegStorage[type][0], ibdsegIndex[type][0], ibdsegStorage[type][1], ibdsegIndex[type][1], false, 2500000, mincons);
         else
            IBDSegOnly(allpairs[type], seg, ibdsegStorage[type][0], ibdsegIndex[type][0], ibdsegStorage[type][1], ibdsegIndex[type][1]);
         pi[type].Dimension(ndim);
         pi[type].Zero();
         for(int k = 0; k < 2; k++){
            int localvalue = k+1;
            int tempcount = (ibdsegStorage[type][k].Length()>>1);
            for(int i = 0; i < tempcount; i++){
               int startword = ((ibdsegStorage[type][k][i*2]-1)>>6)+1;
               int stopword = ((ibdsegStorage[type][k][i*2+1]+1)>>6)-1;
               int w = startword-chrsegMin-1;
               if(w >= 0 && ibdsegStorage[type][k][i*2] <= (startword<<6)-33)
                  pi[type][w] += localvalue;
               int localcount = stopword-chrsegMin+1;
               for(w++; w < localcount; w++)
                  pi[type][w] += localvalue;
               if(w < ndim && (ibdsegStorage[type][k][i*2+1]&0x3F) >= 32)
                  pi[type][w] += localvalue;
            }
         }
      }
      delta.Dimension(ndim);
      for(int w = 0; w < ndim; w++)
         delta[w] = pi[2][w] * 0.5 / pairCount[2] - pi[1][w] * 0.5 / pairCount[1];
      for(int p = 0; p < FIXEDPERMCOUNT; p++){
         allnewdelta[p].Dimension(ndim);
         allnewdelta[p].Set(_NAN_);
      }

#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount)
{
#endif
      Random localrand0;
      IntArray perm0(idCount);
#ifdef _OPENMP
   #pragma omp for
#endif
      for(int p = 0; p < FIXEDPERMCOUNT; p++){
         localrand0.Reset(seeds[p]);
         for(int j = 0; j < validCount; j++)// validCount
            perm0[j] = j;
         for(int j = validCount; j > 1; j--){
            int pos = localrand0.NextInt()%j;
            int swap = perm0[j-1];
            perm0[j-1] = perm0[pos];
            perm0[pos] = swap;
         }
         IntArray newpi0[2];
         for(int newtype = 0; newtype < 2; newtype++){
            newpi0[newtype].Dimension(ndim);
            newpi0[newtype].Zero();
         }
         for(int type = 0; type < 3; type++){
            for(int k = 0; k < 2; k++){
               int localvalue = k+1;
               int tempcount = (ibdsegStorage[type][k].Length()>>1);
               for(int i = 0; i < tempcount; i++){
                  int index = ibdsegIndex[type][k][i];
                  int id1 = allpairs[type][index*2];
                  int id2 = allpairs[type][index*2+1];
                  int startword = ((ibdsegStorage[type][k][i*2]-1)>>6)+1;
                  int stopword = ((ibdsegStorage[type][k][i*2+1]+1)>>6)-1;
                  bool negativeDelta=true;
                  int localcount = stopword-chrsegMin+1;
                  int w;
                  for(w = startword-chrsegMin; w < localcount; w++)
                     if(delta[w]>0) {
                        negativeDelta = false;
                        break;
                     }
                  if(negativeDelta) continue;   // skip if delta all negative
                  int newtype = afftype[perm0[inverseid[id1]]] + afftype[perm0[inverseid[id2]]] - 1;
                  if(newtype==-1) continue;
                  w = startword-chrsegMin-1;
                  if(w >= 0 && ibdsegStorage[type][k][i*2] <= (startword<<6)-33)
                     newpi0[newtype][w] += localvalue;
                  for(w++; w < localcount; w++)
                     newpi0[newtype][w] += localvalue;
                  if(w < ndim && (ibdsegStorage[type][k][i*2+1]&0x3F) >= 32)
                     newpi0[newtype][w] += localvalue;
               }  // end of interval i
            }  // end of loop k for IBD1 or IBD2
         }  // end of loop type for URP, DRP, or ARP
         for(int w = 0; w < ndim; w++)
            if(delta[w] > 0 && (newpi0[0][w] || newpi0[1][w]))
               allnewdelta[p][w] =  newpi0[1][w] * 0.5 / pairCount[2] - newpi0[0][w] * 0.5 / pairCount[1];
      }  // end of permute p
#ifdef _OPENMP
}  // extra bracket for omp
#endif

      int thread = 0;
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) private(thread)
{
   thread = omp_get_thread_num();
#endif
   IntArray newpi[FREEPERMCOUNT][2];
   extremeCount[thread].Dimension(ndim);
   extremeCount[thread].Zero();
   IntArray perm[FREEPERMCOUNT];
   for(int p = 0; p < FIXEDPERMCOUNT; p++){
      allextremeCount[p][thread].Dimension(ndim);
      allextremeCount[p][thread].Zero();
   }
   for(int p = 0; p < FREEPERMCOUNT; p++)
      perm[p].Dimension(idCount); // validCount instead
   Random localrand;
#ifdef _OPENMP
   #pragma omp for
#endif
      for(int p = 0; p < permCount; p++){
         localrand.Reset(seeds[p]);
         for(int p2 = 0; p2 < FREEPERMCOUNT; p2++){
            for(int j = 0; j < validCount; j++)// validCount
               perm[p2][j] = j;
            for(int j = validCount; j > 1; j--){
               int pos = localrand.NextInt()%j;
               int swap = perm[p2][j-1];
               perm[p2][j-1] = perm[p2][pos];
               perm[p2][pos] = swap;
            }
            for(int newtype = 0; newtype < 2; newtype++){
               newpi[p2][newtype].Dimension(ndim);
               newpi[p2][newtype].Zero();
            }
         }
         for(int type = 0; type < 3; type++){
            for(int k = 0; k < 2; k++){
               int localvalue = k+1;
               int tempcount = (ibdsegStorage[type][k].Length()>>1);
               for(int i = 0; i < tempcount; i++){
                  int index = ibdsegIndex[type][k][i];
                  int id1 = allpairs[type][index*2];
                  int id2 = allpairs[type][index*2+1];
                  int startword = ((ibdsegStorage[type][k][i*2]-1)>>6)+1;
                  int stopword = ((ibdsegStorage[type][k][i*2+1]+1)>>6)-1;
                  bool negativeDelta=true;
                  int localcount = stopword-chrsegMin+1;
                  int w;
                  for(w = startword-chrsegMin; w < localcount; w++)
                     if(delta[w]>0) {
                        negativeDelta = false;
                        break;
                     }
                  if(negativeDelta) continue;   // skip if delta all negative
                  for(int p2 = 0; p2 < FREEPERMCOUNT; p2++){
                     int newtype = afftype[perm[p2][inverseid[id1]]] + afftype[perm[p2][inverseid[id2]]] - 1;
                     if(newtype==-1) continue;
                     w = startword-chrsegMin-1;
                     if(w >= 0 && ibdsegStorage[type][k][i*2] <= (startword<<6)-33)
                        newpi[p2][newtype][w] += localvalue;
                     for(w++; w < localcount; w++)
                        newpi[p2][newtype][w] += localvalue;
                     if(w < ndim && (ibdsegStorage[type][k][i*2+1]&0x3F) >= 32)
                        newpi[p2][newtype][w] += localvalue;
                  }  // endof permute p2
               }  // end of interval i
            }  // end of loop k for IBD1 or IBD2
         }  // end of loop type for URP, DRP, or ARP
         for(int p2 = 0; p2 < FREEPERMCOUNT; p2++)
            for(int w = 0; w < ndim; w++)
               if(delta[w] > 0 && (newpi[p2][0][w] || newpi[p2][1][w])){
                  double newdelta =  newpi[p2][1][w] * 0.5 / pairCount[2] - newpi[p2][0][w] * 0.5 / pairCount[1];
                  if(newdelta >= delta[w]) extremeCount[thread][w] ++;
                  for(int i = 0; i < FIXEDPERMCOUNT; i++)
                     if(allnewdelta[i][w]!=_NAN_ && newdelta >= allnewdelta[i][w])
                        allextremeCount[i][thread][w] ++;
               }
      }  // end of permute p
#ifdef _OPENMP
}  // extra bracket for omp
#endif
      int chr = chromosomes[chrsegMin<<6];
      double localPi[3];
      for(int w = 0; w < ndim; w++){
         for(int type = 0; type < 3; type++)
            localPi[type] = pi[type][w] * 0.5 / pairCount[type];
         double pvalue;
         for(int p = 0; p < FIXEDPERMCOUNT; p++){
            if(allnewdelta[p][w] == _NAN_) continue;
            if(allnewdelta[p][w] > 0){
               pvalue = 0.0;
               for(int t = 0; t < defaultMaxCoreCount; t++)
                  pvalue += allextremeCount[p][t][w];
               pvalue /= (permCount*FREEPERMCOUNT);
               if(pvalue < smallestPV[p]) smallestPV[p] = pvalue;
            }
         }
         pvalue = 0.0;
         if(delta[w] <= 0) pvalue = 1.0;
         else{
            for(int t = 0; t < defaultMaxCoreCount; t++)
               pvalue += extremeCount[t][w];
            pvalue /= (permCount*FREEPERMCOUNT);
         }
         double delta = localPi[2] - localPi[1];
         int base = ((w+chrsegMin)<<6);
         fprintf(fp, "%d\t%.3lf\t%s\t%s\t%.2G\t%.2G\t%.2G\t%.2G\t%.2G\n",
            chr, (bp[base+31]+bp[base+32])*0.0000005,
            (const char*)snpName[base+31], (const char*)snpName[base+32],
            localPi[0], localPi[1], localPi[2], delta, pvalue);
      }  // end of word loop
   }  // end of seg loop
   fclose(fp);

   QuickIndex idx;
   idx.Index(smallestPV);
//   smallestPV.Print();
   printf("P value range in %d permutations is from %.2G to %.2G\n", FIXEDPERMCOUNT, smallestPV[idx[0]], smallestPV[idx[199]]);
   printf("1%% genome-wide significance P value cutoff is %.2G\n", smallestPV[idx[1]]);
   printf("5%% genome-wide significance P value cutoff is %.2G\n", smallestPV[idx[9]]);
   printf("10%% genome-wide significance P value cutoff is %.2G\n", smallestPV[idx[19]]);
   printf("IBD mapping ends at %s", currentTime());
   printf("Genome-wide IBD-mapping scan results saved in file %s\n", (const char*)outfile);
}

void Engine::LocalH2()
{
   printf("\nOptions in effect:\n");
   printf("\t--ibdH2\n");
   if(mincons)
      printf("\t--mincons %d\n", mincons);
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
   int traitCount = traits.Length();
   if(traitCount == 0){printf("No traits are specified.\n"); return;}
   printf("Local heritability inference starts at %s", currentTime());

   const int BLOCKSIZE = 8;
   IntArray allpairs(0);
   for(int s = 0; s < idCount; s += BLOCKSIZE)
      for(int s2 = s; s2 < s+BLOCKSIZE && s2 < idCount; s2++)
         for(int t = s; t < idCount; t += BLOCKSIZE)
            for(int t2 = t; t2 < t+BLOCKSIZE && t2 < idCount; t2++){
               if(t==s && t2 <= s2) continue;
               allpairs.Push(s2);
               allpairs.Push(t2);
            }
   int pairCount = allpairs.Length()>>1;
   if(pairCount == 0) {printf("No pairs are found.\n"); return;}
   String outfile(prefix);
   outfile.Add(".ih2");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "Chr\tPos\tFlankSNP1\tFlankSNP2\tTrait\tN_Obs\tN_Pred\tR2\tH2\n");
   fprintf(fp, "\n");
   IntArray ibdsegStorage[2], ibdsegIndex[2], pi[3];
   IntArray **connection = new IntArray *[idCount];
   Matrix observed(idCount, traitCount);
   observed.Set(_NAN_);
   Vector observedMean(traitCount);
   observedMean.Zero();
   IntArray samplesize(traitCount);
   samplesize.Zero();
   for(int i = 0; i < idCount; i++){
      int id = phenoid[i];
      for(int t = 0; t < traitCount; t++)
         if(ped[id].traits[traits[t]] != _NAN_){
            observed[i][t] = ped[id].traits[traits[t]];
            observedMean[t] += observed[i][t];
            samplesize[t]++;
         }
   }
   for(int t = 0; t < traitCount; t++)
      observedMean[t] /= samplesize[t];

   for(int t = 0; t < traitCount; t++)
      printf("%s\t%.3lf\n", (const char*)ped.traitNames[traits[t]], observedMean[t]);

   Matrix *predicted = new Matrix[idCount];
   int segCount = (chrSeg.Length()>>2);
   printf("Scanning genome...\n");
   for(int seg = 0; seg < segCount; seg++){
      int chrsegMin = chrSeg[seg<<2];
      int chrsegMax = chrSeg[(seg<<2)|1];
      int ndim = chrsegMax - chrsegMin + 1;
      for(int i = 0; i < idCount; i++){
         predicted[i].Dimension(ndim, traitCount);
         predicted[i].Zero();
         connection[i] = new IntArray[ndim];
         for(int w = 0; w < ndim; w++){
            connection[i][w].Dimension(traitCount);
            connection[i][w].Zero();
         }
      }
      if(mincons)
         IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], false, 2500000, mincons);
      else
         IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1]);
      for(int k = 0; k < 2; k++){
         int tempcount = (ibdsegStorage[k].Length()>>1);
         for(int i = 0; i < tempcount; i++){
            int index = ibdsegIndex[k][i];
            int id1 = allpairs[index*2];
            int id2 = allpairs[index*2+1];
            int startword = ((ibdsegStorage[k][i*2]-1)>>6)+1;
            int stopword = ((ibdsegStorage[k][i*2+1]+1)>>6)-1;
            int localcount = stopword-chrsegMin;
            int w = startword-chrsegMin-1;
            if(w >= 0 && ibdsegStorage[k][i*2] <= (startword<<6)-33)
               for(int t = 0; t < traitCount; t++)
                  if(observed[id2][t] != _NAN_ && observed[id1][t] != _NAN_){
                     predicted[id1][w][t] += observed[id2][t];
                     predicted[id2][w][t] += observed[id1][t];
                     connection[id1][w][t] ++;
                     connection[id2][w][t] ++;
                  }
            for(w++; w <= localcount; w++)
               for(int t = 0; t < traitCount; t++)
                  if(observed[id2][t] != _NAN_ && observed[id1][t] != _NAN_){
                     predicted[id1][w][t] += observed[id2][t];
                     predicted[id2][w][t] += observed[id1][t];
                     connection[id1][w][t] ++;
                     connection[id2][w][t] ++;
                  }
            if(w < ndim && (ibdsegStorage[k][i*2+1]&0x3F) >= 32)
               for(int t = 0; t < traitCount; t++)
                  if(observed[id2][t] != _NAN_ && observed[id1][t] != _NAN_){
                     predicted[id1][w][t] += observed[id2][t];
                     predicted[id2][w][t] += observed[id1][t];
                     connection[id1][w][t] ++;
                     connection[id2][w][t] ++;
                  }
         }  // end of segment i loop
      }  // end of IBD1 or IBD2 k loop
      int chr = chromosomes[chrsegMin<<6];
      Vector mean(traitCount);
      IntArray predictedSS(traitCount);
      double lxx, lxy, lyy, meanX;
      for(int w = 0; w < ndim; w++){
         mean.Zero();
         predictedSS.Zero();
         for(int i = 0; i < idCount; i++)
            for(int t = 0; t < traitCount; t++)
               if(connection[i][w][t]){
                  predicted[i][w][t] /= connection[i][w][t];
                  mean[t] += predicted[i][w][t];
                  predictedSS[t]++;
               }else
                  predicted[i][w][t] = _NAN_;

         int base = ((w+chrsegMin)<<6);
         for(int t = 0; t < traitCount; t++){
            mean[t] /= predictedSS[t];
            lxy = lxx = lyy = meanX = 0;
            for(int i = 0; i < idCount; i++)
               if(observed[i][t] != _NAN_ && predicted[i][w][t] != _NAN_){
                  lxy += predicted[i][w][t] * observed[i][t];
                  lyy += predicted[i][w][t] * predicted[i][w][t];
                  lxx += observed[i][t] * observed[i][t];
                  meanX += observed[i][t];
               }
            meanX /= predictedSS[t];
            lxy -= mean[t] * meanX * predictedSS[t];
            lxx -= meanX * meanX * predictedSS[t];
            lyy -= mean[t] * mean[t] * predictedSS[t];
            double h2 = lxy / lxx * 2;
            double r2 = lxy*lxy / (lxx*lyy);
            fprintf(fp, "%d\t%.3lf\t%s\t%s\t%s\t%d\t%d\t%.3lf\t%.3lf\n",
               chr, (bp[base+31]+bp[base+32])*0.0000005,
               (const char*)snpName[base+31], (const char*)snpName[base+32],
               (const char*)ped.traitNames[traits[t]],
               samplesize[t], predictedSS[t], r2, h2);
         }
      }
      for(int i = 0; i < idCount; i++)
         delete []connection[i];
   }  // end of seg loop
   fclose(fp);
   printf("Local heritability inference ends at %s", currentTime());
   printf("Local heritability inference results saved in file %s\n", (const char*)outfile);
}

void Engine::IBDMI()
{
   int nperm;
   printf("\nOptions in effect:\n");
   printf("\t--ibdMI\n");
   if(mincons)
      printf("\t--mincons %d\n", mincons);
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

   printf("IBD-based MI starts at %s", currentTime());
   IntArray allpairs(0);
   const int BLOCKSIZE = 8;
   for(int s = 0; s < idCount; s += BLOCKSIZE)
      for(int s2 = s; s2 < s+BLOCKSIZE && s2 < idCount; s2++)
         for(int t = s; t < idCount; t += BLOCKSIZE)
            for(int t2 = t; t2 < t+BLOCKSIZE && t2 < idCount; t2++){
               if(t==s && t2 <= s2) continue;
               allpairs.Push(s2);
               allpairs.Push(t2);
            }
   char header[256];
   sprintf(header, "SNP\tChr\tPos\tN_IBD\tN_Inf\tN_MI\tR_InfMI\tRate_MI\n");
   String outfile(prefix);
   outfile.Add(".mi");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "%s", header);

   char revbase[256];
   for(int i = 0; i < 8; i++)
      revbase[shortbase[i]] = i;
   char rightmost[256];
   for(int i = 0; i < 256; i++)
      rightmost[i] = revbase[i&(-i)];
   unsigned long long IBS0, informative, word;
   IntArray ibdsegStorage[2], ibdsegIndex[2];
   IntArray *MI, *ibdCount, *infoCount;
   MI = new IntArray[defaultMaxCoreCount];
   ibdCount = new IntArray[defaultMaxCoreCount];
   infoCount = new IntArray[defaultMaxCoreCount];
   printf("Scanning genome...\n");
   int segCount = (chrSeg.Length()>>2);
   for(int seg = 0; seg < segCount; seg++){
      int chrsegMin = chrSeg[seg<<2];
      int chrsegMax = chrSeg[(seg<<2)|1];
      if(seg == 0 || chromosomes[chrsegMin<<6]!=chromosomes[chrSeg[(seg-1)<<2]<<6])
         printf("  Chr %d...\n", chromosomes[chrsegMin<<6]);
      int ndim = chrsegMax - chrsegMin + 1;
      if(mincons)
         IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], false, 2500000, mincons);
      else
         IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1]);
      int tempcount = (ibdsegStorage[0].Length()>>1);
      int thread = 0;
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) private(thread, IBS0, informative, word)
{
      thread = omp_get_thread_num();
#endif
      MI[thread].Dimension(ndim<<6);
      MI[thread].Zero();
      infoCount[thread].Dimension(ndim<<6);
      infoCount[thread].Zero();
      ibdCount[thread].Dimension(ndim);
      ibdCount[thread].Zero();
#ifdef _OPENMP
   #pragma omp for
#endif
      for(int i = 0; i < tempcount; i++){
         int index = ibdsegIndex[0][i];
         int id1 = allpairs[index*2];     // unaff
         int id2 = allpairs[index*2+1];   // aff
         int startword = ((ibdsegStorage[0][i*2]-1)>>6)+1;
         int stopword = ((ibdsegStorage[0][i*2+1]+1)>>6)-1;
         for(int w = startword; w <= stopword; w++){
            ibdCount[thread][w-chrsegMin]++;
            IBS0 = LG[0][id1][w] & LG[0][id2][w] & (LG[1][id1][w] ^ LG[1][id2][w]);
            informative = (LG[0][id1][w] | LG[0][id2][w]) & (LG[1][id1][w] & LG[1][id2][w]);
            for(int k = 0; k < 8; k++){
               int base = (((w-chrsegMin)<<6)+(k<<3));
               for(word = (IBS0 & 0xFF); word; word &= (word-1))
                  MI[thread][base+rightmost[word]]++;
               for(word = (informative & 0xFF); word; word &= (word-1))
                  infoCount[thread][base+rightmost[word]]++;
               IBS0 >>= 8;
               informative >>= 8;
            }
         }
      }
#ifdef _OPENMP
}
#endif
      int chr = chromosomes[chrsegMin<<6];
      for(int w = 0; w < ndim; w++){
         int base = (chrsegMin+w)*64;
         int local_ibdCount = ibdCount[0][w];
         for(int c = 1; c < defaultMaxCoreCount; c++)
            local_ibdCount += ibdCount[c][w];
         for(int k = 0; k < 64; k++){
            int m = w*64+k;
            int local_MI = MI[0][m];
            int local_infoCount = infoCount[0][m];
            for(int c = 1; c < defaultMaxCoreCount; c++){
               local_MI += MI[c][m];
               local_infoCount += infoCount[c][m];
            }
            if(!local_ibdCount || !local_MI) continue;
            double error = local_MI * 1.0 / local_ibdCount;
            double infErr = local_MI * 1.0 / (local_infoCount+local_MI);
            if(error > 0.001)
               fprintf(fp, "%s\t%d\t%.3lf\t%d\t%d\t%d\t%.3lf\t%.3lf\n",
                  (const char*)snpName[base+k], chr, bp[base+k]*0.000001,
                  local_ibdCount, local_infoCount+local_MI, local_MI, infErr, error);
         }
      }  // end of word loop
   }  // end of seg loop
   fclose(fp);
   printf("IBD-based MI ends at %s", currentTime());
   printf("IBD-based MI results saved in file %s\n", (const char*)outfile);
}

void Engine::PopulationIBD()
{
   printf("\nOptions in effect:\n");
   printf("\t--popibd\n");
   if(mincons)
      printf("\t--mincons %d\n", mincons);
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

   IntArray ibdsegStorage[2], ibdsegIndex[2];
   IntArray allpairs;
   IntArray *pi = new IntArray [popCount];
   printf("Scanning genome...\n");
   String outfile(prefix);
   outfile.Add(".popibd");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "Chr\tPos\tFlankSNP1\tFlankSNP2");
   for(int pop = 0; pop < popCount; pop++)
      fprintf(fp, "\tPi_%d", pop+1);
   fprintf(fp, "\n");
   const int BLOCKSIZE = 8;
   int segCount = (chrSeg.Length()>>2);
   for(int seg = 0; seg < segCount; seg++){
      int chrsegMin = chrSeg[seg<<2];
      int chrsegMax = chrSeg[(seg<<2)|1];
      if(seg == 0 || chromosomes[chrsegMin<<6]!=chromosomes[chrSeg[(seg-1)<<2]<<6])
         printf("  Chr %d...\n", chromosomes[chrsegMin<<6]);
      int ndim = chrsegMax - chrsegMin + 1;
      int chr = chromosomes[chrsegMin<<6];
      for(int pop = 0; pop < popCount; pop++){
         pi[pop].Dimension(ndim);
         pi[pop].Zero();
         allpairs.Dimension(0);
         for(int s = 0; s < affCount[pop]; s += BLOCKSIZE)
            for(int s2 = s; s2 < s+BLOCKSIZE && s2 < affCount[pop]; s2++)
               for(int t = s; t < affCount[pop]; t += BLOCKSIZE)
                  for(int t2 = t; t2 < t+BLOCKSIZE && t2 < affCount[pop]; t2++){
                     if(t==s && t2 <= s2) continue;
                     allpairs.Push(aff[pop][s2]);
                     allpairs.Push(aff[pop][t2]);
                  }
         if(mincons)
            IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], false, 2500000, mincons);
         else
            IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1]);
         for(int k = 0; k < 2; k++){
            int localvalue = k+1;
            int tempcount = (ibdsegStorage[k].Length()>>1);
            for(int i = 0; i < tempcount; i++){
               int startword = ((ibdsegStorage[k][i*2]-1)>>6)+1;
               int stopword = ((ibdsegStorage[k][i*2+1]+1)>>6)-1;
               int w = startword-chrsegMin-1;
               if(w >= 0 && ibdsegStorage[k][i*2] <= (startword<<6)-33)
                  pi[pop][w] += localvalue;
               int localcount = stopword-chrsegMin+1;
               for(w++; w < localcount; w++)
                  pi[pop][w] += localvalue;
               if(w < ndim && (ibdsegStorage[k][i*2+1]&0x3F) >= 32)
                  pi[pop][w] += localvalue;
            }
         }
      }  // end of pop loop
      for(int w = 0; w < ndim; w++){
         int base = ((w+chrsegMin)<<6);
         fprintf(fp, "%d\t%.3lf\t%s\t%s",
            chr, (bp[base+31]+bp[base+32])*0.0000005,
            (const char*)snpName[base+31], (const char*)snpName[base+32]);
         for(int pop = 0; pop < popCount; pop++){
            double localpi = pi[pop][w] * 0.5 / (affCount[pop]*(affCount[pop]-1)/2);
            fprintf(fp, "\t%.3lf", localpi);
         }
         fprintf(fp, "\n");
      }  // end of word w loop
   }  // end of seg loop
   fclose(fp);
   printf("Population IBD analysis ends at %s", currentTime());
   printf("Population IBD saved in file %s\n", (const char*)outfile);
}

void Engine::IBDGDT()
{
   int nperm;
   printf("\nOptions in effect:\n");
   printf("\t--ibdgdt\n");
   if(mincons)
      printf("\t--mincons %d\n", mincons);
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

   printf("IBD-based GDT mapping starts at %s", currentTime());

   IntArray allpairs(0), aff[2];
   int pairCount, affcount[2];
   for(int k = 0; k < 2; k++)
      aff[k].Dimension(0);
   for(int i = 0; i < idCount; i++)
      if(ped[phenoid[i]].affections[0] == 2)
         aff[1].Push(i);
      else if(ped[phenoid[i]].affections[0] == 1)
         aff[0].Push(i);
   for(int k = 0; k < 2; k++)
      affcount[k] = aff[k].Length();

   const int BLOCKSIZE = 8;
   for(int s = 0; s < affcount[0]; s += BLOCKSIZE)
      for(int s2 = s; s2 < s+BLOCKSIZE && s2 < affcount[0]; s2++)
         for(int t = 0; t < affcount[1]; t += BLOCKSIZE)
            for(int t2 = t; t2 < t+BLOCKSIZE && t2 < affcount[1]; t2++){
               allpairs.Push(aff[0][s2]);
               allpairs.Push(aff[1][t2]);
            }
   pairCount = allpairs.Length()>>1;
   printf("%d affected, %d unaffected, #%d DRPs are used for IBD-based GDT mapping.\n",
      affcount[1], affcount[0], pairCount);
   char header[256];
   sprintf(header, "SNP\tChr\tPos\tT\tNT\tChisq\n");
   String outfile(prefix);
   outfile.Add(".ibdgdt");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "%s", header);

   unsigned long long informative, word;
   unsigned long long *T, *NT;
   IntArray ibdsegStorage[2], ibdsegIndex[2];
   int TCount[4], NTCount[4];
   printf("Scanning genome...\n");
   int segCount = (chrSeg.Length()>>2);
   for(int seg = 0; seg < segCount; seg++){
      int chrsegMin = chrSeg[seg<<2];
      int chrsegMax = chrSeg[(seg<<2)|1];
      if(seg == 0 || chromosomes[chrsegMin<<6]!=chromosomes[chrSeg[(seg-1)<<2]<<6])
         printf("  Chr %d...\n", chromosomes[chrsegMin<<6]);
      int ndim = chrsegMax - chrsegMin + 1;
      T = new unsigned long long [ndim<<5];
      NT = new unsigned long long [ndim<<5];
      for(int m = 0; m < (ndim<<5); m++)
         T[m]=NT[m]=0;
      if(mincons)
         IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], false, 2500000, mincons);
      else
         IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1]);
      int tempcount = (ibdsegStorage[0].Length()>>1);
      for(int i = 0; i < tempcount; i++){
         int index = ibdsegIndex[0][i];
         int id1 = allpairs[index*2];     // unaff
         int id2 = allpairs[index*2+1];   // aff
         int startword = ((ibdsegStorage[0][i*2]-1)>>6)+1;
         int stopword = ((ibdsegStorage[0][i*2+1]+1)>>6)-1;
         for(int w = startword; w <= stopword; w++){
            informative = (LG[0][id1][w] | LG[1][id1][w]) &  // id1 not missing
                     (LG[0][id2][w] | LG[1][id2][w]) &      // id2 not missing
                     (LG[0][id1][w] ^ LG[0][id2][w]);       // Aa vs AA or aa
            int bbb = (w-chrsegMin)*32;

            word = informative & ((LG[0][id2][w] & LG[1][id2][w]) | (~LG[1][id1][w])); // Aff=AA OR Unaff=aa
            T[bbb] += (word&0x0000000100000001);
            for(int k = 1; k < 32; k++)
               T[bbb|k] += ((word>>k)&0x0000000100000001);

            word = informative & ((LG[0][id1][w] & LG[1][id1][w]) | (~LG[1][id2][w])); // Aff=aa OR Unaff=AA
            NT[bbb] += (word&0x0000000100000001);
            for(int k = 1; k < 32; k++)
               NT[bbb|k] += ((word>>k)&0x0000000100000001);
         }
      }
      int chr = chromosomes[chrsegMin<<6];
      for(int w = 0; w < ndim; w++){
         int base = (chrsegMin+w)*64;
         for(int k = 0; k < 2; k++)
            for(int j = 0; j < 32; j++){
               TCount[k] = (T[w*32+j]>>(k*32))&0xFFFFFFFF;
               NTCount[k] = (NT[w*32+j]>>(k*32))&0xFFFFFFFF;
               int pos = base+k*32+j;
               double tdt = (TCount[k]-NTCount[k]);
               if(TCount[k] || NTCount[k])
                  tdt = tdt*tdt/(TCount[k]+NTCount[k]);
               else
                  tdt = 0;
               fprintf(fp, "%s\t%d\t%.3lf\t%d\t%d\t%.2lf\n",
                  (const char*)snpName[pos], chr, bp[pos]*0.000001,
                  TCount[k], NTCount[k], tdt);
            }
         }  // end of word loop
      delete []T;
      delete []NT;
   }  // end of seg loop
   fclose(fp);
   printf("IBD mapping ends at %s", currentTime());
   printf("Genome-wide IBD-mapping scan results saved in file %s\n", (const char*)outfile);
}

void Engine::PopulationDistance()
{
   printf("\nOptions in effect:\n");
   printf("\t--popdist\n");
   if(mincons)
      printf("\t--mincons %d\n", mincons);
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
   printf("Population distance analysis starts at %s", currentTime());

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

   int distCount = popCount*(popCount+1)/2;
   Vector popdist(distCount);
   popdist.Zero();
   IntArray ibdsegStorage[2], ibdsegIndex[2];
   IntArray allpairs;
   const int BLOCKSIZE = 8;
   printf("Scanning genome...\n");

   int segCount = (chrSeg.Length()>>2);
   for(int seg = 0; seg < segCount; seg++){
      int index = 0;
      for(int p1 = 0; p1 < popCount; p1++)
         for(int p2 = p1; p2 < popCount; p2++){
            allpairs.Dimension(0);
            if(p1==p2)
            for(int s = 0; s < affCount[p1]; s += BLOCKSIZE)
               for(int s2 = s; s2 < s+BLOCKSIZE && s2 < affCount[p1]; s2++)
                  for(int t = s; t < affCount[p2]; t += BLOCKSIZE)
                     for(int t2 = t; t2 < t+BLOCKSIZE && t2 < affCount[p2]; t2++){
                        if(t==s && t2 <= s2) continue;
                        allpairs.Push(aff[p1][s2]);
                        allpairs.Push(aff[p2][t2]);
               }
            else
            for(int s = 0; s < affCount[p1]; s += BLOCKSIZE)
               for(int s2 = s; s2 < s+BLOCKSIZE && s2 < affCount[p1]; s2++)
                  for(int t = 0; t < affCount[p2]; t += BLOCKSIZE)
                     for(int t2 = t; t2 < t+BLOCKSIZE && t2 < affCount[p2]; t2++){
                        allpairs.Push(aff[p1][s2]);
                        allpairs.Push(aff[p2][t2]);
               }
            if(mincons)
               IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], true, 2500000, mincons);
            else
               IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], true);
            int pairCount = allpairs.Length()/2;
            double temp = 0.0;
            for(int pair = 0; pair < pairCount; pair++)
               temp += (ibdsegStorage[0][pair]*0.5 + ibdsegStorage[1][pair])*0.000001;
            popdist[index] += temp;
            index++;
         }  // end of p2 loop
   }  // end of seg loop
   int index = 0;
   for(int p1 = 0; p1 < popCount; p1++)
      for(int p2 = p1; p2 < popCount; p2++){
         popdist[index] /= totalLength;
         if(p1==p2)
            popdist[index] /= (affCount[p1]*(affCount[p2]-1)/2);
         else
            popdist[index] /= (affCount[p1]*affCount[p2]);
         popdist[index] = 1 - popdist[index];
         index ++;
      }

   String outfile(prefix);
   outfile.Add(".dst");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "POP1\tPOP2\tN\tPropIBD\tDistIBD\n");
   printf("PropIBD");
   for(int p = 0; p < popCount; p++)
      printf("\t%d", p+1);
   printf("\n");
   index = 0;
   for(int p1 = 0; p1 < popCount; p1++){
      printf("%d", p1+1);
      for(int p2 = 0; p2 < p1; p2++)
         printf("\t");
      for(int p2 = p1; p2 < popCount; p2++){
         fprintf(fp, "%d\t%d\t%d\t%.4lf\t%.4lf\n",
            p1+1, p2+1,
            p1==p2? (affCount[p1]*(affCount[p2]-1)/2):(affCount[p1]*affCount[p2]),
            1-popdist[index], popdist[index]);
         printf("\t%.4lf", 1-popdist[index]);
         index ++;
      }
      printf("\n");
   }
   fclose(fp);
   printf("Population distance analysis ends at %s", currentTime());
   printf("Population distances saved in file %s\n", (const char*)outfile);
}

void Engine::AncestryInference()
{
   printf("\nOptions in effect:\n");
   printf("\t--ancestry\n");
   if(mincons)
      printf("\t--mincons %d\n", mincons);
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

   if(ped.affectionCount==0){
      printf("Reference populations need to be specified as affection status.\n");
      return;
   }
   printf("Ancestry inference starts at %s", currentTime());
   IntArray aff[3], ancestry[4];
   for(int k = 0; k < 3; k++)
      aff[k].Dimension(0);
   for(int i = 0; i < idCount; i++){
      if(ped[phenoid[i]].affections[0] == 2)
         aff[2].Push(i);
      else if(ped[phenoid[i]].affections[0] == 1)
         aff[1].Push(i);
      else
         aff[0].Push(i);
   }
   int affCount[3];
   for(int a = 0; a < 3; a++)
      affCount[a] = aff[a].Length();
   if(!affCount[0]) {printf("No ancestry to be inferred. Note these samples should have affection status -9.\n"); return;}
   if(!affCount[1]) {printf("No ancestry 1 reference samples. Note these samples should have affection status 1.\n"); return;}
   if(!affCount[2]) {printf("No ancestry 2 reference samples. Note these samples should have affection status 2.\n"); return;}
   printf("Reference 1 sample count: %d\n", affCount[1]);
   printf("Reference 2 sample count: %d\n", affCount[2]);
   printf("To-be-inferred sample count: %d\n", affCount[0]);

   for(int a = 0; a < 4; a++){
      ancestry[a].Dimension(affCount[0]);
      ancestry[a].Zero();
   }
   IntArray *ibdsegStorage[2][2], *ibdsegIndex[2][2];
   for(int i = 0; i < 2; i++)
      for(int j = 0; j < 2; j++){
         ibdsegStorage[i][j] = new IntArray [affCount[0]];
         ibdsegIndex[i][j] = new IntArray [affCount[0]];
      }
   IntArray allpairs;
   int segCount = (chrSeg.Length()>>2);
   printf("Scanning genome...\n");
   for(int seg = 0; seg < segCount; seg++){
      int chrsegMin = chrSeg[seg<<2];
      int chrsegMax = chrSeg[(seg<<2)|1];
      int ndim = chrsegMax - chrsegMin + 1;
      for(int a = 0; a < 2; a++)   // affection status
         for(int person = 0; person < affCount[0]; person++){
            allpairs.Dimension(0);
            int id = aff[0][person];
            for(int j = 0; j < affCount[a+1]; j++){
               allpairs.Push(id);
               allpairs.Push(aff[a+1][j]);
            }
            if(mincons)
               IBDSegOnly(allpairs, seg, ibdsegStorage[a][0][person], ibdsegIndex[a][0][person], ibdsegStorage[a][1][person], ibdsegIndex[a][1][person], false, 2500000, mincons);
            else
               IBDSegOnly(allpairs, seg, ibdsegStorage[a][0][person], ibdsegIndex[a][0][person], ibdsegStorage[a][1][person], ibdsegIndex[a][1][person]);
         }
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount)
#endif
      for(int person = 0; person < affCount[0]; person++){
         IntArray *connections[2], connection[2];
         for(int a = 0; a < 2; a++){   // affection status
            IntArray allpairs(0);
            int id = aff[0][person];
            for(int j = 0; j < affCount[a+1]; j++){
               allpairs.Push(id);
               allpairs.Push(aff[a+1][j]);
            }
            connection[a].Dimension(ndim);
            connection[a].Zero();
            connections[a] = new IntArray[ndim];
            for(int w = 0; w < ndim; w++)
               connections[a][w].Dimension(0);
            for(int k = 0; k < 2; k++){
               int tempcount = (ibdsegStorage[a][k][person].Length()>>1);
               for(int i = 0; i < tempcount; i++){
                  int startword = ((ibdsegStorage[a][k][person][i*2]-1)>>6)+1;
                  int stopword = ((ibdsegStorage[a][k][person][i*2+1]+1)>>6)-1;
                  int index = ibdsegIndex[a][k][person][i]*2;
                  int id2 = allpairs[index+1];
                  int localcount = stopword-chrsegMin;
                  for(int w = startword-chrsegMin; w <= localcount; w++){
                     if(k==0){
                        connection[a][w] ++;
                        connections[a][w].Push(id2);
                     }else
                        connection[a][w] += 1000000;
                  }
                  int w = startword-chrsegMin-1;
                  if(w >=0 && ibdsegStorage[a][k][person][i*2] <= (startword<<6)-33){
                     if(k==0){
                        connection[a][w] ++;
                        connections[a][w].Push(id2);
                     }else
                        connection[a][w] += 1000000;
                  }
                  w = stopword-chrsegMin+1;
                  if(w < ndim && ibdsegStorage[a][k][person][i*2+1] >= ((stopword<<6)|0x3F)+33){
                     if(k==0){
                        connection[a][w] ++;
                        connections[a][w].Push(id2);
                     }else
                        connection[a][w] += 1000000;
                  }
               }  // end of ibdseg i loop
            }  // end of k loop
         }  // end of affection a loop
         for(int w = 0; w < ndim; w++){
            if((!connection[0][w] && !connection[1][w]) || (connection[0][w]>=1000000 && connection[1][w]>=1000000))
               ancestry[0][person] ++; // NA
            else if(!connection[0][w] || connection[1][w]>=1000000)
               ancestry[3][person] ++; // 22
            else if(!connection[1][w] || connection[0][w]>=1000000)
               ancestry[1][person] ++; // 11
            else{
               int m = chrsegMin+w;
               bool IBS0Flag12 = false;
               unsigned long long IBS0;
               for(int i = 0; i < connection[0][w]; i++){
                  int id1 = connections[0][w][i];
                  for(int j = 0; j < connection[1][w]; j++){
                     int id2 = connections[1][w][j];
                     IBS0 = LG[0][id1][m] & LG[0][id2][m] & (LG[1][id1][m] ^ LG[1][id2][m]);
                     if(!IBS0 && !(IBS0 & (IBS0-1))) {
                        IBS0Flag12 = true;
                        break;
                     }
                  }
                  if(IBS0Flag12) break;
               }
               if(IBS0Flag12)  // Ref1 and Ref2 not IBD
                  ancestry[2][person] ++; // 12
               else{
                  bool IBS0Flag11 = false;
                  for(int i = 0; i < connection[0][w]; i++){
                     int id1 = connections[0][w][i];
                     for(int j = i+1; j < connection[0][w]; j++){
                        int id2 = connections[0][w][j];
                        IBS0 = LG[0][id1][m] & LG[0][id2][m] & (LG[1][id1][m] ^ LG[1][id2][m]);
                        if(!IBS0 && !(IBS0 & (IBS0-1))) {
                           IBS0Flag11 = true;
                           break;
                        }
                     }
                     if(IBS0Flag11) break;
                  }
                  bool IBS0Flag22 = false;
                  for(int i = 0; i < connection[1][w]; i++){
                     int id1 = connections[1][w][i];
                     for(int j = i+1; j < connection[1][w]; j++){
                        int id2 = connections[1][w][j];
                        IBS0 = LG[0][id1][m] & LG[0][id2][m] & (LG[1][id1][m] ^ LG[1][id2][m]);
                        if(!IBS0 && !(IBS0 & (IBS0-1))) {
                           IBS0Flag22 = true;
                           break;
                        }
                     }
                     if(IBS0Flag22) break;
                  }
                  if(IBS0Flag11 && IBS0Flag22)
                     ancestry[0][person] ++; // NA
                  else if(IBS0Flag11)
                     ancestry[1][person] ++; // 11
                  else if(IBS0Flag22)
                     ancestry[3][person] ++; // 22
                  else  // 1 IBD in Ref1 and Ref2
                     ancestry[0][person] ++; // NA
               }
            }
         }  // end word w loop
         for(int a = 0; a < 2; a++)
            delete []connections[a];
      }  // end of person loop
   }  // end of seg loop
   String outfile(prefix);
   outfile.Add(".anc");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "FID\tIID\tNMISS\tNMIX\tN1\tN2\tAnc_P1\tAnc_P2\tAdmix\tAncestry\n");
   double parent[2];
   for(int person = 0; person < affCount[0]; person++){
      int id = phenoid[aff[0][person]];
      int total = ancestry[1][person] + ancestry[2][person] + ancestry[3][person];
      double Anc = (ancestry[3][person] + ancestry[2][person]*0.5) / total;
      double Adm = ancestry[2][person] * 1.0 / total;
      double temp = Anc*Anc - Anc + Adm*0.5;
      temp = temp >= 0? sqrt(temp): 0;
      parent[0] = Anc - temp;
      parent[1] = Anc + temp;
      for(int k = 0; k < 2; k++)
         if(parent[k]<0) {
            parent[k] = 0;
            parent[1-k] = Anc * 2;
         }else if(parent[k]>1) {
            parent[k] = 1;
            parent[1-k] = Anc * 2 - 1;
         }
      fprintf(fp, "%s\t%s\t%d\t%d\t%d\t%d\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n",
         (const char*)ped[id].famid, (const char*)ped[id].pid,
         ancestry[0][person], ancestry[2][person], ancestry[1][person], ancestry[3][person],
         parent[0], parent[1], Adm, Anc);
   }
   fclose(fp);
   printf("Ancestry inference ends at %s", currentTime());
   printf("Ancestry results saved in file %s\n", (const char*)outfile);
}

void Engine::AUCpredicting(IntArray & allchr, IntArray & allpos)
{
   printf("\nOptions in effect:\n");
   printf("\t--aucmap\n");
   if(mincons)
      printf("\t--mincons %d\n", mincons);
   printf("\t--position");
   for(int i = 0; i < allchr.Length(); i++)
      printf(" %d:%d", allchr[i], allpos[i]);
   printf("\n");
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

   printf("AUC-based IBD mapping starts at %s", currentTime());
   IntArray allpairs[3], aff[2];
   int pairCount[3], affcount[2];
   for(int k = 0; k < 3; k++)
      allpairs[k].Dimension(0);
   for(int k = 0; k < 2; k++)
      aff[k].Dimension(0);
   for(int i = 0; i < idCount; i++){
      if(ped[phenoid[i]].affections[0] == 2)
         aff[1].Push(i);
      else if(ped[phenoid[i]].affections[0] == 1)
         aff[0].Push(i);
   }
   const int BLOCKSIZE = 8;
   for(int k = 0; k < 2; k++){
      affcount[k] = aff[k].Length();
      for(int s = 0; s < affcount[k]; s += BLOCKSIZE)
         for(int s2 = s; s2 < s+BLOCKSIZE && s2 < affcount[k]; s2++)
            for(int t = s; t < affcount[k]; t += BLOCKSIZE)
               for(int t2 = t; t2 < t+BLOCKSIZE && t2 < affcount[k]; t2++){
                  if(t==s && t2 <= s2) continue;
                  allpairs[k*2].Push(aff[k][s2]);
                  allpairs[k*2].Push(aff[k][t2]);
               }
   }
   for(int s = 0; s < affcount[0]; s += BLOCKSIZE)
      for(int s2 = s; s2 < s+BLOCKSIZE && s2 < affcount[0]; s2++)
         for(int t = 0; t < affcount[1]; t += BLOCKSIZE)
            for(int t2 = t; t2 < t+BLOCKSIZE && t2 < affcount[1]; t2++){
               allpairs[1].Push(aff[0][s2]);
               allpairs[1].Push(aff[1][t2]);
            }
   for(int k = 0; k < 3; k++)
      pairCount[k] = allpairs[k].Length()>>1;
   if(pairCount[2] == 0) {printf("No ARPs are found.\n"); return;}
   else printf("#%d URPs, #%d DRPs, #%d ARPs are used for AUC mapping.\n",
      pairCount[0], pairCount[1], pairCount[2]);
   IntArray ibdsegStorage[2], ibdsegIndex[2];
   int posCount = allchr.Length();
   IntArray *connections[2];
   IntArray *connectionCount[2];
   for(int a = 0; a < 2; a++){
      connections[a] = new IntArray [idCount];
      connectionCount[a] = new IntArray[idCount];
      for(int i = 0; i < idCount; i++){
         connections[a][i].Dimension(0);
         connectionCount[a][i].Dimension(posCount);
         connectionCount[a][i].Zero();
      }
   }
   IntArray affections(idCount);
   for(int i = 0; i < idCount; i++)
      affections[i] = ped[phenoid[i]].affections[0]-1;
   int segCount = (chrSeg.Length()>>2);
   bool stopFlag;
   int chrIndex;
   for(int seg = 0; seg < segCount; seg++){
      int chrsegMin = chrSeg[seg<<2];
      stopFlag = true;
      for(int i = 0; i < posCount; i++)
         if(chromosomes[chrsegMin<<6] == allchr[i]){
            stopFlag = false;
            chrIndex = i;
            break;
         }
      if(stopFlag) continue;
      int chrsegMax = chrSeg[(seg<<2)|1];
      if(bp[chrsegMin<<6] >= allpos[chrIndex] || bp[(chrsegMax<<6)|0x3F] <= allpos[chrIndex])
         continue;
      for(int type = 0; type < 3; type++){
         int afftype[2];
         if(type == 0){
            afftype[0] = afftype[1] = 0;
         }else if(type == 1){
            afftype[0] = 0; afftype[1] = 1;
         }else{
            afftype[0] = afftype[1] = 1;
         }
         if(mincons)
            IBDSegOnly(allpairs[type], seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], false, 2500000, mincons);
         else
            IBDSegOnly(allpairs[type], seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1]);
         for(int k = 0; k < 2; k++){
            int tempcount = (ibdsegStorage[k].Length()>>1);
            for(int i = 0; i < tempcount; i++){
               if(bp[ibdsegStorage[k][i*2]] >= allpos[chrIndex]
                  || bp[ibdsegStorage[k][i*2+1]] <= allpos[chrIndex])
                  continue;
               int index = ibdsegIndex[k][i];
               int id1 = allpairs[type][index*2];
               int id2 = allpairs[type][index*2+1];
               connectionCount[afftype[1]][id1][chrIndex] ++;
               connectionCount[afftype[0]][id2][chrIndex] ++;
               bool exist = false;
               for(int j = 0; j < connections[afftype[1]][id1].Length(); j++)
                  if(connections[afftype[1]][id1][j] == id2)
                     exist = true;
               if(!exist)
                  connections[afftype[1]][id1].Push(id2);
               exist = false;
               for(int j = 0; j < connections[afftype[0]][id2].Length(); j++)
                  if(connections[afftype[0]][id2][j] == id1)
                     exist = true;
               if(!exist)
                  connections[afftype[0]][id2].Push(id1);
            }  // end of segment i loop
         }  // end of IBD1 or IBD2 k loop
      }  // end of pair type loop
   }  // end of seg loop
   double AUC;
   int NAA, NUU, NUA, NMISS;
   NAA=NUU=NUA=NMISS=0;
   Vector risks(idCount);
   int connection[2];
   for(int i = 0; i < idCount; i++){
      for(int a = 0; a < 2; a++)
         connection[a] = connections[a][i].Length();
      if(connection[0] + connection[1] > 0)  // 10 or more relatives
         risks[i] = connection[1]*1.0 / (connection[0] + connection[1]);
      else {
         risks[i] = -9;
         NMISS ++;
      }
      if(affections[i]==0) {
         NUU += connection[0];
         NUA += connection[1];
      }else if(affections[i] == 1)
         NAA += connection[1];
   }  // AUC of Aff and connection[1]/(connection[0]+connection[1])
   NUU /= 2;
   NAA /= 2;
   AUC = ComputeAUC(risks, affections, 0);
   printf("Chr");
   for(int i = 0; i < posCount; i++)
      printf("\t%d:%d", allchr[i], allpos[i]);
   printf("\n");
   printf("Unaffected IBD relative pair count: %d\n", NUU);
   printf("Discordant IBD relative pair count: %d\n", NUA);
   printf("Affected IBD relative pair count: %d\n", NAA);
   printf("No relative count: %d\n", NMISS);
   printf("AUC: %.3lf\n", AUC);

   if(posCount>1){
      for(int p = 0; p < posCount; p++){
         for(int i = 0; i < idCount; i++){
            for(int a = 0; a < 2; a++)
               connection[a] = connectionCount[a][i][p];
            if(connection[0] + connection[1] > 0)
               risks[i] = connection[1]*1.0 / (connection[0] + connection[1]);
            else
               risks[i] = -9;
         }
         AUC = ComputeAUC(risks, affections, 0);
         printf("Chr %d:%d\tAUC=%.3lf\n", allchr[p], allpos[p], AUC);
      }
      for(int i = 0; i < idCount; i++){
         connection[0]=connection[1]=0;
         for(int a = 0; a < 2; a++)
            for(int p = 0; p < posCount; p++)
               connection[a] += connectionCount[a][i][p];
         if(connection[0] + connection[1] > 0)
            risks[i] = connection[1]*1.0 / (connection[0] + connection[1]);
         else
            risks[i] = -9;
      }
      AUC = ComputeAUC(risks, affections, 0);
      printf("Additive model: AUC=%.3lf\n", AUC);

      for(int i = 0; i < idCount; i++){
         risks[i] = 1;
         bool validFlag = false;
         for(int p = 0; p < posCount; p++){
            for(int a = 0; a < 2; a++)
               connection[a] = connectionCount[a][i][p];
            if(connection[0] || connection[1]){
               risks[i] *= connection[1]*1.0 / (connection[0] + connection[1]);
               validFlag = true;
            }
         }
         if(!validFlag) risks[i] = -9;
      }
      AUC = ComputeAUC(risks, affections, 0);
      printf("Multiplicative model: AUC=%.3lf\n", AUC);

      for(int i = 0; i < idCount; i++){
         risks[i] = -9;
         for(int p = 0; p < posCount; p++){
            for(int a = 0; a < 2; a++)
               connection[a] = connectionCount[a][i][p];
            if(connection[0] || connection[1]){
               double temp = connection[1]*1.0 / (connection[0] + connection[1]);
               if(temp > risks[i]) risks[i] = temp;
            }
         }
      }
      AUC = ComputeAUC(risks, affections, 0);
      printf("Max model: AUC=%.3lf\n", AUC);

      double T = aff[1].Length()*1.0 / (aff[0].Length()+aff[1].Length());
      for(int i = 0; i < idCount; i++){
         risks[i] = -9;
         for(int a = 0; a < 2; a++)
            connection[a] = connectionCount[a][i][0];
         if(connection[0] || connection[1])
            risks[i] = connection[1]*1.0 / (connection[0] + connection[1]);
         if(risks[i] < T)
         for(int p = 1; p < posCount; p++){
            for(int a = 0; a < 2; a++)
               connection[a] = connectionCount[a][i][p];
            if(connection[0] || connection[1]){
               double temp = connection[1]*1.0 / (connection[0] + connection[1]);
               if(temp > risks[i]) risks[i] = temp;
            }
         }
      }
      AUC = ComputeAUC(risks, affections, 0);
      printf("Position1-preferred model: AUC=%.3lf\n", AUC);
   }
   printf("AUC-based IBD mapping ends at %s", currentTime());
}

void Engine::AUCmapping()
{
   printf("\nOptions in effect:\n");
   printf("\t--aucmap\n");
   if(mincons)
      printf("\t--mincons %d\n", mincons);
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
   printf("IBD-based AUC mapping starts at %s", currentTime());

   IntArray affections(idCount);
   for(int i = 0; i < idCount; i++)
      affections[i] = ped[phenoid[i]].affections[0]-1;
   const int BLOCKSIZE = 256;
   IntArray allpairs[3];
   int pairCount[3], affCount[2];
   affCount[0] = affCount[1] = 0;
   for(int i = 0; i < idCount; i++)
      if(affections[i] == 1)
         affCount[1] ++;
      else if(affections[i] == 0)
         affCount[0] ++;
   pairCount[0] = affCount[0] * (affCount[0]-1) / 2;
   pairCount[1] = affCount[0] * affCount[1];
   pairCount[2] = affCount[1] * (affCount[1]-1) / 2;
   if(pairCount[2] == 0) {printf("No ARPs are found.\n"); return;}
   else printf("#%d URPs, #%d DRPs, #%d ARPs are used for AUC mapping.\n",
      pairCount[0], pairCount[1], pairCount[2]);
   char header[256];
   sprintf(header, "Chr\tPos\tFlankSNP1\tFlankSNP2\tPI_URP\tPI_DRP\tPI_ARP\tNMISS\tN_URP\tN_DRP\tN_ARP\tSuccess\tAUC\n");
   bool headerPrinted=false;
   String outfile(prefix);
   outfile.Add(".aucmap");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "%s", header);
   IntArray ibdsegStorage[2], ibdsegIndex[2], pi[3];
   IntArray *connection[2];
   for(int a = 0; a < 2; a++)
      connection[a] = new IntArray [idCount];
   IntArray subpairs;
   int segCount = (chrSeg.Length()>>2);
   printf("Scanning genome...\n");
   int afftype[2];
   for(int seg = 0; seg < segCount; seg++){
      int chrsegMin = chrSeg[seg<<2];
      int chrsegMax = chrSeg[(seg<<2)|1];
      int ndim = chrsegMax - chrsegMin + 1;
      for(int a = 0; a < 2; a++)
         for(int i = 0; i < idCount; i++){
            connection[a][i].Dimension(ndim);
            connection[a][i].Zero();
         }
      for(int type = 0; type < 3; type++){
         pi[type].Dimension(ndim);
         pi[type].Zero();
      }
      for(int s = 0; s < idCount; s += BLOCKSIZE){
         int sMax = s+BLOCKSIZE < idCount? s+BLOCKSIZE: idCount;
         for(int t = s; t < idCount; t += BLOCKSIZE){
            int tMax = t+BLOCKSIZE < idCount? t+BLOCKSIZE: idCount;
            for(int type = 0; type < 3; type ++)
               allpairs[type].Dimension(0);
            for(int s2 = s; s2 < sMax; s2++)
               for(int t2 = t; t2 < tMax; t2++){
                  if(t==s && t2 <= s2) continue;
                  if(affections[s2]<0 || affections[t2]<0) continue;
                  int type = affections[s2]+affections[t2];
                  if(type == 1 && affections[s2]){
                     allpairs[type].Push(t2);
                     allpairs[type].Push(s2);
                  }else{
                     allpairs[type].Push(s2);
                     allpairs[type].Push(t2);
                  }
               }
            for(int type = 0; type < 3; type++){
               if(type == 0){
                  afftype[0] = afftype[1] = 0;
               }else if(type == 1){
                  afftype[0] = 0; afftype[1] = 1;
               }else{
                  afftype[0] = afftype[1] = 1;
               }
               if(mincons)
                  IBDSegOnly(allpairs[type], seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], false, 2500000, mincons);
               else
                  IBDSegOnly(allpairs[type], seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1]);
               for(int k = 0; k < 2; k++){
                  int localvalue = k+1;
                  int tempcount = (ibdsegStorage[k].Length()>>1);
                  for(int i = 0; i < tempcount; i++){
                     int index = ibdsegIndex[k][i];
                     int id1 = allpairs[type][index*2];
                     int id2 = allpairs[type][index*2+1];
                     int startword = ((ibdsegStorage[k][i*2]-1)>>6)+1;
                     int stopword = ((ibdsegStorage[k][i*2+1]+1)>>6)-1;
                     int startw = startword-chrsegMin;
                     if(startw >= 1 && ibdsegStorage[k][i*2] <= (startword<<6)-33)
                        startw --;
                     int stopw = stopword-chrsegMin;
                     if(stopw+1 < ndim && (ibdsegStorage[k][i*2+1]&0x3F) >= 32)
                        stopw ++;
                     for(int w = startw; w <= stopw; w++){
                        pi[type][w] += localvalue;
                        connection[afftype[1]][id1][w] += localvalue;
                        connection[afftype[0]][id2][w] += localvalue;
                     }  // end of w loop for word
                  }  // end of segment i loop
               }  // end of IBD1 or IBD2 k loop
            }  // end of type loop for IBD type
         }  // end of t loop for pair2
      }  // end of s type loop for pair1
      int chr = chromosomes[chrsegMin<<6];
      double localPi[3];
      Vector AUC(ndim);
      IntArray NAA(ndim);
      IntArray NUU(ndim);
      IntArray NUA(ndim);
      IntArray NMISS(ndim);
      NAA.Zero();
      NUU.Zero();
      NUA.Zero();
      NMISS.Zero();
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount)
{
#endif
   Vector risks(idCount);
#ifdef _OPENMP
   #pragma omp for
#endif
      for(int w = 0; w < ndim; w++){
         for(int i = 0; i < idCount; i++){
            if(connection[0][i][w] || connection[1][i][w])
               risks[i] = connection[1][i][w]*1.0 / (connection[0][i][w] + connection[1][i][w]);
            else {
               risks[i] = -9;
               if(affections[i]>=0) NMISS[w]++;
            }
            if(affections[i]==0) {
               NUU[w] += connection[0][i][w];
               NUA[w] += connection[1][i][w];
            }else if(affections[i] == 1)
               NAA[w] += connection[1][i][w];
         }  // AUC of Aff and connection[1]/(connection[0]+connection[1])
         NUU[w] = int(NUU[w]*0.25+0.5);
         NUA[w] = int(NUA[w]*0.5+0.5);
         NAA[w] = int(NAA[w]*0.25+0.5);
         AUC[w] = ComputeAUC(risks, affections, 0);
      }
#ifdef _OPENMP
}
#endif
      double maxAUC = 0.5999;
      char buffer[256];
      for(int w = 0; w < ndim; w++){
         for(int type = 0; type < 3; type++)
            localPi[type] = pi[type][w] * 0.25 / pairCount[type];
         double success = 1 - NMISS[w]*1.0/(affCount[0]+affCount[1]);
         int base = ((w+chrsegMin)<<6);
         fprintf(fp, "%d\t%.3lf\t%s\t%s\t%.2G\t%.2G\t%.2G\t%d\t%d\t%d\t%d\t%.3lf\t%.3lf\n",
            chr, (bp[base+31]+bp[base+32])*0.0000005,
            (const char*)snpName[base+31], (const char*)snpName[base+32],
            localPi[0], localPi[1], localPi[2], NMISS[w], NUU[w], NUA[w], NAA[w], success, AUC[w]);
         if(AUC[w] > maxAUC && success > 0.9){
            maxAUC = AUC[w];
            sprintf(buffer, "%d\t%.3lf\t%s\t%s\t%.2G\t%.2G\t%.2G\t%d\t%d\t%d\t%d\t%.3lf\t%.3lf\n",
            chr, (bp[base+31]+bp[base+32])*0.0000005,
            (const char*)snpName[base+31], (const char*)snpName[base+32],
            localPi[0], localPi[1], localPi[2], NMISS[w], NUU[w], NUA[w], NAA[w], success, AUC[w]);
         }
      }  // end of word loop
      if(maxAUC >= 0.6) {
         if(!headerPrinted){
            printf("%s", header);
            headerPrinted = true;
         }
         printf("%s", buffer);
      }
   }  // end of seg loop
   fclose(fp);
   printf("IBD-based AUC mapping ends at %s", currentTime());
   printf("Genome-wide AUC-mapping scan results saved in file %s\n", (const char*)outfile);
}
               /*
               int w = startword-chrsegMin-1;
               if(w >= 0 && ibdsegStorage[k][i*2] <= (startword<<6)-33){
                  pi[type][w] += localvalue;
                  connection[afftype[1]][id1][w] += localvalue;
                  connection[afftype[0]][id2][w] += localvalue;
               }
               for(w++; w <= localcount; w++){
                  pi[type][w] += localvalue;
                  connection[afftype[1]][id1][w] += localvalue;
                  connection[afftype[0]][id2][w] += localvalue;
               }
               if(w < ndim && (ibdsegStorage[k][i*2+1]&0x3F) >= 32){
                  pi[type][w] += localvalue;
                  connection[afftype[1]][id1][w] += localvalue;
                  connection[afftype[0]][id2][w] += localvalue;
               } */

                     /*
   IntArray allpairs[3], aff[2];
   int pairCount[3], affcount[2];
   for(int k = 0; k < 3; k++)
      allpairs[k].Dimension(0);
   for(int k = 0; k < 2; k++)
      aff[k].Dimension(0);
   for(int i = 0; i < idCount; i++){
      if(ped[phenoid[i]].affections[0] == 2)
         aff[1].Push(i);
      else if(ped[phenoid[i]].affections[0] == 1)
         aff[0].Push(i);
   }
   const int BLOCKSIZE = 8;
   for(int k = 0; k < 2; k++){
      affcount[k] = aff[k].Length();
      for(int s = 0; s < affcount[k]; s += BLOCKSIZE)
         for(int s2 = s; s2 < s+BLOCKSIZE && s2 < affcount[k]; s2++)
            for(int t = s; t < affcount[k]; t += BLOCKSIZE)
               for(int t2 = t; t2 < t+BLOCKSIZE && t2 < affcount[k]; t2++){
                  if(t==s && t2 <= s2) continue;
                  allpairs[k*2].Push(aff[k][s2]);
                  allpairs[k*2].Push(aff[k][t2]);
               }
   }
   for(int s = 0; s < affcount[0]; s += BLOCKSIZE)
      for(int s2 = s; s2 < s+BLOCKSIZE && s2 < affcount[0]; s2++)
         for(int t = 0; t < affcount[1]; t += BLOCKSIZE)
            for(int t2 = t; t2 < t+BLOCKSIZE && t2 < affcount[1]; t2++){
               allpairs[1].Push(aff[0][s2]);
               allpairs[1].Push(aff[1][t2]);
            }
   for(int k = 0; k < 3; k++)
      pairCount[k] = allpairs[k].Length()>>1;
   if(pairCount[2] == 0) {printf("No ARPs are found.\n"); return;}
   else printf("#%d URPs, #%d DRPs, #%d ARPs are used for AUC mapping.\n",
      pairCount[0], pairCount[1], pairCount[2]);
*/


void Engine::NPL()
{
   printf("\nOptions in effect:\n");
   printf("\t--npl\n");
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

   IntArray allpairs[3], sibs[2], allsibs(0), unrelatedpairs(0);
   int pairCount[3], sibcount[2];
   for(int k = 0; k < 3; k++)
      allpairs[k].Dimension(0);
   for(int f = 0; f < ped.familyCount; f++){
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(ped[i].sibCount > 1 && ped[i].sibs[0]->serial == i){
            for(int k = 0; k < 2; k++) sibs[k].Dimension(0);
            for(int s = 0; s < ped[i].sibCount; s++)
               if(geno[ped[i].sibs[s]->serial] != -1){
                  if(ped[ped[i].sibs[s]->serial].affections[0]==1)
                     sibs[0].Push(geno[ped[i].sibs[s]->serial]);
                  else if(ped[ped[i].sibs[s]->serial].affections[0]==2)
                     sibs[1].Push(geno[ped[i].sibs[s]->serial]);
               }
            for(int k = 0; k < 2; k++){
               sibcount[k] = sibs[k].Length();
               for(int s = 0; s < sibcount[k]; s++)
                  for(int t = s+1; t < sibcount[k]; t++){
                     allpairs[k*2].Push(sibs[k][s]);
                     allpairs[k*2].Push(sibs[k][t]);
                  }
            }
            if(sibcount[1] > 1)
               for(int s = 0; s < sibcount[1]; s++)
                  allsibs.Push(phenoid[sibs[1][s]]);
            for(int s = 0; s < sibcount[0]; s++)
               for(int t = 0; t < sibcount[1]; t++){
                  allpairs[1].Push(sibs[0][s]);
                  allpairs[1].Push(sibs[1][t]);
               }
         }
   }
   int allsibcount = allsibs.Length();
   for(int i = 0; i < allsibcount; i++)
      for(int j = i+1; j < allsibcount; j++)
         if(ped[allsibs[i]].famid != ped[allsibs[j]].famid){
            unrelatedpairs.Push(geno[allsibs[i]]);
            unrelatedpairs.Push(geno[allsibs[j]]);
         }
   int unrelatedpairsCount = unrelatedpairs.Length()/2;

   for(int k = 0; k < 3; k++)
      pairCount[k] = allpairs[k].Length()>>1;
   if(pairCount[2] == 0) {printf("No ASPs are found.\n"); return;}
   else if(pairCount[1] > pairCount[2]*0.1 && pairCount[1] >= 20)
      printf("%d ASPs, %d DSPs and %d unrelated pairs are used for ASP/DSP NPL scan.\n",
         pairCount[2], pairCount[1], unrelatedpairsCount);
   else
      printf("%d ASPs and %d unrelated pairs are used for ASP NPL scan.\n",
         pairCount[2], unrelatedpairsCount);
   char header[256], buffer[256];
   if(pairCount[1] > pairCount[2]*0.1 && pairCount[1] >= 20) // enough DSP
      sprintf(header, "Chr\tPos\tFlankSNP1\tFlankSNP2\tPI_USP\tPI_DSP\tPI_ASP\tPI_AP\tLOD_Raw\tLOD_ASP\tLODwDSP\n");
   else
      sprintf(header, "Chr\tPos\tFlankSNP1\tFlankSNP2\tPI_ASP\tPI_AP\tLOD_Raw\tLOD_ASP\n");
   String outfile(prefix);
   outfile.Add(".npl");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "%s", header);
   bool headerPrinted=false;
   IntArray ibdsegStorage[2], ibdsegIndex[2];
   Vector pi[3], pi0;
   int segCount = (chrSeg.Length()>>2);
   const double LODfactor = 0.5 / log(10.0);
   for(int seg = 0; seg < segCount; seg++){
      int chrsegMin = chrSeg[seg<<2];
      int chrsegMax = chrSeg[(seg<<2)|1];
      int ndim = chrsegMax - chrsegMin + 1;
      IBDSegOnly(unrelatedpairs, seg, ibdsegStorage[0], ibdsegIndex[0],
         ibdsegStorage[1], ibdsegIndex[1], false, 2000000, 10000);
      pi0.Dimension(ndim);
      pi0.Zero();
      for(int k = 0; k < 2; k++){
         int localvalue = k+1;
         int tempcount = (ibdsegStorage[k].Length()>>1);
         for(int i = 0; i < tempcount; i++){
            int startword = ((ibdsegStorage[k][i*2]-1)>>6)+1;
            int stopword = ((ibdsegStorage[k][i*2+1]+1)>>6)-1;
            int localcount = stopword-chrsegMin;
            int w = startword-chrsegMin-1;
            if(w >= 0 && ibdsegStorage[k][i*2] <= (startword<<6)-33)
               pi0[w] += localvalue;
            for(w++; w <= localcount; w++)
               pi0[w] += localvalue;
            if(w < ndim && (ibdsegStorage[k][i*2+1]&0x3F) >= 32)
               pi0[w] += localvalue;
         }
      }  // end of loop k for IBD1 or IBD2
      for(int w = 0; w <= (chrsegMax-chrsegMin); w++)
         pi0[w] *= (0.5 / unrelatedpairsCount);
      for(int type = 0; type < 3; type++){
         if(pairCount[type]==0) continue;
         IBDSegOnly(allpairs[type], seg, ibdsegStorage[0], ibdsegIndex[0],
            ibdsegStorage[1], ibdsegIndex[1], false, 2000000, 10000);
         pi[type].Dimension(ndim);
         pi[type].Zero();
         for(int k = 0; k < 2; k++){
               int localvalue = k+1;
               int tempcount = (ibdsegStorage[k].Length()>>1);
               for(int i = 0; i < tempcount; i++){
                  int startword = ((ibdsegStorage[k][i*2]-1)>>6)+1;
                  int stopword = ((ibdsegStorage[k][i*2+1]+1)>>6)-1;
                  int localcount = stopword-chrsegMin;
                  int w = startword-chrsegMin-1;
                  if(w >= 0 && ibdsegStorage[k][i*2] <= (startword<<6)-33)
                     pi[type][w] += localvalue;
                  for(w++; w <= localcount; w++)
                     pi[type][w] += localvalue;
                  if(w < ndim && (ibdsegStorage[k][i*2+1]&0x3F) >= 32)
                     pi[type][w] += localvalue;
               }
         }  // end of loop k for IBD1 or IBD2
         for(int w = 0; w <= (chrsegMax-chrsegMin); w++)
            pi[type][w] *= (0.5 / pairCount[type]);
      }  // end of loop type for USP, DSP, or ASP
      int chr = chromosomes[chrsegMin<<6];
      double maxLOD = 0.0;
      for(int w = 0; w <= chrsegMax - chrsegMin; w++){
         int base = ((w+chrsegMin)<<6);
         double diff = pi[2][w]-0.5;
         double LOD2 = pairCount[2]*8.0*diff*diff*LODfactor;
         if(pi[2][w] < 0.5)
            LOD2 = -LOD2;
         diff = pi[2][w] - pi0[w] - 0.5;
         double LOD_ASP = pairCount[2] * 8.0 * diff * diff * LODfactor;// Z = sqrt(8N)(pi-pi0-1/2)
         if(diff < 0) LOD_ASP = -LOD_ASP;
         if(pairCount[1] > pairCount[2]*0.1 && pairCount[1]>=20){ // DSP exists
            diff = pi[2][w]-pi[1][w];
            double LOD_DSP = 8.0*pairCount[1]*pairCount[2]/(pairCount[1]+pairCount[2])*
               diff*diff*LODfactor;
            if(pi[2][w] < 0.5) LOD_DSP = 0;
            else if(pi[2][w] < pi[1][w]) LOD_DSP = -LOD_DSP;
            fprintf(fp, "%d\t%.3lf\t%s\t%s\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.2lf\t%.2lf\t%.2lf\n",
               chr, (bp[base+31]+bp[base+32])*0.0000005,
               (const char*)snpName[base+31], (const char*)snpName[base+32],
               pi[0][w], pi[1][w], pi[2][w], pi0[w], LOD2, LOD_ASP, LOD_DSP);
            if(LOD_DSP > 3 && LOD_DSP > maxLOD){
               maxLOD = LOD_DSP;
               sprintf(buffer, "%d\t%.3lf\t%s\t%s\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.2lf\t%.2lf\t%.2lf\n",
               chr, (bp[base+31]+bp[base+32])*0.0000005,
               (const char*)snpName[base+31], (const char*)snpName[base+32],
               pi[0][w], pi[1][w], pi[2][w], pi0[w], LOD2, LOD_ASP, LOD_DSP);
            }
         }else{   // ASP only
            fprintf(fp, "%d\t%.3lf\t%s\t%s\t%.3lf\t%.3lf\t%.2lf\t%.2lf\n",
               chr, (bp[base+31]+bp[base+32])*0.0000005,
               (const char*)snpName[base+31], (const char*)snpName[base+32],
               pi[2][w], pi0[w], LOD2, LOD_ASP);
            if(LOD_ASP > 3 && LOD_ASP > maxLOD){
               maxLOD = LOD_ASP;
               sprintf(buffer, "%d\t%.3lf\t%s\t%s\t%.3lf\t%.3lf\t%.2lf\t%.2lf\n",
               chr, (bp[base+31]+bp[base+32])*0.0000005,
               (const char*)snpName[base+31], (const char*)snpName[base+32],
               pi[2][w], pi0[w], LOD2, LOD_ASP);
            }
         }
      }
      if(maxLOD>3) {
         if(!headerPrinted){
            printf("%s", header);
            headerPrinted = true;
         }
         printf("%s", buffer);
      }
   }
   fclose(fp);
   printf("Genome-wide NPL scan results saved in file %s\n", (const char*)outfile);
}

void Engine::IBDSegOnly(IntArray & pairList, int segment,
   IntArray & ibdsegStorage1, IntArray & ibdsegIndex1,
   IntArray & ibdsegStorage2, IntArray & ibdsegIndex2, bool LengthOnly, int MINSEGLENGTH, int MINCCOUNT)
{
   double ibdprop, maxLength;
   int pairCount = pairList.Length()/2;
   unsigned long long int word;
   int id1, id2, cCount, icCount, segstart, localMin, localMax, minExtraBit, maxExtraBit;
   bool skipFlag;
   int thread = 0;
   IntArray *allsegStorage[2], *allsegIndex[2];
   if(LengthOnly){
      ibdsegStorage1.Dimension(pairCount);
      ibdsegStorage1.Zero();
      ibdsegStorage2.Dimension(pairCount);
      ibdsegStorage2.Zero();
   }else  // detailed segments to be returned
      for(int k = 0; k < 2; k++){
         allsegStorage[k] = new IntArray [defaultMaxCoreCount];
         allsegIndex[k] = new IntArray [defaultMaxCoreCount];
         for(int c = 0; c < defaultMaxCoreCount; c++){
            allsegStorage[k][c].Dimension(0);
            allsegIndex[k][c].Dimension(0);
         }
      }

#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(id1, id2, ibdprop, maxLength, word, segstart, thread, \
      skipFlag, localMin, localMax, minExtraBit, maxExtraBit, cCount, icCount)
{
#endif
   IntArray startPos[2], stopPos[2], startExtraBit[2], stopExtraBit[2];
   IntArray tempStart, tempStop, mergedStart, mergedStop, mergedBit, newchrSeg, newExtraBit;
   IntArray cCounts;
#ifdef _OPENMP
   thread = omp_get_thread_num();
   #pragma omp for
#endif
   for(int pair = 0; pair < pairCount; pair++){
      id1 = pairList[pair<<1];
      id2 = pairList[(pair<<1)|1];
      for(int k = 0; k < 2; k++){
         startPos[k].Dimension(0);
         stopPos[k].Dimension(0);
         startExtraBit[k].Dimension(0);
         stopExtraBit[k].Dimension(0);
      }
//      cCounts.Dimension(0);
      int chrsegMin = localMin = chrSeg[segment<<2];
      int chrsegMax = localMax = chrSeg[(segment<<2)|1];
      minExtraBit = chrSeg[(segment<<2)|2];
      maxExtraBit = chrSeg[(segment<<2)|3];
      for(; localMin <= chrsegMax; localMin++){
         for(; (localMin <= chrsegMax) &&
            (((LG[0][id1][localMin] & LG[0][id2][localMin] &
            (LG[1][id1][localMin] ^ LG[1][id2][localMin])) || //AA x aa
            ((LG[1][id1][localMin] | LG[1][id2][localMin])==0))); // Neither is AA or Aa
            localMin++); // keep passing if localMin includes AA x aa
         if((localMin < chrsegMax) &&
            (LG[0][id1][localMin+1] & LG[0][id2][localMin+1] &
            (LG[1][id1][localMin+1] ^ LG[1][id2][localMin+1]))==0)
            break;   // get out of the loop only when two 0 AA x aa words in a row
      }
      if(localMin >= chrsegMax ||
         bp[(chrsegMax<<6)|0x3F] - bp[localMin<<6] < MINSEGLENGTH)
         continue;// pass if the segment is shorter than 2.5MB
      if(localMin > chrsegMin){
         word = LG[0][id1][localMin-1] & LG[0][id2][localMin-1] & (LG[1][id1][localMin-1] ^ LG[1][id2][localMin-1]);
         if(word == 0){
            localMin --;
            minExtraBit = 0;
         }else
            for(minExtraBit = 0; (word & (1<<(63-minExtraBit))) == 0; minExtraBit++);
      }
      for(; localMax >= localMin; localMax--){
         for(; (LG[0][id1][localMax] & LG[0][id2][localMax] &
            (LG[1][id1][localMax] ^ LG[1][id2][localMax])) ||  //AA x aa
            ((LG[1][id1][localMax] | LG[1][id2][localMax])==0);//None is AA or Aa
            localMax--); // keep passing if localMax includes AA x aa
         if((localMax > localMin) &&
            (LG[0][id1][localMax-1] & LG[0][id2][localMax-1] &
            (LG[1][id1][localMax-1] ^ LG[1][id2][localMax-1]))==0)
            break;   // get out of the loop only when two 0 AA x aa words in a row
      }
      if(localMax < localMin + 2 ||   // 1 or 2 words are not enough for an IBD segment
         bp[(localMax<<6)|0x3F] - bp[localMin<<6] < MINSEGLENGTH)
         continue;// pass if the segment is shorter than 2.5MB
      if(localMax < chrsegMax){
         word = LG[0][id1][localMax+1] & LG[0][id2][localMax+1] & (LG[1][id1][localMax+1] ^ LG[1][id2][localMax+1]);
         if(word == 0){
            localMax ++;
            maxExtraBit = 0;
         }else
            for(maxExtraBit = 0; (word & (1<<maxExtraBit)) == 0; maxExtraBit++);
      }
      // IBD2 segment analysis starts here
      tempStart.Dimension(0); tempStop.Dimension(0);
      for(int m = localMin; m <= localMax; m++){
         for(; (m <= localMax) &&  // keep passing if m word includes IC or does not include C
         ((((LG[0][id1][m] ^ LG[0][id2][m]) | (LG[1][id1][m] ^ LG[1][id2][m])) &
         (LG[0][id1][m] | LG[1][id1][m]) & (LG[0][id2][m] | LG[1][id2][m])) ||
         (LG[1][id1][m] & LG[1][id2][m] & (~(LG[0][id1][m]^LG[0][id2][m])))==0);
         m++); //includes IC (Het x Hom or AA x aa) OR does not include any C (AA x AA or Het x Het)
         for(cCount = 0, segstart = m; m <= localMax &&  // no IC
         ((((LG[0][id1][m] ^ LG[0][id2][m]) | (LG[1][id1][m] ^ LG[1][id2][m])) &
         (LG[0][id1][m] | LG[1][id1][m]) & (LG[0][id2][m] | LG[1][id2][m]))) == 0; m++)
            cCount += popcount(LG[1][id1][m] & LG[1][id2][m] & (~(LG[0][id1][m]^LG[0][id2][m])));
         if( (cCount < 10) || // for sparse array
            (cCount < 20 && bp[((m-1)<<6)|0x3F] - bp[segstart<<6] < 50000) ) // for dense array or WGS
            continue;
         for(word=0; (segstart >= localMin) && (word==0 || (word & (word-1)) == 0); segstart--)
            word = (((LG[0][id1][segstart] ^ LG[0][id2][segstart]) |
               (LG[1][id1][segstart] ^ LG[1][id2][segstart])) &
               (LG[0][id1][segstart] | LG[1][id1][segstart]) &
               (LG[0][id2][segstart] | LG[1][id2][segstart]));// IC
         if(segstart < localMin && (word==0 || (word & (word-1))==0))
            segstart ++;
         else
            segstart += 2; // now 0 or 1 inconsistency
         tempStart.Push(segstart);
//         cCounts.Push(cCount);
         for(word = 0, m--; (m <= localMax) && (word == 0 || (word & (word-1)) == 0); m++)
            word = (((LG[0][id1][m] ^ LG[0][id2][m]) |
            (LG[1][id1][m] ^ LG[1][id2][m])) &
            (LG[0][id1][m] | LG[1][id1][m]) &
            (LG[0][id2][m] | LG[1][id2][m]));// IC
         if(m > localMax && (word == 0 || (word & (word-1)) == 0))
            m--;
         else
            m -= 2; // 0 or 1 inconsistency
         tempStop.Push(m);
      }  // end of scan indexed by m
      skipFlag = false;
      newchrSeg.Dimension(0);
      newExtraBit.Dimension(0);
      int tempcount = tempStart.Length();
      if(tempcount == 0){  // no IBD2 segments
         newchrSeg.Push(localMin);
         newchrSeg.Push(localMax);
         newExtraBit.Push(minExtraBit);
         newExtraBit.Push(maxExtraBit);
         skipFlag = true;
      }
      if(!skipFlag){
         mergedStart.Dimension(0);
         mergedStop.Dimension(0);
         mergedBit.Dimension(0);
         for(int t = 0; t < tempcount-1; t++){
            int gap = bp[tempStart[t+1]<<6]-bp[(tempStop[t]<<6)|0x3F];
            if(tempStart[t+1] - tempStop[t] < 3){
               tempStart[t+1] = tempStart[t];// merge if 1 word in-between
            }else if((gap < 2500000 && tempStart[t+1] - tempStop[t] < 20) || (gap < 0.5)){
               cCount = 0; // consistency C (AA x AA or Het x Het) count
               icCount = -2;   // inconsistency IC (AA x aa or Het x Hom) count
               for(int m = tempStop[t]+1; m < tempStart[t+1]; m++){
                  icCount += popcount((((LG[0][id1][m] ^ LG[0][id2][m]) | (LG[1][id1][m] ^ LG[1][id2][m])) &
                     (LG[0][id1][m] | LG[1][id1][m]) & (LG[0][id2][m] | LG[1][id2][m])));
                  cCount += popcount(LG[1][id1][m] & LG[1][id2][m] & (~(LG[0][id1][m]^LG[0][id2][m])));
               }
               if(cCount > icCount*3) // merge if C_IBD2 >
                  tempStart[t+1] = tempStart[t];
               else if(bp[(tempStop[t]<<6)|0x3F] - bp[tempStart[t]<<6] > MINSEGLENGTH){
                  mergedStart.Push(tempStart[t]);  // IBD2 segments need to be > 2.5MB
                  mergedStop.Push(tempStop[t]);
               } // else discard the left interval
            }else if(bp[(tempStop[t]<<6)|0x3F] - bp[tempStart[t]<<6] > MINSEGLENGTH){
               mergedStart.Push(tempStart[t]);  // IBD2 segments need to be > 2.5MB
               mergedStop.Push(tempStop[t]);    // No gap to consider
            } // else discard the left interval
         }
         if(bp[(tempStop[tempcount-1]<<6)|0x3F] - bp[tempStart[tempcount-1]<<6] > MINSEGLENGTH){
            mergedStart.Push(tempStart[tempcount-1]);  // IBD2 segments need to be > 2.5MB
            mergedStop.Push(tempStop[tempcount-1]);
         }
         tempcount = mergedStart.Length();
         if(tempcount == 0){  // no IBD2 segments
            newchrSeg.Push(localMin);
            newchrSeg.Push(localMax);
            newExtraBit.Push(minExtraBit);
            newExtraBit.Push(maxExtraBit);
            skipFlag = true;
         }else
            for(int t = 0; t < tempcount; t++){
               if(mergedStart[t] == localMin)
                  mergedBit.Push(minExtraBit);
               else{
                  int m = mergedStart[t] - 1;
                  word = ((LG[0][id1][m] ^ LG[0][id2][m]) | (LG[1][id1][m] ^ LG[1][id2][m]))
                     & (LG[0][id1][m] | LG[1][id1][m]) & (LG[0][id2][m] | LG[1][id2][m]);
                  int bit = 0;
                  for(; bit < 63 && (word & (1<<(63-bit))) == 0; bit++);
                  mergedBit.Push(bit);
               }
               if(mergedStop[t] == localMax)
                  mergedBit.Push(maxExtraBit);
               else{
                  int m = mergedStop[t] + 1;
                  word = ((LG[0][id1][m] ^ LG[0][id2][m]) | (LG[1][id1][m] ^ LG[1][id2][m]))
                     & (LG[0][id1][m] | LG[1][id1][m]) & (LG[0][id2][m] | LG[1][id2][m]);
                  int bit = 0;
                  for(; bit < 63 && (word & (1<<bit)) == 0; bit++);
                  mergedBit.Push(bit);
               }
            }
      }
      if(!skipFlag){
         for(int t = 0; t < tempcount; t++){
            startPos[1].Push(mergedStart[t]);
            stopPos[1].Push(mergedStop[t]);
            startExtraBit[1].Push(mergedStart[t]>chrsegMin? mergedBit[t*2]+1: mergedBit[t*2]);
            stopExtraBit[1].Push(mergedStop[t]<chrsegMax? mergedBit[t*2+1]+1: mergedBit[t*2+1]);
         }
         if(mergedStart[0] > localMin){
            newchrSeg.Push(localMin);
            newExtraBit.Push(minExtraBit);
            if(mergedBit[0]){
               newchrSeg.Push(mergedStart[0]-2);
               newExtraBit.Push(64-mergedBit[0]);
            }else{
               newchrSeg.Push(mergedStart[0]-1);
               newExtraBit.Push(0);
            }
         }
         for(int t = 0; t < tempcount-1; t++){
            if(mergedBit[t*2+1]){
               newchrSeg.Push(mergedStop[t]+2);
               newExtraBit.Push(64 - mergedBit[t*2+1]);
            }else{
               newchrSeg.Push(mergedStop[t]+1);
               newExtraBit.Push(0);
            }
            if(mergedBit[(t+1)*2]){
               newchrSeg.Push(mergedStart[t+1]-2);
               newExtraBit.Push(64 - mergedBit[(t+1)*2]);
            }else{
               newchrSeg.Push(mergedStart[t+1]-1);
               newExtraBit.Push(0);
            }
         }
         if(mergedStop[tempcount-1] < localMax){
            if(mergedBit[tempcount*2-1]){
               newchrSeg.Push(mergedStop[tempcount-1]+2);
               newExtraBit.Push(64 - mergedBit[tempcount*2-1]);
            }else{
               newchrSeg.Push(mergedStop[tempcount-1]+1);
               newExtraBit.Push(0);
            }
            newchrSeg.Push(localMax);
            newExtraBit.Push(maxExtraBit);
         }
      }

      // IBD1 segment analysis starts here
      int newcount = newchrSeg.Length()/2;
      for(int s = 0; s < newcount; s++){
         int chrsegMin = newchrSeg[s<<1];
         int chrsegMax = newchrSeg[(s<<1)|1];
         tempStart.Dimension(0);
         tempStop.Dimension(0);
         cCounts.Dimension(0);
         for(int m = chrsegMin; m <= chrsegMax; m++){
            for(; (m <= chrsegMax) &&  // keep passing if m word
            ((LG[0][id1][m] & LG[0][id2][m] & (LG[1][id1][m] ^ LG[1][id2][m])) ||
            ((LG[1][id1][m] & LG[1][id2][m] & (LG[0][id1][m] | LG[0][id2][m]))==0));
            m++); //includes IC (AA x aa) OR does not include any C (AA x AA or AA x Aa)
            for(cCount = 0, segstart = m; m <= chrsegMax &&
               (LG[0][id1][m] & LG[0][id2][m] & (LG[1][id1][m] ^ LG[1][id2][m]))==0; m++)
               cCount += popcount(LG[1][id1][m] & LG[1][id2][m] & (LG[0][id1][m] | LG[0][id2][m]));
            if( (cCount < 10) || // for sparse array
               (cCount < 20 && bp[((m-1)<<6)|0x3F] - bp[segstart<<6] < 20000) ) // for dense array or WGS
               continue;
            for(segstart--; segstart >= chrsegMin &&  // icCount==0
               (LG[0][id1][segstart] & LG[0][id2][segstart] & (LG[1][id1][segstart] ^ LG[1][id2][segstart]))==0;
               segstart--); // keep passing if free of IC (AA x aa)
            tempStart.Push(segstart+1);
            tempStop.Push(m-1);
            cCounts.Push(cCount);
         }  // end of scan indexed by m
         tempcount = tempStart.Length();
         if(tempcount == 0) continue;
         mergedStart.Dimension(0);
         mergedStop.Dimension(0);
         for(int t = 0; t < tempcount-1; t++){
            int gap = bp[tempStart[t+1]<<6]-bp[(tempStop[t]<<6)|0x3F];
            if(tempStart[t+1] - tempStop[t] < 3){
               tempStart[t+1] = tempStart[t];// merge if 1 word in-between
            }else if((gap < 2500000 && tempStart[t+1] - tempStop[t] < 20) || (gap < 500000)){
               cCount = 0; // consistency C (AA x AA or AA x Aa) count
               icCount = -2;   // inconsistency IC (AA x aa) count
               for(int m = tempStop[t]+1; m < tempStart[t+1]; m++){
                  icCount += popcount(LG[0][id1][m] & LG[0][id2][m] & (LG[1][id1][m] ^ LG[1][id2][m]));
                  cCount += popcount(LG[1][id1][m] & LG[1][id2][m] & (LG[0][id1][m] | LG[0][id2][m]));
               }
               if(cCount > icCount*6){ // merge if C_IBD1 > 85.7%
                  tempStart[t+1] = tempStart[t];
                  cCounts[t+1] += cCounts[t] + cCount;
               }else if((bp[(tempStop[t]<<6)|0x3F] - bp[tempStart[t]<<6] > MINSEGLENGTH)||
                  (cCounts[t] >= MINCCOUNT) ){
                  mergedStart.Push(tempStart[t]);  // IBD1 segments need to be > 2.5MB
                  mergedStop.Push(tempStop[t]);
               } // else discard the left interval
            }else if((bp[(tempStop[t]<<6)|0x3F] - bp[tempStart[t]<<6] > MINSEGLENGTH)||
               (cCounts[t] >= MINCCOUNT) ){
               mergedStart.Push(tempStart[t]);  // IBD1 segments need to be > 2.5MB
               mergedStop.Push(tempStop[t]);    // No gap to consider
            } // else discard the left interval
         }
         if((bp[(tempStop[tempcount-1]<<6)|0x3F] - bp[tempStart[tempcount-1]<<6] > MINSEGLENGTH) ||
            (cCounts[tempcount-1] >= MINCCOUNT) ){
            mergedStart.Push(tempStart[tempcount-1]);  // IBD1 segments need to be > 2.5MB
            mergedStop.Push(tempStop[tempcount-1]);
         }
         tempcount = mergedStart.Length();
         for(int t = 0; t < tempcount; t++){
            startPos[0].Push(mergedStart[t]);
            stopPos[0].Push(mergedStop[t]);
         }
         for(int t = 0; t < tempcount; t++){
            int m = mergedStart[t];
            if(m > chrsegMin){
               m --; // AA x aa
               word = LG[0][id1][m] & LG[0][id2][m] & (LG[1][id1][m] ^ LG[1][id2][m]);
               int bit = 0;
               for(; bit < 63 && (word & (1<<(63-bit))) == 0; bit++);
               startExtraBit[0].Push(bit);
            }else
               startExtraBit[0].Push(newExtraBit[s*2]);
            m = mergedStop[t];
            if(m < chrsegMax){
               m ++; // AA x aa
               word = LG[0][id1][m] & LG[0][id2][m] & (LG[1][id1][m] ^ LG[1][id2][m]);
               int bit = 0;
               for(; bit < 63 && (word & (1<<bit)) == 0; bit++);
               stopExtraBit[0].Push(bit);
            }else
               stopExtraBit[0].Push(newExtraBit[s*2+1]);
         }// end of t
      }// end of s

      if(startPos[0].Length() || startPos[1].Length())
         for(int k = 0; k < 2; k++){
            int newsegCount = startPos[k].Length();
            if(LengthOnly){
               double length = 0.0;
               for(int s = 0; s < newsegCount; s++){
                  localMin = (startPos[k][s]<<6)-startExtraBit[k][s];
                  if(stopPos[k][s] == longCount-1)
                     localMax = markerCount - 1;
                  else
                     localMax = ((stopPos[k][s]<<6)|0x3F)+stopExtraBit[k][s];
                  length += bp[localMax] - bp[localMin];
               }
               if(k==0)
                  ibdsegStorage1[pair] = int(length+0.5);
               else
                  ibdsegStorage2[pair] = int(length+0.5);
            }else // Detailed segments to be returned
               for(int s = 0; s < newsegCount; s++){
                  localMin = (startPos[k][s]<<6)-startExtraBit[k][s];
                  if(stopPos[k][s] == longCount-1)
                     localMax = markerCount - 1;
                  else
                     localMax = ((stopPos[k][s]<<6)|0x3F)+stopExtraBit[k][s];
                  allsegStorage[k][thread].Push(localMin);
                  allsegStorage[k][thread].Push(localMax);
                  allsegIndex[k][thread].Push(pair);
               }  // end of s
         }  // end of k
   }  // end of all pairs x segs
#ifdef _OPENMP
}
#endif
   if(!LengthOnly){
      ibdsegStorage1 = allsegStorage[0][0];
      ibdsegStorage2 = allsegStorage[1][0];
      if(ibdsegIndex1) ibdsegIndex1 = allsegIndex[0][0];
      if(ibdsegIndex2) ibdsegIndex2 = allsegIndex[1][0];
#ifdef _OPENMP
      for(thread = 1; thread < defaultMaxCoreCount; thread++){
         ibdsegStorage1.Append(allsegStorage[0][thread]);
         ibdsegStorage2.Append(allsegStorage[1][thread]);
         if(ibdsegIndex1)
            ibdsegIndex1.Append(allsegIndex[0][thread]);
         if(ibdsegIndex2)
            ibdsegIndex2.Append(allsegIndex[1][thread]);
      }
#endif
      for(int k = 0; k < 2; k++){
         delete []allsegStorage[k];
         delete []allsegIndex[k];
      }
   }
}




/*
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
   IntArray pi[2];
   printf("Scanning genome...\n");
   String outfile(prefix);
   outfile.Add(".homomap");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "Chr\tPos\tFlankSNP1\tFlankSNP2\tN_Con\tHomoCon\tN_Cas\tHomoCas\tOR\tSE\tLOD\n");
   int segCount = (chrSeg.Length()>>2);
   for(int seg = 0; seg < segCount; seg++){
      int chrsegMin = chrSeg[seg<<2];
      int chrsegMax = chrSeg[(seg<<2)|1];
      int ndim = chrsegMax - chrsegMin + 1;
      int chr = chromosomes[chrsegMin<<6];
      for(int pop = 0; pop < 2; pop++){
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
            double length = positions[stopPos] - positions[startPos];
            int id = rohIndex[i];
            FROH[pop][id] += length;
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
      for(int w = 0; w < ndim; w++){
         int base = ((w+chrsegMin)<<6);
         if(!pi[0][w] || !pi[1][w]) continue;
         fprintf(fp, "%d\t%.3lf\t%s\t%s",
            chr, (positions[base+31]+positions[base+32])/2,
            (const char*)snpName[base+31], (const char*)snpName[base+32]);
         int a = affCount[0] - pi[0][w];
         int b = affCount[1] - pi[1][w];
         if(!a || !b) continue;
         double OR = pi[1][w] * a * 1.0 / (pi[0][w] * b);
         double se = sqrt(1.0/pi[0][w] + 1.0/pi[1][w] + 1.0/a + 1.0/b);
         double LOD = log(OR) / se;
         LOD = LOD*LOD*LODfactor;
         if(OR < 1) LOD = -LOD;
         fprintf(fp, "\t%d\t%d\t%d\t%d\t%.3lf\t%.3lf\t%.2lf\n",
            affCount[0], pi[0][w], affCount[1], pi[1][w], OR, se, LOD);
      }  // end of word w loop
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
   printf("Homozygosity mapping scan results saved in file %s\n", (const char*)outfile);
}
*/




