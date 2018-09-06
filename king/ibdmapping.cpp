//////////////////////////////////////////////////////////////////////
// ibdmapping.cpp
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
// August 22, 2018

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
#ifdef __ZLIB_AVAILABLE__
  #include <zlib.h>
#endif

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
   IntArray ibdsegCount, ibdsegStorage[2], ibdsegIndex[2];
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
         IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], false, 2.5, mincons);
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
            chromosomes[segstart], positions[segstart], positions[segstop],
            (const char*)snpName[segstart], (const char*)snpName[segstop],
            segstop-segstart+1, positions[segstop]-positions[segstart]);
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
            chromosomes[segstart], positions[segstart], positions[segstop],
            (const char*)snpName[segstart], (const char*)snpName[segstop],
            segstop-segstart+1, positions[segstop]-positions[segstart]);
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
         IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], true, 2.5, mincons);
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

void Engine::IBDVC()
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
   if(N_Trait>1)
      printf("Individuals with less than %d phenotypes are removed.\n", N_Trait);
   if(normalization)printf("Inverse normal transformation is applied to phenotypes.\n");

   int N_ID = ID.Length();
   printf("Calculating similarity matrix starts at %s", currentTime());
   Matrix D(N_ID, N_ID);
   D.Zero();
   for(int i = 0; i < N_ID; i++) D[i][i] = 1.0;
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
      if(mincons)
         IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], true, 2.5, mincons);
      else
         IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], true);
      int pairCount = pairIndex[0].Length();
      for(int p = 0; p < pairCount; p++)
         D[pairIndex[0][p]][pairIndex[1][p]] += (ibdsegStorage[0][p]*0.5 + ibdsegStorage[1][p])*0.000001;
   }  // end of seg loop
   for(int i = 0; i < N_ID; i++)
      for(int j = i+1; j < N_ID; j++){
         D[i][j] /= totalLength;
         D[j][i] = D[i][j];
      }
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
   }

   double variances[2];
   Vector coef(N_Covariate+1);
   double IVY, PVYY, IVI, PVP;
   Vector XVY(N_Covariate+1);
   Matrix XVX(N_Covariate+1, N_Covariate+1);
   printf("\nPolygenic parameter estimates\n");
   printf("%-20s %7s %7s %9s", "TraitName", "Herit", "N", "Mu");
   for(int i = 0; i < N_Covariate; i++)
      printf(" %9s", (const char*)ped.covariateNames[covariates[i]]);
   printf("\n");
   for(int t = 0; t < N_Trait; t++){
      coef.Zero();
      IVY = PVYY = IVI = PVP = 0;
      XVY.Zero();
      XVX.Zero();
      for(int i = 0; i < N_ID; i++){
         for(int k = 0; k < N_Covariate+1; k++){
            XVY[k] += X[i][k] * Y[i][t];
            for(int l = 0; l < N_Covariate+1; l++)
               XVX[k][l] += X[i][k] * X[i][l];
         }
      }
      if(N_Covariate == 1) coef[0] = XVY[0] / XVX[0][0];
      else{
         Matrix tempM;
         tempM.CholeskyInvert(XVX);
         for(int i = 0; i < N_Covariate+1; i++){
            coef[i] = 0;
            for(int u = 0; u < N_Covariate+1; u++)
               coef[i] += tempM[i][u] * XVY[u];
         }
      }
      for(int i = 0; i < N_ID; i++)
         for(int k = 0; k < N_Covariate+1; k++)
            Y[i][t] -= coef[k] * X[i][k];
      for(int i = 0; i < N_ID; i++){
         IVY += Y[i][t] * Y[i][t];
         IVI ++;
      }
      for(int u = 0; u < N_ID; u++)
         for(int v = u+1; v < N_ID; v++) {
            PVYY += D[u][v] * Y[u][t] * Y[v][t];
            PVP += D[u][v] * D[u][v];
         }
      variances[0] = (IVY - IVI * PVYY / PVP) / IVI;
      if(variances[0]<0) variances[0] =  0.000001;
      variances[1] = PVYY / PVP;
      if(variances[1]<0) variances[1] =  0.000001;
      double h2 = variances[1] / (variances[0]+variances[1]);
      printf("%-20s %7.5lf %7d", (const char*)ped.traitNames[traits[t]], h2, N_ID);
      for(int i = 0; i < N_Covariate+1; i++)
         printf(" %9.5lf", coef[i]);
      printf("\n");
   }
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
   IntArray ibdsegCount, ibdsegStorage[3][2], ibdsegIndex[3][2], pi[3];
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
            IBDSegOnly(allpairs[type], seg, ibdsegStorage[type][0], ibdsegIndex[type][0], ibdsegStorage[type][1], ibdsegIndex[type][1], false, 2.5, mincons);
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
            chr, (positions[base+31]+positions[base+32])/2,
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
         IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], false, 2.5, mincons);
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
               chr, (positions[base+31]+positions[base+32])/2,
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
   int pairCount = allpairs.Length()>>1;
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
         IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], false, 2.5, mincons);
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
                  (const char*)snpName[base+k], chr, positions[base+k],
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
            IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], false, 2.5, mincons);
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
            chr, (positions[base+31]+positions[base+32])/2,
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
         IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], false, 2.5, mincons);
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
                  (const char*)snpName[pos], chr, positions[pos],
                  TCount[k], NTCount[k], tdt);
            }
         }  // end of word loop
      delete T;
      delete NT;
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
               IBDSegOnly(allpairs, seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], true, 2.5, mincons);
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
               IBDSegOnly(allpairs, seg, ibdsegStorage[a][0][person], ibdsegIndex[a][0][person], ibdsegStorage[a][1][person], ibdsegIndex[a][1][person], false, 2.5, mincons);
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
      if(positions[chrsegMin<<6] >= allpos[chrIndex]*0.000001
         || positions[(chrsegMax<<6)|0x3F] <= allpos[chrIndex]*0.000001)
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
            IBDSegOnly(allpairs[type], seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], false, 2.5, mincons);
         else
            IBDSegOnly(allpairs[type], seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1]);
         for(int k = 0; k < 2; k++){
            int tempcount = (ibdsegStorage[k].Length()>>1);
            for(int i = 0; i < tempcount; i++){
               if(positions[ibdsegStorage[k][i*2]] >= allpos[chrIndex]*0.000001
                  || positions[ibdsegStorage[k][i*2+1]] <= allpos[chrIndex]*0.000001)
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
   AUC = ComputeAUC(risks, affections, false, 1000);
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
         AUC = ComputeAUC(risks, affections, false, 1000);
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
      AUC = ComputeAUC(risks, affections, false, 1000);
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
      AUC = ComputeAUC(risks, affections, false, 1000);
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
      AUC = ComputeAUC(risks, affections, false, 1000);
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
      AUC = ComputeAUC(risks, affections, false, 1000);
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
   IntArray affections(idCount);
   for(int i = 0; i < idCount; i++)
      affections[i] = ped[phenoid[i]].affections[0]-1;
   int segCount = (chrSeg.Length()>>2);
   printf("Scanning genome...\n");
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
         int afftype[2];
         if(type == 0){
            afftype[0] = afftype[1] = 0;
         }else if(type == 1){
            afftype[0] = 0; afftype[1] = 1;
         }else{
            afftype[0] = afftype[1] = 1;
         }
         if(mincons)
            IBDSegOnly(allpairs[type], seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1], false, 2.5, mincons);
         else
            IBDSegOnly(allpairs[type], seg, ibdsegStorage[0], ibdsegIndex[0], ibdsegStorage[1], ibdsegIndex[1]);
         pi[type].Dimension(ndim);
         pi[type].Zero();
         for(int k = 0; k < 2; k++){
            int localvalue = k+1;
            int tempcount = (ibdsegStorage[k].Length()>>1);
            for(int i = 0; i < tempcount; i++){
               int index = ibdsegIndex[k][i];
               int id1 = allpairs[type][index*2];
               int id2 = allpairs[type][index*2+1];
               int startword = ((ibdsegStorage[k][i*2]-1)>>6)+1;
               int stopword = ((ibdsegStorage[k][i*2+1]+1)>>6)-1;
               int localcount = stopword-chrsegMin;
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
               }
            }  // end of segment i loop
         }  // end of IBD1 or IBD2 k loop
      }  // end of pair type loop
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
         AUC[w] = ComputeAUC(risks, affections, false, 100);
         if(AUC[w] > 0.6)
            AUC[w] = ComputeAUC(risks, affections, false, 1000);
      }
#ifdef _OPENMP
}
#endif
      double maxAUC = 0.5999;
      char buffer[256];
      for(int w = 0; w < ndim; w++){
         for(int type = 0; type < 3; type++)
            localPi[type] = pi[type][w] * 0.25 / pairCount[type];
         double success = 1 - NMISS[w]*1.0/(affcount[0]+affcount[1]);
         int base = ((w+chrsegMin)<<6);
         fprintf(fp, "%d\t%.3lf\t%s\t%s\t%.2G\t%.2G\t%.2G\t%d\t%d\t%d\t%d\t%.3lf\t%.3lf\n",
            chr, (positions[base+31]+positions[base+32])/2,
            (const char*)snpName[base+31], (const char*)snpName[base+32],
            localPi[0], localPi[1], localPi[2], NMISS[w], NUU[w], NUA[w], NAA[w], success, AUC[w]);
         if(AUC[w] > maxAUC && success > 0.9){
            maxAUC = AUC[w];
            sprintf(buffer, "%d\t%.3lf\t%s\t%s\t%.2G\t%.2G\t%.2G\t%d\t%d\t%d\t%d\t%.3lf\t%.3lf\n",
            chr, (positions[base+31]+positions[base+32])/2,
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

void Engine::NPL()
{
   printf("\nOptions in effect:\n");
   printf("\t--npl\n");
   if(mincons)
      printf("\t--mincons %d\n", mincons);
   if(relativedegree)
      printf("\t--degree %d\n", relativedegree);
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

   IntArray allpairs[3], sibs[2];
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
            for(int s = 0; s < sibcount[0]; s++)
               for(int t = 0; t < sibcount[1]; t++){
                  allpairs[1].Push(sibs[0][s]);
                  allpairs[1].Push(sibs[1][t]);
               }
         }
   }
   for(int k = 0; k < 3; k++)
      pairCount[k] = allpairs[k].Length()>>1;
   if(pairCount[2] == 0) {printf("No ASPs are found.\n"); return;}
   else printf("#%d ASPs (and #%d DSPs) are used for NPL scan.\n",
      pairCount[2], pairCount[1]);
   char header[256], buffer[256];
   if(pairCount[1])
      sprintf(header, "Chr\tPos\tFlankSNP1\tFlankSNP2\tPI_USP\tPI_DSP\tPI_ASP\tLOD10\tLOD3\tLOD_ASP\tLOD_DSP\n");
   else
      sprintf(header, "Chr\tPos\tFlankSNP1\tFlankSNP2\tPI_ASP\tLOD10\tLOD3\tLOD_ASP\n");
   String outfile(prefix);
   outfile.Add(".npl");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "%s", header);
   bool headerPrinted=false;
   IntArray ibdsegCount, ibdsegStorage[2], ibdsegIndex[2];
   Vector pi[2][3];
   int segCount = (chrSeg.Length()>>2);
   const double LODfactor = 0.5 / log(10.0);
   for(int seg = 0; seg < segCount; seg++){
      int chrsegMin = chrSeg[seg<<2];
      int chrsegMax = chrSeg[(seg<<2)|1];
      int ndim = chrsegMax - chrsegMin + 1;
      for(int longseg = 0; longseg < 2; longseg++)
         for(int type = 0; type < 3; type++){
            if(pairCount[type]==0) continue;
            if(mincons)
               IBDSegOnly(allpairs[type], seg, ibdsegStorage[0], ibdsegIndex[0],
               ibdsegStorage[1], ibdsegIndex[1], false, longseg==1?10:3, mincons);
            else
               IBDSegOnly(allpairs[type], seg, ibdsegStorage[0], ibdsegIndex[0],
               ibdsegStorage[1], ibdsegIndex[1], false, longseg==1?10:3, 10000);
            pi[longseg][type].Dimension(ndim);
            pi[longseg][type].Zero();
            for(int k = 0; k < 2; k++){
               int localvalue = k+1;
               int tempcount = (ibdsegStorage[k].Length()>>1);
               for(int i = 0; i < tempcount; i++){
                  int startword = ((ibdsegStorage[k][i*2]-1)>>6)+1;
                  int stopword = ((ibdsegStorage[k][i*2+1]+1)>>6)-1;
                  int localcount = stopword-chrsegMin;
                  int w = startword-chrsegMin-1;
                  if(w >= 0 && ibdsegStorage[k][i*2] <= (startword<<6)-33)
                     pi[longseg][type][w] += localvalue;
                  for(w++; w <= localcount; w++)
                     pi[longseg][type][w] += localvalue;
                  if(w < ndim && (ibdsegStorage[k][i*2+1]&0x3F) >= 32)
                     pi[longseg][type][w] += localvalue;
               }
            }  // end of loop k for IBD1 or IBD2
            for(int w = 0; w <= (chrsegMax-chrsegMin); w++)
               pi[longseg][type][w] *= (0.5 / pairCount[type]);
         }  // end of loop type for USP, DSP, or ASP
      int chr = chromosomes[chrsegMin<<6];
      double maxLOD = 0.0;
      for(int w = 0; w <= chrsegMax - chrsegMin; w++){
         int base = ((w+chrsegMin)<<6);
         double LOD3 = pairCount[2]*8.0*(pi[0][2][w]-0.5)*(pi[0][2][w]-0.5)*LODfactor;
         double LOD10 = pairCount[2]*8.0*(pi[1][2][w]-0.5)*(pi[1][2][w]-0.5)*LODfactor;
         if(pi[0][2][w] < 0.5)
            LOD3 = -LOD3;
         if(pi[1][2][w] < 0.5)
            LOD10 = -LOD10;
         double LOD_ASP = (LOD3 > LOD10+1)? LOD10: LOD3;   // Z = sqrt(8N)(pi-1/2)
         if(pairCount[1]){ // DSP exists
            double LOD_DSP = 8.0*pairCount[1]*pairCount[2]/(pairCount[1]+pairCount[2])*
               (pi[0][2][w]-pi[0][1][w])*(pi[0][2][w]-pi[0][1][w])*LODfactor;
            if(pi[0][2][w] < pi[0][1][w]) LOD_DSP = -LOD_DSP;
            fprintf(fp, "%d\t%.3lf\t%s\t%s\t%.3lf\t%.3lf\t%.3lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n",
               chr, (positions[base+31]+positions[base+32])/2,
               (const char*)snpName[base+31], (const char*)snpName[base+32],
               pi[0][0][w], pi[0][1][w], pi[0][2][w], LOD10, LOD3, LOD_ASP, LOD_DSP);
            if(LOD_DSP > 3 && LOD_DSP > maxLOD){
               maxLOD = LOD_DSP;
               sprintf(buffer, "%d\t%.3lf\t%s\t%s\t%.3lf\t%.3lf\t%.3lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n",
               chr, (positions[base+31]+positions[base+32])/2,
               (const char*)snpName[base+31], (const char*)snpName[base+32],
               pi[0][0][w], pi[0][1][w], pi[0][2][w], LOD10, LOD3, LOD_ASP, LOD_DSP);
            }
         }else{   // ASP only
            fprintf(fp, "%d\t%.3lf\t%s\t%s\t%.3lf\t%.2lf\t%.2lf\t%.2lf\n",
               chr, (positions[base+31]+positions[base+32])/2,
               (const char*)snpName[base+31], (const char*)snpName[base+32],
               pi[0][2][w], LOD10, LOD3, LOD_ASP);
            if(LOD_ASP > 3 && LOD_ASP > maxLOD){
               maxLOD = LOD_ASP;
               sprintf(buffer, "%d\t%.3lf\t%s\t%s\t%.3lf\t%.2lf\t%.2lf\t%.2lf\n",
               chr, (positions[base+31]+positions[base+32])/2,
               (const char*)snpName[base+31], (const char*)snpName[base+32],
               pi[0][2][w], LOD10, LOD3, LOD_ASP);
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
   IntArray & ibdsegStorage2, IntArray & ibdsegIndex2, bool LengthOnly, double MINSEGLENGTH, int MINCCOUNT)
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
         positions[(chrsegMax<<6)|0x3F] - positions[localMin<<6] < MINSEGLENGTH)
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
         positions[(localMax<<6)|0x3F] - positions[localMin<<6] < MINSEGLENGTH)
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
            (cCount < 20 && positions[((m-1)<<6)|0x3F] - positions[segstart<<6] < 0.05) ) // for dense array or WGS
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
            double gap = positions[tempStart[t+1]<<6]-positions[(tempStop[t]<<6)|0x3F];
            if(tempStart[t+1] - tempStop[t] < 3){
               tempStart[t+1] = tempStart[t];// merge if 1 word in-between
            }else if((gap < 2.5 && tempStart[t+1] - tempStop[t] < 20) || (gap < 0.5)){
               cCount = 0; // consistency C (AA x AA or Het x Het) count
               icCount = -2;   // inconsistency IC (AA x aa or Het x Hom) count
               for(int m = tempStop[t]+1; m < tempStart[t+1]; m++){
                  icCount += popcount((((LG[0][id1][m] ^ LG[0][id2][m]) | (LG[1][id1][m] ^ LG[1][id2][m])) &
                     (LG[0][id1][m] | LG[1][id1][m]) & (LG[0][id2][m] | LG[1][id2][m])));
                  cCount += popcount(LG[1][id1][m] & LG[1][id2][m] & (~(LG[0][id1][m]^LG[0][id2][m])));
               }
               if(cCount > icCount*3) // merge if C_IBD2 >
                  tempStart[t+1] = tempStart[t];
               else if(positions[(tempStop[t]<<6)|0x3F] - positions[tempStart[t]<<6] > MINSEGLENGTH){
                  mergedStart.Push(tempStart[t]);  // IBD2 segments need to be > 2.5MB
                  mergedStop.Push(tempStop[t]);
               } // else discard the left interval
            }else if(positions[(tempStop[t]<<6)|0x3F] - positions[tempStart[t]<<6] > MINSEGLENGTH){
               mergedStart.Push(tempStart[t]);  // IBD2 segments need to be > 2.5MB
               mergedStop.Push(tempStop[t]);    // No gap to consider
            } // else discard the left interval
         }
         if(positions[(tempStop[tempcount-1]<<6)|0x3F] - positions[tempStart[tempcount-1]<<6] > MINSEGLENGTH){
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
               (cCount < 20 && positions[((m-1)<<6)|0x3F] - positions[segstart<<6] < 0.02) ) // for dense array or WGS
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
            double gap = positions[tempStart[t+1]<<6]-positions[(tempStop[t]<<6)|0x3F];
            if(tempStart[t+1] - tempStop[t] < 3){
               tempStart[t+1] = tempStart[t];// merge if 1 word in-between
            }else if((gap < 2.5 && tempStart[t+1] - tempStop[t] < 20) || (gap < 0.5)){
               cCount = 0; // consistency C (AA x AA or AA x Aa) count
               icCount = -2;   // inconsistency IC (AA x aa) count
               for(int m = tempStop[t]+1; m < tempStart[t+1]; m++){
                  icCount += popcount(LG[0][id1][m] & LG[0][id2][m] & (LG[1][id1][m] ^ LG[1][id2][m]));
                  cCount += popcount(LG[1][id1][m] & LG[1][id2][m] & (LG[0][id1][m] | LG[0][id2][m]));
               }
               if(cCount > icCount*6){ // merge if C_IBD1 > 85.7%
                  tempStart[t+1] = tempStart[t];
                  cCounts[t+1] += cCounts[t] + cCount;
               }else if((positions[(tempStop[t]<<6)|0x3F] - positions[tempStart[t]<<6] > MINSEGLENGTH)||
                  (cCounts[t] >= MINCCOUNT) ){
                  mergedStart.Push(tempStart[t]);  // IBD1 segments need to be > 2.5MB
                  mergedStop.Push(tempStop[t]);
               } // else discard the left interval
            }else if((positions[(tempStop[t]<<6)|0x3F] - positions[tempStart[t]<<6] > MINSEGLENGTH)||
               (cCounts[t] >= MINCCOUNT) ){
               mergedStart.Push(tempStart[t]);  // IBD1 segments need to be > 2.5MB
               mergedStop.Push(tempStop[t]);    // No gap to consider
            } // else discard the left interval
         }
         if((positions[(tempStop[tempcount-1]<<6)|0x3F] - positions[tempStart[tempcount-1]<<6] > MINSEGLENGTH) ||
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
                  length += positions[localMax] - positions[localMin];
               }
               if(k==0)
                  ibdsegStorage1[pair] = int(length*1000000+0.5);
               else
                  ibdsegStorage2[pair] = int(length*1000000+0.5);
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




