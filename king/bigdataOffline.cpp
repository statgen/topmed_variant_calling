//////////////////////////////////////////////////////////////////////
// bigdataOffline.cpp
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
// March 20, 2019

#include <math.h>
#include "analysis.h"
#include "KingCore.h"
#include "Kinship.h"
#include "KinshipX.h"
#include "MathStats.h"
#include "MathSVD.h"
#include "QuickIndex.h"
#include "MathCholesky.h"
#include "Intervals.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

void Engine::ComputeBigDataSecondDegree()
{
   if(Bit64 != 64) {printf("Cannot run --exact analysis without 64-bit system.\n"); return;}
   if(longCount==0) {printf("No genotype data.\n"); return;}
   printf("Autosome genotypes stored in %d words for each of %d individuals.\n", longCount, idCount);
   printf("\nOptions in effect:\n");
   printf("\t--exact\n");
   if(CoreCount)
      printf("\t--cpus %d\n", CoreCount);
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");
   bool IBDvalidFlag = PreSegment();
   if(!IBDvalidFlag) {
      printf("%s\n", (const char*)segmessage);
      printf("Cannot run --exact analysis without IBD segments.\n");
      return;
   }
   int notMissingCount, afterCount[6], i1, i2, degree, beforeCount[6];
   double ibs0, kinship, smaller, CHet, ibdprop, ibd2prop;
   IntArray pairs, pairI, HetHetCounts, IBS0Counts, het1Counts, het2Counts, HomHomCounts, IBSCounts;
   Vector ibdprops, maxLengths, ibd2props, maxLengths2, pis;
   Matrix Phi;
   String type, outfile;
   outfile.Copy(prefix);
   outfile.Add(".kin");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "FID\tID1\tID2\tN_SNP\tZ0\tPhi\tHetHet\tIBS0\tHetConc\tHomIBS0\tKinship\tIBD1Seg\tIBD2Seg\tPropIBD\tInfType\tExact\n");
   for(int i = 0; i < 6; i++) beforeCount[i] = afterCount[i] = 0;
   Kinship kin;
   IntArray *pair2Index;
   IntArray tempArray, tempArray2, ibd1seg;
   String exactType;
   IntArray *allrelatives, *PR, *ibd1segs;
   for(int f = 0; f < ped.familyCount; f++){
      int idfCount = id[f].Length();
      if(idfCount < 2) continue;   // no pairs in family f
      kin.Setup(*ped.families[f]);
      pairs.Dimension(0);
      pairI.Dimension(0);
      Phi.Dimension(idfCount, idfCount);
      for(int i = 0; i < idfCount; i++)
         for(int j = i+1; j < idfCount; j++){
            pairs.Push(geno[id[f][i]]);
            pairs.Push(geno[id[f][j]]);
            pairI.Push(i);
            pairI.Push(j);
            double phi = kin(ped[id[f][i]], ped[id[f][j]]);
            Phi[i][j] = phi;
            double pi0 = 0.0;
            if(phi < 0.2)
               pi0 = 1-4*phi;
            else if(phi < 0.3 && ped[id[f][i]].isSib(ped[id[f][j]]))
               pi0 = 0.25;
            Phi[j][i] = pi0;
         }
      int pairCount = pairs.Length()/2;
      KinshipInSubset64Bit(pairs, HetHetCounts, IBS0Counts, het1Counts, het2Counts, HomHomCounts, IBSCounts);
      IntArray *ibd1segs = new IntArray [pairCount];
      IBDSegInSubset64Bit(pairs, ibdprops, maxLengths, ibd2props, maxLengths2, ibd1segs);

      allrelatives = new IntArray[idfCount];
      PR = new IntArray[idfCount];
      for(int i = 0; i < idfCount; i++) {
         allrelatives[i].Dimension(0);
         PR[i].Dimension(0);
      }
      pis.Dimension(pairCount);
      pair2Index = new IntArray[idfCount];
      for(int i = 0; i < idfCount; i++)
         pair2Index[i].Dimension(idfCount);
      for(int p = 0; p < pairCount; p++){
         if(!het1Counts[p] && !het2Counts[p] && !HetHetCounts[p]) continue;
         int i1 = pairI[p*2]; int i2 = pairI[p*2+1];
         double CHet = HetHetCounts[p] * 1.0 / (HetHetCounts[p]+het1Counts[p]+het2Counts[p]);
         double ibdprop = ibdprops[p];
         double ibd2prop = ibd2props[p];
         double pi = ibd2prop + ibdprop * 0.5;
         pair2Index[i1][i2] = pair2Index[i2][i1] = p;
         pis[p] = pi;
         if(CHet<0.8){
            bool isFS = false;
            if(pi > 0.3535534 && (ibdprop + ibd2prop < 0.96) &&
               (ibdprop + ibd2prop < 0.9 || ibd2prop > 0.08) && (ibd2prop > 0.15 || pi > 0.4) )
                  isFS = true;
            if(pi > 0.04419417 && !isFS){// closer than 4th-degree
               allrelatives[i1].Push(i2);
               allrelatives[i2].Push(i1);
            }
         }
      }
      for(int i = 0; i < idfCount; i++){
         int allrelativeCount = allrelatives[i].Length();
         for(int j = 0; j < allrelativeCount; j++){
            int R1 = allrelatives[i][j];
            for(int k = j+1; k < allrelativeCount; k++){
               int R2 = allrelatives[i][k];
               if(RoRP(ibd1segs[pair2Index[i][R1]], ibd1segs[pair2Index[i][R2]],
                  ibd1segs[pair2Index[R1][R2]], bp) < 0.125){
                  if(PR[i].Find(R1)==-1) PR[i].Push(R1);
                  if(PR[i].Find(R2)==-1) PR[i].Push(R2);
               }
            }
         }
      }
      for(int p = 0; p < pairCount; p++){
         if(!het1Counts[p] && !het2Counts[p] && !HetHetCounts[p]) continue;
         i1 = pairI[p*2]; i2 = pairI[p*2+1];
         kinship = (HetHetCounts[p] - IBS0Counts[p]*2.0) / (HetHetCounts[p]*2+het1Counts[p]+het2Counts[p]);
         notMissingCount = HetHetCounts[p]+het1Counts[p]+het2Counts[p]+HomHomCounts[p];
         ibs0 = IBS0Counts[p]*1.0/notMissingCount;
         CHet = HetHetCounts[p] * 1.0 / (HetHetCounts[p]+het1Counts[p]+het2Counts[p]);
         degree = 4;
         if(Phi[i1][i2] > 0.0442)
            degree = int(-log(Phi[i1][i2])/log(2.0) - 0.5);
         if(degree < 4){
            beforeCount[degree] ++;
            if(degree == 1)
            if(Phi[i2][i1]==0) beforeCount[5] ++;
         }else
            beforeCount[4] ++;
         ibdprop = ibdprops[p];
         ibd2prop = ibd2props[p];
         double pi = ibd2prop + ibdprop * 0.5;
         if(CHet<0.8){  // not MZ/Dup
            if(pi > 0.3535534){  // 1st-degree
               if(ibdprop + ibd2prop > 0.96 || (ibdprop + ibd2prop > 0.9 && ibd2prop <= 0.08))
                  type = "PO";
               else if(ibd2prop > 0.08)
                  type = "FS";
               else
                  type = "2nd";
            }else if(pi > 0.1767767){  // 2nd-degree
               if(pi > 0.32 && ibd2prop > 0.15)
                  type = "FS";
               else
                  type = "2nd";
            }else if(pi > 0.08838835)
               type = "3rd";
            else if(pi > 0.04419417)
               type = "4th";
            else
               type = "UN";
         }else // Duplicate
            type="Dup/MZ";
         if(type=="Dup/MZ"){
            afterCount[0]++;
         }else if(type=="PO"){
            afterCount[1]++;
            afterCount[5]++;
         }else if(type=="FS"){
            afterCount[1]++;
         }else if(type=="2nd"){
            afterCount[2]++;
         }else if(type=="3rd"){
            afterCount[3]++;
         }else if(type=="4th"){
            afterCount[4]++;
         }else{   // type="UN"
            afterCount[4]++;
         }
         exactType.Clear();
         if(type == "2nd"){
            int pr1Count = PR[i1].Length();
            int pr2Count = PR[i2].Length();
            int relative = -1;
            double maxPI = 0;
            for(int i = 0; i < pr1Count; i++)
               for(int j = 0; j < pr2Count; j++)
                  if(PR[i1][i] == PR[i2][j]){
                     int R = PR[i1][i];
                     double pi1 = pis[pair2Index[i1][R]];
                     double pi2 = pis[pair2Index[i2][R]];
                     if(pi1 < 0.355 && pi2 < 0.355 && (pi1 > 0.0884 || pi2 > 0.0884)){
                        double joinRatio = RoRP(ibd1segs[pair2Index[R][i1]],
                           ibd1segs[pair2Index[R][i2]], ibd1segs[pair2Index[i1][i2]], bp);
                        if(joinRatio > 0.8){ // relevant common relative
                           if(pi1 > maxPI) {maxPI = pi1; relative = R;}
                           if(pi2 > maxPI) {maxPI = pi2; relative = R;}
                        }
                     }
                  }
            if(relative != -1){
               printf("(%s, %s)'s common parental relative is %s, ",
                  (const char*)ped[id[f][i1]].pid,
                  (const char*)ped[id[f][i2]].pid,
                  (const char*)ped[id[f][relative]].pid);
               tempArray = ibd1segs[pair2Index[i1][relative]];  // PR
               tempArray2 = ibd1segs[pair2Index[i2][relative]];  // OR
               SegmentIntersect(tempArray, tempArray2, ibd1seg);
               double intersectLength = SegmentLength(ibd1seg, bp);
               double temp2Length = SegmentLength(tempArray2, bp);   // OR length in Mb
               double tempLength = SegmentLength(tempArray, bp);   // OR length in Mb
               if(temp2Length > 100000000){  // OR Length > 100Mb
                  double ratio1 = intersectLength * 1.0 / temp2Length;
                  int type1 = int(0.5-log(ratio1)/log(2));
                  double ratio2 = intersectLength * 1.0 / tempLength;
                  int type2 = int(0.5-log(ratio2)/log(2));
                  int diff = -9;
                  if(type1 <= 3 && type2 <= 3){
                     if(type1 < 0) type1 = 0;
                     else if(type1 > 2) type1 = 2;
                     if(type2 < 0) type2 = 0;
                     else if(type2 > 2) type2 = 2;
                     diff = ratio2 > ratio1? int(log(ratio2/ratio1)/log(2)+0.5):int(log(ratio2/ratio1)/log(2)-0.5) ;
                     switch(diff){
                        case 0: exactType = "HS"; break;
                        case 1: exactType = "NA/NU"; break;
                        case -1: exactType = "AN/UN"; break;
                        case 2: exactType = "GcGp"; break;
                        case -2: exactType = "GpGc"; break;
                     }
                  }
                  printf("G1=%d (%.2lf), G2=%d (%.2lf), Diff=%d (%.2lf), ExactType: %s\n",
                     type1, ratio1,
                     type2, ratio2,
                     diff, ratio1/ratio2,
                     exactType.IsEmpty()?"--":(const char*)exactType);
               }
            }else{   // no common relatives
               int IsHS1 = -1;
               for(int i = 0; i < pr1Count; i++){
                  int R1 = PR[i1][i];
                  if(R1==i2) continue;
                  if(RoRP(ibd1segs[pair2Index[i1][R1]], ibd1segs[pair2Index[i1][i2]],
                     ibd1segs[pair2Index[R1][i2]], bp) < 0.125){
                        IsHS1 = R1;
                        break;
                     }
               }
               int IsHS2 = -1;
               for(int i = 0; i < pr2Count; i++){
                  int R2 = PR[i2][i];
                  if(R2==i1) continue;
                  if(RoRP(ibd1segs[pair2Index[i2][R2]], ibd1segs[pair2Index[i2][i1]],
                     ibd1segs[pair2Index[R2][i1]], bp) < 0.125){
                        IsHS2 = R2;
                        break;
                     }
               }
               if(IsHS1 > -1 && IsHS2 > -1){
                  exactType = "HS";
                  printf("(%s, %s) must be HS for having unrelated relatives %s and %s\n",
                     (const char*)ped[id[f][i1]].pid,
                     (const char*)ped[id[f][i2]].pid,
                     (const char*)ped[id[f][IsHS1]].pid,
                     (const char*)ped[id[f][IsHS2]].pid);
               }

            }
         }
         fprintf(fp, "%s\t%s\t%s\t%d\t%.3lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%s\t%s\n",
            (const char*)ped[id[f][i1]].famid,
            (const char*)ped[id[f][i1]].pid,
            (const char*)ped[id[f][i2]].pid,
            notMissingCount, Phi[i2][i1], Phi[i1][i2],
            HetHetCounts[p]*1.0/notMissingCount, ibs0,
            CHet, // CHet
            ibs0/(ibs0+(IBSCounts[p]-HetHetCounts[p])*1.0/notMissingCount), kinship,
            ibdprops[p], ibd2props[p], ibd2props[p] + ibdprops[p]*0.5,
            (const char*)type, exactType.IsEmpty()? "--": (const char*)exactType);
      }  // end of pairs
      delete []pair2Index;
      delete []allrelatives;
      delete []PR;
      delete []ibd1segs;
   }  // End of f loop for families
   fclose(fp);

   bool pedigreeFlag = false;
   for(int i = 0; i < 6; i++)
      if(beforeCount[i]) pedigreeFlag = true;
   if(pedigreeFlag){
      printf("Within-family kinship data saved in file %s\n", (const char*)outfile);
      printRelationship(beforeCount, afterCount);
   }
}

void Engine::ComputeBigDataKinshipAdjustPC()
{
   if(shortCount==0) error("No genotype data");
   printf("Autosome genotypes stored in %d", shortCount);
   printf(" words for each of %d individuals.\n", idCount);

   printf("\nOptions in effect:\n");
   if(errorrateCutoff!=_NAN_)
      printf("\t--errorrate %G\n", errorrateCutoff);
   printf("\t--related\n");
   if(relativedegree>1)
      printf("\t--degree %d\n", relativedegree);
   printf("\t--adjustPC %d\n", BetaSum.Length()-1);
   if(bigdataFlag){
      if(faster)
         printf("\t--faster %d\n", faster);
      if(slower)
         printf("\t--slower %d\n", slower);
   }
   if(CoreCount)
      printf("\t--cpus %d\n", CoreCount);
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);

   printf("\n");

   int prescanCount[10];
   for(int i = 0; i < 10; i++)
      prescanCount[i] = 128 * (1<<(i*3));
   int short_prescan = relativedegree < 10? prescanCount[relativedegree-1]: prescanCount[9];
   if(short_prescan > shortCount || (!bigdataFlag)) short_prescan=shortCount;

   if(bigdataFlag){
      if(faster > 1){
         short_prescan /= faster;
         if(short_prescan < 10) short_prescan = 10;
      }
      if(slower > 1){
         short_prescan *= slower;
         if(short_prescan > shortCount) short_prescan = shortCount;
      }
   }

   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   unsigned long int onebit[65536];
   for(long int i = 0; i < 65536; i++){
      if(oneoneCount[i] > 8)
         onebit[i] = 0xFFFFFFFF;
      else{
         onebit[i] = 0;
         for(int j = 15; j >= 0; j--){
            if(i & shortbase[j]){
               onebit[i] <<= 4;
               onebit[i] |= j;
            }
         }
      }
   }
   unsigned long int ulong;
   double F1, F2;// inbreeding coefficient
   double Fst, var1, var2;
   double kinshipUnadjusted, kinshipAdjusted;
   int nCov = BetaSum.Length();
   if(nCov<2) printf("Dimension of beta is %d.\n", nCov);
   Vector *BetaSumPerID = new Vector[idCount];
   Matrix *BetaSquareSumPerID = new Matrix[idCount];
   for(int i = 0; i < idCount; i++){
      BetaSumPerID[i] = BetaSum;
      BetaSquareSumPerID[i] = BetaSquareSum;
   }
   Matrix pcX(idCount, nCov);
   for(int i = 0; i < idCount; i++){
      pcX[i][0] = 1.0;
      for(int j = 1; j < nCov; j++)
         pcX[i][j] = ped[phenoid[i]].covariates[covariatePC[j-1]];
   }
   Vector varP_ID(idCount);
   varP_ID.Zero();
   Vector varP(idCount);
   varP.Zero();

   double myBetaSum[100];
   double myBetaSquareSum[100][100];
   int start1, stop1, start2, stop2, start3, stop3;
   start1 = 0;
   stop1 = short_prescan;
   if(stop1 > shortCount) stop1 = shortCount;
   start2 = stop1;
   stop2 = stop1 + (short_prescan<<3);
   if(stop2 > shortCount) stop2 = shortCount;
   start3 = stop2;
   stop3 = shortCount;

   int relativeCount = 0;
   int rawrelativeCount = 0;
   int midrelativeCount = 0;
   double lowerbound = pow(2.0, -(relativedegree+1.5));
   double kincutoff1 = 0.125; // 1-degree
   if(short_prescan==shortCount)
      kincutoff1 = lowerbound;
   else if(relativedegree == 2)
      kincutoff1 = 0.04419417;
   else if(relativedegree == 3)
      kincutoff1 = 0.015625;
   else if(relativedegree > 3)
      kincutoff1 = 0;
   double kincutoff2 = sqrt(lowerbound * kincutoff1);
   int m1, m2, m3, m4;
   unsigned int **missingWordInOnePerson = new unsigned int *[idCount];
   int *missingWordInOnePersonCount[3];
   int *missingInOnePersonCount[3];
   int *hetInOnePersonCount[3];
   for(int i = 0; i < 3; i++){
      missingWordInOnePersonCount[i] = new int[idCount];
      missingInOnePersonCount[i] = new int[idCount];
      hetInOnePersonCount[i] = new int[idCount];
   }
   for(int i = 0; i < idCount; i++)
      for(int j = 0; j < 3; j++)
         missingWordInOnePersonCount[j][i] = missingInOnePersonCount[j][i] =
            hetInOnePersonCount[j][i] = 0;
   unsigned short int *masks = new unsigned short int[shortCount];
   for(int i = 0; i < shortCount;  i++)
      masks[i] = 0xFFFF;
   if(stop2 == shortCount)
      masks[stop2-1] = markerCount%16? ((1 << (markerCount % 16))-1): 0xFFFF;
   masks[shortCount-1] = markerCount%16? ((1 << (markerCount % 16))-1): 0xFFFF;

#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultEfficientCoreCount)
#endif
   for(int i = 0; i < idCount; i++){
      for(int m1 = 0; m1 < shortCount; m1++){
            int m3 = ~GG[0][i][m1] & ~GG[1][i][m1] & 0xFFFF;
            if(m3==0) continue;  // missing somewhere
            for(int m2 = 0; m2 < 16; m2++)
               if(m3 & shortbase[m2]){
                  int m4 = m1*16 + m2;
                  for(int u = 0; u < nCov; u++)
                     BetaSumPerID[i][u] -= freqBeta[m4][u];
                  for(int u = 0; u < nCov; u++)
                     for(int v = 0; v < nCov; v++)
                        BetaSquareSumPerID[i][u][v] -= freqBeta[m4][u] * freqBeta[m4][v];
            }
         }
      for(int u = 0; u < nCov; u++)
         varP_ID[i] += BetaSumPerID[i][u] * pcX[i][u];
      for(int u = 0; u < nCov; u++)
         for(int v = 0; v < nCov; v++)
            varP_ID[i] -= BetaSquareSumPerID[i][u][v] * pcX[i][u] * pcX[i][v];
      varP_ID[i] *= 2.0;
      for(int u = 0; u < nCov; u++)
         varP[i] += BetaSum[u] * pcX[i][u];
      for(int u = 0; u < nCov; u++)
         for(int v = 0; v < nCov; v++)
            varP[i] -= BetaSquareSum[u][v] * pcX[i][u] * pcX[i][v];
      varP[i] *= 2.0;
      IntArray tArray(0);
      for(int m = 0; m < stop1; m++)
         if(((~GG[0][i][m]) & (~GG[1][i][m]) & masks[m])!=0){
            missingWordInOnePersonCount[0][i] ++;
            missingInOnePersonCount[0][i] += oneoneCount[(~GG[0][i][m]) & (~GG[1][i][m]) & masks[m]];
            tArray.Push(m);
          }
      missingWordInOnePersonCount[1][i] = missingWordInOnePersonCount[0][i];
      for(int m = start2; m < stop2; m++)
         if(((~GG[0][i][m]) & (~GG[1][i][m]) & masks[m])!=0){
            missingWordInOnePersonCount[1][i] ++;
            missingInOnePersonCount[1][i] += oneoneCount[(~GG[0][i][m]) & (~GG[1][i][m]) & masks[m]];
            tArray.Push(m);
          }
      missingWordInOnePersonCount[2][i] = missingWordInOnePersonCount[1][i];
      for(int m = stop2; m < shortCount; m++)
         if(((~GG[0][i][m]) & (~GG[1][i][m]) & masks[m])!=0){
            missingWordInOnePersonCount[2][i] ++;
            missingInOnePersonCount[2][i] += oneoneCount[(~GG[0][i][m]) & (~GG[1][i][m]) & masks[m]];
            tArray.Push(m);
          }
      if(tArray.Length()){
         missingWordInOnePerson[i] = new unsigned int [tArray.Length()];
         for(int m = 0; m < tArray.Length(); m++)
            missingWordInOnePerson[i][m] = tArray[m];
      }else
         missingWordInOnePerson[i] = NULL;
      for(int m = 0; m < stop1; m++)
         hetInOnePersonCount[0][i] += oneoneCount[(~GG[0][i][m]) & GG[1][i][m]];
      for(int m = start2; m < stop2; m++)
         hetInOnePersonCount[1][i] += oneoneCount[(~GG[0][i][m]) & GG[1][i][m]];
      for(int m = stop2; m < shortCount; m++)
         hetInOnePersonCount[2][i] += oneoneCount[(~GG[0][i][m]) & GG[1][i][m]];
   }
   delete []masks;
   const int cutoffMissingCount = stop1*16-MINSNPCOUNT/10;

   String outfile;
   outfile.Copy(prefix);
   outfile.Add(".kin");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "FID\tID1\tID2\tN_SNP\tZ0\tPhi\tHetHet\tIBS0\tF1\tF2\tFst\tKinship\tError\n");
   double kinship, inflation, smaller, larger;
   int id1, id2;
   Kinship kin;
   int het1, het2, homhom, unequal;
   int HetHetCount, IBS0Count, het1Count, het2Count, notMissingCount, HetHomCount;
   double phi, pi0, errorFlag;
   int beforeCount[6], afterCount[6];
   Vector IBS0PO(0), IBS0FS(0), IBS0L1(0);
   int degree; double ibs0;
   for(int i = 0; i < 6; i++) beforeCount[i] = afterCount[i] = 0;
   for(int f = 0; f < ped.familyCount; f++){
      if(id[f].Length()<2) continue;
      kin.Setup(*ped.families[f]);
      for(int i = 0; i < id[f].Length(); i++)
         for(int j = i+1; j < id[f].Length(); j++){
            id1 = geno[id[f][i]]; id2 = geno[id[f][j]];
            HetHetCount = IBS0Count = het1Count = het2Count = notMissingCount = 0;
            for(int m = 0; m < shortCount; m++){
               HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
               IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
               notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
               het1Count += oneoneCount[GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m]]; // Het1Hom2
               het2Count += oneoneCount[GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]]; // Hom1Het2
            }
            HetHomCount = het1Count + het2Count;
            for(int u = 0; u < nCov; u++)
               myBetaSum[u] = BetaSumPerID[id1][u] + BetaSumPerID[id2][u] - BetaSum[u];
            for(int u = 0; u < nCov; u++)
               for(int v = 0; v < nCov; v++)
                  myBetaSquareSum[u][v] = BetaSquareSumPerID[id1][u][v] + BetaSquareSumPerID[id2][u][v] - BetaSquareSum[u][v];
            for(m1 = m2 = 0; m1 < missingWordInOnePersonCount[2][id1]; m1++){
               m3 = missingWordInOnePerson[id1][m1];
               for(; m2 < missingWordInOnePersonCount[2][id2] &&
                  missingWordInOnePerson[id2][m2] < m3; m2++);
               if( m2 == missingWordInOnePersonCount[2][id2] ) break;
               m4 = missingWordInOnePerson[id2][m2];
               if(m3 == m4){  // both missing
                  ulong = ~GG[0][id1][m4] & ~GG[1][id1][m4] & ~GG[0][id2][m4] & ~GG[1][id2][m4] & 0xFFFF;
                  if(ulong==0) continue;
                  ulong = onebit[ulong];
                  if(ulong==0xFFFFFFFF){  // more than 8 1's
                     ulong = ~GG[0][id1][m4] & ~GG[1][id1][m4] & ~GG[0][id2][m4] & ~GG[1][id2][m4] & 0xFFFF;
                     for(int k = 0; k < 16; k++)
                        if(ulong & shortbase[k]){
                           m3 = m4*16 + k;
                           for(int u = 0; u < nCov; u++)
                              myBetaSum[u] += freqBeta[m3][u];
                           for(int u = 0; u < nCov; u++)
                              for(int v = 0; v < nCov; v++)
                                 myBetaSquareSum[u][v] += freqBeta[m3][u] * freqBeta[m3][v];
                        }
                  }else{   // 1-8 1's
                     do{
                        m3 = m4*16 + (ulong&0xF);
                        for(int u = 0; u < nCov; u++)
                           myBetaSum[u] += freqBeta[m3][u];
                        for(int u = 0; u < nCov; u++)
                           for(int v = 0; v < nCov; v++)
                              myBetaSquareSum[u][v] += freqBeta[m3][u] * freqBeta[m3][v];
                        ulong >>= 4;
                     }while(ulong);
                  }
               }
            }
            var1 = 0;
            for(int u = 0; u < nCov; u++)
               var1 += myBetaSum[u] * pcX[id1][u];
            for(int u = 0; u < nCov; u++)
               for(int v = 0; v < nCov; v++)
                  var1 -= myBetaSquareSum[u][v] * pcX[id1][u] * pcX[id1][v];
            var1 *= 2.0;
            var2 = 0;
            for(int u = 0; u < nCov; u++)
               var2 += myBetaSum[u] * pcX[id2][u];
            for(int u = 0; u < nCov; u++)
               for(int v = 0; v < nCov; v++)
                  var2 -= myBetaSquareSum[u][v] * pcX[id2][u] * pcX[id2][v];
            var2 *= 2.0;
            Fst = 0.0;
            for(int u = 0; u < nCov; u++)
               for(int v = 0; v < nCov; v++)
                  Fst += myBetaSquareSum[u][v] * (pcX[id1][u]-pcX[id2][u]) * (pcX[id1][v]-pcX[id2][v]);
            het1Count += HetHetCount;
            het2Count += HetHetCount;
            F1 = 1.0 - het1Count/var1;
            F2 = 1.0 - het2Count/var2;
            var1 *= (1+F1);
            var2 *= (1+F2);
            kinship = (Fst + (var1+var2-HetHomCount)*0.25 - IBS0Count) / sqrt(var1*var2);
            Fst /= sqrt(var1*var2);

            phi = kin(ped[id[f][i]], ped[id[f][j]]);
            pi0 = 0.0;
            if(phi < 0.2)
               pi0 = 1-4*phi;
            else if(phi < 0.3 && ped[id[f][i]].isSib(ped[id[f][j]]))
               pi0 = 0.25;
            errorFlag = 0;
            inflation = (phi > 0)? kinship / phi: -1;
            if(phi > 0.03){  // up to 4th-degree relative
               if(inflation > 2 || inflation < 0.5)
                  errorFlag = 1;
               else if(inflation > 1.4142 || inflation < 0.70711)
                  errorFlag = 0.5;
            }else if(phi < 0.005){  // unrelated pair
               if(kinship > 0.0442) errorFlag = 1;
               else if(kinship > 0.0221) errorFlag = 0.5;
            }else{   // distant relatives
               if(kinship < -0.0221) errorFlag = 1;
               else errorFlag = 0.5;
            }
            ibs0 = IBS0Count*1.0/notMissingCount;
            degree = 4;
            if(phi > 0.0442)
               degree = int(-log(phi)/log(2.0) - 0.5);
            if(degree < 4){
               beforeCount[degree] ++;
               if(degree == 1)
                  if(pi0==0) beforeCount[5] ++;
            }else
               beforeCount[4] ++;
            degree = 4;
            if(kinship > 0.0442)
               degree = int(-log(kinship)/log(2.0) - 0.5);
            if(degree < 4){
               afterCount[degree] ++;
               if(degree == 1 && errorrateCutoff != _NAN_ && ibs0 < errorrateCutoff)
                  afterCount[5]++;
            }else
               afterCount[4] ++;

            if(degree==1 && phi==0.25){
               if(pi0==0)
                  IBS0PO.Push(ibs0);
               else
                  IBS0FS.Push(ibs0);
            }
            if(degree==1)
               IBS0L1.Push(ibs0);
            fprintf(fp, "%s\t%s\t%s\t%d\t%.3lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%G\n",
               (const char*)ped[id[f][i]].famid, (const char*)ped[id[f][i]].pid,
               (const char*)ped[id[f][j]].pid, notMissingCount, pi0, phi,
               HetHetCount*1.0/notMissingCount, ibs0, F1, F2, Fst, kinship, errorFlag);
         }
   }
   fclose(fp);

   bool pedigreeFlag = false;
   for(int i = 0; i < 6; i++)
      if(beforeCount[i]) pedigreeFlag = true;
   if(pedigreeFlag){
      if(errorrateCutoff==_NAN_){
         if(IBS0PO.Length() && IBS0FS.Length())
            errorrateCutoff = (IBS0PO.Sum()/IBS0PO.Length() + IBS0FS.Sum()/IBS0FS.Length())*0.5;
         else if(IBS0FS.Length())
            errorrateCutoff = IBS0FS.Sum()/IBS0FS.Length()*0.5;
         else if(IBS0PO.Length())
            errorrateCutoff = 0.005;
         for(int i = 0; i < IBS0L1.Length(); i++)
            if(IBS0L1[i] < errorrateCutoff)
               afterCount[5]++;
      }

      printf("Within-family kinship data saved in file %s\n", (const char*)outfile);
      printRelationship(beforeCount, afterCount);
   }else
      printf("Each family consists of one individual.\n");

   if(ped.familyCount < 2) {
      if(ped.familyCount==1 && ped.families[0]->famid=="0")
         warning("All individuals with family ID 0 are considered as relatives.\n");
      printf("There is only one family.\n");
      return;
   }
   printf("Relationship inference across families starts at %s", currentTime());

   char buffer[1024];
   StringArray buffers;
   buffers.Dimension(defaultMaxCoreCount);
   for(int i = 0; i < buffers.Length(); i++)
      buffers[i].Clear();

#ifdef _OPENMP
   printf("%d CPU cores are used.\n", defaultMaxCoreCount);
#endif

   const int LOOPBLOCKINGSIZE=100;
   IntArray loopIndex[2];
   loopIndex[0].Dimension(0); loopIndex[1].Dimension(0);
   for(int i = 0; i < idCount; i += LOOPBLOCKINGSIZE)
      for(int j = i; j < idCount; j += LOOPBLOCKINGSIZE){
         loopIndex[0].Push(i);
         loopIndex[1].Push(j);
      }
   int thread = 0;
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(HetHetCount, IBS0Count, het1Count, het2Count, notMissingCount, \
      HetHomCount, id1, id2, kinship, buffer, m1, m2, m3, m4, var1, var2, \
      myBetaSquareSum, myBetaSum, F1, F2, Fst, thread) \
      reduction(+:rawrelativeCount, midrelativeCount, relativeCount)
{
   thread = omp_get_thread_num();
   #pragma omp for
#endif
   for(int k = 0; k < loopIndex[0].Length(); k++){
      int i = loopIndex[0][k];
      int j = loopIndex[1][k];
      for(id1 = i; (id1 < i + LOOPBLOCKINGSIZE) && (id1 < idCount); id1++){
         if(missingInOnePersonCount[0][id1]>=cutoffMissingCount) continue;
         for(id2 = j; (id2 < j + LOOPBLOCKINGSIZE) && (id2 < idCount); id2++){
            if(i==j && id1 >= id2) continue;
            if(missingInOnePersonCount[0][id2]>=cutoffMissingCount) continue;
            if(ped[phenoid[id1]].famid == ped[phenoid[id2]].famid) continue;
            het1Count = hetInOnePersonCount[0][id1] + hetInOnePersonCount[1][id1] + hetInOnePersonCount[2][id1];
            het2Count = hetInOnePersonCount[0][id2] + hetInOnePersonCount[1][id2] + hetInOnePersonCount[2][id2];
            F1 = 1.0 - het1Count/varP_ID[id1];
            F2 = 1.0 - het2Count/varP_ID[id2];
            var1 = (1+F1)*varP[id1];
            var2 = (1+F2)*varP[id2];
            Fst = 0.0;
            for(int u = 0; u < nCov; u++)
               for(int v = 0; v < nCov; v++)
                  Fst += BetaSquareSum[u][v] * (pcX[id1][u]-pcX[id2][u]) * (pcX[id1][v]-pcX[id2][v]);
            if(Fst > 0.177 * sqrt(var1*var2)) continue;  // from two different populations
            // Stage 1: all pairs
            double cutoffDist = ((Fst + (var1+var2)*0.25)/sqrt(var1*var2) - kincutoff1) *
               4 * sqrt(1.0*hetInOnePersonCount[0][id1]*hetInOnePersonCount[0][id2]/((1-F1)*(1-F2)));
            for(HetHomCount = -missingInOnePersonCount[0][id1]-missingInOnePersonCount[0][id2],
               m1 = 0; m1 < stop1; m1++) // Lower bound of HetHom Count, minus HetMiss and MissMiss
                  HetHomCount += oneoneCount[GG[0][id1][m1]^GG[0][id2][m1]];
            if(HetHomCount > cutoffDist) continue;
            for(IBS0Count=0, m1 = 0; m1 < stop1; m1++)
               IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
            if(HetHomCount+IBS0Count*4 > cutoffDist) continue;
            for(m1 = m2 = 0; m1 < missingWordInOnePersonCount[0][id1]; m1++){
               m3 = missingWordInOnePerson[id1][m1]; // m for id1
               for(; m2 < missingWordInOnePersonCount[0][id2] &&
                  missingWordInOnePerson[id2][m2] < m3; m2++);
               if( m2 == missingWordInOnePersonCount[0][id2] ) break;
               m4 = missingWordInOnePerson[id2][m2];  // m for id2
               if(m3 == m4)   // id1 and id2 are both missing
                  HetHomCount += 2 * oneoneCount[(~GG[0][id1][m3])&(~GG[1][id1][m3])&(~GG[0][id2][m4])&(~GG[1][id2][m4])&0xFFFF];
            }
            for(het1Count = hetInOnePersonCount[0][id1], m1 = 0;
               m1 < missingWordInOnePersonCount[0][id2]; m1++){
                  m2 = missingWordInOnePerson[id2][m1];
                  m3 = oneoneCount[(~GG[0][id2][m2])&(~GG[1][id2][m2])&(~GG[0][id1][m2])&GG[1][id1][m2]];
                  het1Count -= m3;
                  HetHomCount += m3;
            }
            for(het2Count = hetInOnePersonCount[0][id2], m1 = 0;
               m1 < missingWordInOnePersonCount[0][id1]; m1++){
                  m2 = missingWordInOnePerson[id1][m1];
                  m3 = oneoneCount[(~GG[0][id1][m2])&(~GG[1][id1][m2])&(~GG[0][id2][m2])&GG[1][id2][m2]];
                  het2Count -= m3;
                  HetHomCount += m3;
            }
            if(HetHomCount+IBS0Count*4 >= cutoffDist) continue;
            rawrelativeCount ++;
            // Stage 2: few pairs
            for(HetHomCount -= (missingInOnePersonCount[1][id1]+missingInOnePersonCount[1][id2]),  // different non-missing genotypes
               m1 = start2; m1 < stop2; m1++) // different-genotpe Count
                  HetHomCount += oneoneCount[GG[0][id1][m1]^GG[0][id2][m1]];
            het1Count += hetInOnePersonCount[1][id1];
            het2Count += hetInOnePersonCount[1][id2];
            for(m1 = start2; m1 < stop2; m1++)
               IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
            // now add the MissMissCount back
            for(m1 = missingWordInOnePersonCount[0][id1], m2 = missingWordInOnePersonCount[0][id2];
               m1 < missingWordInOnePersonCount[1][id1]; m1++){
               m3 = missingWordInOnePerson[id1][m1]; // m for id1
               for(; m2 < missingWordInOnePersonCount[1][id2] && missingWordInOnePerson[id2][m2] < m3; m2++);
               if(m2 == missingWordInOnePersonCount[1][id2] ) break;
               m4 = missingWordInOnePerson[id2][m2];  // m for id2
               if(m3 == m4)   // id1 and id2 are both missing
                  HetHomCount += 2 * oneoneCount[(~GG[0][id1][m3])&(~GG[1][id1][m3])&(~GG[0][id2][m4])&(~GG[1][id2][m4])&0xFFFF];
            }
            for(m1 = missingWordInOnePersonCount[0][id2]; m1 < missingWordInOnePersonCount[1][id2]; m1++){
               m2 = missingWordInOnePerson[id2][m1];
               m3 = oneoneCount[(~GG[0][id2][m2])&(~GG[1][id2][m2])&(~GG[0][id1][m2])&GG[1][id1][m2]];
               het1Count -= m3;
               HetHomCount += m3;
            }
            for(m1 = missingWordInOnePersonCount[0][id1]; m1 < missingWordInOnePersonCount[1][id1]; m1++){
               m2 = missingWordInOnePerson[id1][m1];
               m3 = oneoneCount[(~GG[0][id1][m2])&(~GG[1][id1][m2])&(~GG[0][id2][m2])&GG[1][id2][m2]];
               het2Count -= m3;
               HetHomCount += m3;
            }
            HetHetCount = (het1Count + het2Count - HetHomCount)/2;
            kinship = (Fst+(var1+var2)*0.25)/sqrt(var1*var2) -
               (HetHomCount*0.25+IBS0Count)*sqrt((1-F1)*(1-F2)/het1Count/het2Count);
            if(kinship < kincutoff2) continue;
            midrelativeCount ++;

            // Stage 3: exact
            for(HetHomCount -= (missingInOnePersonCount[2][id1]+missingInOnePersonCount[2][id2]),  // different non-missing genotypes
               m1 = start3; m1 < stop3; m1++) // different-genotpe Count
                  HetHomCount += oneoneCount[GG[0][id1][m1]^GG[0][id2][m1]];
            het1Count += hetInOnePersonCount[2][id1];
            het2Count += hetInOnePersonCount[2][id2];
            for(m1 = start3; m1 < stop3; m1++)
               IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
            // now add the MissMissCount back
            for(m1 = missingWordInOnePersonCount[1][id1], m2 = missingWordInOnePersonCount[1][id2];
               m1 < missingWordInOnePersonCount[2][id1]; m1++){
               m3 = missingWordInOnePerson[id1][m1]; // m for id1
               for(; m2 < missingWordInOnePersonCount[2][id2] && missingWordInOnePerson[id2][m2] < m3; m2++);
               if(m2 == missingWordInOnePersonCount[2][id2]) break;
               m4 = missingWordInOnePerson[id2][m2];  // m for id2
               if(m3 == m4)   // id1 and id2 are both missing
                  HetHomCount += 2 * oneoneCount[(~GG[0][id1][m3])&(~GG[1][id1][m3])&(~GG[0][id2][m4])&(~GG[1][id2][m4])&0xFFFF];
            }
            for(m1 = missingWordInOnePersonCount[1][id2]; m1 < missingWordInOnePersonCount[2][id2]; m1++){
               m2 = missingWordInOnePerson[id2][m1];
               m3 = oneoneCount[(~GG[0][id2][m2])&(~GG[1][id2][m2])&(~GG[0][id1][m2])&GG[1][id1][m2]];
               het1Count -= m3;
               HetHomCount += m3;
            }
            for(m1 = missingWordInOnePersonCount[1][id1]; m1 < missingWordInOnePersonCount[2][id1]; m1++){
               m2 = missingWordInOnePerson[id1][m1];
               m3 = oneoneCount[(~GG[0][id1][m2])&(~GG[1][id1][m2])&(~GG[0][id2][m2])&GG[1][id2][m2]];
               het2Count -= m3;
               HetHomCount += m3;
            }
            HetHetCount = (het1Count + het2Count - HetHomCount)/2;
            for(int u = 0; u < nCov; u++)
               myBetaSum[u] = BetaSumPerID[id1][u] + BetaSumPerID[id2][u] - BetaSum[u];
            for(int u = 0; u < nCov; u++)
               for(int v = 0; v < nCov; v++)
                  myBetaSquareSum[u][v] = BetaSquareSumPerID[id1][u][v] + BetaSquareSumPerID[id2][u][v] - BetaSquareSum[u][v];
            for(m1 = m2 = 0; m1 < missingWordInOnePersonCount[2][id1]; m1++){
               m3 = missingWordInOnePerson[id1][m1];
               for(; m2 < missingWordInOnePersonCount[2][id2] &&
                  missingWordInOnePerson[id2][m2] < m3; m2++);
               if( m2 == missingWordInOnePersonCount[2][id2] ) break;
               m4 = missingWordInOnePerson[id2][m2];
               if(m3 == m4){  // both missing
                  ulong = ~GG[0][id1][m4] & ~GG[1][id1][m4] & ~GG[0][id2][m4] & ~GG[1][id2][m4] & 0xFFFF;
                  if(ulong==0) continue;  // no MissMiss
                  ulong = onebit[ulong];  // 1's positions
                  if(ulong==0xFFFFFFFF){  // more than 8 1's
                     ulong = ~GG[0][id1][m4] & ~GG[1][id1][m4] & ~GG[0][id2][m4] & ~GG[1][id2][m4] & 0xFFFF;
                     for(int kk = 0; kk < 16; kk++)
                        if(ulong & shortbase[kk]){
                           m3 = m4*16 + kk;
                           for(int u = 0; u < nCov; u++)
                              myBetaSum[u] += freqBeta[m3][u];
                           for(int u = 0; u < nCov; u++)
                              for(int v = 0; v < nCov; v++)
                                 myBetaSquareSum[u][v] += freqBeta[m3][u] * freqBeta[m3][v];
                        }
                  }else{   // 1-8 1's
                     do{
                        m3 = m4*16 + (ulong&0xF);
                        for(int u = 0; u < nCov; u++)
                           myBetaSum[u] += freqBeta[m3][u];
                        for(int u = 0; u < nCov; u++)
                           for(int v = 0; v < nCov; v++)
                              myBetaSquareSum[u][v] += freqBeta[m3][u] * freqBeta[m3][v];
                        ulong >>= 4;
                     }while(ulong);
                  }
               }
            }
            var1 = 0;
            for(int u = 0; u < nCov; u++)
               var1 += myBetaSum[u] * pcX[id1][u];
            for(int u = 0; u < nCov; u++)
               for(int v = 0; v < nCov; v++)
                  var1 -= myBetaSquareSum[u][v] * pcX[id1][u] * pcX[id1][v];
            var1 *= 2.0;
            var2 = 0;
            for(int u = 0; u < nCov; u++)
               var2 += myBetaSum[u] * pcX[id2][u];
            for(int u = 0; u < nCov; u++)
               for(int v = 0; v < nCov; v++)
                  var2 -= myBetaSquareSum[u][v] * pcX[id2][u] * pcX[id2][v];
            var2 *= 2.0;
            Fst = 0.0;
            for(int u = 0; u < nCov; u++)
               for(int v = 0; v < nCov; v++)
                  Fst += myBetaSquareSum[u][v] * (pcX[id1][u]-pcX[id2][u]) * (pcX[id1][v]-pcX[id2][v]);
            F1 = 1.0 - het1Count/var1;
            F2 = 1.0 - het2Count/var2;
            var1 *= (1+F1);
            var2 *= (1+F2);
            kinship = (Fst + (var1+var2-HetHomCount)*0.25 - IBS0Count)/sqrt(var1*var2);
            Fst /= sqrt(var1*var2);
            if(kinship < lowerbound) continue;
            relativeCount ++;
            notMissingCount = 0;
            for(int m = 0; m < shortCount; m++)
               notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
            sprintf(buffer, "%s\t%s\t%s\t%s\t%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\n",
               (const char*)ped[phenoid[id1]].famid, (const char*)ped[phenoid[id1]].pid,
               (const char*)ped[phenoid[id2]].famid, (const char*)ped[phenoid[id2]].pid,
               notMissingCount, HetHetCount*1.0/notMissingCount, IBS0Count*1.0/notMissingCount, F1, F2, Fst, kinship);
            buffers[thread].Add(buffer);
         }
      }
   }
#ifdef _OPENMP
}
#endif

   for(int i = 0; i < idCount; i++)
      if(missingWordInOnePerson[i]) delete []missingWordInOnePerson[i];
   for(int i = 0; i < 3; i++){
      if(missingWordInOnePersonCount[i]) delete []missingWordInOnePersonCount[i];
      if(missingInOnePersonCount[i]) delete []missingInOnePersonCount[i];
      if(hetInOnePersonCount[i]) delete []hetInOnePersonCount[i];
   }

   outfile.Copy(prefix);
   outfile.Add(".kin0");
   fp = fopen(outfile, "wt");
   fprintf(fp, "FID1\tID1\tFID2\tID2\tN_SNP\tHetHet\tIBS0\tF1\tF2\tFst\tKinship\n");
   for(int b = 0; b < buffers.Length(); b++)
      buffers[b].Write(fp);
   fclose(fp);
   printf("                                         ends at %s", currentTime());
   if(relativedegree){
      printf("Relationship inference speed-up by screening top informative SNPs:\n");
      printf("  Stage 1 (with %d SNPs): %d relative pairs detected (with Kinship > %.5lf)\n",
         stop1*16>markerCount? markerCount: stop1*16,
         rawrelativeCount, kincutoff1);
      printf("  Stage 2 (with %d SNPs): %d out of %d remain related (with Kinship > %.5lf)\n",
         stop2*16>markerCount? markerCount: stop2*16,
         midrelativeCount, rawrelativeCount, kincutoff2);
      printf("  Final Stage (with %d SNPs): %d pairs of relatives up to %d-degree are identified\n",
         markerCount, relativeCount, relativedegree);
      printf("Between-family relatives (kinship >= %.5lf) saved in file %s\n",
         lowerbound, (const char*)outfile);
      if(relativedegree == 1){
         printf("\nNote only duplicates and 1st-degree relative pairs are included here.\n");
         printf("A higher degree of relationship can be also included by specifying ");
         printf("'--degree 2' (for up to 2nd-degree) or '--degree 3' (for up to 3rd-degree).\n");
      }
      printf("\n");
   }
}

/*
void Engine::ComputeBigDataSecondDegree()
{
   printf("Autosome genotypes stored in %d", shortCount);
   printf(" words for each of %d individuals.\n", idCount);

   printf("\nOptions in effect:\n");
   printf("\t--close\n");
   if(minMAF!=_NAN_)
      printf("\t--minMAF %lf\n", minMAF);
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");

   int id1, id2, id3;
   int HetHetCount, Het1Het3Count, Het2Het3Count, nonMissingCount, HetHetHetCount;
   int het1, het2, het3, nonmissing1, nonmissing2, nonmissing3;

   String outfile;
   outfile.Copy(prefix);
   outfile.Add(".d23");
   printf("There are %d rare SNPs with MAF < %lf.\n", veryrareCount, minMAF==_NAN_? 0.01: minMAF);
   FILE *fp = fopen(outfile, "wt");
   printf("Inferring exact second/thrid-degree relationship starts at %s", currentTime());

   fprintf(fp, "FID1\tID1\tFID2\tID2\tKinship\tFID3\tID3\tN_HHH\tH_3|1\tH_3|2\tH_3|12\tH_1|23\tH_2|13\tN_Gen1\tN_Gen2\tN_Ancestor\tGenDiff\n");
   int m1, m2;
   double kinship, smaller;
   int IBS0Count, het1Count, het2Count, notMissingCount;
   double H31, H32, H312, H123, H213;
   double R12; // H123 / H213;
   double increase;

   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   IntArray SNPCount(idCount);
   SNPCount.Zero();
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = 0; i < id[f].Length(); i++){
         notMissingCount = 0;
         for(int m = 0; m < shortCount; m++)
            notMissingCount += oneoneCount[GG[0][geno[id[f][i]]][m] | GG[1][geno[id[f][i]]][m]];
         SNPCount[geno[id[f][i]]] = notMissingCount;
         if(notMissingCount < 1000) {
            printf("Fam %s ID %s excluded for having %d non-missing SNPs\n",
               (const char*)ped[id[f][i]].famid, (const char*)ped[id[f][i]].pid, notMissingCount);
            continue;
         }
      }

   int tdegree=3;
   int prescanCount[10];
   for(int i = 0; i < 10; i++)
      prescanCount[i] = 128 * (1<<(i*3));
   int short_prescan = tdegree < 10? prescanCount[tdegree-1]: prescanCount[9];
   if(short_prescan > shortCount || (!bigdataFlag)) short_prescan=shortCount;
   double lowerbound = pow(2.0, -(tdegree+1.5));
   double kincutoff1 = 0.125; // 1-degree
   if(short_prescan==shortCount)
      kincutoff1 = lowerbound;
   else if(tdegree == 2)
      kincutoff1 = 0.04419417;
   else if(tdegree == 3)
      kincutoff1 = 0.015625;
   else if(tdegree > 3)
      kincutoff1 = 0;
   double kincutoff2 = sqrt(lowerbound * kincutoff1);
   int rel=0;
   bool controversial;
   int genDiff, genMax;
   int type_d2[][4]={{0, 1, 0, 0}, {0, 2, 0, 0}, {0, 0, 3, 0}, {0,0,0,0}};
   int type_d3[][4]={{0, 1, 0, 0}, {0, 0, 2, 0}, {0, 0, 3, 0}, {0,0,0,4}};
   StringArray label_d2(0);
   label_d2.Push("N/V"); label_d2.Push("HS"); label_d2.Push("AV"); label_d2.Push("GG");
   StringArray label_d3(0);
   label_d3.Push("N/V"); label_d3.Push("FC"); label_d3.Push("UN3D");
   label_d3.Push("GAV"); label_d3.Push("GGG");
   int gen1Count, gen1bCount, gen2Count, ancestorCount;

   for(int f1 = 0; f1 < ped.familyCount; f1++)
      for(int i = 0; i < id[f1].Length(); i++){
         id1 = geno[id[f1][i]];
         if(SNPCount[id1]<1000) continue;
         for(int f2 = f1+1; f2 < ped.familyCount; f2++)
            for(int j = 0; j < id[f2].Length(); j++){
               id2 = geno[id[f2][j]];
               if(SNPCount[id2]<1000) continue;
               // Stage 1: many pairs
               HetHetCount = IBS0Count = het1Count = het2Count = 0;
               for(int m = 0; m < short_prescan; m++){
                  HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
                  IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
               }
               if(HetHetCount + 2 * IBS0Count > 500){
                  if(HetHetCount <= 2*IBS0Count) continue;
                  for(int m = 0; m < short_prescan; m++){
                     het1Count += oneoneCount[GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m]]; // Het1Hom2
                     het2Count += oneoneCount[GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]]; // Hom1Het2
                  }
                  smaller = HetHetCount + (het1Count < het2Count? het1Count: het2Count);
                  kinship = 0.5 - ((het1Count+het2Count)*0.25+IBS0Count)/smaller;
                  if(kinship < kincutoff1) continue;
               }else{ // insufficient informative SNPs
                  for(int m = 0; m < short_prescan; m++){
                     het1Count += oneoneCount[GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m]]; // Het1Hom2
                     het2Count += oneoneCount[GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]]; // Hom1Het2
                  }
               }
               m1 = short_prescan;
               // Stage 2: few pairs
               if(tdegree < 4){
                  m2 = m1 + (short_prescan<<3);
                  if(m2 > shortCount) m2 = shortCount;
                  for(int m = m1; m < m2; m++){
                     HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
                     IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
                  }
                  if(HetHetCount <= 2*IBS0Count) continue;
                  for(int m = m1; m < m2; m++){
                     het1Count += oneoneCount[GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m]]; // Het1Hom2
                     het2Count += oneoneCount[GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]]; // Hom1Het2
                  }
                  smaller = HetHetCount + (het1Count < het2Count? het1Count: het2Count);
                  kinship = 0.5 - ((het1Count+het2Count)*0.25+IBS0Count)/smaller;
                  if(kinship < kincutoff2) continue;
                  m1 = m2;
               }
               // Stage 3: exact
               for(int m = m1; m < shortCount; m++){
                  HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
                  IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
               }
               if(HetHetCount <= 2*IBS0Count) continue;
               for(int m = m1; m < shortCount; m++){
                  het1Count += oneoneCount[GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m]]; // Het1Hom2
                  het2Count += oneoneCount[GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]]; // Hom1Het2
               }
               smaller = HetHetCount + (het1Count < het2Count? het1Count: het2Count);
               kinship = 0.5 - ((het1Count+het2Count)*0.25+IBS0Count)/smaller;
               if(kinship < lowerbound) continue;
               if(kinship > 0.1767767) continue;
               rel = 0;
               controversial=false;
               for(int f3 = 0; f3 < ped.familyCount; f3++)
               for(int k = 0; k < id[f3].Length(); k++){
                  id3 = geno[id[f3][k]];
                  if(id3==id1) continue;
                  if(id3==id2) continue;
                  if(SNPCount[id3]<1000) continue;
                  HetHetCount = HetHetHetCount = het1Count
                     = het2Count = Het1Het3Count = Het2Het3Count = 0;
                  for(int m = shortCount-1; m > shortCount - (veryrareCount / 16); m--){
                     het1 = (~GG[0][id1][m]) & GG[1][id1][m];
                     nonmissing1 = GG[0][id1][m] | GG[1][id1][m];
                     het2 = (~GG[0][id2][m]) & GG[1][id2][m];
                     nonmissing2 = GG[0][id2][m] | GG[1][id2][m];
                     het3 = (~GG[0][id3][m]) & GG[1][id3][m];
                     nonmissing3 = GG[0][id3][m] | GG[1][id3][m];
                     het1 &= nonmissing2 & nonmissing3;
                     het2 &= nonmissing1 & nonmissing3;
                     het3 &= nonmissing1 & nonmissing2;
                     het1Count += oneoneCount[het1];
                     het2Count += oneoneCount[het2];
                     HetHetCount += oneoneCount[het1 & het2];
                     Het1Het3Count += oneoneCount[het1 & het3];
                     Het2Het3Count += oneoneCount[het2 & het3];
                     HetHetHetCount += oneoneCount[het1 & het2 & het3];
                  }  // end of marker
                  if(HetHetHetCount < 10) continue;
                  H123 = (double)HetHetHetCount / Het2Het3Count;
                  H213 = (double)HetHetHetCount / Het1Het3Count;
                  if(kinship > 0.08838835 && H123 < 0.375 && H213 < 0.375) continue;   // no information for 2nd
                  if(kinship < 0.08838835 && H123 < 0.1875 && H213 < 0.1875) continue;  // no information for 3rd
                  if(kinship > 0.08838835 && het1Count < het2Count * 0.1767767){
                     printf("Relationship between %s:%s and %s:%s (kinship=%.3lf) may be third-degree\n",
                        (const char*)ped[id[f1][i]].famid,
                        (const char*)ped[id[f1][i]].pid,
                        (const char*)ped[id[f2][j]].famid,
                        (const char*)ped[id[f2][j]].pid,
                        kinship);
                     continue;
                  }
                  if(kinship < 0.08838835 && het1Count < het2Count * 0.08838835){
                     printf("Relationship between %s:%s and %s:%s (kinship=%.3lf) may be fourth-degree\n",
                        (const char*)ped[id[f1][i]].famid,
                        (const char*)ped[id[f1][i]].pid,
                        (const char*)ped[id[f2][j]].famid,
                        (const char*)ped[id[f2][j]].pid,
                        kinship);
                     continue;
                  }
                  H312 = (double)HetHetHetCount / HetHetCount;
                  H31 = (double)Het1Het3Count / het1Count;
                  H32 = (double)Het2Het3Count / het2Count;
                  if(H312<0.001) continue; // too distant to be useful
                  if(H312>0.3535534) continue; // existence of 1st-degree relative, use better info
                  gen1Count = int(-log((double)HetHetHetCount / Het2Het3Count) / 0.6931472 + 0.5);
                  gen1bCount = int(-log(H31 / H312) / 0.6931472 + 0.5);
                  gen2Count = int(-log((double)HetHetHetCount / Het1Het3Count) / 0.6931472 + 0.5);
                  if(gen1Count == gen1bCount) ancestorCount = 1;
                  else if(gen1Count == gen1bCount+1) ancestorCount = 2;
                  else continue;
                  genDiff = int(fabs(log((double)Het1Het3Count / Het2Het3Count) / 0.6931472) + 0.5);

                  fprintf(fp, "%s\t%s\t%s\t%s\t%.3lf\t%s\t%s\t%d\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%d\t%d\t%d\t%d",
                     (const char*)ped[id[f1][i]].famid,
                     (const char*)ped[id[f1][i]].pid,
                     (const char*)ped[id[f2][j]].famid,
                     (const char*)ped[id[f2][j]].pid,
                     kinship,
                     (const char*)ped[id[f3][k]].famid,
                     (const char*)ped[id[f3][k]].pid,
                     HetHetHetCount, H31, H32, H312,
                     H123, H213, gen1Count, gen2Count, ancestorCount, genDiff);
                     fprintf(fp, "\n");
               } // end of k
         }  // end of j
      }  // end of i

   fclose(fp);
   printf("                              ends at %s", currentTime());
   printf("Exact inference of second-degree relationships is saved in file %s\n",
         (const char*)outfile);
}
*/
void Engine::ComputeBigDataPO()
{
   printf("Autosome genotypes stored in %d", shortCount);
   printf(" words for each of %d individuals.\n", idCount);

   printf("\nOptions in effect:\n");
   printf("\t--porel\n");
   if(errorrateCutoff!=_NAN_)
      printf("\t--errorrate %lf\n", errorrateCutoff);
   if(minMAF!=_NAN_)
      printf("\t--minMAF %lf\n", minMAF);
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");

   int id1, id2, id3;
   int HetHetCount, nonMissingCount, HetCount;
   int het1, het2, het3, nonmissing1, nonmissing2, nonmissing3;
   int het1Count, het2Count, het3Count, het12Count, het13Count, het23Count, het123Count;
   double H21, H23, H213, H12, H13, H123;

   String outfile;
   outfile.Copy(prefix);
   outfile.Add(".por");
   printf("There are %d rare SNPs with MAF < %lf.\n", veryrareCount, minMAF==_NAN_? 0.01: minMAF);
   FILE *fp = fopen(outfile, "wt");
   printf("Inferring relatives of parent/offspring starts at %s", currentTime());

   fprintf(fp, "FID_OFF\tID_OFF\tFID_PAR\tID_PAR\tFID_REL\tID_REL\tN_OR\tH_P|O\tH_P|R\tH_P|OR\n");
   int m1, m2;
   double kinship, smaller;
   int IBS0Count, notMissingCount;

   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   IntArray SNPCount(idCount);
   SNPCount.Zero();
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = 0; i < id[f].Length(); i++){
         notMissingCount = 0;
         for(int m = 0; m < shortCount; m++)
            notMissingCount += oneoneCount[GG[0][geno[id[f][i]]][m] | GG[1][geno[id[f][i]]][m]];
         SNPCount[geno[id[f][i]]] = notMissingCount;
         if(notMissingCount < 1000) {
            printf("Fam %s ID %s excluded for having %d non-missing SNPs\n",
               (const char*)ped[id[f][i]].famid, (const char*)ped[id[f][i]].pid, notMissingCount);
            continue;
         }
      }

   int tdegree=1;
   if(errorrateCutoff==_NAN_) errorrateCutoff=0.001;

   int prescanCount[10];
   for(int i = 0; i < 10; i++)
      prescanCount[i] = 128 * (1<<(i*3));
   int short_prescan = tdegree < 10? prescanCount[tdegree-1]: prescanCount[9];
   if(short_prescan > shortCount || (!bigdataFlag)) short_prescan=shortCount;
   double lowerbound = pow(2.0, -(tdegree+1.5));
   double kincutoff1 = 0.125; // 1-degree
   if(short_prescan==shortCount)
      kincutoff1 = lowerbound;
   else if(tdegree == 2)
      kincutoff1 = 0.04419417;
   else if(tdegree == 3)
      kincutoff1 = 0.015625;
   else if(tdegree > 3)
      kincutoff1 = 0;
   double kincutoff2 = sqrt(lowerbound * kincutoff1);

   for(int f1 = 0; f1 < ped.familyCount; f1++)
      for(int i = 0; i < id[f1].Length(); i++){
         id1 = geno[id[f1][i]];
         if(SNPCount[id1]<1000) continue;
         for(int f2 = f1; f2 < ped.familyCount; f2++)
            for(int j = 0; j < id[f2].Length(); j++){
               id2 = geno[id[f2][j]];
               if(id2==id1) continue;
               if(SNPCount[id2]<1000) continue;
               // Stage 1: many pairs
               HetHetCount = IBS0Count = het1Count = het2Count = 0;
               for(int m = 0; m < short_prescan; m++){
                  HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
                  IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
               }
               if(HetHetCount + 2 * IBS0Count > 500){
                  if(HetHetCount <= 2*IBS0Count) continue;
                  for(int m = 0; m < short_prescan; m++){
                     het1Count += oneoneCount[GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m]]; // Het1Hom2
                     het2Count += oneoneCount[GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]]; // Hom1Het2
                  }
                  smaller = HetHetCount + (het1Count < het2Count? het1Count: het2Count);
                  kinship = 0.5 - ((het1Count+het2Count)*0.25+IBS0Count)/smaller;
                  if(kinship < kincutoff1) continue;
               }else{ // insufficient informative SNPs
                  for(int m = 0; m < short_prescan; m++){
                     het1Count += oneoneCount[GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m]]; // Het1Hom2
                     het2Count += oneoneCount[GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]]; // Hom1Het2
                  }
               }
               m1 = short_prescan;
               // Stage 2: few pairs
               if(tdegree < 4){
                  m2 = m1 + (short_prescan<<3);
                  if(m2 > shortCount) m2 = shortCount;
                  for(int m = m1; m < m2; m++){
                     HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
                     IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
                  }
                  if(HetHetCount <= 2*IBS0Count) continue;
                  for(int m = m1; m < m2; m++){
                     het1Count += oneoneCount[GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m]]; // Het1Hom2
                     het2Count += oneoneCount[GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]]; // Hom1Het2
                  }
                  smaller = HetHetCount + (het1Count < het2Count? het1Count: het2Count);
                  kinship = 0.5 - ((het1Count+het2Count)*0.25+IBS0Count)/smaller;
                  if(kinship < kincutoff2) continue;
                  m1 = m2;
               }
               // Stage 3: exact
               for(int m = m1; m < shortCount; m++){
                  HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
                  IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
               }
               if(HetHetCount <= 2*IBS0Count) continue;
               for(int m = m1; m < shortCount; m++){
                  het1Count += oneoneCount[GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m]]; // Het1Hom2
                  het2Count += oneoneCount[GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]]; // Hom1Het2
               }
               smaller = HetHetCount + (het1Count < het2Count? het1Count: het2Count);
               kinship = 0.5 - ((het1Count+het2Count)*0.25+IBS0Count)/smaller;
               if(kinship < lowerbound) continue;
               if(kinship > 0.375) continue;
               notMissingCount = 0;
               for(int m = 0; m < shortCount; m++)
                  notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
               if(IBS0Count > errorrateCutoff * notMissingCount) continue;
 /*
               printf("PO pair %s:%s and %s:%s (kinship=%.3lf)...\n",
                        (const char*)ped[id[f1][i]].famid,
                        (const char*)ped[id[f1][i]].pid,
                        (const char*)ped[id[f2][j]].famid,
                        (const char*)ped[id[f2][j]].pid,
                        kinship);
   */
               for(int f3 = 0; f3 < ped.familyCount; f3++)
               for(int k = 0; k < id[f3].Length(); k++){
                  id3 = geno[id[f3][k]];
                  if(id3==id1) continue;
                  if(id3==id2) continue;
                  if(SNPCount[id3]<1000) continue;

                  het1Count = het2Count = het3Count = het12Count =
                     het13Count = het23Count = het123Count = 0;
                  for(int m = shortCount-1; m > shortCount - (veryrareCount / 16); m--){
                     het1 = (~GG[0][id1][m]) & GG[1][id1][m];
                     nonmissing1 = GG[0][id1][m] | GG[1][id1][m];
                     het2 = (~GG[0][id2][m]) & GG[1][id2][m];
                     nonmissing2 = GG[0][id2][m] | GG[1][id2][m];
                     het3 = (~GG[0][id3][m]) & GG[1][id3][m];
                     nonmissing3 = GG[0][id3][m] | GG[1][id3][m];
                     het1 &= nonmissing2 & nonmissing3;
                     het2 &= nonmissing1 & nonmissing3;
                     het3 &= nonmissing1 & nonmissing2;

                     het1Count += oneoneCount[het1];
                     het2Count += oneoneCount[het2];
                     het3Count += oneoneCount[het3];
                     het12Count += oneoneCount[het1&het2];
                     het13Count += oneoneCount[het1&het3];
                     het23Count += oneoneCount[het2&het3];
                     het123Count += oneoneCount[het1&het2&het3];
                  }  // end of marker
                  if(het13Count > 100){ // 2 is the parent
                     H213 = (double)het123Count / het13Count;
                     H21 = (double)het12Count / het1Count;
                     H23 = (double)het23Count / het3Count;
                     if(H213 > 0.75 || H213 < 0.25)
                        fprintf(fp, "%s\t%s\t%s\t%s\t%s\t%s\t%d\t%.3lf\t%.3lf\t%.3lf\n",
                     (const char*)ped[id[f1][i]].famid,
                     (const char*)ped[id[f1][i]].pid,
                     (const char*)ped[id[f2][j]].famid,
                     (const char*)ped[id[f2][j]].pid,
                     (const char*)ped[id[f3][k]].famid,
                     (const char*)ped[id[f3][k]].pid,
                     het13Count, H21, H23, H213);
                  }
                  if(het23Count > 100){ // 1 is the parent
                     H123 = (double)het123Count / het23Count;
                     H12 = (double)het12Count / het2Count;
                     H13 = (double)het13Count / het3Count;
                     if(H123 > 0.75 || H123 < 0.25)
                        fprintf(fp, "%s\t%s\t%s\t%s\t%s\t%s\t%d\t%.3lf\t%.3lf\t%.3lf\n",
                     (const char*)ped[id[f2][j]].famid,
                     (const char*)ped[id[f2][j]].pid,
                     (const char*)ped[id[f1][i]].famid,
                     (const char*)ped[id[f1][i]].pid,
                     (const char*)ped[id[f3][k]].famid,
                     (const char*)ped[id[f3][k]].pid,
                     het23Count, H12, H13, H123);
                  }
               } // end of k
         }  // end of j
      }  // end of i

   fclose(fp);
   printf("                              ends at %s", currentTime());
   printf("Distant relationships are saved in file %s\n",
         (const char*)outfile);
}

void Engine::ComputeBigDataDistant()
{
   printf("Autosome genotypes stored in %d", shortCount);
   printf(" words for each of %d individuals.\n", idCount);

   printf("\nOptions in effect:\n");
   printf("\t--distant\n");
   if(relativedegree>1)
      printf("\t--degree %d\n", relativedegree);
   if(minMAF!=_NAN_)
      printf("\t--minMAF %lf\n", minMAF);
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");

   int id1, id2, id3;
   int HetHetCount, Het1Het3Count, Het2Het3Count, nonMissingCount, HetCount, HetHetHetCount;
   int het1, het2, het3, nonmissing1, nonmissing2, nonmissing3, het, hethet, hethethet;

   String outfile;
   outfile.Copy(prefix);
   outfile.Add(".dis");
   printf("There are %d rare SNPs with MAF < %lf.\n", veryrareCount, minMAF==_NAN_? 0.01: minMAF);
   FILE *fp = fopen(outfile, "wt");
   printf("Inferring distant relatives starts at %s", currentTime());

   fprintf(fp, "FID1\tID1\tFID2\tID2\tKinship\tN_H\tN_HH\tFID3\tID3\tN_HHH\tH_3|1\tH_3|2\tH_3|12\n");
   int m1, m2;
   double kinship, smaller;
   int IBS0Count, het1Count, het2Count, notMissingCount;
   double H31, H32, H312;

   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   IntArray SNPCount(idCount);
   SNPCount.Zero();
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = 0; i < id[f].Length(); i++){
         notMissingCount = 0;
         for(int m = 0; m < shortCount; m++)
            notMissingCount += oneoneCount[GG[0][geno[id[f][i]]][m] | GG[1][geno[id[f][i]]][m]];
         SNPCount[geno[id[f][i]]] = notMissingCount;
         if(notMissingCount < 1000) {
            printf("Fam %s ID %s excluded for having %d non-missing SNPs\n",
               (const char*)ped[id[f][i]].famid, (const char*)ped[id[f][i]].pid, notMissingCount);
            continue;
         }
      }

   int tdegree=3;
   if(relativedegree>1) tdegree = relativedegree;

   int prescanCount[10];
   for(int i = 0; i < 10; i++)
      prescanCount[i] = 128 * (1<<(i*3));
   int short_prescan = tdegree < 10? prescanCount[tdegree-1]: prescanCount[9];
   if(short_prescan > shortCount || (!bigdataFlag)) short_prescan=shortCount;
   double lowerbound = pow(2.0, -(tdegree+1.5));
   double kincutoff1 = 0.125; // 1-degree
   if(short_prescan==shortCount)
      kincutoff1 = lowerbound;
   else if(tdegree == 2)
      kincutoff1 = 0.04419417;
   else if(tdegree == 3)
      kincutoff1 = 0.015625;
   else if(tdegree > 3)
      kincutoff1 = 0;
   double kincutoff2 = sqrt(lowerbound * kincutoff1);

   for(int f1 = 0; f1 < ped.familyCount; f1++)
      for(int i = 0; i < id[f1].Length(); i++){
         id1 = geno[id[f1][i]];
         if(SNPCount[id1]<1000) continue;
         for(int f2 = f1+1; f2 < ped.familyCount; f2++)
            for(int j = 0; j < id[f2].Length(); j++){
               id2 = geno[id[f2][j]];
               if(SNPCount[id2]<1000) continue;
               // Stage 1: many pairs
               HetHetCount = IBS0Count = het1Count = het2Count = 0;
               for(int m = 0; m < short_prescan; m++){
                  HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
                  IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
               }
               if(HetHetCount + 2 * IBS0Count > 500){
                  if(HetHetCount <= 2*IBS0Count) continue;
                  for(int m = 0; m < short_prescan; m++){
                     het1Count += oneoneCount[GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m]]; // Het1Hom2
                     het2Count += oneoneCount[GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]]; // Hom1Het2
                  }
                  smaller = HetHetCount + (het1Count < het2Count? het1Count: het2Count);
                  kinship = 0.5 - ((het1Count+het2Count)*0.25+IBS0Count)/smaller;
                  if(kinship < kincutoff1) continue;
               }else{ // insufficient informative SNPs
                  for(int m = 0; m < short_prescan; m++){
                     het1Count += oneoneCount[GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m]]; // Het1Hom2
                     het2Count += oneoneCount[GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]]; // Hom1Het2
                  }
               }
               m1 = short_prescan;
               // Stage 2: few pairs
               if(tdegree < 4){
                  m2 = m1 + (short_prescan<<3);
                  if(m2 > shortCount) m2 = shortCount;
                  for(int m = m1; m < m2; m++){
                     HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
                     IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
                  }
                  if(HetHetCount <= 2*IBS0Count) continue;
                  for(int m = m1; m < m2; m++){
                     het1Count += oneoneCount[GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m]]; // Het1Hom2
                     het2Count += oneoneCount[GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]]; // Hom1Het2
                  }
                  smaller = HetHetCount + (het1Count < het2Count? het1Count: het2Count);
                  kinship = 0.5 - ((het1Count+het2Count)*0.25+IBS0Count)/smaller;
                  if(kinship < kincutoff2) continue;
                  m1 = m2;
               }
               // Stage 3: exact
               for(int m = m1; m < shortCount; m++){
                  HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
                  IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
               }
               if(HetHetCount <= 2*IBS0Count) continue;
               for(int m = m1; m < shortCount; m++){
                  het1Count += oneoneCount[GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m]]; // Het1Hom2
                  het2Count += oneoneCount[GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]]; // Hom1Het2
               }
               smaller = HetHetCount + (het1Count < het2Count? het1Count: het2Count);
               kinship = 0.5 - ((het1Count+het2Count)*0.25+IBS0Count)/smaller;
               if(kinship < lowerbound) continue;
               if(kinship > 0.1767767) continue; 
               for(int f3 = 0; f3 < ped.familyCount; f3++)
               for(int k = 0; k < id[f3].Length(); k++){
                  id3 = geno[id[f3][k]];
                  if(id3==id1) continue;
                  if(id3==id2) continue;
                  if(SNPCount[id3]<1000) continue;

                  bool firsttimeFlag=true;
                  HetCount = HetHetCount = HetHetHetCount = het1Count
                     = het2Count = Het1Het3Count = Het2Het3Count = 0;
                  for(int m = shortCount-1; m > shortCount - (veryrareCount / 16); m--){
                     het1 = (~GG[0][id1][m]) & GG[1][id1][m];
                     nonmissing1 = GG[0][id1][m] | GG[1][id1][m];
                     het2 = (~GG[0][id2][m]) & GG[1][id2][m];
                     nonmissing2 = GG[0][id2][m] | GG[1][id2][m];
                     het3 = (~GG[0][id3][m]) & GG[1][id3][m];
                     nonmissing3 = GG[0][id3][m] | GG[1][id3][m];

                     het1 &= nonmissing2 & nonmissing3;
                     het2 &= nonmissing1 & nonmissing3;
                     het3 &= nonmissing1 & nonmissing2;
                     het = het1 | het2;
                     hethet = het1 & het2;
                     hethethet = het1 & het2 & het3 & 65535;
                     HetHetCount += oneoneCount[hethet];
                     HetHetHetCount += oneoneCount[hethethet];
                     HetCount += oneoneCount[het];
                     het1Count += oneoneCount[het1];
                     het2Count += oneoneCount[het2];
                     Het1Het3Count += oneoneCount[het1 & het3 & nonmissing2];
                     Het2Het3Count += oneoneCount[het2 & het3 & nonmissing1];
                  }  // end of marker
                  if(HetHetCount < 10 || het1Count==0 || het2Count==0) continue;
                  H312 = (double)HetHetHetCount / HetHetCount;
                  H31 = (double)Het1Het3Count / het1Count;
                  H32 = (double)Het2Het3Count / het2Count;
                  if(H312 < 0.1767767 && H31 < 0.08838835 && H32 < 0.08838835) continue;
                  fprintf(fp, "%s\t%s\t%s\t%s\t%.3lf\t%d\t%d\t%s\t%s\t%d\t%.3lf\t%.3lf\t%.3lf\n",
                     (const char*)ped[id[f1][i]].famid,
                     (const char*)ped[id[f1][i]].pid,
                     (const char*)ped[id[f2][j]].famid,
                     (const char*)ped[id[f2][j]].pid,
                     kinship, HetCount, HetHetCount,
                     (const char*)ped[id[f3][k]].famid,
                     (const char*)ped[id[f3][k]].pid,
                     HetHetHetCount, H31, H32, H312);
                  if(firsttimeFlag){
                     printf("Relatedness between %s:%s and %s:%s (kinship=%.3lf) is confirmed\n",
                        (const char*)ped[id[f1][i]].famid,
                        (const char*)ped[id[f1][i]].pid,
                        (const char*)ped[id[f2][j]].famid,
                        (const char*)ped[id[f2][j]].pid,
                        kinship);
                     firsttimeFlag=false;
                  }
               } // end of k
         }  // end of j
      }  // end of i

   fclose(fp);
   printf("                              ends at %s", currentTime());
   printf("Distant relationships are saved in file %s\n",
         (const char*)outfile);
}


