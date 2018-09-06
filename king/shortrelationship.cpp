//////////////////////////////////////////////////////////////////////
// shortrelationship.cpp
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

#include <math.h>
#include "analysis.h"
#include "Kinship.h"
#include "KinshipX.h"
#include "MathStats.h"
#include "MathSVD.h"
#include "QuickIndex.h"
#include "MathCholesky.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

void Engine::ComputeShortRobustKinshipAdjustPC()
{
   if(diagFlag) countGenotype();
   printf("Autosome genotypes stored in %d", shortCount);
   printf(" words for each of %d individuals.\n", idCount);
   printf("\nOptions in effect:\n");
   printf("\t--kinship\n");
   printf("\t--adjustPC %d\n", BetaSum.Length()-1);
   if(minMAF > 0.05)
      printf("\t--minMAF %lf\n", minMAF);
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");

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
   unsigned int **missingWordInOnePerson = new unsigned int *[idCount];
   int *missingWordInOnePersonCount = new int[idCount];
   int *missingInOnePersonCount = new int[idCount];
   int *hetInOnePersonCount = new int[idCount];
   int m1, m2, m3, m4;
   unsigned long int ulong;
   for(int i = 0; i < idCount; i++)
      missingWordInOnePersonCount[i] = missingInOnePersonCount[i] = hetInOnePersonCount[i] = 0;
   double F1, F2;// inbreeding coefficient
   double Fst;
   double kinshipAdjusted;
   int nCov = BetaSum.Length();
   if(nCov<2) printf("Dimension of beta is %d.\n", nCov);
   Vector *BetaSumPerID = new Vector[idCount];
   Matrix *BetaSquareSumPerID = new Matrix[idCount];
   for(int i = 0; i < idCount; i++){
      BetaSumPerID[i] = BetaSum;
      BetaSquareSumPerID[i] = BetaSquareSum;
   }
   int AACount, AaCount, MissingCount;
   const double RAREFREQ = minMAF > 0.05? minMAF: 0.05;
   int wordCount = 1;
   for(; wordCount < shortCount; wordCount <<= 1){
      AACount = AaCount = MissingCount = 0;
      for(int i = 0; i < idCount; i++){
         AACount += oneoneCount[GG[0][i][wordCount] & GG[1][i][wordCount]];
         AaCount += oneoneCount[(~GG[0][i][wordCount]) & GG[1][i][wordCount]];
         MissingCount += oneoneCount[~(GG[0][i][wordCount] | GG[1][i][wordCount]) & 65535];
      }
      if( (AACount + AaCount * 0.5) / (idCount*16 - MissingCount) < RAREFREQ) break;
   }
   wordCount >>= 1;
   for(; wordCount < shortCount; wordCount++){
      AACount = AaCount = MissingCount = 0;
      for(int i = 0; i < idCount; i++){
         AACount += oneoneCount[GG[0][i][wordCount] & GG[1][i][wordCount]];
         AaCount += oneoneCount[(~GG[0][i][wordCount]) & GG[1][i][wordCount]];
         MissingCount += oneoneCount[~(GG[0][i][wordCount] | GG[1][i][wordCount]) & 65535];
      }
      if((AACount + AaCount * 0.5) / (idCount*16 - MissingCount) < RAREFREQ) break;
   }
   wordCount--;
   if((wordCount == shortCount) && (markerCount%16!=0)) wordCount--;
   int commonCount = wordCount*16;
   printf("%d common variants (MAF>%.2lf) are used for AdjustPC relationship inference.\n",
      commonCount, RAREFREQ);

#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount)
#endif
   for(int i = 0; i < idCount; i++){
      for(int m1 = 0; m1 < wordCount; m1++){
         int m3 = ~(GG[0][i][m1] | GG[1][i][m1]) & 0xFFFF;
         if(m3==0) continue;  // no missingness
         for(int m2 = 0; m2 < 16; m2++)
            if(m3 & shortbase[m2]){
               int m4 = m1*16 + m2; // missing at m4
               for(int u = 0; u < nCov; u++)
                  BetaSumPerID[i][u] -= freqBeta[m4][u];
               for(int u = 0; u < nCov; u++)
                  for(int v = 0; v < nCov; v++)
                     BetaSquareSumPerID[i][u][v] -= freqBeta[m4][u] * freqBeta[m4][v];
         }
      }
      IntArray tArray(0);
      for(int m = 0; m < wordCount; m++)
         if((~(GG[0][i][m] | GG[1][i][m]) & 0xFFFF) != 0){
            missingWordInOnePersonCount[i] ++;
            missingInOnePersonCount[i] += oneoneCount[~(GG[0][i][m] | GG[1][i][m]) & 65535];
            tArray.Push(m);
         }
      if(missingWordInOnePersonCount[i]){
         missingWordInOnePerson[i] = new unsigned int [missingWordInOnePersonCount[i]];
         for(int m = 0; m < tArray.Length(); m++)
            missingWordInOnePerson[i][m] = tArray[m];
      }else
         missingWordInOnePerson[i] = NULL;
      for(int m = 0; m < wordCount; m++)
         hetInOnePersonCount[i] += oneoneCount[(~GG[0][i][m]) & GG[1][i][m]];
   }  // end of idCount loop
   const int cutoffMissingCount=commonCount-MINSNPCOUNT/10; // still needs a few

   Matrix pcX(idCount, nCov);
   for(int i = 0; i < idCount; i++){
      pcX[i][0] = 1.0;
      for(int j = 1; j < nCov; j++)
         pcX[i][j] = ped[phenoid[i]].covariates[covariatePC[j-1]];
   }
   Vector var(idCount);
   var.Zero();
   Matrix VarAtWord(idCount, wordCount);
   VarAtWord.Zero();
   Matrix VarAtSNP(idCount, commonCount);
   VarAtSNP.Zero();
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount)
#endif
   for(int i = 0; i < idCount; i++){
      for(int m = 0; m < commonCount; m++){
         double mu = 0.0;
         for(int j = 0; j < nCov; j++)
            mu += freqBeta[m][j] * pcX[i][j];
         if(mu > 0 && mu < 1)
            VarAtSNP[i][m] = mu * (1-mu) * 2.0;
      }
      for(int w = 0; w < wordCount; w++){
         int mMax = (w+1)*16;
         double sum = 0.0;
         for(int m = w*16; m < mMax; m++)
            sum += VarAtSNP[i][m];
         VarAtWord[i][w] = sum;
         var[i] += sum;
      }
   }

   double myBetaSum[100];
   double myBetaSquareSum[100][100];
   double var1, var2;
   double lambda;

   String outfile;
   outfile.Copy(prefix);
   outfile.Add(".kin");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "FID\tID1\tID2\tN_SNP\tZ0\tPhi\tHetHet\tIBS0\tKING-R\tF1\tF2\tFst\tLambda\tKinship\tError\n");
   double kinship, inflation;
   int id1, id2;
   Kinship kin;
   int HetHomCount, HetHetCount, IBS0Count, het1Count, het2Count, notMissingCount;
   double phi, pi0, errorFlag;
   int beforeCount[6], afterCount[6];
   int degree; double ibs0;
   for(int i = 0; i < 6; i++) beforeCount[i] = afterCount[i] = 0;
   for(int f = 0; f < ped.familyCount; f++){
      kin.Setup(*ped.families[f]);
      for(int i = 0; i < id[f].Length(); i++)
      for(int j = i+1; j < id[f].Length(); j++){
         id1 = geno[id[f][i]]; id2 = geno[id[f][j]];
         if(missingInOnePersonCount[id1]>=cutoffMissingCount) continue;
         if(missingInOnePersonCount[id2]>=cutoffMissingCount) continue;
         for(IBS0Count=0, m1 = 0; m1 < wordCount; m1++)
            IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
         for(HetHomCount = -missingInOnePersonCount[id1]-missingInOnePersonCount[id2],
            m1 = 0; m1 < wordCount; m1++) // Lower bound of HetHom Count, minus HetMiss and MissMiss
               HetHomCount += oneoneCount[GG[0][id1][m1]^GG[0][id2][m1]];
         int m1Max = missingWordInOnePersonCount[id1];
         int m2Max = missingWordInOnePersonCount[id2];
         for(m1 = m2 = notMissingCount = 0; m1 < m1Max; m1++){
            m3 = missingWordInOnePerson[id1][m1]; // m for id1
            for(; m2 < m2Max && missingWordInOnePerson[id2][m2] < m3; m2++);
            if( m2 == m2Max ) break;
            if(m3 == missingWordInOnePerson[id2][m2])   // id1 and id2 are both missing
               notMissingCount += oneoneCount[~(GG[0][id1][m3] | GG[0][id2][m3] | GG[1][id1][m3] | GG[1][id2][m3])&0xFFFF];
         }  // MissMissCount
         HetHomCount += 2 * notMissingCount;
         notMissingCount += commonCount - missingInOnePersonCount[id1] - missingInOnePersonCount[id2];
         for(het1Count = hetInOnePersonCount[id1], m1 = 0;
            m1 < missingWordInOnePersonCount[id2]; m1++){
            m2 = missingWordInOnePerson[id2][m1];
            m3 = oneoneCount[~(GG[0][id1][m2] | GG[0][id2][m2] | GG[1][id2][m2])&GG[1][id1][m2]];
            het1Count -= m3;
            HetHomCount += m3;
         }  // Add HetMiss
         for(het2Count = hetInOnePersonCount[id2], m1 = 0; m1 < m1Max; m1++){
            m2 = missingWordInOnePerson[id1][m1];
            m3 = oneoneCount[~(GG[0][id1][m2] | GG[0][id2][m2] | GG[1][id1][m2]) & GG[1][id2][m2]];
            het2Count -= m3;
            HetHomCount += m3;
         }  // Add MissHet
         m3 = het1Count < het2Count? het1Count: het2Count;
         if(m3 == 0)
            kinship = 0.0;
         else
            kinship = 0.5 - (HetHomCount*0.25+IBS0Count) / m3;
         for(int u = 0; u < nCov; u++)
            myBetaSum[u] = BetaSumPerID[id1][u] + BetaSumPerID[id2][u] - BetaSum[u];
         for(int u = 0; u < nCov; u++)
            for(int v = 0; v < nCov; v++)
               myBetaSquareSum[u][v] = BetaSquareSumPerID[id1][u][v] + BetaSquareSumPerID[id2][u][v] - BetaSquareSum[u][v];
         for(m1 = m2 = 0; m1 < m1Max; m1++){
            m3 = missingWordInOnePerson[id1][m1];
            for(; m2 < m2Max && missingWordInOnePerson[id2][m2] < m3; m2++);
            if( m2 == m2Max ) break;
            m4 = missingWordInOnePerson[id2][m2];
            if(m3 == m4){  // both missing
               ulong = ~(GG[0][id1][m4] | GG[0][id2][m4] | GG[1][id1][m4] | GG[1][id2][m4]) & 0xFFFF;
               if(ulong==0) continue;  // no MissMiss
               ulong = onebit[ulong];  // 1's positions
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
         F1 = 1.0 - het1Count/var1; // inbreeding coefficient for id1
         F2 = 1.0 - het2Count/var2; // inbreeding coefficient for id2
         var1 *= (1+F1); // new variance for individual 1
         var2 *= (1+F2); // new variance for individual 2
         Fst = 0.0;
         for(int u = 0; u < nCov; u++)
            for(int v = 0; v < nCov; v++)
               Fst += myBetaSquareSum[u][v] * (pcX[id1][u]-pcX[id2][u]) * (pcX[id1][v]-pcX[id2][v]);
         kinshipAdjusted = (Fst+(var1+var2-HetHomCount)*0.25-IBS0Count)/sqrt(var1*var2);
         lambda = 1.0;  // by default
         if(kinshipAdjusted > 0){   // possibly related
            lambda = 0;
            for(int w = 0; w < wordCount; w++)
               lambda += sqrt(VarAtWord[id1][w]*VarAtWord[id2][w]);
            lambda = sqrt(var[id1]*var[id2]) / lambda;
            kinshipAdjusted *= lambda;
            if(kinshipAdjusted > 0.03125){   // most likely related
               kinshipAdjusted /= lambda;
               lambda = 0;
               for(int m = 0; m < commonCount; m++)
                  lambda += sqrt(VarAtSNP[id1][m]*VarAtSNP[id2][m]);
               lambda = sqrt(var[id1]*var[id2]) / lambda;
               kinshipAdjusted *= lambda;
            }
         }
         Fst /= sqrt(var1*var2);

         HetHetCount = (het1Count+het2Count-HetHomCount)/2;
            phi = kin(ped[id[f][i]], ped[id[f][j]]);
            pi0 = 0.0;
            if(phi < 0.2)
               pi0 = 1-4*phi;
            else if(phi < 0.3 && ped[id[f][i]].isSib(ped[id[f][j]]))
               pi0 = 0.25;
            errorFlag = 0;
            inflation = (phi > 0)? kinshipAdjusted / phi: -1;
            if(phi > 0.03){  // up to 4th-degree relative
               if(inflation > 2 || inflation < 0.5)
                  errorFlag = 1;
               else if(inflation > 1.4142 || inflation < 0.70711)
                  errorFlag = 0.5;
            }else if(phi < 0.005){  // unrelated pair
               if(kinshipAdjusted > 0.0442) errorFlag = 1;
               else if(kinshipAdjusted > 0.0221) errorFlag = 0.5;
            }else{   // distant relatives
               if(kinshipAdjusted < -0.0221) errorFlag = 1;
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
            if(kinshipAdjusted > 0.0442)
               degree = int(-log(kinshipAdjusted)/log(2.0) - 0.5);
            if(degree < 4){
               afterCount[degree] ++;
               if(degree == 1){
                  if(errorrateCutoff == _NAN_){
                    if(ibs0 < 0.008)
                        afterCount[5] ++;
                  }else{
                     if(ibs0 < errorrateCutoff)
                        afterCount[5]++;
                  }
               }
            }else
               afterCount[4] ++;

            fprintf(fp, "%s\t%s\t%s\t%d\t%.3lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%G\n",
               (const char*)ped[id[f][i]].famid, (const char*)ped[id[f][i]].pid,
               (const char*)ped[id[f][j]].pid, notMissingCount,
               pi0, phi, HetHetCount*1.0/notMissingCount, ibs0, kinship, F1, F2, Fst, lambda, kinshipAdjusted, errorFlag);
         }
   }
   fclose(fp);

   bool pedigreeFlag = false;
   for(int i = 0; i < 6; i++)
      if(beforeCount[i]) pedigreeFlag = true;
   if(pedigreeFlag){
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


   int thread = 0;
   FILE **fps;
   outfile.Copy(prefix);
   outfile.Add(".kin0");

  IntArray ompindex(idCount-1);
#ifdef _OPENMP
   for(int i = 0; i < (idCount-1)/2; i++){
      ompindex[2*i] = i;
      ompindex[2*i+1] = idCount-2-i;
   }
   if(idCount%2==0)
      ompindex[idCount-2] = (idCount-1)/2;
#else
   for(int i = 0; i < idCount-1; i++)
      ompindex[i] = i;
#endif
#ifdef _OPENMP
   printf("%d CPU cores are used.\n", defaultMaxCoreCount);
   fps = new FILE *[defaultMaxCoreCount];
   StringArray outfiles(defaultMaxCoreCount);
   for(int c = 1; c < defaultMaxCoreCount; c++){
      outfiles[c].Copy(prefix);
      outfiles[c] += (c+1);
      outfiles[c].Add("$$$.kin0");
      fps[c] = fopen(outfiles[c], "wt");
   }
#else
   fps = new FILE *[1];
#endif
   fps[0] = fopen(outfile, "wt");
   fprintf(fps[0], "FID1\tID1\tFID2\tID2\tN_SNP\tHetHet\tIBS0\tKING-R\tF1\tF2\tFst\tLambda\tKinship\n");
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
   private(HetHomCount, IBS0Count, het1Count, het2Count, notMissingCount, ulong, \
      id1, id2, kinship, thread, m1, m2, m3, m4, myBetaSquareSum, myBetaSum, \
      F1, F2, Fst, var1, var2, kinshipAdjusted, lambda)
#endif
   for(int i = 0; i < idCount-1; i++){
      id1 = ompindex[i];
      if(missingInOnePersonCount[id1]>=cutoffMissingCount) continue;
      for(id2 = id1 + 1; id2 < idCount; id2++){
         if(missingInOnePersonCount[id2]>=cutoffMissingCount) continue;
         if(ped[phenoid[id1]].famid == ped[phenoid[id2]].famid) continue;
         for(IBS0Count=0, m1 = 0; m1 < wordCount; m1++)
            IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
         for(HetHomCount = -missingInOnePersonCount[id1]-missingInOnePersonCount[id2],
            m1 = 0; m1 < wordCount; m1++) // Lower bound of HetHom Count, minus HetMiss and MissMiss
               HetHomCount += oneoneCount[GG[0][id1][m1]^GG[0][id2][m1]];
         int m1Max = missingWordInOnePersonCount[id1];
         int m2Max = missingWordInOnePersonCount[id2];
         for(m1 = m2 = notMissingCount = 0; m1 < m1Max; m1++){
            m3 = missingWordInOnePerson[id1][m1]; // m for id1
            for(; m2 < m2Max && missingWordInOnePerson[id2][m2] < m3; m2++);
            if( m2 == m2Max ) break;
            if(m3 == missingWordInOnePerson[id2][m2])   // id1 and id2 are both missing
               notMissingCount += oneoneCount[~(GG[0][id1][m3] | GG[0][id2][m3] | GG[1][id1][m3] | GG[1][id2][m3])&0xFFFF];
         }  // MissMissCount
         HetHomCount += 2 * notMissingCount;
         notMissingCount += commonCount - missingInOnePersonCount[id1] - missingInOnePersonCount[id2];
         for(m1 = m3 = 0; m1 < m2Max; m1++){
            m2 = missingWordInOnePerson[id2][m1];
            m3 += oneoneCount[~(GG[0][id1][m2] | GG[0][id2][m2] | GG[1][id2][m2])&GG[1][id1][m2]];
         }  // Add HetMiss
         het1Count = hetInOnePersonCount[id1] - m3;
         HetHomCount += m3;
         for(m1 = m3 = 0; m1 < m1Max; m1++){
            m2 = missingWordInOnePerson[id1][m1];
            m3 += oneoneCount[~(GG[0][id1][m2] | GG[0][id2][m2] | GG[1][id1][m2]) & GG[1][id2][m2]];
         }  // Add MissHet
         het2Count = hetInOnePersonCount[id2] - m3;
         HetHomCount += m3;
         m3 = het1Count < het2Count? het1Count: het2Count;
         if(m3 == 0)
            kinship = 0.0;
         else
            kinship = 0.5 - (HetHomCount*0.25+IBS0Count) / m3;

         for(int u = 0; u < nCov; u++)
            myBetaSum[u] = BetaSumPerID[id1][u] + BetaSumPerID[id2][u] - BetaSum[u];
         for(int u = 0; u < nCov; u++)
            for(int v = 0; v < nCov; v++)
               myBetaSquareSum[u][v] = BetaSquareSumPerID[id1][u][v] + BetaSquareSumPerID[id2][u][v] - BetaSquareSum[u][v];
         for(m1 = m2 = 0; m1 < m1Max; m1++){
            m3 = missingWordInOnePerson[id1][m1];
            for(; m2 < m2Max && missingWordInOnePerson[id2][m2] < m3; m2++);
            if( m2 == m2Max ) break;
            m4 = missingWordInOnePerson[id2][m2];
            if(m3 == m4){  // both missing
               ulong = ~(GG[0][id1][m4] | GG[0][id2][m4] | GG[1][id1][m4] | GG[1][id2][m4]) & 0xFFFF;
               if(ulong==0) continue;  // no MissMiss
               ulong = onebit[ulong];  // 1's positions
               if(ulong==0xFFFFFFFF){  // more than 8 1's
                  ulong = ~(GG[0][id1][m4] | GG[0][id2][m4] | GG[1][id1][m4] | GG[1][id2][m4]) & 0xFFFF;
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
         F1 = 1.0 - het1Count/var1; // inbreeding coefficient for id1
         F2 = 1.0 - het2Count/var2; // inbreeding coefficient for id2
         var1 *= (1+F1); // new variance for individual 1
         var2 *= (1+F2); // new variance for individual 2
         Fst = 0.0;
         for(int u = 0; u < nCov; u++)
            for(int v = 0; v < nCov; v++)
               Fst += myBetaSquareSum[u][v] * (pcX[id1][u]-pcX[id2][u]) * (pcX[id1][v]-pcX[id2][v]);
         kinshipAdjusted = (Fst+(var1+var2-HetHomCount)*0.25-IBS0Count)/sqrt(var1*var2);
         lambda = 1.0;  // by default
         if(kinshipAdjusted > 0){   // possibly related
            lambda = 0.0;
            for(int w = 0; w < wordCount; w++)
               lambda += sqrt(VarAtWord[id1][w]*VarAtWord[id2][w]);
            lambda = sqrt(var[id1]*var[id2]) / lambda;
            kinshipAdjusted *= lambda;
            if(kinshipAdjusted > 0.03125){   // most likely related
               kinshipAdjusted /= lambda;
               lambda = 0.0;
               for(int m = 0; m < commonCount; m++)
                  lambda += sqrt(VarAtSNP[id1][m]*VarAtSNP[id2][m]);
               lambda = sqrt(var[id1]*var[id2]) / lambda;
               kinshipAdjusted *= lambda;
            }
         }
         Fst /= sqrt(var1*var2);
#ifdef _OPENMP
         thread = omp_get_thread_num();
#endif
         fprintf(fps[thread], "%s\t%s\t%s\t%s\t%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\n",
            (const char*)ped[phenoid[id1]].famid, (const char*)ped[phenoid[id1]].pid,
            (const char*)ped[phenoid[id2]].famid, (const char*)ped[phenoid[id2]].pid,
            notMissingCount, (het1Count+het2Count-HetHomCount)*0.5/notMissingCount,
            IBS0Count*1.0/notMissingCount, kinship, F1, F2, Fst, lambda, kinshipAdjusted);
      }
   }
#ifdef _OPENMP
   for(int c = 0; c < defaultMaxCoreCount; c++)
      fclose(fps[c]);
   delete []BetaSumPerID;
   delete []BetaSquareSumPerID;
   for(int i = 0; i < idCount; i++)
      delete []missingWordInOnePerson[i];
   delete []missingWordInOnePersonCount;
   delete []missingInOnePersonCount;
   delete []hetInOnePersonCount;

   fps[0] = fopen(outfile, "at");
   char buffer[1024];
   for(int c = 1; c < defaultMaxCoreCount; c++){
      fps[c] = fopen(outfiles[c], "rt");
      while(fgets(buffer, 1024, fps[c]) != NULL){
         fputs(buffer, fps[0]);
      }
      fclose(fps[c]);
      remove(outfiles[c]);
   }
#endif
   fclose(fps[0]);
   printf("                                         ends at %s", currentTime());
   if(relativedegree){
      printf("Between-family relatives (kinship >= %.4lf) saved in file %s\n",
         pow(2, -relativedegree-1.5), (const char*)outfile);
      if(relativedegree == 1){
         printf("\nNote only duplicates and 1st-degree relative pairs are included here.\n");
         printf("A higher degree of relationship can be also included by specifying ");
         printf("'--degree 2' (for up to 2nd-degree) or '--degree 3' (for up to 3rd-degree).\n\n");
      }
   }else
      printf("Between-family kinship data saved in file %s\n", (const char*)outfile);
}

void Engine::ComputeFilteredRobustKinship()
{
/*
   if(geno.Length()==0) BuildShortBinary();
   if(QCbySNP || QCbySample || QCwipe){
      if(QCbySample) QC_By_Sample();
      if(QCbySNP) QC_By_SNP();
      if(QCwipe) QC_WipeMI();
      return;
   }
   */
   if(diagFlag) countGenotype();
   printf("Autosome genotypes stored in %d", shortCount);
   printf(" words for each of %d individuals.\n", idCount);
/*
   if(individualInfo) {
      OutputIndividualInfo();
      if(!kinFlag) return;
   }
   */
   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];

   IntArray tArray;
   unsigned int **missingWordInOnePerson = new unsigned int *[idCount];
   int *missingWordInOnePersonCount = new int[idCount];
   int *missingInOnePersonCount = new int[idCount];
   int *hetInOnePersonCount = new int[idCount];
   int m1, m2, m3, m4;
   for(int i = 0; i < idCount; i++)
      missingWordInOnePersonCount[i] = missingInOnePersonCount[i] = hetInOnePersonCount[i] = 0;
   int Mask = (1 << (markerCount % 16))-1;
   if(markerCount%16==0) Mask = 0xFFFF;
   for(int i = 0; i < idCount; i++){
      tArray.Dimension(0);
      for(int m = 0; m < shortCount-1; m++)
         if(((~GG[0][i][m]) & (~GG[1][i][m]) & 0xFFFF)!=0){
            missingWordInOnePersonCount[i] ++;
            missingInOnePersonCount[i] += oneoneCount[(~GG[0][i][m]) & (~GG[1][i][m]) & 65535];
            tArray.Push(m);
         }
      int m = shortCount-1;
      if(((~GG[0][i][m]) & (~GG[1][i][m]) & Mask)!=0){
         missingWordInOnePersonCount[i] ++;
         missingInOnePersonCount[i] += oneoneCount[(~GG[0][i][m]) & (~GG[1][i][m]) & Mask];
         tArray.Push(m);
      }
      if(missingWordInOnePersonCount[i]){
         missingWordInOnePerson[i] = new unsigned  int [missingWordInOnePersonCount[i]];
         for(int m = 0; m < tArray.Length(); m++)
            missingWordInOnePerson[i][m] = tArray[m];
      }else
         missingWordInOnePerson[i] = NULL;
      for(int m = 0; m < shortCount; m++)
         hetInOnePersonCount[i] += oneoneCount[(~GG[0][i][m]) & GG[1][i][m]];
   }

   String outfile;
   outfile.Copy(prefix);
   outfile.Add(".kin");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "FID\tID1\tID2\tN_SNP\tZ0\tPhi\tHetHet\tIBS0\tKinship\tError\n");
   double kinship, inflation;
   int id1, id2;
   Kinship kin;
   int HetHetCount, IBS0Count, het1Count, het2Count, notMissingCount;
   double phi, pi0, errorFlag;
   int beforeCount[6], afterCount[6];
   int degree; double ibs0;
   for(int i = 0; i < 6; i++) beforeCount[i] = afterCount[i] = 0;
   for(int f = 0; f < ped.familyCount; f++){
      kin.Setup(*ped.families[f]);
      for(int i = 0; i < id[f].Length(); i++)
         for(int j = i+1; j < id[f].Length(); j++){
            id1 = geno[id[f][i]]; id2 = geno[id[f][j]];
            HetHetCount = IBS0Count = 0;
            for(int m = 0; m < shortCount; m++){
               HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
               IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
            }
            if(HetHetCount==0 && IBS0Count==0) continue;
            het1Count = hetInOnePersonCount[id1] + hetInOnePersonCount[id2];
            for(m1 = 0; m1 < missingWordInOnePersonCount[id1]; m1++){
               m2 = missingWordInOnePerson[id1][m1];
               het1Count -= oneoneCount[(~GG[0][id1][m2])&(~GG[1][id1][m2])&(~GG[0][id2][m2])&GG[1][id2][m2]];
            }
            for(m1 = 0; m1 < missingWordInOnePersonCount[id2]; m1++){
               m2 = missingWordInOnePerson[id2][m1];
               het1Count -= oneoneCount[(~GG[0][id2][m2])&(~GG[1][id2][m2])&(~GG[0][id1][m2])&GG[1][id1][m2]];
            }
            notMissingCount = markerCount - missingInOnePersonCount[id1] - missingInOnePersonCount[id2];
            for(m1 = m2 = 0; m1 < missingWordInOnePersonCount[id1]; m1++){
               m3 = missingWordInOnePerson[id1][m1];
               for(; m2 < missingWordInOnePersonCount[id2] &&
                  missingWordInOnePerson[id2][m2] < m3; m2++);
               if( m2 == missingWordInOnePersonCount[id2] ) break;
               m4 = missingWordInOnePerson[id2][m2];
               if(m3 == m4)
                  notMissingCount += oneoneCount[(~GG[0][id1][m3])&(~GG[1][id1][m3])&(~GG[0][id2][m4])&(~GG[1][id2][m4])&0xFFFF];
            }
            kinship = (HetHetCount - IBS0Count*2.0) / het1Count;

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
               if(degree == 1){
                  if(errorrateCutoff == _NAN_){
                    if(ibs0 < 0.008)
                        afterCount[5] ++;
                  }else{
                     if(ibs0 < errorrateCutoff)
                        afterCount[5]++;
                  }
               }
            }else
               afterCount[4] ++;

            fprintf(fp, "%s\t%s\t%s\t%d\t%.3lf\t%.4lf\t%.3lf\t%.4lf\t%.4lf\t%G\n",
               (const char*)ped[id[f][i]].famid, (const char*)ped[id[f][i]].pid,
               (const char*)ped[id[f][j]].pid, notMissingCount, pi0, phi,
               HetHetCount*1.0/notMissingCount, ibs0, kinship, errorFlag);
         }
   }
   fclose(fp);

   bool pedigreeFlag = false;
   for(int i = 0; i < 6; i++)
      if(beforeCount[i]) pedigreeFlag = true;
   if(pedigreeFlag){
      printf("Within-family kinship data saved in file %s\n", (const char*)outfile);
      printRelationship(beforeCount, afterCount);
   }else
      printf("Each family consists of one individual.\n");

   if(ped.familyCount < 2) {
      if(ped.familyCount==1 && ped.families[0]->famid=="0")
         warning("All individuals with family ID 0 are considered as relatives.\n");
      if(xmarkerCount >= MINSNPCOUNT){
         printf("\nX-chromosome analysis...\n");
         ComputeShortRobustXKinship();
      }
      printf("There is only one family.\n");
      return;
   }

   printf("Relationship inference across families starts at %s", currentTime());
   int thread = 0;
   printf("Only individuals with kinship coefficent >  %.4lf will be printed\n", kinFilter);
   const double cutoff0 = 1-2*kinFilter;
   FILE **fps;
   outfile.Copy(prefix);
   outfile.Add(".kin0");

  IntArray ompindex(idCount-1);
#ifdef _OPENMP
   for(int i = 0; i < (idCount-1)/2; i++){
      ompindex[2*i] = i;
      ompindex[2*i+1] = idCount-2-i;
   }
   if(idCount%2==0)
      ompindex[idCount-2] = (idCount-1)/2;
#else
   for(int i = 0; i < idCount-1; i++)
      ompindex[i] = i;
#endif
#ifdef _OPENMP
   printf("%d CPU cores are used.\n", defaultMaxCoreCount);
   fps = new FILE *[defaultMaxCoreCount];
   StringArray outfiles(defaultMaxCoreCount);
   for(int c = 1; c < defaultMaxCoreCount; c++){
      outfiles[c].Copy(prefix);
      outfiles[c] += (c+1);
      outfiles[c].Add("$$$.kin0");
      fps[c] = fopen(outfiles[c], "wt");
   }
#else
   fps = new FILE *[1];
#endif
   fps[0] = fopen(outfile, "wt");
   fprintf(fps[0], "FID1\tID1\tFID2\tID2\tN_SNP\tHetHet\tIBS0\tKinship\n");
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
   private(HetHetCount, IBS0Count, het1Count, het2Count, notMissingCount, \
      id1, id2, kinship, thread, m1, m2, m3, m4)
#endif
   for(int i = 0; i < idCount-1; i++){
      id1 = ompindex[i];
      if(ped[phenoid[id1]].ngeno==0) continue;
      for(id2 = id1 + 1; id2 < idCount; id2++){
         if(ped[phenoid[id2]].ngeno==0) continue;
         if(ped[phenoid[id1]].famid == ped[phenoid[id2]].famid) continue;
         for(het1Count = hetInOnePersonCount[id1], m1 = 0; m1 < missingWordInOnePersonCount[id2]; m1++){
            m2 = missingWordInOnePerson[id2][m1];
            het1Count -= oneoneCount[(~GG[0][id2][m2])&(~GG[1][id2][m2])&(~GG[0][id1][m2])&GG[1][id1][m2]];
         }
         for(het2Count = hetInOnePersonCount[id2], m1 = 0; m1 < missingWordInOnePersonCount[id1]; m1++){
            m2 = missingWordInOnePerson[id1][m1];
            het2Count -= oneoneCount[(~GG[0][id1][m2])&(~GG[1][id1][m2])&(~GG[0][id2][m2])&GG[1][id2][m2]];
         }
         double cutoff = cutoff0 * 2 * (het1Count < het2Count? het1Count: het2Count);
         int HetHomCount = IBS0Count = 0;
         for(int m = 0; (m < shortCount) && (HetHomCount+IBS0Count*4 <= cutoff); m++){
            HetHomCount += oneoneCount[
               ( (~GG[0][id1][m]) & GG[1][id1][m] & GG[0][id2][m] ) | // Het x Hom
               ( (~GG[0][id2][m]) & GG[1][id2][m] & GG[0][id1][m] ) ];// Hom x Het
            IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
         }
         if(HetHomCount+IBS0Count*4 > cutoff) continue;
         for(m1 = m2 = 0, notMissingCount = markerCount - missingInOnePersonCount[id1] - missingInOnePersonCount[id2];
            m1 < missingWordInOnePersonCount[id1]; m1++){
            m3 = missingWordInOnePerson[id1][m1]; // m for id1
            for(; m2 < missingWordInOnePersonCount[id2] &&
               missingWordInOnePerson[id2][m2] < m3; m2++);
            if( m2 == missingWordInOnePersonCount[id2] ) break;
            m4 = missingWordInOnePerson[id2][m2];  // m for id2
            if(m3 == m4)   // id1 and id2 are both missing
               notMissingCount += oneoneCount[(~GG[0][id1][m3])&(~GG[1][id1][m3])&(~GG[0][id2][m4])&(~GG[1][id2][m4])&0xFFFF];
         }
         kinship = 0.5 - (HetHomCount*0.25 + IBS0Count)/(het1Count < het2Count? het1Count: het2Count);
         HetHetCount = (het1Count + het2Count - HetHomCount)/2;
#ifdef _OPENMP
         thread = omp_get_thread_num();
#endif
         fprintf(fps[thread], "%s\t%s\t%s\t%s\t%d\t%.4lf\t%.4lf\t%.4lf\n",
            (const char*)ped[phenoid[id1]].famid, (const char*)ped[phenoid[id1]].pid,
            (const char*)ped[phenoid[id2]].famid, (const char*)ped[phenoid[id2]].pid,
            notMissingCount, HetHetCount*1.0/notMissingCount,
            IBS0Count*1.0/notMissingCount, kinship);
      }
   }
#ifdef _OPENMP
   for(int c = 0; c < defaultMaxCoreCount; c++)
      fclose(fps[c]);
   for(int i = 0; i < idCount; i++)
      delete []missingWordInOnePerson[i];
   delete []missingWordInOnePersonCount;
   delete []missingInOnePersonCount;
   delete []hetInOnePersonCount;

   fps[0] = fopen(outfile, "at");
   char buffer[1024];
   for(int c = 1; c < defaultMaxCoreCount; c++){
      fps[c] = fopen(outfiles[c], "rt");
      while(fgets(buffer, 1024, fps[c]) != NULL){
         fputs(buffer, fps[0]);
      }
      fclose(fps[c]);
      remove(outfiles[c]);
   }
#endif
   fclose(fps[0]);
   printf("                                         ends at %s", currentTime());
   if(relativedegree){
      printf("Between-family relatives (kinship >= %.4lf) saved in file %s\n",
         pow(2, -relativedegree-1.5), (const char*)outfile);
      if(relativedegree == 1){
         printf("\nNote only duplicates and 1st-degree relative pairs are included here.\n");
         printf("A higher degree of relationship can be also included by specifying ");
         printf("'--degree 2' (for up to 2nd-degree) or '--degree 3' (for up to 3rd-degree).\n\n");
      }
   }else
      printf("Between-family kinship data saved in file %s\n", (const char*)outfile);
}

void Engine::IBDLength_trioMI(IntArray &trio, int order, Vector &lengths)
{
   const double GAPMIN = 0.5;
   const double SEGMIN = 1.0;
   const double SEGOFINTEREST = 20.0;

   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   int HetHetCount, IBS0Count, het1Count, het2Count, notMissingCount;
   int informativeCount;
   int informative, nonmissing;
   int id1, id2, id3;
   double totalLength;
   double runlength, start, stop;
   bool runFlag;

   double *startpos = new double[shortCount];
   double *stoppos = new double[shortCount];
   int *chrpos = new int[shortCount];
   bool *shortseg = new bool[shortCount];
   bool *shortgap = new bool[shortCount];
   for(int m = 0; m < positions.Length(); m+=16){
      if(chromosomes[m]>0 && chromosomes[m] < SEXCHR)
         chrpos[m/16] = chromosomes[m];
      else
         chrpos[m/16] = -1;
      startpos[m/16] = (chromosomes[m]-1)*1000 + positions[m];
   }
   for(int m = 15; m < positions.Length(); m+=16){
      if(chrpos[m/16] != chromosomes[m])
         chrpos[m/16] = -1;
      stoppos[m/16] = (chromosomes[m]-1)*1000 + positions[m];
   }
   stoppos[(positions.Length()-1)/16] = (chromosomes[positions.Length()-1]-1)*1000 + positions[positions.Length()-1];
   for(int m = 0; m < shortCount; m++){
      shortseg[m] = true;
      if(chrpos[m]==-1) shortseg[m] = false;
      if(stoppos[m] - startpos[m] > SEGMIN) shortseg[m] = false;
      if(stoppos[m] - startpos[m] < 0) shortseg[m] = false;
      shortgap[m] = true;
      if(m > 0){
         if(chrpos[m] != chrpos[m-1]) shortgap[m] = false;
         else if(startpos[m] - stoppos[m-1] > GAPMIN) shortgap[m] = false;
         else if(startpos[m] - stoppos[m-1] < 0) shortgap[m] = false;
      }
   }

   for(int i = 0; i < trio.Length()/4; i++){
      int f = trio[i*4];
      id1 = geno[id[f][trio[i*4+1]]];
      id2 = geno[id[f][trio[i*4+2]]];
      id3 = geno[id[f][trio[i*4+3]]];
      totalLength = 0;
      runFlag = false;
      start = 0; stop = -9;
      for(int m = 0; m < shortCount; m++){
         if(!shortseg[m]){  // long segment
            if(runFlag){   // count the previous run length
               runlength = stop - start;
               if(runlength > SEGOFINTEREST)
                  totalLength += runlength;
               runFlag = false;
               start = 0; stop = -9;
            }
            continue;
         }// now focus on short segments only

         if(order==1) // MI: AA1 -> aa2, or AA3 x AA1 -> Aa2
         IBS0Count = oneoneCount[
            (GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m]))
            | ((~GG[0][id2][m]) & GG[1][id2][m] & GG[0][id1][m]
            & GG[0][id3][m] & (~(GG[1][id1][m] ^ GG[1][id3][m])) )];
         else if(order==0)  // MI: AA1 -> aa2
            IBS0Count = oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
         else error("Order undefined in IBDLength_trioMI");
         // informative: at least one person is Aa
         informative = (~GG[0][id1][m] & GG[1][id1][m]) |
            (~GG[0][id2][m] & GG[1][id2][m]) | (~GG[0][id3][m] & GG[1][id3][m]);
         nonmissing = (GG[0][id1][m] | GG[1][id1][m]) &
            (GG[0][id2][m] | GG[1][id2][m]) & (GG[0][id3][m] | GG[1][id3][m]) &
            ~(GG[0][id3][m] & GG[0][id2][m] & (GG[1][id3][m] ^ GG[1][id2][m]));
         informativeCount = oneoneCount[informative & nonmissing];
         if(runFlag)
            if((!shortgap[m]) || ((IBS0Count > 0) && shortgap[m])){
               runlength = stop - start;
               if(runlength > SEGOFINTEREST)
                  totalLength += runlength;
               runFlag = false;
               start = 0; stop = -9;
            }
         if(IBS0Count==0 && informativeCount>1){// share segment
            if((!shortgap[m])||(!runFlag))
               start = startpos[m];
            stop = stoppos[m];
            runFlag = true;
         }
      }
      lengths[i] = totalLength;
   }
   if(startpos) delete startpos;
   if(stoppos) delete stoppos;
   if(chrpos) delete chrpos;
   if(shortseg) delete shortseg;
   if(shortgap) delete shortgap;
}

void Engine::IBDLength_3Pairs(IntArray &trio, int order, Vector &lengths)
{
   const double GAPMIN = 0.5;
   const double SEGMIN = 1.0;
   const double SEGOFINTEREST = 5.0;

   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   int HetHetCount, IBS0Count, het1Count, het2Count, notMissingCount;
   int informativeCount;
   int informative, nonmissing;
   int id1, id2, id3;
   double totalLength;
   double runlength, start, stop;
   bool runFlag;

   double *startpos = new double[shortCount];
   double *stoppos = new double[shortCount];
   int *chrpos = new int[shortCount];
   bool *shortseg = new bool[shortCount];
   bool *shortgap = new bool[shortCount];
   for(int m = 0; m < positions.Length(); m+=16){
      if(chromosomes[m]>0 && chromosomes[m] < SEXCHR)
         chrpos[m/16] = chromosomes[m];
      else
         chrpos[m/16] = -1;
      startpos[m/16] = (chromosomes[m]-1)*1000 + positions[m];
   }
   for(int m = 15; m < positions.Length(); m+=16){
      if(chrpos[m/16] != chromosomes[m])
         chrpos[m/16] = -1;
      stoppos[m/16] = (chromosomes[m]-1)*1000 + positions[m];
   }
   stoppos[(positions.Length()-1)/16] = (chromosomes[positions.Length()-1]-1)*1000 + positions[positions.Length()-1];
   for(int m = 0; m < shortCount; m++){
      shortseg[m] = true;
      if(chrpos[m]==-1) shortseg[m] = false;
      if(stoppos[m] - startpos[m] > SEGMIN) shortseg[m] = false;
      if(stoppos[m] - startpos[m] < 0) shortseg[m] = false;
      shortgap[m] = true;
      if(m > 0){
         if(chrpos[m] != chrpos[m-1]) shortgap[m] = false;
         else if(startpos[m] - stoppos[m-1] > GAPMIN) shortgap[m] = false;
         else if(startpos[m] - stoppos[m-1] < 0) shortgap[m] = false;
      }
   }

   for(int i = 0; i < trio.Length()/4; i++){
      int f = trio[i*4];
      id1 = geno[id[f][trio[i*4+1]]];
      id2 = geno[id[f][trio[i*4+2]]];
      id3 = geno[id[f][trio[i*4+3]]];
      totalLength = 0;
      runFlag = false;
      start = 0; stop = -9;
      for(int m = 0; m < shortCount; m++){
         if(!shortseg[m]){  // long segment
            if(runFlag){   // count the previous run length
               runlength = stop - start;
               if(runlength > SEGOFINTEREST)
                  totalLength += runlength;
               runFlag = false;
               start = 0; stop = -9;
            }
            continue;
         }// now focus on short segments only

         if(order==0) // IBD12 & IBD13 & IBD23
            IBS0Count = oneoneCount[
            (GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m]))
            | (GG[0][id1][m] & GG[0][id3][m] & (GG[1][id1][m] ^ GG[1][id3][m]))
            | (GG[0][id3][m] & GG[0][id2][m] & (GG[1][id3][m] ^ GG[1][id2][m])) ];
         else if(order==1) // IBD12 & IBD13
            IBS0Count = oneoneCount[
            (GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m]))
            | (GG[0][id1][m] & GG[0][id3][m] & (GG[1][id1][m] ^ GG[1][id3][m]))];
         else if(order==2) // IBD12 & IBD23
            IBS0Count = oneoneCount[
            (GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m]))
            | (GG[0][id2][m] & GG[0][id3][m] & (GG[1][id2][m] ^ GG[1][id3][m]))];
         else if(order==3) // IBD13 & IBD23
            IBS0Count = oneoneCount[
            (GG[0][id3][m] & GG[0][id2][m] & (GG[1][id3][m] ^ GG[1][id2][m]))
            | (GG[0][id1][m] & GG[0][id3][m] & (GG[1][id1][m] ^ GG[1][id3][m]))];
         else
            error("IBD order undefined");
         // informative: at least one person is Aa
         informative = (~GG[0][id1][m] & GG[1][id1][m]) |
            (~GG[0][id2][m] & GG[1][id2][m]) | (~GG[0][id3][m] & GG[1][id3][m]);
         nonmissing = (GG[0][id1][m] | GG[1][id1][m]) &
            (GG[0][id2][m] | GG[1][id2][m]) & (GG[0][id3][m] | GG[1][id3][m]);
         informativeCount = oneoneCount[informative & nonmissing];
         if(runFlag)
            if((!shortgap[m]) || ((IBS0Count > 0) && shortgap[m])){
               runlength = stop - start;
               if(runlength > SEGOFINTEREST)
                  totalLength += runlength;
               runFlag = false;
               start = 0; stop = -9;
            }
         if(IBS0Count==0 && informativeCount>1){// share segment
            if((!shortgap[m])||(!runFlag))
               start = startpos[m];
            stop = stoppos[m];
            runFlag = true;
         }
      }
      lengths[i] = totalLength;
   }
   if(startpos) delete startpos;
   if(stoppos) delete stoppos;
   if(chrpos) delete chrpos;
   if(shortseg) delete shortseg;
   if(shortgap) delete shortgap;
}

void Engine::ComputeRunLengthInOneTrio(String &trio_string)
{
   if(diagFlag) countGenotype();
   printf("Autosome genotypes stored in %d", shortCount);
   printf(" words for each of %d individuals.\n", idCount);
   StringArray trioNames;
   trioNames.AddTokens(trio_string, ",");
   IntArray trioID(1);

   trioID[0] = -1;
   for(int f = 0; f < ped.familyCount; f++){
      for(int i = 0; i < id[f].Length(); i++)
         if(ped[id[f][i]].pid==trioNames[0])
            trioID.Push(i);
      for(int i = 0; i < id[f].Length(); i++)
         if(ped[id[f][i]].pid==trioNames[1])
            trioID.Push(i);
      for(int i = 0; i < id[f].Length(); i++)
         if(ped[id[f][i]].pid==trioNames[2])
            trioID.Push(i);
      if(trioID.Length()>1){
         trioID[0] = f;
         break;
      }
   }
   if(trioID.Length() < 4) error("Cannot find %s in a family", (const char*)trio_string);
   Vector lengthMI(1);
   IBDLength_trioMI(trioID, 1, lengthMI);
   Vector lengthPair(1);
   IBDLength_trioMI(trioID, 0, lengthPair);
   printf("The trio includes (%s, %s, %s).\n",
      (const char*)trioNames[0], (const char*)trioNames[1], (const char*)trioNames[2]);
   printf("The length of segments IBD between %s and %s is %.1lf.\n",
      (const char*)trioNames[0], (const char*)trioNames[1], lengthPair[0]);
   printf("The length of segments in which %s and %s are parents of %s is %.1lf.\n",
      (const char*)trioNames[0], (const char*)trioNames[2], (const char*)trioNames[1],
      lengthMI[0]);
}

void Engine::ComputeRunLength()
{
   if(diagFlag) countGenotype();
   printf("Autosome genotypes stored in %d", shortCount);
   printf(" words for each of %d individuals.\n", idCount);
   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   // pairwise relationship inference
   double kinship;
   int id1, id2, id3;
   Kinship kin;
   int HetHetCount, IBS0Count, het1Count;
   double phi, pi0;
   int Ntrio=0;
   IntArray *related1 = new IntArray[ped.familyCount];
   IntArray *related2 = new IntArray[ped.familyCount];
   Vector *allkinship = new Vector[ped.familyCount];

   IntArray trioID(0);
   IntArray trioPair(0);

   if((related1 == NULL) || (related2==NULL))
      error("Cannot allocate memory in function trioIBD");
   for(int f = 0; f < ped.familyCount; f++){
      related1[f].Dimension(0);
      related2[f].Dimension(0);
      allkinship[f].Dimension(0);
      kin.Setup(*ped.families[f]);
      for(int i = 0; i < id[f].Length(); i++)
         for(int j = i+1; j < id[f].Length(); j++){
            id1 = geno[id[f][i]]; id2 = geno[id[f][j]];
            HetHetCount = IBS0Count = het1Count = 0;
            for(int m = 0; m < shortCount; m++){
               HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
               IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
               het1Count += oneoneCount[(GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m]) |
                  (GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m])]; // HomHet
            }
            if(het1Count+HetHetCount==0) continue;
            kinship = (HetHetCount - IBS0Count*2.0) / (HetHetCount*2+het1Count);
            phi = kin(ped[id[f][i]], ped[id[f][j]]);
            if(phi > 0 || kinship > 0.011){ // up to 5th-degree relationship
               related1[f].Push(i);
               related2[f].Push(j);
               allkinship[f].Push(kinship);
            }
         }

      for(int i = 0; i < id[f].Length(); i++){
         for(int j = 0; j < related1[f].Length(); j++){
            id1 = related1[f][j]; id2 = related2[f][j];
            if(id1 <= i) continue; // ordered trio: i < id1 < id2
            int k1;
            for(k1 = 0; k1 < related1[f].Length(); k1++)
               if((related1[f][k1]==i && related2[f][k1]==id1) ||
                  (related1[f][k1]==id1 && related2[f][k1]==i)) break;
            if(k1 == related1[f].Length()) continue;
            int k2;
            for(k2 = 0; k2 < related1[f].Length(); k2++)
               if((related1[f][k2]==i && related2[f][k2]==id2) ||
                  (related1[f][k2]==id2 && related2[f][k2]==i)) break;
            if(k2 == related1[f].Length()) continue;
            trioID.Push(f); trioID.Push(i); trioID.Push(id1); trioID.Push(id2);
            trioPair.Push(f); trioPair.Push(k1); trioPair.Push(k2); trioPair.Push(j);
            Ntrio ++;
         }
      }
   }

   Vector all3pairs(Ntrio);
   Vector join1(Ntrio);
   Vector join2(Ntrio);
   Vector join3(Ntrio);
   IBDLength_3Pairs(trioID, 0, all3pairs);
   IBDLength_3Pairs(trioID, 1, join1);
   IBDLength_3Pairs(trioID, 2, join2);
   IBDLength_3Pairs(trioID, 3, join3);
   String outfile;
   outfile.Copy(prefix);
   outfile.Add(".run");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "FID\tID1\tID2\tID3\tKinship12\tKinship13\tKinship23\tIBDtrio\tIBDjoin1\tIBDjoin2\tIBDjoin3\n");
   int index = 0;
   for(int i = 0; i < trioID.Length()/4; i++, index++){
      int f = trioID[i*4];
      fprintf(fp, "%s\t%s\t%s\t%s\t%.3lf\t%.3lf\t%.3lf\t%.1lf\t%.1lf\t%.1lf\t%.1lf\n",
         (const char*)ped[id[f][trioID[i*4+1]]].famid,
         (const char*)ped[id[f][trioID[i*4+1]]].pid,
         (const char*)ped[id[f][trioID[i*4+2]]].pid,
         (const char*)ped[id[f][trioID[i*4+3]]].pid,
         allkinship[f][trioPair[i*4+1]], allkinship[f][trioPair[i*4+2]],
         allkinship[f][trioPair[i*4+3]], all3pairs[i], join1[i], join2[i], join3[i]);
   }
   fclose(fp);
   printf("Within-family run length data saved in file %s\n", (const char*)outfile);
   if(related1) delete []related1;
   if(related2) delete []related2;
   if(allkinship) delete []allkinship;
}


void Engine::ComputeShortDuplicate()
{
   int id1, id2, notMissingCount;
   int HetHetCount, HetHomCount, HomHetCount, SameHomCount, DiffHomCount;
   double E_IBS;
   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   if(geno.Length()==0) BuildShortBinary();
   printf("Genotypes stored in %d", shortCount);
   if(xshortCount) printf(" + %d X", xshortCount);
   if(yshortCount) printf(" + %d Y", yshortCount);
   if(mtshortCount) printf(" + %d MT", mtshortCount);
   printf(" words for each of %d individuals.\n", idCount);

   double het21, het12, hom21, hom12;
   int nondupCount;
   int rpCount = 0;
   String outfile;
   outfile.Copy(prefix);
   outfile.Add(".con");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "FID1\tID1\tFID2\tID2\tN\tN_IBS0\tN_IBS1\tN_IBS2\tIBS\tConcord\tHetConc");
   fprintf(fp, "\tHet2|1\tHet1|2\tHom2|1\tHom1|2\n");
   printf("Concordance statistics start at %s", currentTime());

   IntArray list(0);
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = 0; i < id[f].Length(); i++)
         list.Push(id[f][i]);
   for(int i = 0; i < list.Length(); i++){
      id1 = geno[list[i]];
      for(int j = i+1; j < list.Length(); j++){
         id2 = geno[list[j]];
         // screen if > 100,000 SNPs
         if(shortCount > 6250){
            nondupCount = 0;
            for(int m = 0; m < shortCount; m+=100)
               nondupCount += (oneoneCount[
                  ((~GG[0][id1][m]) & (GG[1][id1][m]) & GG[0][id2][m]) |
                  (GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m])]
                  - oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]]);
            if(nondupCount > 0)
               continue; // not duplicates, heterozygote concordance rate < 50%
         }

         HetHetCount = HetHomCount = HomHetCount = SameHomCount = DiffHomCount = 0;
         for(int m = 0; m < shortCount; m++){
            HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
            HetHomCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & GG[0][id2][m]];
            HomHetCount += oneoneCount[GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]];
         }
         if(HetHomCount+HomHetCount >= HetHetCount*0.4285714) // 0.3/0.7
              continue; // not duplicates, heterozygote concordance rate <= 70%
         for(int m = 0; m < shortCount; m++){
            DiffHomCount += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
            SameHomCount += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & ~(GG[1][id1][m]^GG[1][id2][m])];
         }
         notMissingCount = HetHetCount + HetHomCount + HomHetCount + SameHomCount + DiffHomCount;
         E_IBS = 1+(HetHetCount+SameHomCount-DiffHomCount)*1.0/notMissingCount;

         //HetHet|het1, HetHet|het2, SameHom|hom1, SameHom|hom2
         het21 = HetHetCount*1.0/(HetHomCount+HetHetCount);
         het12 = HetHetCount*1.0/(HomHetCount+HetHetCount);
         hom21 = SameHomCount*1.0/(HomHetCount+SameHomCount+DiffHomCount);
         hom12 = SameHomCount*1.0/(HetHomCount+SameHomCount+DiffHomCount);

         fprintf(fp, "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.3lf\t%.5lf\t%.5lf",
               (const char*)ped[list[i]].famid, (const char*)ped[list[i]].pid,
               (const char*)ped[list[j]].famid, (const char*)ped[list[j]].pid,
               notMissingCount, DiffHomCount, HetHomCount+HomHetCount,
               HetHetCount+SameHomCount, E_IBS,
               (HetHetCount+SameHomCount)*1.0/notMissingCount,
               HetHetCount * 1.0 /(HetHomCount+HomHetCount+HetHetCount));
         fprintf(fp, "\t%.4lf\t%.4lf\t%.4lf\t%.4lf", het21, het12, hom21, hom12);
         fprintf(fp, "\n");
         rpCount ++;
      }
   }
   fclose(fp);
   printf("                        ends at %s", currentTime());
   if(rpCount)
      printf("%d pairs of duplicates with heterozygote concordance rate>70%% are saved in file %s\n",
         rpCount, (const char*)outfile);
   else
      printf("No duplicates are detected.\n");
}


void Engine::ComputeLongRobustKinship()
{
   printf("Autosome genotypes stored in %d", longCount);
   printf(" words for each of %d individuals.\n", idCount);
   unsigned long long int longword;
   unsigned int **missingWordInOnePerson = new unsigned int *[idCount];
   int *missingWordInOnePersonCount = new int[idCount];
   int *missingInOnePersonCount = new int[idCount];
   int *hetInOnePersonCount = new int[idCount];
   int m1, m2, m3, m4;
   IntArray tArray;
   const unsigned long long int maskbit16 = 0xFFFF;
   const unsigned long long int maskbit32 = 0xFFFFFFFF;
   const unsigned long long int maskbit64 = 0xFFFFFFFFFFFFFFFF;
   // for debugging purpose
   unsigned long long int tmask = maskbit64;
   unsigned long long int Mask = (unsigned long long int)(1 << (markerCount % Bit64))-1;
   if(Bit64==64 && (markerCount%Bit64==0)) Mask = tmask;
   for(int i = 0; i < idCount; i++){
      missingWordInOnePersonCount[i] = missingInOnePersonCount[i] = hetInOnePersonCount[i] = 0;
      tArray.Dimension(0);
      for(int m = 0; m < longCount-1; m++)
         if((~LG[0][i][m]) & (~LG[1][i][m]) & tmask){
            missingWordInOnePersonCount[i] ++;
            for(longword = (~LG[0][i][m]) & (~LG[1][i][m]) & tmask; longword;
               missingInOnePersonCount[i]++) longword &= longword - 1;
            tArray.Push(m);
         }
      int m = longCount-1;
      if((~LG[0][i][m]) & (~LG[1][i][m] & Mask)){
         missingWordInOnePersonCount[i] ++;
         for(longword = (~LG[0][i][m]) & (~LG[1][i][m]) & Mask; longword;
            missingInOnePersonCount[i]++) longword &= longword - 1;
         tArray.Push(m);
      }
      if(missingWordInOnePersonCount[i]){
         missingWordInOnePerson[i] = new unsigned int [missingWordInOnePersonCount[i]];
         for(int m = 0; m < tArray.Length(); m++)
            missingWordInOnePerson[i][m] = tArray[m];
      }else
         missingWordInOnePerson[i] = NULL;
      for(int m = 0; m < longCount; m++)
         for(longword = (~LG[0][i][m]) & LG[1][i][m]; longword;
            hetInOnePersonCount[i]++) longword &= longword - 1;
   }

   String outfile;
   outfile.Copy(prefix);
   outfile.Add(".kin");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "FID\tID1\tID2\tN_SNP\tZ0\tPhi\tHetHet\tIBS0\tKinship\tError\n");
   double kinship, inflation;
   int id1, id2;
   Kinship kin;
   int HetHetCount, IBS0Count, het1Count, het2Count, notMissingCount;
   double phi, pi0, errorFlag;
   int beforeCount[6], afterCount[6];
   int degree; double ibs0;
   for(int i = 0; i < 6; i++) beforeCount[i] = afterCount[i] = 0;
   for(int f = 0; f < ped.familyCount; f++){
      kin.Setup(*ped.families[f]);
      for(int i = 0; i < id[f].Length(); i++)
      for(int j = i+1; j < id[f].Length(); j++){
         id1 = geno[id[f][i]]; id2 = geno[id[f][j]];
         HetHetCount = IBS0Count = 0;
         for(int m = 0; m < longCount; m++){ // most computationally intensive
            for(longword = LG[0][id1][m] & LG[0][id2][m] & (LG[1][id1][m] ^ LG[1][id2][m]);
               longword; IBS0Count++) longword &= (longword-1);
            for(longword = (~LG[0][id1][m]) & (LG[1][id1][m]) & (~LG[0][id2][m]) & LG[1][id2][m];
               longword; HetHetCount++) longword &= (longword-1);
         }
         het1Count = hetInOnePersonCount[id1]+hetInOnePersonCount[id2];
         for(m1 = 0; m1 < missingWordInOnePersonCount[id2]; m1++){
            m2 = missingWordInOnePerson[id2][m1];
            for(longword = (~LG[0][id2][m2])&(~LG[1][id2][m2])&(~LG[0][id1][m2])&LG[1][id1][m2];
               longword; het1Count--) longword &= (longword-1);
         }
         for(m1 = 0; m1 < missingWordInOnePersonCount[id1]; m1++){
            m2 = missingWordInOnePerson[id1][m1];
            for(longword = (~LG[0][id1][m2])&(~LG[1][id1][m2])&(~LG[0][id2][m2])&LG[1][id2][m2];
               longword; het1Count--) longword &= (longword-1);
         }
         notMissingCount = markerCount - missingInOnePersonCount[id1] - missingInOnePersonCount[id2];
         for(m1 = m2 = 0; m1 < missingWordInOnePersonCount[id1]; m1++){
            m3 = missingWordInOnePerson[id1][m1]; // m for id1
            for(; m2 < missingWordInOnePersonCount[id2] &&
               missingWordInOnePerson[id2][m2] < m3; m2++);
            if( m2 == missingWordInOnePersonCount[id2] ) break;
            m4 = missingWordInOnePerson[id2][m2];  // m for id2
            if(m3 == m4)   // id1 and id2 are both missing
               for(longword = (~LG[0][id1][m3])&(~LG[1][id1][m3])&(~LG[0][id2][m4])&(~LG[1][id2][m4]) & tmask;
                  longword; notMissingCount++) longword &= (longword-1);
         }
            kinship = (HetHetCount - IBS0Count*2.0) / het1Count;
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
               if(degree == 1){
                  if(errorrateCutoff == _NAN_){
                    if(ibs0 < 0.008)
                        afterCount[5] ++;
                  }else{
                     if(ibs0 < errorrateCutoff)
                        afterCount[5]++;
                  }
               }
            }else
               afterCount[4] ++;

            fprintf(fp, "%s\t%s\t%s\t%d\t%.3lf\t%.4lf\t%.3lf\t%.4lf\t%.4lf\t%G\n",
               (const char*)ped[id[f][i]].famid, (const char*)ped[id[f][i]].pid,
               (const char*)ped[id[f][j]].pid, notMissingCount, pi0, phi,
               HetHetCount*1.0/notMissingCount, ibs0, kinship, errorFlag);
         }
   }
   fclose(fp);

   bool pedigreeFlag = false;
   for(int i = 0; i < 6; i++)
      if(beforeCount[i]) pedigreeFlag = true;
   if(pedigreeFlag){
      printf("Within-family kinship data saved in file %s\n", (const char*)outfile);
      printRelationship(beforeCount, afterCount);
   }else
      printf("Each family consists of one individual.\n");

   if(ped.familyCount < 2) {
      if(ped.familyCount==1 && ped.families[0]->famid=="0")
         warning("All individuals with family ID 0 are considered as relatives.\n");
      if(xmarkerCount >= MINSNPCOUNT){
         printf("\nX-chromosome analysis...\n");
         ComputeShortRobustXKinship();
      }
      printf("There is only one family.\n");
      return;
   }

   printf("Relationship inference across families starts at %s", currentTime());
   int thread = 0;
   FILE **fps;
   outfile.Copy(prefix);
   outfile.Add(".kin0");

  IntArray ompindex(idCount-1);
#ifdef _OPENMP
   for(int i = 0; i < (idCount-1)/2; i++){
      ompindex[2*i] = i;
      ompindex[2*i+1] = idCount-2-i;
   }
   if(idCount%2==0)
      ompindex[idCount-2] = (idCount-1)/2;
#else
   for(int i = 0; i < idCount-1; i++)
      ompindex[i] = i;
#endif
#ifdef _OPENMP
   printf("%d CPU cores are used.\n", defaultMaxCoreCount);
   fps = new FILE *[defaultMaxCoreCount];
   StringArray outfiles(defaultMaxCoreCount);
   for(int c = 1; c < defaultMaxCoreCount; c++){
      outfiles[c].Copy(prefix);
      outfiles[c] += (c+1);
      outfiles[c].Add("$$$.kin0");
      fps[c] = fopen(outfiles[c], "wt");
   }
#else
   fps = new FILE *[1];
#endif
   fps[0] = fopen(outfile, "wt");
   fprintf(fps[0], "FID1\tID1\tFID2\tID2\tN_SNP\tHetHet\tIBS0\tKinship\n");

#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
   private(HetHetCount, IBS0Count, het1Count, het2Count, notMissingCount, longword, \
      id1, id2, kinship, thread, m1, m2, m3, m4)
#endif
   for(int i = 0; i < idCount-1; i++){
      id1 = ompindex[i];
      if(ped[phenoid[id1]].ngeno==0) continue;
      for(id2 = id1 + 1; id2 < idCount; id2++){
         if(ped[phenoid[id2]].ngeno==0) continue;
         if(ped[phenoid[id1]].famid == ped[phenoid[id2]].famid) continue;
         HetHetCount = IBS0Count = 0;
         for(int m = 0; m < longCount; m++){ // most computationally intensive
            for(longword = LG[0][id1][m] & LG[0][id2][m] & (LG[1][id1][m] ^ LG[1][id2][m]);
               longword; IBS0Count++) longword &= (longword-1);
            for(longword = (~LG[0][id1][m]) & (LG[1][id1][m]) & (~LG[0][id2][m]) & LG[1][id2][m];
               longword; HetHetCount++) longword &= (longword-1);
         }
         for(het1Count = hetInOnePersonCount[id1], m1 = 0; m1 < missingWordInOnePersonCount[id2]; m1++){
            m2 = missingWordInOnePerson[id2][m1];
            for(longword = (~LG[0][id2][m2])&(~LG[1][id2][m2])&(~LG[0][id1][m2])&LG[1][id1][m2];
               longword; het1Count--) longword &= (longword-1);
         }
         for(het2Count = hetInOnePersonCount[id2], m1 = 0; m1 < missingWordInOnePersonCount[id1]; m1++){
            m2 = missingWordInOnePerson[id1][m1];
            for(longword = (~LG[0][id1][m2])&(~LG[1][id1][m2])&(~LG[0][id2][m2])&LG[1][id2][m2];
               longword; het2Count--) longword &= (longword-1);
         }
         for(m1 = m2 = 0, notMissingCount = markerCount - missingInOnePersonCount[id1] - missingInOnePersonCount[id2];
            m1 < missingWordInOnePersonCount[id1]; m1++){
            m3 = missingWordInOnePerson[id1][m1]; // m for id1
            for(; m2 < missingWordInOnePersonCount[id2] &&
               missingWordInOnePerson[id2][m2] < m3; m2++);
            if( m2 == missingWordInOnePersonCount[id2] ) break;
            m4 = missingWordInOnePerson[id2][m2];  // m for id2
            if(m3 == m4)  // id1 and id2 are both missing
               for(longword = (~LG[0][id1][m3])&(~LG[1][id1][m3])&(~LG[0][id2][m4])&(~LG[1][id2][m4])&tmask;
                  longword; notMissingCount++) longword &= (longword-1);
         }
         m3 = (het1Count < het2Count? het1Count: het2Count);
         if(m3 == 0)
            kinship = 0.0;
         else
            kinship = 0.5 - ((het1Count+het2Count-HetHetCount*2)*0.25+IBS0Count) / m3;
#ifdef _OPENMP
         thread = omp_get_thread_num();
#endif
         fprintf(fps[thread], "%s\t%s\t%s\t%s\t%d\t%.4lf\t%.4lf\t%.4lf\n",
            (const char*)ped[phenoid[id1]].famid, (const char*)ped[phenoid[id1]].pid,
            (const char*)ped[phenoid[id2]].famid, (const char*)ped[phenoid[id2]].pid,
            notMissingCount, HetHetCount*1.0/notMissingCount,
            IBS0Count*1.0/notMissingCount, kinship);
      }
   }
#ifdef _OPENMP
   for(int c = 0; c < defaultMaxCoreCount; c++)
      fclose(fps[c]);
   for(int i = 0; i < idCount; i++)
      delete []missingWordInOnePerson[i];
   delete []missingWordInOnePersonCount;
   delete []missingInOnePersonCount;
   delete []hetInOnePersonCount;

   fps[0] = fopen(outfile, "at");
   char buffer[1024];
   for(int c = 1; c < defaultMaxCoreCount; c++){
      fps[c] = fopen(outfiles[c], "rt");
      while(fgets(buffer, 1024, fps[c]) != NULL){
         fputs(buffer, fps[0]);
      }
      fclose(fps[c]);
      remove(outfiles[c]);
   }
#endif
   fclose(fps[0]);
   printf("                                         ends at %s", currentTime());
   if(relativedegree){
      printf("Between-family relatives (kinship >= %.4lf) saved in file %s\n",
         pow(2, -relativedegree-1.5), (const char*)outfile);
      if(relativedegree == 1){
         printf("\nNote only duplicates and 1st-degree relative pairs are included here.\n");
         printf("A higher degree of relationship can be also included by specifying ");
         printf("'--degree 2' (for up to 2nd-degree) or '--degree 3' (for up to 3rd-degree).\n\n");
      }
   }else
      printf("Between-family kinship data saved in file %s\n", (const char*)outfile);
//   if(xmarkerCount >= MINSNPCOUNT && relativedegree==0){
//      printf("\nX-chromosome analysis...\n");
//      ComputeShortRobustXKinship();
//   }
}



void Engine::ComputeShortFastXHomoKinship()
{
   if(xshortCount == 0) {
      printf("There are no X-chromosome SNP data.\n");
      return;
   }
   printf("X-chromosome genotypes stored in %d words for each of %d individuals.\n",
      xshortCount, idCount);

   Vector xfrequencies;
   ComputeAlleleFrequency(xfrequencies);

   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];

   double *pq = new double[xshortCount*16];
   double *p2q2 = new double[xshortCount*16];
   for(int i = 0; i < xshortCount*16; i++)
      pq[i] = p2q2[i] = 0;
   double pqTotal = 0;
   double p2q2Total = 0;
   for(int m = 0; m < xmarkerCount; m++){
      double freq = xfrequencies[m];
      pq[m] = freq * (1-freq);
      p2q2[m] = pq[m]*pq[m];
      pqTotal += pq[m];
      p2q2Total += p2q2[m];
   }
   int id1, id2, pos;
   unsigned short int notMissing, missing;
   double *pqTotals = new double[idCount];
   double *p2q2Totals = new double[idCount];
   for(int i = 0; i < idCount; i++)
      pqTotals[i] = p2q2Totals[i] = 0;

   for(int f = 0; f < ped.familyCount; f++)
      for(int i = 0; i < id[f].Length(); i++){
         id1 = geno[id[f][i]];
         for(int m = 0; m < xshortCount; m++){
            missing = (~XG[0][id1][m]) & (~XG[1][id1][m]);
            for(int k = 0; k < 16; k++)
               if(missing & shortbase[k]){
                  pos = m*16+k;
                  pqTotals[id1] += pq[pos];
                  p2q2Totals[id1] += p2q2[pos];
               }
         }
         pqTotals[id1] = pqTotal - pqTotals[id1];
         p2q2Totals[id1] = p2q2Total - p2q2Totals[id1];
      }
   String outfile;
   outfile.Copy(prefix);
   outfile.Add("X.kin");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "FID\tID1\tID2\tN_SNP\tPhiX\tPhi");
//   if(ibdFlag) fprintf(fp, "\tIBD0X\tIBD1X\tIBD2X");
   fprintf(fp, "\tKinshipX\n");
   double kinship, inflation;

   Kinship kin;
   KinshipX kinx;
   int HetHetCount, IBS0Count, het1Count, notMissingCount;
   double phi, phix, errorFlag;
   double pqSum, p2q2Sum, ibd0, ibd1, ibd2;
   int rpCount = 0;
   for(int f = 0; f < ped.familyCount; f++){
      kin.Setup(*ped.families[f]);
      kinx.Setup(*ped.families[f]);
      for(int i = 0; i < id[f].Length(); i++)
         for(int j = i+1; j < id[f].Length(); j++){
            id1 = geno[id[f][i]]; id2 = geno[id[f][j]];
            IBS0Count = notMissingCount = het1Count = 0;
            pqSum = p2q2Sum = 0;
            for(int m = 0; m < xshortCount; m++){
               notMissing = (XG[0][id1][m] | XG[1][id1][m]) & (XG[0][id2][m] | XG[1][id2][m]);
               IBS0Count += oneoneCount[XG[0][id1][m] & XG[0][id2][m] & (XG[1][id1][m] ^ XG[1][id2][m])];
               notMissingCount += oneoneCount[notMissing];
               het1Count += oneoneCount[(XG[0][id2][m] & (~XG[0][id1][m]) & XG[1][id1][m]) |
                  (XG[0][id1][m] & (~XG[0][id2][m]) & XG[1][id2][m]) ];
               missing = ~notMissing & (XG[0][id1][m] | XG[1][id1][m] | XG[0][id2][m] | XG[1][id2][m]);
               if(missing)
                  for(int k = 0; k < 16; k++)
                     if(missing & shortbase[k]){
                        pqSum += pq[m*16+k];
                        p2q2Sum += p2q2[m*16+k];
                     }
            }
            pqSum = (pqTotals[id1] + pqTotals[id2] - pqSum) * 0.5;
            p2q2Sum = (p2q2Totals[id1] + p2q2Totals[id2] - p2q2Sum) * 0.5;
            if(ped[id[f][i]].sex==2 && ped[id[f][j]].sex==2){
               kinship = 0.5 - (het1Count+4*IBS0Count)*0.125/pqSum;
               // kinship = 0.5 - het1Count*0.125/pqSum + (HetHetCount - IBS0Count*2.0)*0.25/pqSum;
            }else if(ped[id[f][i]].sex==1 && ped[id[f][j]].sex==1){
               kinship = 1 - (het1Count+4*IBS0Count)*0.125/pqSum;
               //kinship = 1 - het1Count*0.125/pqSum + (HetHetCount - IBS0Count*2.0)*0.25/pqSum;
            }else
               kinship = 0.75 - (het1Count+4*IBS0Count)*0.125/pqSum;
               //kinship = 0.75 - het1Count*0.125/pqSum + (HetHetCount - IBS0Count*2.0)*0.25/pqSum;
            phi = kin(ped[id[f][i]], ped[id[f][j]]);
            phix = kinx(ped[id[f][i]], ped[id[f][j]]);
            fprintf(fp, "%s\t%s\t%s\t%d\t%.3lf\t%.4lf",
               (const char*)ped[id[f][i]].famid, (const char*)ped[id[f][i]].pid,
               (const char*)ped[id[f][j]].pid, notMissingCount, phix, phi);
/*            if(ibdFlag){
               ibd0=ibd1=ibd2=0;
               if(ped[id[f][i]].sex==2 && ped[id[f][j]].sex==2){  // F,F
                  ibd0 = IBS0Count * 0.5 / p2q2Sum;
                  ibd1 = 2-2*ibd0-4*kinship;
                  ibd2 = 4*kinship+ibd0-1;
               }else if(ped[id[f][i]].sex==1 && ped[id[f][j]].sex==1){//M,M
                  ibd2 = kinship;
                  ibd0 = 1-kinship;
               }else if(ped[id[f][i]].sex!=0 && ped[id[f][j]].sex!=0){ // M,F
                  ibd1=2*kinship;
                  ibd0 = 1-ibd1;
               }
               fprintf(fp, "\t%.4lf\t%.4lf\t%.4lf", ibd0, ibd1, ibd2);
            }
*/
            fprintf(fp, "\t%.4lf", kinship);
            fprintf(fp, "\n");
            rpCount ++;
         }
   }
   fclose(fp);
   if(rpCount)
      printf("Within-family kinship data saved in file %s\n", (const char*)outfile);
   else
      printf("Each family consists of one individual.\n");

   if(ped.familyCount < 2) {
      if(ped.familyCount==1 && ped.families[0]->famid=="0")
         warning("All individuals with family ID 0 are considered as relatives.\n");
      return;
   }

   printf("Relationship inference across families starts at %s", currentTime());

   outfile.Copy(prefix);
   outfile.Add("X.kin0");
   fp = fopen(outfile, "wt");
   fprintf(fp, "FID1\tID1\tFID2\tID2\tN_SNP");
//   if(ibdFlag) fprintf(fp, "\tIBD0X\tIBD1X\tIBD2X");
   fprintf(fp, "\tKinshipX\n");

   for(int f1 = 0; f1 < ped.familyCount; f1++)
   for(int i = 0; i < id[f1].Length(); i++){
      id1 = geno[id[f1][i]];
      for(int f2 = f1+1; f2 < ped.familyCount; f2++)
      for(int j = 0; j < id[f2].Length(); j++){
            id2 = geno[id[f2][j]];
            IBS0Count = notMissingCount = het1Count = 0;
            pqSum = p2q2Sum = 0;
            for(int m = 0; m < xshortCount; m++){
               notMissing = (XG[0][id1][m] | XG[1][id1][m]) & (XG[0][id2][m] | XG[1][id2][m]);
               IBS0Count += oneoneCount[XG[0][id1][m] & XG[0][id2][m] & (XG[1][id1][m] ^ XG[1][id2][m])];
               notMissingCount += oneoneCount[notMissing];
               het1Count += oneoneCount[(XG[0][id2][m] & (~XG[0][id1][m]) & XG[1][id1][m]) |
                  (XG[0][id1][m] & (~XG[0][id2][m]) & XG[1][id2][m]) ];
               missing = ~notMissing & (XG[0][id1][m] | XG[1][id1][m] | XG[0][id2][m] | XG[1][id2][m]);
               if(missing)
                  for(int k = 0; k < 16; k++)
                     if(missing & shortbase[k]){
                        pqSum += pq[m*16+k];
                        p2q2Sum += p2q2[m*16+k];
                     }
            }
            pqSum = (pqTotals[id1] + pqTotals[id2] - pqSum) * 0.5;
            p2q2Sum = (p2q2Totals[id1] + p2q2Totals[id2] - p2q2Sum) * 0.5;
            if(ped[id[f1][i]].sex==2 && ped[id[f2][j]].sex==2){
               kinship = 0.5 - (het1Count+4*IBS0Count)*0.125/pqSum;
               //kinship = 0.5 - het1Count*0.125/pqSum + (HetHetCount - IBS0Count*2.0)*0.25/pqSum;
            }else if(ped[id[f1][i]].sex==1 && ped[id[f2][j]].sex==1){
               kinship = 1.0 - (het1Count+4*IBS0Count)*0.125/pqSum;
               //kinship = 1 - het1Count*0.125/pqSum + (HetHetCount - IBS0Count*2.0)*0.25/pqSum;
            }else
               kinship = 0.75 - (het1Count+4*IBS0Count)*0.125/pqSum;
               //kinship = 0.75 - het1Count*0.125/pqSum + (HetHetCount - IBS0Count*2.0)*0.25/pqSum;
            fprintf(fp, "%s\t%s\t%s\t%s\t%d",
               (const char*)ped[id[f1][i]].famid, (const char*)ped[id[f1][i]].pid,
               (const char*)ped[id[f2][j]].famid, (const char*)ped[id[f2][j]].pid,
               notMissingCount);
/*            if(ibdFlag){
               ibd0=ibd1=ibd2=0;
               if(ped[id[f1][i]].sex==2 && ped[id[f2][j]].sex==2){  // F,F
                  ibd0 = IBS0Count * 0.5 / p2q2Sum;
                  ibd1 = 2-2*ibd0-4*kinship;
                  ibd2 = 4*kinship+ibd0-1;
               }else if(ped[id[f1][i]].sex==1 && ped[id[f2][j]].sex==1){//M,M
                  ibd2 = kinship;
                  ibd0 = 1-kinship;
               }else if(ped[id[f1][i]].sex!=0 && ped[id[f2][j]].sex!=0){ // M,F
                  ibd1=2*kinship;
                  ibd0 = 1-ibd1;
               }
               fprintf(fp, "\t%.4lf\t%.4lf\t%.4lf", ibd0, ibd1, ibd2);
            }
*/
            fprintf(fp, "\t%.4lf", kinship);
            fprintf(fp, "\n");
         }
   }
   fclose(fp);
   delete []pq;
   delete []p2q2;
   delete []pqTotals;
   delete []p2q2Totals;
   printf("                                         ends at %s", currentTime());
   printf("Between-family kinship data saved in file %s\n", (const char*)outfile);
}





void Engine::ComputeShortSimpleIBS()
{
   Kinship kin;
   int id1, id2, IBS0Count, IBS2Count, notMissingCount;
   int xIBS0Count, xIBS2Count, xnotMissingCount;
   double phi, pi0, E_IBS, xE_IBS;
   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++){
      oneoneCount[i] = 0;
      for(int j = 0; j < 16; j++)
         if(i & shortbase[j]) oneoneCount[i]++;
   }
   if(geno.Length()==0) BuildShortBinary();
   printf("Genotypes stored in %d", shortCount);
   if(xshortCount) printf(" + %d X", xshortCount);
   if(yshortCount) printf(" + %d Y", yshortCount);
   if(mtshortCount) printf(" + %d MT", mtshortCount);
   printf(" words for each of %d individuals.\n", idCount);

   int rpCount = 0;
   String outfile;
   outfile.Copy(prefix);
   outfile.Add(".ibs");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "FID\tID1\tID2\tZ0\tPhi\tN\tN_IBS0\tN_IBS1\tN_IBS2\tIBS");
   if(xshortCount)
      fprintf(fp, "\tNX\tN_XIBS0\tN_XIBS1\tN_XIBS2\tXIBS");
   if(yshortCount)
      fprintf(fp, "\tNY\tN_YIBS0\tYIBS");
   if(mtshortCount)
      fprintf(fp, "\tNMT\tN_MTIBS0\tMTIBS");
   fprintf(fp, "\n");
   for(int f = 0; f < ped.familyCount; f++){
      kin.Setup(*ped.families[f]);
      for(int i = 0; i < id[f].Length(); i++)
         for(int j = i+1; j < id[f].Length(); j++){
            id1 = geno[id[f][i]]; id2 = geno[id[f][j]];
            IBS0Count = IBS2Count = notMissingCount = 0;
            for(int m = 0; m < shortCount; m++){
               notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
               IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
               IBS2Count += oneoneCount[~(GG[0][id1][m]^GG[0][id2][m]) & ~(GG[1][id1][m]^GG[1][id2][m]) & (GG[0][id1][m] | GG[1][id1][m])];
            }
            E_IBS = notMissingCount? 1+(IBS2Count-IBS0Count)*1.0/notMissingCount:0;
            phi = kin(ped[id[f][i]], ped[id[f][j]]);
            pi0 = 0.0;
            if(phi < 0.2)
               pi0 = 1-4*phi;
            else if(phi < 0.3 && ped[id[f][i]].isSib(ped[id[f][j]]))
               pi0 = 0.25;

            xIBS0Count = xIBS2Count = xnotMissingCount = 0;
            for(int m = 0; m < xshortCount; m++){
               xIBS0Count += oneoneCount[XG[0][id1][m] & XG[0][id2][m] & (XG[1][id1][m] ^ XG[1][id2][m])];
               xIBS2Count += oneoneCount[~(XG[0][id1][m]^XG[0][id2][m]) & ~(XG[1][id1][m]^XG[1][id2][m]) & (XG[0][id1][m] | XG[1][id1][m])];
               xnotMissingCount += oneoneCount[(XG[0][id1][m] | XG[1][id1][m]) & (XG[0][id2][m] | XG[1][id2][m])];
            }
            xE_IBS = xnotMissingCount? 1+(xIBS2Count-xIBS0Count)*1.0/xnotMissingCount: 0;

            if((!notMissingCount) && (!xnotMissingCount)) continue;
            fprintf(fp, "%s\t%s\t%s\t%.3lf\t%.4lf\t%d\t%d\t%d\t%d\t%.3lf",
               (const char*)ped[id[f][i]].famid, (const char*)ped[id[f][i]].pid,
               (const char*)ped[id[f][j]].pid,
               pi0, phi, notMissingCount,
               IBS0Count, notMissingCount - IBS0Count - IBS2Count, IBS2Count,
               E_IBS);
            if(xshortCount)
               fprintf(fp, "\t%d\t%d\t%d\t%d\t%.3lf",
               xnotMissingCount, xIBS0Count, xnotMissingCount - xIBS0Count - xIBS2Count,
               xIBS2Count, xE_IBS);

            xIBS0Count = xIBS2Count = xnotMissingCount = 0;
            for(int m = 0; m < yshortCount; m++){
               xIBS0Count += oneoneCount[YG[0][id1][m] & YG[0][id2][m] & (YG[1][id1][m] ^ YG[1][id2][m])];
               xIBS2Count += oneoneCount[~(YG[0][id1][m]^YG[0][id2][m]) & ~(YG[1][id1][m]^YG[1][id2][m]) & (YG[0][id1][m] | YG[1][id1][m])];
               xnotMissingCount += oneoneCount[(YG[0][id1][m] | YG[1][id1][m]) & (YG[0][id2][m] | YG[1][id2][m])];
            }
            xE_IBS = xnotMissingCount? 1+(xIBS2Count-xIBS0Count)*1.0/xnotMissingCount: 0;
            if(yshortCount)
               fprintf(fp, "\t%d\t%d\t%.3lf",
               xnotMissingCount, xIBS0Count, xE_IBS);

            xIBS0Count = xIBS2Count = xnotMissingCount = 0;
            for(int m = 0; m < mtshortCount; m++){
               xIBS0Count += oneoneCount[MG[0][id1][m] & MG[0][id2][m] & (MG[1][id1][m] ^ MG[1][id2][m])];
               xIBS2Count += oneoneCount[~(MG[0][id1][m]^MG[0][id2][m]) & ~(MG[1][id1][m]^MG[1][id2][m]) & (MG[0][id1][m] | MG[1][id1][m])];
               xnotMissingCount += oneoneCount[(MG[0][id1][m] | MG[1][id1][m]) & (MG[0][id2][m] | MG[1][id2][m])];
            }
            xE_IBS = xnotMissingCount? 1+(xIBS2Count-xIBS0Count)*1.0/xnotMissingCount: 0;
            if(mtshortCount)
               fprintf(fp, "\t%d\t%d\t%.3lf",
               xnotMissingCount, xIBS0Count, xE_IBS);

            fprintf(fp, "\n");
            rpCount ++;
         }
   }
   fclose(fp);
   if(rpCount)
      printf("Within-family IBS data saved in file %s\n", (const char*)outfile);
   else
      printf("Each family consists of one individual.\n");

   if(ped.familyCount < 2) {
      if(ped.familyCount==1 && ped.families[0]->famid=="0")
         warning("All individuals with family ID 0 are considered as relatives.\n");
      printf("There is only one family.\n");
      return;
   }

   printf("IBS statistics across families starts at %s", currentTime());

   int thread = 0;
   FILE **fps;
   outfile.Copy(prefix);
   outfile.Add(".ibs0");
  IntArray ompindex(idCount-1);
#ifdef _OPENMP
   for(int i = 0; i < (idCount-1)/2; i++){
      ompindex[2*i] = i;
      ompindex[2*i+1] = idCount-2-i;
   }
   if(idCount%2==0)
      ompindex[idCount-2] = (idCount-1)/2;
#else
   for(int i = 0; i < idCount-1; i++)
      ompindex[i] = i;
#endif
#ifdef _OPENMP
   printf("%d CPU cores are used.\n", defaultMaxCoreCount);
   fps = new FILE *[defaultMaxCoreCount];
   StringArray outfiles(defaultMaxCoreCount);
   for(int c = 1; c < defaultMaxCoreCount; c++){
      outfiles[c].Copy(prefix);
      outfiles[c] += (c+1);
      outfiles[c].Add("$$$.ibs0");
      fps[c] = fopen(outfiles[c], "wt");
   }
#else
   fps = new FILE *[1];
#endif
   fps[0] = fopen(outfile, "wt");
   fprintf(fps[0], "FID1\tID1\tFID2\tID2\tN\tN_IBS0\tN_IBS1\tN_IBS2\tIBS\n");
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
   private(IBS0Count, IBS2Count, notMissingCount, id1, id2, thread)
#endif
   for(int i = 0; i < idCount-1; i++){
      id1 = ompindex[i];
      if(ped[phenoid[id1]].ngeno==0) continue;
      for(id2 = id1 + 1; id2 < idCount; id2++){
         if(ped[phenoid[id2]].ngeno==0) continue;
         if(ped[phenoid[id1]].famid == ped[phenoid[id2]].famid) continue;
         IBS0Count = IBS2Count = notMissingCount = 0;
         for(int m = 0; m < shortCount; m++){
            IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
            IBS2Count += oneoneCount[~(GG[0][id1][m]^GG[0][id2][m]) & ~(GG[1][id1][m]^GG[1][id2][m]) & (GG[0][id1][m] | GG[1][id1][m])];
            notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
         }
#ifdef _OPENMP
         thread = omp_get_thread_num();
#endif
         fprintf(fps[thread], "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.3lf\n",
            (const char*)ped[phenoid[id1]].famid, (const char*)ped[phenoid[id1]].pid,
            (const char*)ped[phenoid[id2]].famid, (const char*)ped[phenoid[id2]].pid,
            notMissingCount, IBS0Count,
            notMissingCount - IBS0Count - IBS2Count, IBS2Count,
            1+(IBS2Count-IBS0Count)*1.0/notMissingCount);
      }
   }
#ifdef _OPENMP
   for(int c = 0; c < defaultMaxCoreCount; c++)
      fclose(fps[c]);
   fps[0] = fopen(outfile, "at");
   char buffer[1024];
   for(int c = 1; c < defaultMaxCoreCount; c++){
      fps[c] = fopen(outfiles[c], "rt");
      while(fgets(buffer, 1024, fps[c]) != NULL){
         fputs(buffer, fps[0]);
      }
      fclose(fps[c]);
      remove(outfiles[c]);
   }
#endif
   fclose(fps[0]);
   printf("                                 ends at %s", currentTime());
   printf("Between-family IBS data saved in file %s\n", (const char*)outfile);
}

void Engine::ComputeShortSibIBS()
{
   if(geno.Length()==0) BuildShortBinary();
   printf("Genotypes stored in %d words for each of %d individuals.\n",
      shortCount, idCount);
   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++){
      oneoneCount[i] = 0;
      for(int j = 0; j < 16; j++)
         if(i & shortbase[j]) oneoneCount[i]++;
   }

   IntArray **sibship;
   IntArray sibCount(ped.familyCount);
   sibship = new IntArray * [ped.familyCount];
   IntArray sib;
   int HetHetCount, IBS0Count, het1Count, het2Count, notMissingCount;
   for(int f = 0; f < ped.familyCount; f++){
      sib.Dimension(0);
      for(int i = 0; i < id[f].Length(); i++)
         if(ped[id[f][i]].sibCount > 1)
            if(sib.Find(ped[id[f][i]].sibs[0]->serial)==-1){
               int count = 0;
               for(int j = 0; j < ped[id[f][i]].sibCount; j++)
                  if(id[f].Find(ped[id[f][i]].sibs[j]->serial)>-1) count++;
               if(count > 1)
                  sib.Push(ped[id[f][i]].sibs[0]->serial);
            }
      sibCount[f] = sib.Length();
      sibship[f] = new IntArray [sibCount[f]];
      for(int i = 0; i < sibCount[f]; i++){
         sibship[f][i].Dimension(0);
         for(int j = 0; j < ped[sib[i]].sibCount; j++)
            if(id[f].Find(ped[sib[i]].sibs[j]->serial)>-1)
               sibship[f][i].Push(ped[sib[i]].sibs[j]->serial);
      }
   }

   String outfile;
   outfile.Copy(prefix);
   outfile.Add(".skin");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "FID\tSibship\tPID\tPhi\tHetHet\tIBS0\tKinship\tError\n");
   double kinship, inflation;
   int id1, id2;
   Kinship kin;
   double phi, errorFlag, smaller;
   bool valid;
   for(int f = 0; f < ped.familyCount; f++){
      kin.Setup(*ped.families[f]);
      for(int s = 0; s < sibCount[f]; s++)
         for(int j = 0; j < id[f].Length(); j++){
            if(sibship[f][s].Find(id[f][j])>-1) continue;
            phi = kin(ped[id[f][j]], ped[sibship[f][s][0]]);
            valid = true;
            for(int i = 1; i < sibship[f][s].Length(); i++)
               if(kin(ped[id[f][j]], ped[sibship[f][s][i]])!=phi){
                  valid = false;
                  break;
               }
            if(!valid) continue;

            HetHetCount = IBS0Count = het1Count = het2Count = notMissingCount = 0;
            for(int i = 0; i < sibship[f][s].Length(); i++){
               id1 = geno[sibship[f][s][i]];
               id2 = geno[id[f][j]];
               for(int m = 0; m < shortCount; m++){
                  HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
                  IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
                  notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
                  het1Count += oneoneCount[(GG[0][id2][m] | GG[1][id2][m]) & (~GG[0][id1][m]) & GG[1][id1][m]];
                  het2Count += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
               }
            }
            if(het1Count == 0 || het2Count == 0) continue;
            kinship = (HetHetCount - IBS0Count*2.0)/(het1Count+het2Count);
            errorFlag = 0;
            inflation = (phi > 0)? kinship / phi: -1;
            if(phi > 0.03){
               if(inflation > 1.4142 || inflation < 0.70711)
                  errorFlag = 0.5;
               else if(inflation > 2 || inflation < 0.5)
                  errorFlag = 1;
            }else if(kinship > 0.022){
               if(kinship > 0.03125) errorFlag = 1;
               else errorFlag = 0.5;
            }
            fprintf(fp, "%s\t%s",
               (const char*)ped[sibship[f][s][0]].famid,
               (const char*)ped[sibship[f][s][0]].pid);
            for(int i = 1; i < sibship[f][s].Length(); i++)
               fprintf(fp, ",%s",
               (const char*)ped[sibship[f][s][i]].pid);
            fprintf(fp, "\t%s\t%.3lf\t%.3lf\t%.3lf\t%.4lf\t%GG",
               (const char*)ped[id[f][j]].pid, phi,
               HetHetCount*1.0/notMissingCount, IBS0Count*1.0/notMissingCount, kinship, errorFlag);
            fprintf(fp, "\n");
         }
   }
   fclose(fp);
   printf("Within-family kinship with sibships saved in file %s\n", (const char*)outfile);

   printf("Relationship inference across families starts at %s", currentTime());

   outfile.Copy(prefix);
   outfile.Add(".skin0");
   fp = fopen(outfile, "wt");
   fprintf(fp, "FID1\tSibship\tFID2\tPID\tHetHet\tIBS0\tKinship\n");

   for(int f1 = 0; f1 < ped.familyCount; f1++)
   for(int s = 0; s < sibCount[f1]; s++)
   for(int f2 = 0; f2 < ped.familyCount; f2++){
      if(f2==f1) continue;
      for(int j = 0; j < id[f2].Length(); j++){
         HetHetCount = IBS0Count = het1Count = het1Count = notMissingCount = 0;
         for(int i = 0; i < sibship[f1][s].Length(); i++){
            id1 = geno[sibship[f1][s][i]];
            id2 = geno[id[f2][j]];
            for(int m = 0; m < shortCount; m++){
               HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
               IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
               notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
               het1Count += oneoneCount[(GG[0][id2][m] | GG[1][id2][m]) & (~GG[0][id1][m]) & GG[1][id1][m]];
               het2Count += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
            }
         }
         if(het1Count == 0 || het2Count == 0) continue;
         smaller = het1Count < het2Count? het1Count: het2Count;
         kinship = 0.5 * (1 - (het1Count+het2Count)*0.5/smaller + (HetHetCount - IBS0Count*2.0)/smaller);
         fprintf(fp, "%s\t%s",
               (const char*)ped[sibship[f1][s][0]].famid,
               (const char*)ped[sibship[f1][s][0]].pid);
         for(int i = 1; i < sibship[f1][s].Length(); i++)
               fprintf(fp, ",%s",
               (const char*)ped[sibship[f1][s][i]].pid);
         fprintf(fp, "\t%s\t%s\t%.3lf\t%.3lf\t%.4lf",
               (const char*)ped[id[f2][j]].famid,
               (const char*)ped[id[f2][j]].pid,
               HetHetCount*1.0/notMissingCount, IBS0Count*1.0/notMissingCount, kinship);
         fprintf(fp, "\n");
      }
   }
   fclose(fp);
   printf("                                         ends at %s", currentTime());
   printf("Between-family kinship with sibships saved in file %s\n", (const char*)outfile);

   delete []sibship;
}
void Engine::WriteMerlin()
{
   String datfile = prefix;
   datfile.Add(".dat");
   String pedfile = prefix;
   pedfile.Add(".ped");
   String mapfile = prefix;
   mapfile.Add(".map");
   if(geno.Length()==0) {// Merlin data yet
      ped.WriteDataFile(datfile);
      ped.WritePedigreeFile(pedfile);
      ped.WriteMapFile(mapfile);
      printf("Genotype data saved as Merlin format in files %s, %s and %s\n",
         (const char*)datfile, (const char*)pedfile, (const char*)mapfile);
      return;
   }
   FILE *fp = fopen(datfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)datfile);
   for(int t = 0; t < ped.affectionCount; t++)
      fprintf(fp, "A %s\n", (const char*)ped.affectionNames[t]);
   for(int t = 0; t < ped.traitCount; t++)
      fprintf(fp, "T %s\n", (const char*)ped.traitNames[t]);
   for(int t = 0; t < ped.covariateCount; t++)
      fprintf(fp, "C %s\n", (const char*)ped.covariateNames[t]);
   if(ZG.Length())
      fprintf(fp, "Z Twin\n");
   if(snpName.Length() == 0)  // Marker info not provided
      for(int i = 0; i < markerCount; i++)
         fprintf(fp, "M SNP%d\n", i);
   else
      for(int i = 0; i < snpName.Length(); i++)
         fprintf(fp, "M %s\n", (const char*)snpName[i]);
   if(xpositions.Length())
      for(int i = 0; i < xsnpName.Length(); i++)
         fprintf(fp, "M %s\n", (const char*)xsnpName[i]);
   if(ypositions.Length())
      for(int i = 0; i < ysnpName.Length(); i++)
         fprintf(fp, "M %s\n", (const char*)ysnpName[i]);
   if(mtpositions.Length())
      for(int i = 0; i < mtsnpName.Length(); i++)
         fprintf(fp, "M %s\n", (const char*)mtsnpName[i]);
   fclose(fp);
   fp = fopen(pedfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
   for(int f = 0; f < ped.familyCount; f++){
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
         fprintf(fp, "%s %s %s %s %d",
            (const char*)ped[i].famid,
            (const char*)ped[i].pid,
            (const char*)ped[i].fatid,
            (const char*)ped[i].motid,
            ped[i].sex);
         for(int t = 0; t < ped.affectionCount; t++)
            fprintf(fp, " %d", ped[i].affections[t]);
         for(int t = 0; t < ped.traitCount; t++)
            if(ped[i].traits[t] == _NAN_)
               fprintf(fp, " X");
            else
               fprintf(fp, " %lf", ped[i].traits[t]);
         for(int t = 0; t < ped.covariateCount; t++)
            if(ped[i].covariates[t] == _NAN_)
               fprintf(fp, " X");
            else
               fprintf(fp, " %lf", ped[i].covariates[t]);
         int k = geno[i];
         if(ZG.Length()){
            if(geno[i] > -1)
               fprintf(fp, " %s", (const char*)ZG[k]);
            else
               fprintf(fp, " X");
         }
         if(geno[i] > -1){   // genotype available
            for(int m = 0; m < markerCount; m++){
               int byte = m/16;
               int offset = m%16;
               if(GG[0][k][byte] & shortbase[offset]){  // homozygote
                  if(GG[1][k][byte] & shortbase[offset])   // AA
                     fprintf(fp, " %s %s",
                     (const char*)alleleLabel[0][m], (const char*)alleleLabel[0][m]);
                  else  // aa
                     fprintf(fp, " %s %s",
                     (const char*)alleleLabel[1][m], (const char*)alleleLabel[1][m]);
               }else if(GG[1][k][byte] & shortbase[offset])   // Aa
                  fprintf(fp, " %s %s",
                  (const char*)alleleLabel[0][m], (const char*)alleleLabel[1][m]);
               else fprintf(fp, " X X");
            }
            if(xpositions.Length())
            for(int m = 0; m < xmarkerCount; m++){
               int byte = m/16;
               int offset = m%16;
               if(XG[0][k][byte] & shortbase[offset]){  // homozygote
                  if(XG[1][k][byte] & shortbase[offset])   // AA
                     fprintf(fp, " %s %s",
                     (const char*)xalleleLabel[0][m], (const char*)xalleleLabel[0][m]);
                  else  // aa
                     fprintf(fp, " %s %s",
                     (const char*)xalleleLabel[1][m], (const char*)xalleleLabel[1][m]);
               }else if(XG[1][k][byte] & shortbase[offset])   // Aa
                  fprintf(fp, " %s %s",
                  (const char*)xalleleLabel[0][m], (const char*)xalleleLabel[1][m]);
               else fprintf(fp, " X X");
            }

            if(ypositions.Length())
            for(int m = 0; m < ymarkerCount; m++){
               int byte = m/16;
               int offset = m%16;
               if(YG[0][k][byte] & shortbase[offset]){  // homozygote
                  if(YG[1][k][byte] & shortbase[offset])   // AA
                     fprintf(fp, " %s %s",
                     (const char*)yalleleLabel[0][m], (const char*)yalleleLabel[0][m]);
                  else  // aa
                     fprintf(fp, " %s %s",
                     (const char*)yalleleLabel[1][m], (const char*)yalleleLabel[1][m]);
               }else if(YG[1][k][byte] & shortbase[offset])   // Aa
                  fprintf(fp, " %s %s",
                  (const char*)yalleleLabel[0][m], (const char*)yalleleLabel[1][m]);
               else fprintf(fp, " X X");
            }
            if(mtpositions.Length())
            for(int m = 0; m < mtmarkerCount; m++){
               int byte = m/16;
               int offset = m%16;
               if(MG[0][k][byte] & shortbase[offset]){  // homozygote
                  if(MG[1][k][byte] & shortbase[offset])   // AA
                     fprintf(fp, " %s %s",
                     (const char*)mtalleleLabel[0][m], (const char*)mtalleleLabel[0][m]);
                  else  // aa
                     fprintf(fp, " %s %s",
                     (const char*)mtalleleLabel[1][m], (const char*)mtalleleLabel[1][m]);
               }else if(MG[1][k][byte] & shortbase[offset])   // Aa
                  fprintf(fp, " %s %s",
                  (const char*)mtalleleLabel[0][m], (const char*)mtalleleLabel[1][m]);
               else fprintf(fp, " X X");
            }
         }else{     // genotype unavailable
            for(int m = 0; m < markerCount; m++)
               fprintf(fp, " X X");
            if(xpositions.Length())
               for(int m = 0; m < xmarkerCount; m++)
                  fprintf(fp, " X X");
            if(ypositions.Length())
               for(int m = 0; m < ymarkerCount; m++)
                  fprintf(fp, " X X");
            if(mtpositions.Length())
               for(int m = 0; m < mtmarkerCount; m++)
                  fprintf(fp, " X X");
         }
         fprintf(fp, "\n");
      }
   }
   fclose(fp);

   if(positions.Length() || xpositions.Length() || ypositions.Length() || mtpositions.Length()){
      fp = fopen(mapfile, "wt");
      if(fp == NULL)
         error("Cannot open %s to write.", (const char*)mapfile);
         for(int m = 0; m < markerCount; m++)
            fprintf(fp, "%d\t%s\t%.6lf\n",
               chromosomes[m], (const char*)snpName[m], positions[m]);
      if(xpositions.Length())
         for(int m = 0; m < xmarkerCount; m++)
            fprintf(fp, "%d\t%s\t%.6lf\n",
               SEXCHR, (const char*)xsnpName[m], xpositions[m]);
      if(ypositions.Length())
         for(int m = 0; m < ymarkerCount; m++)
            fprintf(fp, "%d\t%s\t%.6lf\n",
               SEXCHR+1, (const char*)ysnpName[m], ypositions[m]);
      if(mtpositions.Length())
         for(int m = 0; m < mtmarkerCount; m++)
            fprintf(fp, "%d\t%s\t%.6lf\n",
               SEXCHR+3, (const char*)mtsnpName[m], mtpositions[m]);
      fclose(fp);
      printf("Genotype data saved as MERLIN format in files %s, %s and %s\n",
         (const char*)datfile, (const char*)pedfile, (const char*)mapfile);
   }else
      printf("Genotype data saved as MERLIN format in files %s and %s\n",
         (const char*)datfile, (const char*)pedfile);
}

