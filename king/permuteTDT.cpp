// permuteTDT.cpp
// 6/20/2011 Wei-Min Chen

#include "analysis.h"
#include "MathStats.h"
#include "Kinship.h"
#include "QuickIndex.h"
#include "Random.h"
#include "time.h"

void Engine::permuteRareTDT()
{
   Haplotyping();

   int totalSNP = rmarkers.Length();
   IntArray thetransmission(Ltrio.Length()/3*2);
   thetransmission.Set(-10000);

   // obtain thetransmission
   int ids[3];
   int match[2];
   for(int i = 0; i < Ltrio.Length()/3; i++){
      ids[0] = geno[Ltrio[i*3+1]];
      ids[1] = geno[Ltrio[i*3+2]];
      ids[2] = geno[Ltrio[i*3]]; // offspring
      for(int k = 0; k < 2; k++){
         match[0]=match[1]=1;
         for(int m = 0; m < totalSNP; m++)
            if(hap[k][ids[2]][m] >= 0 && // offspring haplotype
               hap[0][ids[k]][m] >= 0 && // parent haplotype 0
               hap[k][ids[2]][m] != hap[0][ids[k]][m]){
               match[0] = 0;
               break;
            }
         for(int m = 0; m < totalSNP; m++)
            if(hap[k][ids[2]][m] >= 0 &&
               hap[1][ids[k]][m] >= 0 && // parent haplotype 1
               hap[k][ids[2]][m] != hap[1][ids[k]][m]){
               match[1] = 0;
               break;
            }
         if(match[0] == 1 && match[1] == 0)
            thetransmission[i*2+k] = 0;
         else if(match[0] == 0 && match[1] == 1)
            thetransmission[i*2+k] = 1;
         else if(match[0] == 1 && match[1] == 1)
            thetransmission[i*2+k] = -9999;
      }
   }

   if(chisqFilter!=_NAN_)
      printf("Only SNPs with SNP-level TDT statistic > %.2lf are analyzed\n",
         chisqFilter);

   // Get the TDT statistic for the observed data
   double thechisq1 = OnePermuteRareTDT(thetransmission, 1);
   printf("The collapsed (protective effect) TDT P value is %.3G (%d/%d, chisq=%.2lf)\n",
      chidist(thechisq1, 1), TDT_T[1], TDT_T[0], thechisq1);

   double thechisq0 = OnePermuteRareTDT(thetransmission, 0);
   printf("The collapsed (risk effect) TDT P value is %.3G (%d/%d, chisq=%.2lf)\n",
      chidist(thechisq0, 1), TDT_T[1], TDT_T[0], thechisq0);

   // Permute and get TDT distribution
   printf("Permutation start");
   time_t t;
   Random rand((long)time(&t));
   Vector chisq1(0);
   Vector chisq0(0);
   IntArray transmission(Ltrio.Length()/3*2);
   for(int r = 0; r < permuteCount; r++){
      if(r%100==0) printf(".");
      for(int i = 0; i < Ltrio.Length()/3*2; i++)
         transmission[i] = int(rand.Next()+0.5);
      chisq1.Push(OnePermuteRareTDT(transmission, 1));
      chisq0.Push(OnePermuteRareTDT(transmission, 0));
   }
   printf("\n");
   chisq1.Sort();
   chisq0.Sort();
   int rank;
   for(rank = permuteCount-1; rank >=0 && chisq1[rank] > thechisq1; rank--);
   double pvalue = 1-(rank + 0.5)/permuteCount;
   printf("Permuted P value for protective effect is %.4G\n", pvalue);
   for(rank = permuteCount-1; rank >=0 && chisq0[rank] > thechisq0; rank--);
   pvalue = 1-(rank + 0.5)/permuteCount;
   printf("Permuted P value for risk effect is %.4G\n", pvalue);
}

void Engine::permuteTDT()
{
   Haplotyping();

   int totalSNP = rmarkers.Length();
   IntArray thetransmission(Ltrio.Length()/3*2);
   thetransmission.Set(-10000);

   // obtain thetransmission
   int ids[3];
   int match[2];
   for(int i = 0; i < Ltrio.Length()/3; i++){
      ids[0] = geno[Ltrio[i*3+1]];
      ids[1] = geno[Ltrio[i*3+2]];
      ids[2] = geno[Ltrio[i*3]]; // offspring
      for(int k = 0; k < 2; k++){
         match[0]=match[1]=1;
         for(int m = 0; m < totalSNP; m++)
            if(hap[k][ids[2]][m] >= 0 && // offspring haplotype
               hap[0][ids[k]][m] >= 0 && // parent haplotype 0
               hap[k][ids[2]][m] != hap[0][ids[k]][m]){
               match[0] = 0;
               break;
            }
         for(int m = 0; m < totalSNP; m++)
            if(hap[k][ids[2]][m] >= 0 &&
               hap[1][ids[k]][m] >= 0 && // parent haplotype 1
               hap[k][ids[2]][m] != hap[1][ids[k]][m]){
               match[1] = 0;
               break;
            }
         if(match[0] == 1 && match[1] == 0)
            thetransmission[i*2+k] = 0;
         else if(match[0] == 0 && match[1] == 1)
            thetransmission[i*2+k] = 1;
         else if(match[0] == 1 && match[1] == 1)
            thetransmission[i*2+k] = -9999;
      }
   }

   // Get the TDT statistic for the observed data
   double thechisq = OnePermuteTDT(thetransmission);
   printf("The smallest TDT P value is %.3G (%d/%d, chisq=%.2lf) at SNP %s\n",
      chidist(thechisq, 1), TDT_T[1], TDT_T[0], thechisq,
      (const char*)snpName[rmarkers[TDT_PeakPos]]);

   // Permute and get TDT distribution
   printf("Permutation start");
   time_t t;
   Random rand((long)time(&t));
   Vector chisq(0);
   IntArray transmission(Ltrio.Length()/3*2);
   for(int r = 0; r < permuteCount; r++){
      if(r%100==0) printf(".");
      for(int i = 0; i < Ltrio.Length()/3*2; i++)
         transmission[i] = int(rand.Next()+0.5);
      chisq.Push(OnePermuteTDT(transmission));
   }
   printf("\n");
   chisq.Sort();
   printf("P value quantiles: Q0.05=%.2G, Q0.01=%.2G, Q0.001=%.2G, Q0.0001=%.2G\n",
      chidist(chisq[int(permuteCount*0.95-.5)],1),
      chidist(chisq[int(permuteCount*0.99-.5)],1),
      chidist(chisq[int(permuteCount*0.999-.5)],1),
      chidist(chisq[int(permuteCount*0.9999-.5)],1));
   int rank;
   for(rank = permuteCount-1; rank >=0 && chisq[rank] > thechisq; rank--);
   double pvalue = 1-(rank + 0.5)/permuteCount;
   printf("Permuted P value is %.4G\n", pvalue); 
}

double Engine::OnePermuteTDT(IntArray & onetransmission)
{
   double maxchisq = -9999;
   int ids[2];
   IntArray T[2];
   T[0].Dimension(rmarkers.Length());
   T[0].Zero();
   T[1].Dimension(rmarkers.Length());
   T[1].Zero();

   for(int i = 0; i < Ltrio.Length()/3; i++){
      ids[0] = geno[Ltrio[i*3+1]];
      ids[1] = geno[Ltrio[i*3+2]];
      for(int k = 0; k < 2; k++)
         if(onetransmission[i*2+k]>=0)
            for(int m = 0; m < rmarkers.Length(); m++)
               if(hap[0][ids[k]][m] >= 0 && hap[0][ids[k]][m] != hap[1][ids[k]][m])
                  T[hap[onetransmission[i*2+k]][ids[k]][m]][m] ++;
   }
   double chisq, t;
   for(int m = 0; m < rmarkers.Length(); m++){
      t = T[0][m]-T[1][m];
      if(T[0][m]+T[1][m]==0)
         continue;
      chisq = t*t/(T[0][m]+T[1][m]);
      if(chisq > maxchisq){
         maxchisq = chisq;
         TDT_T[0] = T[0][m];
         TDT_T[1] = T[1][m];
         TDT_PeakPos = m;
      }
   }
   if(maxchisq < 0) maxchisq = 0;
   return maxchisq;
}

double Engine::OnePermuteRareTDT(IntArray & onetransmission, int Protective)
{
   int rarestrand, p;
   int ids[2];
   IntArray T[2];
   T[0].Dimension(rmarkers.Length());
   T[0].Zero();
   T[1].Dimension(rmarkers.Length());
   T[1].Zero();

   for(int i = 0; i < Ltrio.Length()/3; i++){
      ids[0] = geno[Ltrio[i*3+1]];
      ids[1] = geno[Ltrio[i*3+2]];
      for(int k = 0; k < 2; k++)
         if(onetransmission[i*2+k]>=0)
            for(int m = 0; m < rmarkers.Length(); m++)
               if(hap[0][ids[k]][m] >= 0 && hap[0][ids[k]][m] != hap[1][ids[k]][m])
                  T[hap[onetransmission[i*2+k]][ids[k]][m]][m] ++;
   }
   IntArray valid(rmarkers.Length());
   valid.Set(1);
   for(int m = 0; m < rmarkers.Length(); m++)
      if( (Protective && T[0][m] < T[1][m]) ||
         (!Protective && T[0][m] > T[1][m]) ||
         (T[0][m]+T[1][m] == 0) )
            valid[m] = 0;

   if(chisqFilter!=_NAN_){
      double t;
      for(int m = 0; m < rmarkers.Length(); m++){
         if(valid[m]==0) continue;
         t = T[0][m] - T[1][m];
         if(t*t < chisqFilter*(T[0][m]+T[1][m]))
            valid[m] = 0;
      }
   }           

   TDT_T[0] = TDT_T[1] = 0;
   for(int i = 0; i < Ltrio.Length()/3; i++){
      ids[0] = geno[Ltrio[i*3+1]];
      ids[1] = geno[Ltrio[i*3+2]];
      for(int k = 0; k < 2; k++){
         p = -1;
         rarestrand = -1;
         for(int m = 0; m < rmarkers.Length(); m++)
            if(valid[m] && hap[0][ids[k]][m] >= 0 && hap[0][ids[k]][m] != hap[1][ids[k]][m]){
               rarestrand = hap[0][ids[k]][m]? 0: 1;
               p = m;
               break;
            }
         if(rarestrand!=-1)
         for(int m = p+1; m < rmarkers.Length(); m++)
            if(valid[m] && hap[rarestrand][ids[k]][m] == 0
               && hap[1-rarestrand][ids[k]][m] == 1){
               rarestrand = -1;
               break;
            }
         if(rarestrand!=-1){
            if(rarestrand == onetransmission[i*2+k]) // rare variant transmitted
               TDT_T[1] ++;
            else  // common variant transmitted
               TDT_T[0] ++;
         }
      }
   }
   if(TDT_T[0] == 0 && TDT_T[1] == 0) return 0;
   double t = TDT_T[0] - TDT_T[1];
   double chisq = t*t/(TDT_T[0]+TDT_T[1]);
   return chisq;
}
