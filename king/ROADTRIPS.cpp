// ROADTRIPS.cpp
// 2/12/2013 Wei-Min Chen

#include "analysis.h"
#include <math.h>
#include "Kinship.h"
#include "KinshipX.h"
#include "MathStats.h"
#include "MathSVD.h"
#include "QuickIndex.h"
#include "MathCholesky.h"

void Engine::ROADTRIPS()
{
   if(ped.affectionCount==0) error("Disease data are not available.");
   if(ped.affectionCount > 1)
      warning("Only the first disease data (%s) is analyzed",
         (const char*)ped.affectionNames[0]);
   if(prevalence == 0)
      printf("Disease prevalence needs to be specified (through --prevalence option)\n");
   else if(prevalence >= 1 || prevalence < 0)
      error("Disease prevalence needs to be specified (through --prevalence option)");
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

   IntArray *gid = new IntArray[ped.familyCount];
   IntArray gidCount(ped.familyCount);
   IntArray *RTid = new IntArray[ped.familyCount];
   IntArray RTidCount(ped.familyCount);
   IntArray vid(idCount);
   vid.Set(-1);
   int vidCount = 0;
   for(int f = 0; f < ped.familyCount; f++){
      gid[f].Dimension(0);
      RTid[f].Dimension(0);
      for(int i = 0; i < id[f].Length(); i++)
         if(ped[id[f][i]].ngeno > MINSNPCOUNT){
            gid[f].Push(vidCount);
            RTid[f].Push(id[f][i]);
            vid[geno[id[f][i]]] = vidCount++;
         }
      gidCount[f] = gid[f].Length();
      for(int i = 0; i < id[f].Length(); i++)
         if(ped[id[f][i]].ngeno <= MINSNPCOUNT && ped[id[f][i]].affections[0])
            RTid[f].Push(id[f][i]);
      RTidCount[f] = RTid[f].Length();
   }
   int id1, id2;
   int IBS0Count, IBS1Count, notMissingCount;
   Matrix D(vidCount, vidCount); // Upper right for correlation, lower left for VDV0
   D.Zero();
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = 0; i < gidCount[f]; i++)
         for(int j = i+1; j < gidCount[f]; j++){
            id1 = geno[RTid[f][i]]; id2 = geno[RTid[f][j]];
            IBS0Count = IBS1Count = notMissingCount = 0;
            for(int m = 0; m < shortCount; m++){
               IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
               notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
               IBS1Count += oneoneCount[ (GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m]) |
                  (GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]) ];
            }
            if(notMissingCount)
               D.data[vid[id1]]->data[vid[id2]] = D.data[vid[id2]]->data[vid[id1]] =
                  (IBS1Count + 4.0*IBS0Count) / notMissingCount;
         }
   for(int f1 = 0; f1 < ped.familyCount; f1++)
      for(int i = 0; i < gidCount[f1]; i++) {
         id1 = geno[RTid[f1][i]]; 
         for(int f2 = f1+1; f2 < ped.familyCount; f2++)
         for(int j = 0; j < gidCount[f2]; j++){
            id2 = geno[RTid[f2][j]];
            IBS0Count = IBS1Count = notMissingCount = 0;
            for(int m = 0; m < shortCount; m++){
               IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
               notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
               IBS1Count += oneoneCount[ (GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m]) |
                  (GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]) ];
            }
            if(notMissingCount)
               D.data[vid[id1]]->data[vid[id2]] = D.data[vid[id2]]->data[vid[id1]]
                  = (IBS1Count + 4.0*IBS0Count) / notMissingCount;
         }
      }

   printf("Genotypes stored in %d integers for each of %d individuals used in analysis.\n",
      shortCount, vidCount);

   Vector tempV(vidCount);
   tempV.Zero();
   for(int j = 0; j < vidCount; j++)
      for(int i = 0; i < vidCount; i++)
         tempV[j] += D[i][j];
   tempV.Multiply(1.0/vidCount);             
   for(int i = 0; i < vidCount; i++)
      for(int j = 0; j < vidCount; j++)
         D[i][j] -= tempV[j];// (I-11'/N) * D
   tempV.Zero();
   for(int i = 0; i < vidCount; i++)
      for(int j = 0; j < vidCount; j++)
         tempV[i] += D[i][j];
   tempV.Multiply(1.0/vidCount);
   for(int i = 0; i < vidCount; i++)
      for(int j = 0; j < vidCount; j++)
         D[i][j] -= tempV[i];// D * (I-11'/N)
   D.Multiply(-0.5);
   for(int i = 0; i < vidCount; i++)
      D[i][i] = sqrt(D[i][i]); // save computation for below
   for(int i = 0; i < vidCount; i++)
      for(int j = i+1; j < vidCount; j++)
         D[i][j] /= (D[i][i]*D[j][j]);
   for(int i = 0; i < vidCount; i++)
      for(int j = 0; j <= i; j++)
         D[i][j] = 0;   // lower left of D is for VDV0
   printf("Similarity matrix consists of sample covariance estimates, standardized by diagonal elements\n");
   int degree[3];
   for(int d = 0; d < 3; d++)
      degree[d] = 0;
   for(int i = 0; i < vidCount; i++)
      for(int j = i+1; j < vidCount; j++)
         if(D[i][j] > 0.8) // MZ twins
            degree[0] ++;
         else if(D[i][j] > 0.354) // 1st-degree relatives
            degree[1] ++;
         else if(D[i][j] > 0.177) // 2nd-degree relatives
            degree[2] ++;
   if(degree[0] || degree[1] || degree[2]){
      printf("Close relatives identified: %d MZ twins/duplicates, %d 1st-degree and %d 2nd-degree relative pairs\n",
         degree[0], degree[1], degree[2]);
      if(degree[0])
      warning("One of the duplicates needs to be removed prior to association analysis.\n");
   }else
      printf("No second-degree or closer relationships were found\n");

   int affCount[3];
   for(int i = 0; i < 3; i++)
      affCount[i] = 0;
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = 0; i < gidCount[f]; i++)
         affCount[ped[RTid[f][i]].affections[0]]++;
   printf("Among %d genotyped individuals, %d affected, and %d unaffected\n",
      vidCount, affCount[2], affCount[1]);

   // Pre-Compute
   Cholesky chol;
   Kinship kin;
   Vector *V0 = new Vector[ped.familyCount];
   Matrix PhiInv0;
   Vector *Phi0One = new Vector[ped.familyCount];
   Matrix *PhiInv = new Matrix[ped.familyCount];
   Matrix *PhiN0 = new Matrix[ped.familyCount];
   Matrix *PhiM0 = new Matrix[ped.familyCount];
   Vector *AN0 = new Vector[ped.familyCount];
   Vector *AM0 = new Vector[ped.familyCount];
   Vector *PhiOne[16];

   Vector *tV = new Vector[ped.familyCount];
   for(int f = 0; f < ped.familyCount; f++)
      tV[f].Dimension(gidCount[f]);
   int vid1, vid2;
   for(int f = 0; f < ped.familyCount; f++){
      if(gidCount[f]==0) continue;
      PhiN0[f].Dimension(gidCount[f], gidCount[f]);
      if(gidCount[f]==1){
         PhiN0[f][0][0] = 1;
         PhiInv0.Dimension(1,1);
         PhiInv0[0][0] = 1;
      }else{
         kin.Setup(*ped.families[f]);
         for(int i = 0; i < gidCount[f]; i++)
            for(int j = i+1; j < gidCount[f]; j++)
               PhiN0[f][i][j] = PhiN0[f][j][i] = kin(ped[RTid[f][i]], ped[RTid[f][j]])*2;
         for(int i = 0; i < gidCount[f]; i++)
            PhiN0[f][i][i] = kin(ped[RTid[f][i]], ped[RTid[f][i]])*2;
         chol.Decompose(PhiN0[f]);
         chol.Invert();
         PhiInv0 = chol.inv;
      }
      AN0[f].Dimension(gidCount[f]);
      for(int i = 0; i < gidCount[f]; i++)
         if(ped[RTid[f][i]].affections[0] == 2) AN0[f][i] = 1;
         else if(ped[RTid[f][i]].affections[0] == 1) AN0[f][i] = -prevalence/(1.0-prevalence);
         else if(ped[RTid[f][i]].affections[0] == 0) AN0[f][i] = 0;

      if(gidCount[f] < RTidCount[f]){
         if(gidCount[f]==1) kin.Setup(*ped.families[f]);
         PhiM0[f].Dimension(gidCount[f], RTidCount[f]-gidCount[f]);
         for(int i = gidCount[f]; i < RTidCount[f]; i++)
            for(int j = 0; j < gidCount[f]; j++)
               PhiM0[f][j][i-gidCount[f]] = kin(ped[RTid[f][i]], ped[RTid[f][j]])*2;
         AM0[f].Dimension(RTidCount[f]-gidCount[f]);
         for(int i = gidCount[f]; i < RTidCount[f]; i++)
            if(ped[RTid[f][i]].affections[0] == 2) AM0[f][i-gidCount[f]] = 1;
            else if(ped[RTid[f][i]].affections[0] == 1) AM0[f][i-gidCount[f]] = -prevalence/(1.0-prevalence);
            else if(ped[RTid[f][i]].affections[0] == 0) AM0[f][i-gidCount[f]] = 0;
      }
      V0[f].Dimension(gidCount[f]);
      V0[f].Zero();
      if(gidCount[f] < RTidCount[f]){   // missing genotypes in the family
         tV[f].Zero();
         for(int i = 0; i < gidCount[f]; i++)
            for(int j = 0; j < RTidCount[f]-gidCount[f]; j++)
               tV[f][i] += PhiM0[f][i][j] * AM0[f][j];
         for(int i = 0; i < gidCount[f]; i++)
            for(int j = 0; j < gidCount[f]; j++)
               V0[f][i] += PhiInv0[i][j] * tV[f][j];
      }
      for(int i = 0; i < gidCount[f]; i++)
         V0[f][i] += AN0[f][i];
      Phi0One[f].Dimension(gidCount[f]);
      Phi0One[f].Zero();
      for(int i = 0; i < gidCount[f]; i++)
         for(int j = 0; j < gidCount[f]; j++)
            Phi0One[f][i] += PhiInv0[i][j];
   }
   double V0One = 0.0;
   double OnePhi0One = 0.0;
   for(int f = 0; f < ped.familyCount; f++){
      if(gidCount[f]==0) continue;
      V0One += V0[f].Sum();
      OnePhi0One += Phi0One[f].Sum();
   }
   double V0Mean = V0One / OnePhi0One;
   for(int f = 0; f < ped.familyCount; f++){
      for(int i = 0; i < gidCount[f]; i++)
         V0[f][i] = V0[f][i] - V0Mean * Phi0One[f][i];
      for(int i = 0; i < gidCount[f]; i++){
         D.data[gid[f][i]]->data[gid[f][i]] = V0[f][i]*V0[f][i];
         for(int j = 0; j < i; j++){
            vid1 = gid[f][i];
            vid2 = gid[f][j];
            D.data[vid1]->data[vid2] = V0[f][i] * D.data[vid2]->data[vid1] * V0[f][j];
         }
      }
   }
   for(int f1 = 0; f1 < ped.familyCount; f1++)
      for(int i = 0; i < gidCount[f1]; i++){
         vid1 = gid[f1][i];
         for(int f2 = 0; f2 < f1; f2++)
            for(int j = 0; j < gidCount[f2]; j++){
               vid2 = gid[f2][j];
               D.data[vid1]->data[vid2] = V0[f1][i] * D.data[vid2]->data[vid1] * V0[f2][j];
            }
      }

   double VDV0total = 0.0;
   for(int i = 0; i < vidCount; i++)
      for(int j = 0; j < i; j++)
         VDV0total += D.data[i]->data[j];
   VDV0total *= 2;
   for(int i = 0; i < vidCount; i++)
      VDV0total += D.data[i]->data[i];

   double VDOne0 = 0.0;
   double OneDOne0 = 0.0;
   for(int f1 = 0; f1 < ped.familyCount; f1++){
      for(int i = 0; i < gidCount[f1]; i++){
         vid1 = gid[f1][i];
         VDOne0 += V0[f1][i] * Phi0One[f1][i];
         OneDOne0 += Phi0One[f1][i] * Phi0One[f1][i];
         for(int j = 0; j < i; j++){
            vid2 = gid[f1][j];
            VDOne0 += V0[f1][i] * D.data[vid2]->data[vid1] * Phi0One[f1][j];
            OneDOne0 += Phi0One[f1][i] * D.data[vid2]->data[vid1] * Phi0One[f1][j];
         }
         for(int j = i+1; j < gidCount[f1]; j++){
            vid2 = gid[f1][j];
            VDOne0 += V0[f1][i] * D.data[vid1]->data[vid2] * Phi0One[f1][j];
            OneDOne0 += Phi0One[f1][i] * D.data[vid1]->data[vid2] * Phi0One[f1][j];
         }
         for(int f2 = 0; f2 < f1; f2++)
            for(int j = 0; j < gidCount[f2]; j++){
               vid2 = gid[f2][j];
               VDOne0 += V0[f1][i] * D.data[vid2]->data[vid1] * Phi0One[f2][j];
               OneDOne0 += Phi0One[f1][i] * D.data[vid2]->data[vid1] * Phi0One[f2][j];
            }
         for(int f2 = f1+1; f2 < ped.familyCount; f2++)
            for(int j = 0; j < gidCount[f2]; j++){
               vid2 = gid[f2][j];
               VDOne0 += V0[f1][i] * D.data[vid1]->data[vid2] * Phi0One[f2][j];
               OneDOne0 += Phi0One[f1][i] * D.data[vid1]->data[vid2] * Phi0One[f2][j];
            }
      }
   }
   
   if(HeritFlag){
      double Corr_AA, Corr_UU, Corr_AU;
      int N_AA, N_UU, N_AU;
      int aff1, aff2;
      Corr_AA = Corr_UU = Corr_AU = 0.0;
      N_AA = N_UU = N_AU = 0;
      for(int f1 = 0; f1 < ped.familyCount; f1++)
      for(int i = 0; i < gidCount[f1]; i++){
         vid1 = gid[f1][i];
         aff1 = ped[RTid[f1][i]].affections[0];
         for(int f2 = f1+1; f2 < ped.familyCount; f2++)
         for(int j = 0; j < gidCount[f2]; j++){
            vid2 = gid[f2][j];
            aff2 = ped[RTid[f2][j]].affections[0];
            if(D[vid1][vid2] > 0.354) continue; // relatives
            if(aff1 == 2 && aff2 == 2){
               N_AA ++;
               Corr_AA += D[vid1][vid2];
            }else if(aff1 == 1 && aff2 == 1){
               N_UU ++;
               Corr_UU += D[vid1][vid2];
            }else if( (aff1 == 1 && aff2 == 2) || (aff1 == 2 && aff2 ==1) ){
               N_AU ++;
               Corr_AU += D[vid1][vid2];
            }
         }
      }
      if(N_AA) Corr_AA /= N_AA;
      if(N_UU) Corr_UU /= N_UU;
      if(N_AU) Corr_AU /= N_AU;
      printf("Average correlation among %d AA, %d UU, and %d AU unrelated pairs is %.3G, %.3G, and %.3G\n",
         N_AA, N_UU, N_AU, Corr_AA, Corr_UU, Corr_AU);    

      printf("Variance of the score is %.3lf\n", VDV0total);
      double VV0 = 0;
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = 0; i < gidCount[f]; i++)
            VV0 += V0[f][i] * V0[f][i];
      printf("Variance assuming no structure is %.3lf\n", VV0);
      return;
   }
   printf("Prevalence of %s: %.3f\n",
      (const char*)ped.affectionNames[0], prevalence);
   printf("ROADTRIPS genome scan starts at %s", currentTime());

   QuickIndex idx;
   double freq[16];
   int AACount[16], AaCount[16], missingCount[16];
   int pos;
   int id, AA, Aa, missing;
   Vector chisq(0);
   String lmmfile = prefix;
   lmmfile.Add("rt.txt");
   FILE *fp = fopen(lmmfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)lmmfile);
   fprintf(fp, "SNP");
   if(chromosomes.Length()) fprintf(fp, "\tChr");
   if(bp.Length()) fprintf(fp, "\tPos");
   fprintf(fp, "\tLabelA\tLabela\tFreqA\tN\tChisq\tPvalue\n");

   IntArray *mIndex[16];
   IntArray *nIndex[16];
   for(int m = 0; m < 16; m++){
      mIndex[m] = new IntArray[ped.familyCount];
      nIndex[m] = new IntArray[ped.familyCount];
      for(int f = 0; f < ped.familyCount; f++){
         mIndex[m][f].Dimension(gidCount[f]+1);
         mIndex[m][f][0] = 0;
      }
      for(int f = 0; f < ped.familyCount; f++){
         nIndex[m][f].Dimension(gidCount[f]+1);
         nIndex[m][f][0] = 0;
      }
   }
   Vector score(16);
   double VDV;
   double VDOne;
   double OneDOne;

   Vector *V[16];
   IntArray NSample(16);
   Vector OnePhiOne(16);
   Vector VOne(16);
   Vector VMean(16);
   Vector VDelta(16);
   for(int m = 0; m < 16; m++){
      PhiOne[m] = new Vector[ped.familyCount];
      for(int f = 0; f < ped.familyCount; f++)
         PhiOne[m][f].Dimension(gidCount[f]);
   }
   double factor;
   for(int m = 0; m < 16; m++){
      V[m] = new Vector[ped.familyCount];
      for(int f = 0; f < ped.familyCount; f++)
         V[m][f].Dimension(gidCount[f]);
   }
   Vector *Y[16];
   for(int m = 0; m < 16; m++){
      Y[m] = new Vector[ped.familyCount];
      for(int f = 0; f < ped.familyCount; f++)
         Y[m][f].Dimension(gidCount[f]);
   }
   for(int b = 0; b < shortCount; b++){ // loop over short integer
      for(int i = 0; i < 16; i++){
         AACount[i] = AaCount[i] = missingCount[i] = 0;
         freq[i] = 0.0;
      }
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = 0; i < gidCount[f]; i++){
            id = geno[RTid[f][i]];
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
      for(int m = 0; m < 16; m++)
         if(missingCount[m] < vidCount)  // frequency of allele A
            freq[m] = (AACount[m] + AaCount[m]*0.5) / (vidCount - missingCount[m]);

      NSample.Zero();
      for(int m = 0; m < 16; m++)
         for(int f = 0; f < ped.familyCount; f++)
            for(int i = 0; i < gidCount[f]; i++)
               V[m][f][i] = V0[f][i] + V0One / OnePhi0One * Phi0One[f][i];
      OnePhiOne.Zero();
      VOne.Zero();
      score.Zero();
      for(int m = 0; m < 16; m++)
         for(int f = 0; f < ped.familyCount; f++)
            Y[m][f].Zero();
      for(int f = 0; f < ped.familyCount; f++){
         if(gidCount[f]==0) continue;
         for(int m = 0; m < 16; m++){
            mIndex[m][f][0] = nIndex[m][f][0] = 0;
            PhiOne[m][f].Zero();
         }
         if(gidCount[f]==1){  // population-based data
            id = geno[RTid[f][0]];
            AA = GG[0][id][b] & GG[1][id][b];
            Aa = (~GG[0][id][b]) & GG[1][id][b];
            missing = (~GG[0][id][b]) & (~GG[1][id][b]) & 65535;
            for(int m = 0; m < 16; m++){
               if(missing & shortbase[m]){
                  mIndex[m][f][0] = 1;
                  mIndex[m][f][1] = 0;
                  V[m][f][0] = 0; // V[m][f].Zero();
               }else{
                  NSample[m]++;
                  nIndex[m][f][0] = 1;
                  nIndex[m][f][1] = 0;
                  if(AA & shortbase[m]){
                     Y[m][f][0] = 1;
                  }else if(Aa & shortbase[m]){
                     Y[m][f][0] = 0.5;
                  }
                  VOne[m] += V[m][f][0];
                  OnePhiOne[m] ++;
                  PhiOne[m][f][0] = 1;
               }
            }
         }else{   // family data
            for(int i = 0; i < gidCount[f]; i++){ // genome scan
               id = geno[RTid[f][i]];
               AA = GG[0][id][b] & GG[1][id][b];
               Aa = (~GG[0][id][b]) & GG[1][id][b];
               missing = (~GG[0][id][b]) & (~GG[1][id][b]) & 65535;
               for(int m = 0; m < 16; m++){
                  if(missing & shortbase[m]){
                     //Y[m][f][i] = freq[m];
                     mIndex[m][f][0]++;
                     mIndex[m][f][mIndex[m][f][0]] = i;
                  }else{
                     NSample[m]++;
                     nIndex[m][f][0]++;
                     nIndex[m][f][nIndex[m][f][0]] = i;
                     if(AA & shortbase[m])
                        Y[m][f][i] = 1;
                     else if(Aa & shortbase[m])
                        Y[m][f][i] = 0.5;
                  }
               }
            }
            for(int m = 0; m < 16; m++){
               if(mIndex[m][f][0] && nIndex[m][f][0]){ // missingness in the family
                  tV[f].Zero();
                  V[m][f].Zero();
                  PhiInv[f].Dimension(nIndex[m][f][0], nIndex[m][f][0]);
                  if(nIndex[m][f][0]==1)  // 1 individual in the family
                     PhiInv[f][0][0] = 1;
                  else{
                     for(int i = 0; i < nIndex[m][f][0]; i++)
                        for(int j = 0; j < nIndex[m][f][0]; j++)
                           PhiInv[f][i][j] = PhiN0[f][nIndex[m][f][i+1]][nIndex[m][f][j+1]];
                     chol.Decompose(PhiInv[f]);
                     chol.Invert();
                     PhiInv[f] = chol.inv;
                  }
                  for(int i = 0; i < nIndex[m][f][0]; i++){
                     id1 = nIndex[m][f][i+1];
                     for(int j = 0; j < mIndex[m][f][0]; j++){
                        id2 = mIndex[m][f][j+1];
                        tV[f][id1] += PhiN0[f][id1][id2] * AN0[f][id2];
                     }
                     for(int j = 0; j < RTidCount[f]-gidCount[f]; j++)
                        tV[f][id1] += PhiM0[f][id1][j] * AM0[f][j];
                  }
                  for(int i = 0; i < nIndex[m][f][0]; i++)
                     for(int j = 0; j < nIndex[m][f][0]; j++)
                        V[m][f][nIndex[m][f][i+1]] +=
                        PhiInv[f][i][j] * tV[f][nIndex[m][f][j+1]];
                  for(int i = 0; i < nIndex[m][f][0]; i++){
                     id1 = nIndex[m][f][i+1];
                     V[m][f][id1] += AN0[f][id1];
                  }
                  VOne[m] += V[m][f].Sum();
                  for(int i = 0; i < nIndex[m][f][0]; i++)
                     for(int j = 0; j < nIndex[m][f][0]; j++)
                        PhiOne[m][f][nIndex[m][f][i+1]] += PhiInv[f][i][j];
                  OnePhiOne[m] += PhiOne[m][f].Sum();
               }else{   // no missingness
                  OnePhiOne[m] += Phi0One[f].Sum();
                  VOne[m] += V[m][f].Sum();
                  PhiOne[m][f] = Phi0One[f];
               }
            }
         }
      }
      for(int m = 0; m < 16; m++){
         VMean[m] = VOne[m] / OnePhiOne[m];
         VDelta[m] = VMean[m] - V0Mean;
      }
      for(int f = 0; f < ped.familyCount; f++){
         if(gidCount[f] < 1) continue;
         for(int m = 0; m < 16; m++){
            if(nIndex[m][f][0] == 0) // all missing in the family
               V[m][f].Zero();
            else
               for(int i = 0; i < gidCount[f]; i++){
                  score[m] += (V[m][f][i]-VMean[m] * PhiOne[m][f][i]) * Y[m][f][i];
                  V[m][f][i] -= V0Mean * PhiOne[m][f][i];
               }
         }
      }

      for(int m = 0; m < 16; m++){
         pos = b*16 + m;
         if(pos >= markerCount) break;
         if(freq[m] < 1E-100 || freq[m] > 1-1E-100) continue;

         VDV = VDV0total;
         VDOne = VDOne0;
         OneDOne = OneDOne0;
         if(NSample[m] != vidCount) // missing at this SNP
         for(int f1 = 0; f1 < ped.familyCount; f1++){
            if(mIndex[m][f1][0]==0 || gidCount[f1]==0) continue; // missingness in f1
            if(gidCount[f1]==1){ // Population-based data, f1 missing
               vid1 = gid[f1][0];
               VDV -= D.data[vid1]->data[vid1];
               VDOne -= V0[f1][0] * Phi0One[f1][0];
               OneDOne -= Phi0One[f1][0] * Phi0One[f1][0];
               for(int f2 = f1+1; f2 < ped.familyCount; f2++){  // f1 < f2
                  factor = mIndex[m][f2][0]? 1: 2;
                  for(int j = 0; j < gidCount[f2]; j++){
                     vid2 = gid[f2][j];
                     VDV -= factor * D.data[vid2]->data[vid1];
                     VDOne -= factor * V0[f1][0] * D.data[vid1]->data[vid2] * Phi0One[f2][j];
                     OneDOne -= factor * Phi0One[f1][0] * D.data[vid1]->data[vid2] * Phi0One[f2][j];
                  }
               }
               for(int f2 = 0; f2 < f1; f2++){ // f1 > f2
                  factor = mIndex[m][f2][0]? 1: 2;
                  for(int j = 0; j < gidCount[f2]; j++){
                     vid2 = gid[f2][j];
                     VDV -= factor * D.data[vid1]->data[vid2];
                     VDOne -= factor * V0[f1][0] * D.data[vid2]->data[vid1] * Phi0One[f2][j];
                     OneDOne -= factor * Phi0One[f1][0] * D.data[vid2]->data[vid1] * Phi0One[f2][j];
                  }
               }
            }else{ // family
               for(int i = 0; i < gidCount[f1]; i++){
                  vid1 = gid[f1][i];
                  VDV -= D.data[vid1]->data[vid1];
                  VDOne -= V0[f1][i] * Phi0One[f1][i];
                  OneDOne -= Phi0One[f1][i] * Phi0One[f1][i];
                  for(int j = 0; j < i; j++){
                     vid2 = gid[f1][j];
                     VDV -= 2*D.data[vid1]->data[vid2];
                     VDOne -= V0[f1][i] * D.data[vid2]->data[vid1] * Phi0One[f1][j];
                     OneDOne -= 2 * Phi0One[f1][i] * D.data[vid2]->data[vid1] * Phi0One[f1][j];
                  }
                  for(int j = i+1; j < gidCount[f1]; j++){
                     vid2 = gid[f1][j];
                     VDOne -= V0[f1][i] * D.data[vid1]->data[vid2] * Phi0One[f1][j];
                  }
               }
               for(int i = 0; i < nIndex[m][f1][0]; i++){
                  id1 = nIndex[m][f1][i+1];
                  vid1 = gid[f1][id1];
                  VDV += V[m][f1][id1] * V[m][f1][id1];
                  VDOne += V[m][f1][id1] * PhiOne[m][f1][id1];
                  OneDOne += PhiOne[m][f1][id1] * PhiOne[m][f1][id1];
                  for(int j = 0; j < i; j++){
                     id2 = nIndex[m][f1][j+1];
                     vid2 = gid[f1][id2];
                     VDV += 2 * V[m][f1][id1] * D.data[vid2]->data[vid1] * V[m][f1][id2];
                     VDOne += V[m][f1][id1] * D.data[vid2]->data[vid1] * PhiOne[m][f1][id2];
                     OneDOne += 2 * PhiOne[m][f1][id1] * D.data[vid2]->data[vid1] * PhiOne[m][f1][id2];
                  }
                  for(int j = i+1; j < nIndex[m][f1][0]; j++){
                     id2 = nIndex[m][f1][j+1];
                     vid2 = gid[f1][id2];
                     VDOne += V[m][f1][id1] * D.data[vid1]->data[vid2] * PhiOne[m][f1][id2];
                  }
               }
               for(int f2 = f1+1; f2 < ped.familyCount; f2++){  // f1 < f2
                  factor = mIndex[m][f2][0]? 1: 2;
                  for(int i = 0; i < gidCount[f1]; i++){
                     vid1 = gid[f1][i];
                     for(int j = 0; j < gidCount[f2]; j++){
                        vid2 = gid[f2][j];
                        VDV -= factor * D.data[vid2]->data[vid1];
                        VDOne -= factor * V0[f1][i] * D.data[vid1]->data[vid2] * Phi0One[f2][j];
                        OneDOne -= factor * Phi0One[f1][i] * D.data[vid1]->data[vid2] * Phi0One[f2][j];
                     }
                  }
                  for(int i = 0; i < nIndex[m][f1][0]; i++){
                     id1 = nIndex[m][f1][i+1];
                     vid1 = gid[f1][id1];
                     for(int j = 0; j < nIndex[m][f2][0]; j++){
                        id2 = nIndex[m][f2][j+1];
                        vid2 = gid[f2][id2];
                        VDV += factor * V[m][f1][id1] * D.data[vid1]->data[vid2] * V[m][f2][id2];
                        VDOne += factor * V[m][f1][id1] * D.data[vid1]->data[vid2] * PhiOne[m][f2][id2];
                        OneDOne += factor * PhiOne[m][f1][id1] * D.data[vid1]->data[vid2] * PhiOne[m][f2][id2];
                     }
                  }
               }
               for(int f2 = 0; f2 < f1; f2++){ // f1 > f2
                  factor = mIndex[m][f2][0]? 1: 2;
                  for(int i = 0; i < gidCount[f1]; i++){
                     vid1 = gid[f1][i];
                     for(int j = 0; j < gidCount[f2]; j++){
                        vid2 = gid[f2][j];
                        VDV -= factor * D.data[vid1]->data[vid2];
                        VDOne -= factor * V0[f1][i] * D.data[vid2]->data[vid1] * Phi0One[f2][j];
                        OneDOne -= factor * Phi0One[f1][i] * D.data[vid2]->data[vid1] * Phi0One[f2][j];
                     }
                  }
                  for(int i = 0; i < nIndex[m][f1][0]; i++){
                     id1 = nIndex[m][f1][i+1];
                     vid1 = gid[f1][id1];
                     for(int j = 0; j < nIndex[m][f2][0]; j++){
                        id2 = nIndex[m][f2][j+1];
                        vid2 = gid[f2][id2];
                        VDV += factor * V[m][f1][id1] * D.data[vid2]->data[vid1] * V[m][f2][id2];
                        VDOne += factor * V[m][f1][id1] * D.data[vid2]->data[vid1] * PhiOne[m][f2][id2];
                        OneDOne += factor * PhiOne[m][f1][id1] * D.data[vid2]->data[vid1] * PhiOne[m][f2][id2];
                     }
                  }
               }
            }
         }
         double newVDV = VDV + 2 * VDOne * VDelta[m] + OneDOne * VDelta[m] * VDelta[m];
         double stat = score[m] * score[m] / (0.5 * freq[m] * (1-freq[m]) * newVDV);

         if(slowFlag){ // debugging purpose
            for(int f = 0; f < ped.familyCount; f++)
               for(int i = 0; i < gidCount[f]; i++)
                  V[m][f][i] += (V0Mean-VMean[m]) * PhiOne[m][f][i];
            VDV = 0;
            for(int f1=0; f1 < ped.familyCount; f1++)
            for(int i = 0; i < nIndex[m][f1][0]; i++){
               id1 = nIndex[m][f1][i+1];
               vid1 = gid[f1][id1];
               for(int f2=0; f2 < ped.familyCount; f2++)
               for(int j = 0; j < nIndex[m][f2][0]; j++){
                  id2 = nIndex[m][f2][j+1];
                  vid2 = gid[f2][id2];
                  if(vid1 == vid2)
                     VDV += V[m][f1][id1] * V[m][f2][id2];
                  else if(vid1 < vid2)
                     VDV += V[m][f1][id1] * D.data[vid1]->data[vid2] * V[m][f2][id2];
                  else
                     VDV += V[m][f1][id1] * D.data[vid2]->data[vid1] * V[m][f2][id2];
               }
            }
            stat = score[m] * score[m] / (0.5 * freq[m] * (1-freq[m]) * VDV);
         }
         double pvalue;
         if(stat==_NAN_ || stat < 0){
            pvalue = 1;
         }else{
            pvalue = chidist(stat, 1);
            chisq.Push(stat);
            if(score[m] < 0)
               stat = -stat;
         }
         if(snpName.Length())
            fprintf(fp, "%s", (const char*)snpName[pos]);
         else
            fprintf(fp, "SNP%d", pos+1);
         if(chromosomes.Length())
            fprintf(fp, "\t%d", chromosomes[pos]);
         if(bp.Length())
            fprintf(fp, "\t%.6lf", bp[pos]*0.000001);
         fprintf(fp, "\t%s\t%s\t%.3lf\t%d",
            (const char*)alleleLabel[0][pos], (const char*)alleleLabel[1][pos],
            freq[m], NSample[m]);
         fprintf(fp, "\t%.3lf\t%.2G\n", stat, pvalue);
      }
   }
   fclose(fp);
   printf("ROADTRIPS genome scan ends at %s", currentTime());
   printf("Genome scan results saved in file %s\n\n", (const char*)lmmfile);
   idx.Index(chisq);
   printf("GC lambda for the genome scan of %s in %d autosome SNPS is %.3lf\n",
         (const char*)ped.affectionNames[0],
         chisq.Length(), chisq[idx[chisq.Length()/2]]/.456);
   printf("\n");

   if(gid) delete []gid;
   if(RTid) delete []RTid;
   if(V0) delete []V0;
   if(PhiInv) delete []PhiInv;
   if(Phi0One) delete []Phi0One;
   if(PhiM0) delete []PhiM0;
   if(PhiN0) delete []PhiN0;
   if(AN0) delete []AN0;
   if(AM0) delete []AM0;
   if(tV) delete []tV;
   for(int m = 0; m < 16; m++){
      delete []mIndex[m];
      delete []nIndex[m];
      delete []V[m];
      delete []Y[m];
      delete []PhiOne[m];
   }
}


