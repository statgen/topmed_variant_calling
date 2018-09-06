////////////////////////////////////////////////////////////////////// 
// famSKAT.cpp
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
// Oct 24, 2017

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

void Engine::SKAT_family_quantitatives(const char* skatfile)
{
   printf("\nOptions in effect:\n");
   printf("\t--skat %s\n", (const char*)skatfile);
   if(KernelShort)
      printf("\t--kernel %c\n", KernelShort);
   if(normalization)
      printf("\t--invnorm\n");
   if(traitList.Length()){
      printf("\t--trait");
      for(int t = 0; t < traits.Length(); t++)
         printf(" %s", (const char*)ped.traitNames[traits[t]]);
      printf("\n");
    }
    if(covariates.Length()){
      printf("\t--covariate");
      for(int c = 0; c < covariates.Length(); c++)
         printf(" %s", (const char*)ped.covariateNames[covariates[c]]);
      printf("\n");
    }
    if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
    printf("\n");

   bool isDisease = false;
   if(traitList.Length()==0 && ped.affectionCount>0)
      isDisease = true;
   else if(traitList.Length()==0 && ped.affectionCount==0) return;

   unrelatedFlag = false;
   int relatedCount = 0;
   Kinship kin;
   for(int f = 0; f < ped.familyCount; f++){
      kin.Setup(*ped.families[f]);
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         for(int j = i+1; j <= ped.families[f]->last; j++)
            if(CheckCovariates(ped[i]) && geno[i]!=-1 &&
               CheckCovariates(ped[j]) && geno[j]!=-1 &&
               kin(ped[i], ped[j])>0)
               relatedCount ++;
   }
   if(relatedCount < 5) {
      unrelatedFlag = true;
      if(relatedCount)
         printf("Samples are treated as unrelated. It's recommanded to remove the %d relative pairs.\n", relatedCount);
      else
         printf("All samples are unrelated\n");
   }else
      printf("There are %d relative pairs in the known pedigrees.\n", relatedCount);

   if(normalization)
      for(int t = 0; t < traits.Length(); t++)
         InvNorm(traits[t]);

   IFILE input = ifopen(skatfile, "rt");
   if(input == NULL)
      error("SKAT parameter file %s cannot be openned.", skatfile);

   StringArray skatGene(0);
   IntArray *skatSNP;
   Vector *skatWeight = NULL;
   bool WeightSpecified = false;
   StringIntHash markerLookup;
   for(int i = 0; i < markerCount; i++)
      markerLookup.SetInteger(snpName[i], i);

   String line;
   StringArray tokens;
   tokens.Clear();
   line.ReadLine(input);
   tokens.AddTokens(line);
   if(tokens.Length() < 2)
      error("Format for SKAT parameter file: GeneName SNP");
   if(tokens.Length() > 2 && (tokens[2] == "Weight" || tokens[2] == "weight" || tokens[2] == "WEIGHT"
      || tokens[2] == "Beta" || tokens[2] == "beta" || tokens[2] == "BETA") )
      WeightSpecified = true;
   else if(tokens[0] == "GeneName" || tokens[0]=="genename" || tokens[0]=="GENENAME"
      || tokens[0]=="Gene" || tokens[0]=="GENE" || tokens[0]=="gene"); // skip
   else ifrewind(input);
   while(!ifeof(input)){
      tokens.Clear();
      line.ReadLine(input);
      tokens.AddTokens(line);
      if(tokens.Length() < 2) continue;
      if(!WeightSpecified)
         if(tokens.Length() > 2 && tokens[2].IsNumber()) WeightSpecified = true;
      int g = skatGene.Find(tokens[0]);
      if(g == -1){
         g = skatGene.Length();
         skatGene.Push(tokens[0]);
      }
   }
   int skatGeneCount = skatGene.Length();
   skatSNP = new IntArray[skatGeneCount];
   for(int i = 0; i < skatGeneCount; i++)
      skatSNP[i].Dimension(0);
   if(WeightSpecified) {
      skatWeight = new Vector[skatGeneCount];
      for(int i = 0; i < skatGeneCount; i++)
         skatWeight[i].Dimension(0);
   }
   ifrewind(input);
   while(!ifeof(input)){
      tokens.Clear();
      line.ReadLine(input);
      tokens.AddTokens(line);
      if(tokens.Length() < 2) continue;
      int g = skatGene.Find(tokens[0]);
      if(g==-1) continue;  // header
//      int k = snpName.Find(tokens[1]);
      int k = markerLookup.Integer(tokens[1]);
      if(k!=-1){
         skatSNP[g].Push(k);
         if(WeightSpecified) skatWeight[g].Push(double(tokens[2]));
      }else
         printf("SNP name %s for gene %s does not exist.\n",
               (const char *)tokens[1], (const char *)tokens[0]);
   }
   ifclose(input);
   switch(KernelShort){
      case 'D':
         printf("SKAT with Euclidean distance kernel.\n");
         break;
      case 'G':
         printf("SKAT with Euclidean distance kernel adjusting for global ancestry.\n");
         break;
      case 'A':
         printf("SKAT with Euclidean distance kernel adjusting for local ancestry.\n");
         break;
      default:
         printf("SKAT with weighted linear kernel.\n");
   };

#ifdef WITH_LAPACK
   printf("  LAPACK is used.\n");
#endif

   String skatout = prefix;
   skatout.Add("skat.tbl");
   FILE *fp = fopen((const char*)skatout, "wt");
   if(fp==NULL) error("Cannot open %s to write.", (const char*)skatout);

   if(!isDisease && (traitList.Length()==0)) return;
   if(isDisease || (traitList.Length()==1) ){
      if(skatWeight){
         if(isDisease)
            fprintf(fp, "Gene\tChr\tStart\tStop\tSNPCount\tD_UU\tD_AU\tD_AA\tH2_PC1\tP_PC1\tQ\tP\n");
         else
            fprintf(fp, "Gene\tChr\tStart\tStop\tSNPCount\tH2_PC1\tP_PC1\tQ\tP\n");
      }else{
         if(isDisease)
            fprintf(fp, "Gene\tChr\tStart\tStop\tSNPCount\tRareCount\tD_UU\tD_AU\tD_AA\tH2_PC1\tP_PC1\tQ\tP\n");
         else
            fprintf(fp, "Gene\tChr\tStart\tStop\tSNPCount\tRareCount\tH2_PC1\tP_PC1\tQ\tP\n");
      }
   }else
      fprintf(fp, "Gene\tChr\tStart\tStop\tTraitName\tN\tSNPCount\tRareCount\tH2_PC1\tP_PC1\tQ\tP\n");

   printf("Genome scan starts at %s", currentTime());

   if(isDisease){
      chisqs = new Vector[1];
      chisqs[0].Dimension(0);
      printf("Model: %s ~ ",(const char*)ped.affectionNames[0]);
      if(covariates.Length()==0) printf("1\n");
      else{
         printf("%s", (const char*)ped.covariateNames[covariates[0]]);
         if(covariates.Length()>1)
            for(int c = 1; c < covariates.Length(); c++)
               printf("+%s", (const char*)ped.covariateNames[covariates[c]]);
         printf("\n");
      }
      SKAT_family_quantitative(-1, skatGene, skatSNP, skatWeight, fp);
   }else{
      chisqs = new Vector[traits.Length()];
      for(int t = 0; t < traits.Length(); t++)
         chisqs[t].Dimension(0);
      for(int t = 0; t < traits.Length(); t++){
         printf("Model: %s ~ ",(const char*)ped.traitNames[traits[t]]);
         if(covariates.Length()==0) printf("1\n");
         else{
            printf("%s", (const char*)ped.covariateNames[covariates[0]]);
            if(covariates.Length()>1)
            for(int c = 1; c < covariates.Length(); c++)
               printf("+%s", (const char*)ped.covariateNames[covariates[c]]);
            printf("\n");
         }
         SKAT_family_quantitative(t, skatGene, skatSNP, skatWeight, fp);
      }
   }

   fclose(fp);
   printf("Genome scan ends at %s", currentTime());
   printf("\nSKAT results saved in file %s\n", (const char*)skatout);
   delete []skatSNP;
   if(skatWeight) delete []skatWeight;
   else {
      if(isDisease) ;// to be implemented for disease
      else PostVC();
   }
}

void Engine::SKAT_family_quantitative(int trait, StringArray & skatGene, IntArray *skatSNP, Vector *skatWeight, FILE *fp)
{
/*
   bool isDisease = false;
   if(trait == -1) isDisease = true;
   if(geno.Length()==0) BuildBinary();

   int testCount = skatGene.Length();
   if(testCount==0) {
      error("Format for SKAT parameter file: GeneName SNPName");
      return;
   }
   Kinship kin;
   SVD svd;
   unsigned char AA, Aa, missing;
   int b, id;
   int nid = 0, NA=0, NU=0;
   if(isDisease){
      for(int i = 0; i < ped.count; i++)
         if(ped[i].affections[0]!=0 && CheckCovariates(ped[i]) && geno[i]!=-1){
            nid++;
            if(ped[i].affections[0]==2) NA++;
            else if(ped[i].affections[0]==1) NU++;
         }
      printf("The dichotomous trait is %s (N=%d, %d affected)\n",
         (const char*)ped.affectionNames[0], nid, NA);
   }else{   // quantitative trait
      for(int i = 0; i < ped.count; i++)
         if(ped[i].traits[traits[trait]]!=_NAN_ && CheckCovariates(ped[i]) && geno[i]!=-1)
            nid++;
      printf("The quantitative trait is %s (N=%d)\n", (const char*)ped.traitNames[traits[trait]], nid);
   }

   int *vid = new int[nid];
   Vector adjResid(nid);
   if(isDisease) ResidualDisease(vid, adjResid);
   else ResidualQuantitative(trait, vid, adjResid);

   double mResid = adjResid.Sum() / nid;
   double numerator, denominator;
   Matrix Omega;
   IntArray ids;
   Cholesky chol;

   double WS[256][255], WM[256][255];
   Vector weight;
   Matrix Geno;
   int pos;
   Matrix gDist;
   if(KernelShort=='G'){ // for global ancestry adjusted
      printf("Calculating pair-wise distances among %d individuals...", nid);
      gDist.Dimension(nid, nid);
      gDist.Zero();
      int HetHetCount, IBS0Count, het1Count, het2Count, notMissingCount, HomHomCount;
      int id1, id2;
      for(int i = 0; i < nid; i++){
         id1 = geno[vid[i]];
         for(int j = i+1; j < nid; j++){
            id2 = geno[vid[j]];
            HetHetCount = IBS0Count = het1Count = het2Count = HomHomCount = 0;
            for(int m = 0; m < byteCount; m++){
               HetHetCount += oneCount[(~G[0][id1][m]) & (G[1][id1][m]) & (~G[0][id2][m]) & G[1][id2][m]];
               IBS0Count += oneCount[G[0][id1][m] & G[0][id2][m] & (G[1][id1][m] ^ G[1][id2][m])];
               HomHomCount += oneCount[G[0][id1][m] & G[0][id2][m] & ~(G[1][id1][m]^G[1][id2][m])];
               het1Count += oneCount[(G[0][id2][m] | G[1][id2][m]) & (~G[0][id1][m]) & G[1][id1][m]];
               het2Count += oneCount[(G[0][id1][m] | G[1][id1][m]) & (~G[0][id2][m]) & G[1][id2][m]];
            }
            notMissingCount = HomHomCount + IBS0Count + het1Count + het2Count - HetHetCount;
            gDist[i][j] = (het1Count + het2Count - 2*HetHetCount + 4*IBS0Count) * 1.0 / notMissingCount;
         }
      }
      printf("Done\n");
   }
   IntArray vbyte, vbit, vNext;
   double pvalue;
   Vector lambda;
   Matrix kernel, kernel2;
   Vector freq(8);
   IntArray fCount(8);
   double D[3];
   Vector GA, GU;
   Vector PC1(nid);
   double H2_PC1, P_PC1, mPC;

   for(int t = 0; t < testCount; t++){ // loop over genes/regions
      int skatSNPCount = skatSNP[t].Length();
      if(skatSNPCount==0) continue;
      if(KernelShort=='L' || KernelShort=='l'){
         Geno.Dimension(nid, skatSNPCount);
         Geno.Zero();
      }
      weight.Dimension(skatSNPCount);
      vbyte.Dimension(0);
      vbit.Dimension(0);
      for(int k = 0; k < skatSNPCount; k++){
         int m = skatSNP[t][k];
         int proposedByte = m/8;
         int proposedBit = m%8;
         if(vbyte.Find(proposedByte)>-1){ // in previous byte
            vbit[vbyte.Find(proposedByte)] |= (1<<proposedBit);
         }else{   // a new byte
            vbyte.Push(proposedByte);
            vbit.Push(1<<proposedBit);
         }
      }
      int cIndex = 0;
      int rareCount = 0;
      for(int v = 0; v < vbyte.Length(); v++){ // loop over vbyte
         b = vbyte[v];
         vNext.Dimension(0);
         for(int j = 0; j < 8; j++)
            if(vbit[v] & base[j]) vNext.Push(j);
         freq.Zero();
         fCount.Zero();
         for(int i = 0; i < nid; i++){
            id = geno[vid[i]];
            AA = G[0][id][b] & G[1][id][b];
            Aa = (~G[0][id][b]) & G[1][id][b];
            missing = (~G[0][id][b]) & (~G[1][id][b]);
            if(skatWeight==NULL) // Need to compute the SKAT weight
               for(int j = 0; j < vNext.Length(); j++){
                  if(!(missing & base[vNext[j]])){
                     fCount[j] += 2;
                     if(AA & base[vNext[j]])
                        freq[j] += 2;
                     else if(Aa & base[vNext[j]])
                        freq[j] ++;
                  }
               }
            if(KernelShort=='L' || KernelShort=='l') // Genotypes are needed
               for(int j = 0; j < vNext.Length(); j++){
                  pos = cIndex+j;
                  if(AA & base[vNext[j]])
                     Geno[i][pos] = 2;
                  else if(Aa & base[vNext[j]])
                     Geno[i][pos] = 1;
                  else if(missing & base[vNext[j]])
                     Geno[i][pos] = -1;
               }
         }
         if(KernelShort=='D' || KernelShort=='A' || KernelShort=='G')
            for(int i = 0; i < 256; i++)
               WS[i][v] = WM[i][v] = 0;
         for(int j = 0; j < vNext.Length(); j++){
            pos = cIndex+j;
            if(skatWeight != NULL) // Weight pre-specified
               weight[pos] = skatWeight[t][pos];
            else{ // Weight to be computed
               if(fCount[j]) {
                  freq[j] /= fCount[j];
                  if(freq[j] < 0.5)
                     weight[pos] = pow(1-freq[j], 24)*25; // sqrt of weight
                  else
                     weight[pos] = pow(freq[j], 24)*25;
               }else // missing
                  freq[j] = 0;
               if(freq[j]<0.05 || freq[j]>0.95) rareCount++;
            }
            if(KernelShort=='D' || KernelShort=='A' || KernelShort=='G'){
               for(int i = 0; i < 256; i++)
                  if(i & base[vNext[j]]){
                     WS[i][v] += weight[pos]*weight[pos];
                     WM[i][v] += weight[pos]*weight[pos]*4*freq[j]*(1-freq[j]);
                  }
            }else{  // 'L' or 'l'
               for(int i = 0; i < nid; i++)
                  if(Geno[i][pos] == -1)
                     Geno[i][pos] = 0;
                  else
                     Geno[i][pos] -= (freq[j]*2);
            }
         } // end of vNext loop
         cIndex += vNext.Length();
      } // loop over vbyte

      if(skatWeight==NULL && rareCount==0) continue; // no rare variants
      double Q = 0;
      if(KernelShort=='D' || KernelShort=='A' || KernelShort=='G'){
         kernel.Dimension(nid, nid);
         kernel.Zero();
         int id1, id2;
         double tD;
         for(int i = 0; i < nid; i++){
            id1 = geno[vid[i]];
            for(int j = i+1; j < nid; j++){
               id2 = geno[vid[j]];
               tD = 0;
               for(int v = 0; v < vbyte.Length(); v++){
                  b = vbyte[v];
                  missing = ~((G[0][id1][b] | G[1][id1][b]) & (G[0][id2][b] | G[1][id2][b]));
                  tD +=
                  4*WS[G[0][id1][b] & G[0][id2][b] & (G[1][id1][b] ^ G[1][id2][b])][v]
                  + WS[(G[0][id2][b] & (~G[0][id1][b]) & G[1][id1][b]) |
                  (G[0][id1][b] & (~G[0][id2][b]) & G[1][id2][b])][v]
                  + WM[missing][v];         // Dist = IBS0*4 + IBS1
               }
               kernel.data[i]->data[j] = kernel.data[j]->data[i] = -tD;
            }
         }
         if(KernelShort=='D' && (!unrelatedFlag)){ // family adjusted
            int baseN=0;
            for(int f = 0; f < ped.familyCount; f++){
               ids.Dimension(0);
               if(isDisease){
                  for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
                     if(ped[i].affections[0] !=_NAN_ && CheckCovariates(ped[i]) && geno[i]!=-1)
                        ids.Push(i);
               }else{
                  for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
                     if(ped[i].traits[traits[trait]]!=_NAN_ && CheckCovariates(ped[i]) && geno[i]!=-1)
                        ids.Push(i);
               }
               int N = ids.Length();
               if(N==0) continue;
               else if(N==1){
                  baseN++;
                  continue;
               }
               kin.Setup(*ped.families[f]);
               for(int i = 0; i < N; i++)
                  for(int j = i+1; j < N; j++)
                     kernel[baseN+i][baseN+j] = kernel[baseN+j][baseN+i]
                     = kernel[baseN+i][baseN+j] / (1-2*kin(ped[ids[i]],ped[ids[j]]));
               baseN += N;
            }
         }else if(KernelShort=='G'){  // Global ancestry adjusted
            for(int i = 0; i < nid; i++)
               for(int j = i+1; j < nid; j++)
                  kernel[i][j] = kernel[j][i] = kernel[i][j] / gDist[i][j];
         }else if(KernelShort=='A'){  // admixture/local ancestry adjusted
         }

         if(isDisease){
            D[0]=D[1]=D[2]=0.0;
            for(int i = 0; i < nid; i++)
               for(int j = i+1; j < nid; j++){
                  if(adjResid[i] < 0 && adjResid[j] < 0) // UU
                     D[0] -= kernel[i][j];
                  else if(adjResid[i] > 0 && adjResid[j] > 0) // DD
                     D[2] -= kernel[i][j];
                  else // DU
                     D[1] -= kernel[i][j];
               }
            D[0] /= (0.5*NU*(NU-1)*skatSNPCount);
            D[2] /= (0.5*NA*(NA-1)*skatSNPCount);
            D[1] /= (1.0*NA*NU*skatSNPCount);
         }
         
         Vector mean(nid);
         mean.Zero();
         for(int i = 0; i < nid; i++){
            for(int j = 0; j < nid; j++)
               mean[i] += kernel[i][j];
            mean[i] /= nid;
         }
         double mean0 = mean.Sum() / nid;
         for(int i = 0; i < nid; i++)
            for(int j = 0; j < nid; j++)
               kernel[i][j] -= (mean[i] + mean[j] - mean0);
         Q = 0.0;
         for(int i = 0; i < nid; i++)
            for(int j = 0; j < nid; j++)
               Q += adjResid[i] * kernel[i][j] * adjResid[j];
         pvalue = ComputeDaviesP(Q, kernel);
      }else{ // Linear kernel: K = Omega^-1/2 G W
         kernel.Dimension(nid, skatSNPCount);
         if(unrelatedFlag || isDisease) // unrelated data
            for(int j = 0; j < skatSNPCount; j++)
               for(int i = 0; i < nid; i++)
                  kernel.data[i]->data[j] = Geno.data[i]->data[j] * weight.data[j];
         else{    // quantitative trait in family 
            int baseN = 0;
            for(int f = 0; f < ped.familyCount; f++){
               ids.Dimension(0);
               if(isDisease){
                  for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
                     if(ped[i].affections[0]!=_NAN_ && CheckCovariates(ped[i]) && geno[i]!=-1)
                        ids.Push(i);
               }else{
                  for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
                     if(ped[i].traits[traits[trait]]!=_NAN_ && CheckCovariates(ped[i]) && geno[i]!=-1)
                        ids.Push(i);
               }
               int N = ids.Length();
               if(N==0) continue;
               kin.Setup(*ped.families[f]);
               Omega.Dimension(N, N);
               if(isDisease || KernelShort=='l') // Kinship matrix
                  for(int i = 0; i < N; i++){
                     Omega[i][i] = 1;
                     for(int j = i+1; j < N; j++)
                        Omega[i][j] = Omega[j][i] = 2*kin(ped[ids[i]], ped[ids[j]]);
                  }
               else  // variance covariance matrix
                  for(int i = 0; i < N; i++){
                     Omega[i][i] = variances[trait];
                     for(int j = i+1; j < N; j++) // Upper triangle matrix
                        Omega[i][j] = Omega[j][i] = 2*kin(ped[ids[i]], ped[ids[j]])*heritabilities[trait]*variances[trait];
                  }
               chol.Decompose(Omega);

               Vector temp(N);
               for(int j = 0; j < skatSNPCount; j++){
                  for(int i = 0; i < N; i++)
                     temp[i] = Geno[baseN+i][j] * weight[j];
                  chol.BackSubst0(temp);
                  for(int i = 0; i < N; i++)
                     kernel[baseN+i][j] = chol.x[i];
               }
               baseN += N;
            } // end of family loop
         } // end of unrelatedFlag
         if(isDisease && (KernelShort=='L' || KernelShort=='l')){
            GA.Dimension(skatSNPCount);
            GU.Dimension(skatSNPCount);
            GA.Zero(); GU.Zero();
            for(int i = 0; i < nid; i++)
               if(adjResid[i] > 0)
                  for(int j = 0; j < skatSNPCount; j++)
                     GA[j] += kernel.data[i]->data[j];
               else
                  for(int j = 0; j < skatSNPCount; j++)
                     GU[j] += kernel.data[i]->data[j];
            D[0]=D[1]=D[2]=0.0;
            for(int j = 0; j < skatSNPCount; j++){ // GG' (=-HDH/2)
               D[0] += GU[j] * GU[j];
               D[1] += GU[j] * GA[j] * 2;
               D[2] += GA[j] * GA[j] ;
            }
            GA.Zero(); GU.Zero();
            for(int i = 0; i < nid; i++)
               if(adjResid.data[i] > 0)
                  for(int j = 0; j < skatSNPCount; j++)
                     GA.data[j] += kernel.data[i]->data[j]*kernel.data[i]->data[j];
               else
                  for(int j = 0; j < skatSNPCount; j++)
                     GU.data[j] += kernel.data[i]->data[j]*kernel.data[i]->data[j];
            for(int j = 0; j < skatSNPCount; j++){ // GG' - offset /2 (=-(HDH+offset)/2=-D/2)
               D[0] -= NU*GU[j];
               D[1] -= (NA*GU[j] + NU*GA[j]);
               D[2] -= NA*GA[j];
            }
            D[0] /= -(0.5 * NU * (NU-1) * skatSNPCount);
            D[2] /= -(0.5 * NA * (NA-1) * skatSNPCount);
            D[1] /= -(1.0 * NA * NU * skatSNPCount);
         }
         Q = 0.0;
         for(int j = 0; j < skatSNPCount; j++){
            double temp = 0;
            for(int i = 0; i < nid; i++)
               temp += kernel[i][j] * adjResid[i];
            Q += temp*temp;
         }
         
//         if(effectFlag || skatSNPCount > nid){  // Dimension: N x M
            lambda.Dimension(0);
#ifdef WITH_LAPACK
//            int dimM = nid;
//            int dimN = skatSNPCount;
            char JOBU = 'S';
            char JOBVT = 'O';
            int LDU = nid;
            int LDVT = 1;
            int dimS = nid>skatSNPCount? skatSNPCount: nid;
            int info;
            double *A = (double *)malloc(nid*skatSNPCount*sizeof(double));
            for(int i = 0; i < nid; i++)
               for(int j = 0; j < skatSNPCount; j++)
                  A[j*nid+i] = kernel[i][j];   
            double *S = (double *)malloc(dimS*sizeof(double));
            double *U = (double *)malloc(LDU*dimS*sizeof(double));
            double *VT = (double *)malloc(LDVT*sizeof(double));
            int LWORK = nid>skatSNPCount? nid*2+8*skatSNPCount: skatSNPCount+8*nid;
            double *WORK = (double *)malloc(LWORK*sizeof(double));
            dgesvd_(&JOBU, &JOBVT, &nid, &skatSNPCount, A, &nid, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &info);
            if(info!=0) error("SVD failed.");
            delete WORK;
            delete A;
            for(int i = 0; i < dimS; i++)
               if(S[i]*S[i] > 1E-100) lambda.Push(S[i]*S[i]);
            for(int i = 0; i < nid; i++)
               PC1[i] = U[i];
            delete VT;
            delete S;
            delete U;
#else
            svd.Decompose(kernel);
            if(svd.n == 0) return;
            for(int i = 0; i < (nid>skatSNPCount? skatSNPCount: nid); i++)
               if(svd.w[i]*svd.w[i] > 1E-100) lambda.Push(svd.w[i]*svd.w[i]);
            QuickIndex idx;
            idx.Index(svd.w);
            for(int i = 0; i < nid; i++)
               PC1[i] = svd.v[i][idx[skatSNPCount-1]];
#endif
            pvalue = lambda.Length()>0? 1-Davies(Q, lambda.data, lambda.Length()): 1;
            if(pvalue > 1) pvalue = 1; //lambda.Print(20);

            mPC = PC1.Sum() / nid;
            numerator = 0;
            for(int i = 0; i < nid; i++)
               numerator += PC1[i] * adjResid[i];
            numerator -= nid*mResid*mPC;
            denominator = 1 - nid*mPC*mPC;
            P_PC1 = ndist(fabs(numerator/sqrt(denominator)))*2;
            H2_PC1 = numerator / denominator;
            H2_PC1 = H2_PC1 * H2_PC1 * (1.0/nid - mPC * mPC);

      }
      if(pvalue<0) pvalue=0;

      fprintf(fp, "%s\t%d\t%d\t%d",
         (const char*)skatGene[t], chromosomes[skatSNP[t][0]],
         int(positions[skatSNP[t][0]]*1000000+0.5),
         int(positions[skatSNP[t][skatSNPCount-1]]*1000000+0.5));
      if((!isDisease) && traits.Length()>1) // multiple traits
            fprintf(fp, "\t%s\t%d", (const char*)ped.traitNames[traits[trait]], nid);
      fprintf(fp, "\t%d", skatSNPCount);
      if(skatWeight==NULL)
         fprintf(fp, "\t%d", rareCount);
      if(isDisease)
         fprintf(fp, "\t%.2lf\t%.2lf\t%.2lf", D[0], D[1], D[2]);
      fprintf(fp, "\t%.3lf\t%.3G", H2_PC1, P_PC1);
      fprintf(fp, "\t%.1lf\t%.3G\n", Q, pvalue);
      
      if(skatWeight==NULL){ // SKAT
         if(isDisease){   // GC function not available yet
         }else{
            if(pvalue < 1E-9) chisqs[trait].Push(99999999); // min cutoff for niv is 2E-10
            else if(pvalue > 0.99999) chisqs[trait].Push(0);
            else{
               double t = ninv(pvalue*0.5);
               chisqs[trait].Push(t*t);
            }
         }
      }
   } // loop over genes/regions

   delete []vid;
*/
}
/*         }else{   // kernel = G'G, Dimension: M x M
            if(skatSNPCount == 0) continue;
            else if(skatSNPCount == 1){
               double s = 0;
               for(int i = 0; i < nid; i++)
                  s += kernel[i][0] * kernel[i][0];
               pvalue = chidist(Q/s, 1);
            }else{
               kernel2.Dimension(skatSNPCount, skatSNPCount);
               kernel2.Zero();
               for(int i = 0; i < skatSNPCount; i++)
                  for(int j = i+1; j < skatSNPCount; j++){
                     for(int k = 0; k < nid; k++)
                        kernel2.data[i]->data[j] += kernel.data[k]->data[i] * kernel.data[k]->data[j];
                     kernel2[j][i] = kernel2[i][j];
                  }
               for(int i = 0; i < skatSNPCount; i++)
                  for(int k = 0; k < nid; k++)
                     kernel2.data[i]->data[i] += kernel.data[k]->data[i] * kernel.data[k]->data[i];
               pvalue = ComputeDaviesP(Q, kernel2);
            }
         } // end of effectFlag
         */
/*            if(pvalue < 1E-9) chisqs[0].Push(99999999); // min cutoff for niv is 2E-10
            else if(pvalue > 0.99999) chisqs[0].Push(0);
            else{
               double t = ninv(pvalue*0.5);
               chisqs[0].Push(t*t);
            }
*/

void Engine::InvNorm(int tr)
{
   IntArray individuals(0);
   Vector myTraits(0);
   QuickIndex idx;
   for (int i = 0; i < ped.count; i++){
      int missing = 0;
      for(int k = 0; k < covariates.Length(); k++)
         if(!ped[i].isControlled(covariates[k])){
            missing = 1;
            break;
         }
      if(!missing && ped[i].isPhenotyped(traits[tr])){
         individuals.Push(i);
         myTraits.Push(ped[i].traits[traits[tr]]);
      }
   }
   idx.Index(myTraits);
   for(int i = 0; i < individuals.Length();){
      int start = i, end = i + 1;
      while(end < individuals.Length() &&
            ped[individuals[idx[end]]].traits[traits[tr]] ==
            ped[individuals[idx[start]]].traits[traits[tr]])
         end++;
      end --;
      double q = ninv((start + (end-start)/2.0 + 0.5)/individuals.Length());
      for(int j = start; j <= end; j++)
         ped[individuals[idx[j]]].traits[traits[tr]] = q;
      i = end + 1;
   }
}


/*
   double *resid = new double[nid];
   Vector qtrait(nid);
   Matrix qcov(nid, covariates.Length()+1);
   qtrait.Zero();
   qcov.Zero();
   int count = 0;
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
      if(ped[i].traits[traits[trait]]!=_NAN_ && CheckCovariates(ped[i]) && geno[i]!=-1){
         vid[count] = i;
         qtrait[count] = ped[i].traits[traits[trait]];
         qcov[count][0] = 1;
         for(int j = 0; j < covariates.Length(); j++)
            qcov[count][j+1] = ped[i].covariates[covariates[j]];
         count++;
       }

   for(int i = 0; i < nid; i++){
      resid[i] = qtrait[i] - means[trait][0];
      for(int j = 0; j < covariates.Length(); j++)
         resid[i] -= means[trait][j+1] * qcov[i][j+1];
   }
   Vector adjResid(nid);
   adjResid.Zero();
   int baseN = 0;
   Kinship kin;
   Matrix Omega;
   IntArray ids;
   Cholesky chol;
   for(int f = 0; f < ped.familyCount; f++){
      ids.Dimension(0);
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(ped[i].traits[traits[trait]]!=_NAN_ && CheckCovariates(ped[i]) && geno[i]!=-1)
            ids.Push(i);
      int N = ids.Length();
      if(N==0) continue;
      kin.Setup(*ped.families[f]);
      Omega.Dimension(N, N);
      for(int i = 0; i < N; i++){
         Omega[i][i] = variances[trait];
         for(int j = i+1; j < N; j++) 
            Omega[i][j] = Omega[j][i] = 2*kin(ped[ids[i]], ped[ids[j]])*heritabilities[trait]*variances[trait];
      }
      chol.Decompose(Omega);
      Vector temp(N);
      for(int i = 0; i < N; i++)
         temp[i] = resid[baseN+i];
      chol.BackSubst0(temp);
      for(int i = 0; i < N; i++)
         adjResid[baseN+i] = chol.x[i];
      baseN += N;
   }

   */
/*
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
               if(frequencies[pos]<0.05 || frequencies[pos]>0.95) rareCount++;

               if(KernelShort=='D' || KernelShort=='A' || KernelShort=='G')
                  for(int i = 0; i < 256; i++)
                     if(i & vbit[v] & base[j]){
                        WS[i][v] += weight[pos]*weight[pos];
                        WM[i][v] += weight[pos]*weight[pos]*4*frequencies[pos]*(1-frequencies[pos]);
                     }
               offset++;
         }
*/


