////////////////////////////////////////////////////////////////////// 
// VCR.cpp
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

#include "analysis.h"
#include "MathStats.h"
#include "Kinship.h"
#include "MathCholesky.h"
#include <math.h>
#include "VCLinear.h"

void Engine::VCR_family_quantitatives(const char* vcrfile)
{
   printf("\nOptions in effect:\n");
   printf("\t--vcr %s\n", (const char*)vcrfile);
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

   if(traitList.Length()==0) return;

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
   }//else printf("There are %d relative pairs in the known pedigrees.\n", relatedCount);

   if(normalization)
      for(int t = 0; t < traits.Length(); t++)
         InvNorm(traits[t]);

   IFILE input = ifopen(vcrfile, "rt");
   if(input == NULL)
      error("VCR parameter file %s cannot be openned.", vcrfile);

   StringArray vcrGene(0);
   IntArray *vcrSNP;
   Vector *vcrWeight = NULL;
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
      error("Format for VCR parameter file: GeneName SNP");
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
      int g = vcrGene.Find(tokens[0]);
      if(g == -1){
         g = vcrGene.Length();
         vcrGene.Push(tokens[0]);
      }
   }
   int vcrGeneCount = vcrGene.Length();
   vcrSNP = new IntArray[vcrGeneCount];
   for(int i = 0; i < vcrGeneCount; i++)
      vcrSNP[i].Dimension(0);
   if(WeightSpecified) {
      vcrWeight = new Vector[vcrGeneCount];
      for(int i = 0; i < vcrGeneCount; i++)
         vcrWeight[i].Dimension(0);
   }
   ifrewind(input);
   while(!ifeof(input)){
      tokens.Clear();
      line.ReadLine(input);
      tokens.AddTokens(line);
      if(tokens.Length() < 2) continue;
      int g = vcrGene.Find(tokens[0]);
      if(g==-1) continue;  // header
      int k = markerLookup.Integer(tokens[1]);
      if(k!=-1){
         vcrSNP[g].Push(k);
         if(WeightSpecified) vcrWeight[g].Push(double(tokens[2]));
      }else
         printf("SNP name %s for gene %s does not exist.\n",
               (const char *)tokens[1], (const char *)tokens[0]);
   }
   ifclose(input);
   String vcrout = prefix;
   vcrout.Add("vcr.tbl");
   FILE *fp = fopen((const char*)vcrout, "wt");
   if(fp==NULL) error("Cannot open %s to write.", (const char*)vcrout);

   if(traitList.Length()==1){
      if(vcrWeight){
         fprintf(fp, "Gene\tChr\tStart\tStop\tSNPCount\tBeta\tChisq\tP\n");
      }else{
         fprintf(fp, "Gene\tChr\tStart\tStop\tSNPCount\tRareCount\tNinfo_W\tBeta_W\tChisq_W\tP_W\tNinfo_B\tBeta_B\tChisq_B\tP_B\tChisq\tP\n");
      }
   }else
      fprintf(fp, "Gene\tChr\tStart\tStop\tTraitName\tN\tSNPCount\tRareCount\tBeta\tChisq\tP\n");

   if(geno.Length()==0) BuildBinary();
   printf("Genome scan starts at %s", currentTime());

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
      VCR_family_quantitative(t, vcrGene, vcrSNP, vcrWeight, fp);
   }
   fclose(fp);
   printf("Genome scan ends at %s", currentTime());
   printf("\nVCR results saved in file %s\n", (const char*)vcrout);
   delete []vcrSNP;
   if(vcrWeight) delete []vcrWeight;
   else
      PostVC();
}

void Engine::VCR_family_quantitative(int trait, StringArray & vcrGene, IntArray *vcrSNP, Vector *vcrWeight, FILE *fp)
{
/*
   int testCount = vcrGene.Length();
   if(testCount==0) {
      error("Format for VCR parameter file: GeneName SNPName");
      return;
   }
   Kinship kin;
   unsigned char AA, Aa, missing;
   int b, id;
   int nid = 0;
   for(int i = 0; i < ped.count; i++)
         if(ped[i].traits[traits[trait]]!=_NAN_ && CheckCovariates(ped[i]) && geno[i]!=-1)
            nid++;
   printf("The quantitative trait is %s (N=%d)\n", (const char*)ped.traitNames[traits[trait]], nid);

   int *vid = new int[nid];
   Matrix Omega;
   Matrix *OI = new Matrix[ped.familyCount];
   IntArray *ids = new IntArray[ped.familyCount];
   Matrix *VI = new Matrix[ped.familyCount];
   Cholesky chol;
   int count = 0;
   Vector resid;
   Vector *adjY = new Vector[ped.familyCount];
   Vector *adjResid = new Vector[ped.familyCount];

   if(unrelatedFlag){printf("Unrelated, to be implemented later.\n"); return;}
   polygenic();

   count = 0;
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(ped[i].traits[traits[trait]]!=_NAN_ && CheckCovariates(ped[i]) && geno[i]!=-1)
            vid[count++] = i;

   int maxN = 0;
   int totalPairCountW = 0;
   int effFamCount = 0;
   for(int f = 0; f < ped.familyCount; f++){
      ids[f].Dimension(0);
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(ped[i].traits[traits[trait]]!=_NAN_ && CheckCovariates(ped[i]) && geno[i]!=-1)
            ids[f].Push(i);
      int N = ids[f].Length();
      if(N<1) continue;
      resid.Dimension(N);
      for(int i = 0; i < N; i++){
         resid[i] = ped[ids[f][i]].traits[traits[trait]] - means[trait][0];
         for(int j = 0; j < covariates.Length(); j++)
            resid[i] -= means[trait][j+1] * ped[ids[f][i]].covariates[covariates[j]];
      }
      kin.Setup(*ped.families[f]);
      Omega.Dimension(N, N);
      for(int i = 0; i < N; i++){
         Omega[i][i] = variances[trait];
         for(int j = i+1; j < N; j++)
            Omega[i][j] = Omega[j][i] = 2*kin(ped[ids[f][i]], ped[ids[f][j]])*heritabilities[trait]*variances[trait];
      }
      chol.Decompose(Omega);
      chol.Invert();
      OI[f] = chol.inv;
      VI[f].Dimension(N*(N-1)/2, N*(N-1)/2);
      adjY[f].Dimension(N*(N-1)/2);
      adjY[f].Zero();
      int index1 = 0, index2 = 0;
      for(int u = 0; u < N; u++)
         for(int v = u+1; v < N; v++){
            index2 = 0;
            for(int l = 0; l < N; l++)
               for(int m = l+1; m < N; m++){
                  VI[f][index1][index2] = OI[f][u][l] * OI[f][v][m] + OI[f][u][m] * OI[f][v][l];
                  adjY[f][index1] += VI[f][index1][index2] * (resid[l]*resid[m]-Omega[l][m]);
                  index2++;
               }
            for(int j = 0; j < N; j++)
               adjY[f][index1] += OI[f][u][j] * OI[f][v][j] * (resid[j]*resid[j]-Omega[j][j]);
            index1++;
         }
      adjResid[f].Dimension(N);
      chol.BackSubst(resid);
      for(int i = 0; i < N; i++)
         adjResid[f][i] = chol.x[i];
      if(N>maxN) maxN = N;
      totalPairCountW += N*(N-1)/2;
      effFamCount ++;
   }  // end of family loop
   printf("Among %d families, the largest family has %d typed individuals.\n", effFamCount, maxN);
   printf("%d pairs of relative pairs are used in the within-family association analysis.\n", totalPairCountW);

   double WS[256][255], WM[256][255];
   Vector weight;
   int pos;
   IntArray vbyte, vbit, vNext;
   Vector kernel;
   Vector freq(8);
   IntArray fCount(8);
   Matrix Dfg(maxN, maxN);
   double Qw, Qb, Q, Tw, Tb, T, Pw, Pb, P, varQw, varQb, wBeta, bBeta;
   int informativePairCountW, informativePairCountB;

   for(int t = 0; t < testCount; t++){ // loop over genes/regions
      int vcrSNPCount = vcrSNP[t].Length();
      if(vcrSNPCount==0) continue;
      weight.Dimension(vcrSNPCount);
      vbyte.Dimension(0);
      vbit.Dimension(0);
      for(int k = 0; k < vcrSNPCount; k++){
         int m = vcrSNP[t][k];
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
            if(vcrWeight==NULL) // Need to compute the VCR weight
               for(int j = 0; j < vNext.Length(); j++){
                  if(!(missing & base[vNext[j]])){
                     fCount[j] += 2;
                     if(AA & base[vNext[j]])
                        freq[j] += 2;
                     else if(Aa & base[vNext[j]])
                        freq[j] ++;
                  }
               }
         }
         for(int i = 0; i < 256; i++)
            WS[i][v] = WM[i][v] = 0;
         for(int j = 0; j < vNext.Length(); j++){
            pos = cIndex+j;
            if(vcrWeight != NULL) // Weight pre-specified
               weight[pos] = vcrWeight[t][pos];
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
            for(int i = 0; i < 256; i++)
               if(i & base[vNext[j]]){
                  WS[i][v] += weight[pos]*weight[pos];
                  WM[i][v] += weight[pos]*weight[pos]*4*freq[j]*(1-freq[j]);
               }
         } // end of vNext loop
         cIndex += vNext.Length();
      } // loop over vbyte

      if(vcrWeight==NULL && rareCount==0) continue; // no rare variants
      int id1, id2;
      double tD;
      informativePairCountW = informativePairCountB = 0;
      Qb = varQb = Qw = varQw = 0;
      for(int f = 0; f < ped.familyCount; f++){
         int N = ids[f].Length();
         if(N < 2) continue;
         kernel.Dimension(0);
         for(int i = 0; i < N; i++){
            id1 = geno[ids[f][i]];
            for(int j = i+1; j < N; j++){
               id2 = geno[ids[f][j]];
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
               kernel.Push(tD);
            }
         }
         int count = 0;
         for(int i = 0; i < N*(N-1)/2; i++)
            if(kernel[i] > 0)
               count ++;
         if(count){
            informativePairCountW += count;
            for(int i = 0; i < N*(N-1)/2; i++)
               Qw -= kernel[i] * adjY[f][i];
            for(int i = 0; i < N*(N-1)/2; i++)
               for(int j = 0; j < N*(N-1)/2; j++)
                  varQw += kernel[i] * VI[f][i][j] * kernel[j];
         }
      }  // end of family loop
      for(int f = 0; f < ped.familyCount; f++){
         int N = ids[f].Length();
         if(N == 0) continue;
         for(int g = f+1; g < ped.familyCount; g++){
            int N2 = ids[g].Length();
            if(N2 == 0) continue;
            int count = 0;
            for(int i = 0; i < N; i++){
               id1 = geno[ids[f][i]];
               for(int j = 0; j < N2; j++){
                  id2 = geno[ids[g][j]];
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
                  Dfg[i][j] = tD;
                  if(tD>0) {
                     Qb -= tD * adjResid[f][i] * adjResid[g][j];
                     count ++;
                  }
               }
            }
            informativePairCountB += count;
            if(count)
               for(int u = 0; u < N; u++)
                  for(int v = 0; v < N2; v++)
                     for(int l = 0; l < N; l++)
                        for(int m = 0; m < N2; m++)
                           varQb += Dfg.data[u]->data[v] * Dfg.data[l]->data[m]
                           * OI[f].data[u]->data[l] * OI[g].data[v]->data[m];
         } // end of family g loop
      }  // end of family f loop
      if(informativePairCountW + informativePairCountB < 100) continue;
      if(informativePairCountW < 10){ // no within-family variation
         wBeta = Qw = Tw = varQw = 0;
         Pw = 1;
      }else{
         wBeta = Qw / varQw;
         Tw = Qw * wBeta;
         Pw = Qw > 0? ndist(sqrt(Tw)): 1;
      }
      bBeta = Qb / varQb;
      Tb = Qb * bBeta;
      Pb = Qb > 0? ndist(sqrt(Tb)): 1;

      T = (Qw + Qb) * (Qw + Qb) / (varQw + varQb);
      P = (Qw + Qb > 0)? ndist(sqrt(T)): 1;

      fprintf(fp, "%s\t%d\t%d\t%d",
         (const char*)vcrGene[t], chromosomes[vcrSNP[t][0]],
         int(positions[vcrSNP[t][0]]*1000000+0.5),
         int(positions[vcrSNP[t][vcrSNPCount-1]]*1000000+0.5));
      if(traits.Length()>1) // multiple traits
         fprintf(fp, "\t%s\t%d", (const char*)ped.traitNames[traits[trait]], nid);
      fprintf(fp, "\t%d", vcrSNPCount);
      if(vcrWeight==NULL) fprintf(fp, "\t%d", rareCount);
      fprintf(fp, "\t%d\t%.4lf\t%.2lf\t%.3G", informativePairCountW, wBeta, Tw, Pw);
      fprintf(fp, "\t%d\t%.4lf\t%.2lf\t%.3G", informativePairCountB, bBeta, Tb, Pb);
      fprintf(fp, "\t%.2lf\t%.3G\n", T, P);

      if(vcrWeight==NULL){ // VCR
         if(P < 1E-9) chisqs[trait].Push(99999999); // min cutoff for niv is 2E-10
         else if(P > 0.99999) chisqs[trait].Push(0);
         else{
            double t = ninv(P*0.5);
            chisqs[trait].Push(t*t);
         }
      }
   } // loop over genes/regions
   delete []vid;
   delete []ids;
   delete []VI;
   delete []adjY;
   delete []adjResid;
   delete []OI;
   */
}

