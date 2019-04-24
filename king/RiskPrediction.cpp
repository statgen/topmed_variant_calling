////////////////////////////////////////////////////////////// ////////
// RiskPrediction.cpp
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
// Nov 21, 2018

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

#define SNPCOL 0
#define EACOL 1
#define AFCOL 2
#define WTCOL 3
#define CHRCOL 4
#define POSCOL 5
#define OACOL 6
#define ONCOL 7

double Engine::ComputeAUC(Vector & risks, IntArray & diseaseStatus, int printFlag)
{
   const int TOTALCUT=100000;
   double sensitivity[TOTALCUT], specificity[TOTALCUT];
   int Pcount[2][TOTALCUT], Ncount[2][TOTALCUT];
   for(int t = 0; t < TOTALCUT; t++)
      Pcount[0][t] = Pcount[1][t] = Ncount[0][t] = Ncount[1][t] = 0;
   Vector temprisk[2];
   for(int a = 0; a < 2; a++)
      temprisk[a].Dimension(0);
   int totalCount = risks.Length();
   for(int i = 0; i < totalCount; i++){
      int aff = diseaseStatus[i];
      if((aff != 1 && aff != 0) || risks[i] < -0.999) continue;
      temprisk[aff].Push(printFlag==2? risks[i]: 1-risks[i]);
   }
   QuickIndex riskIdx;
   for(int a = 0; a < 2; a++){
      riskIdx.Index(temprisk[a]);
      int affCount = temprisk[a].Length();
      int oldT, newT;
      oldT=newT=0;
      for(int i = 0; i < affCount; i++){
         newT = int(temprisk[a][riskIdx[i]]*TOTALCUT);
         for(int t = oldT; t < newT; t++) Pcount[a][t] = i;
         Pcount[a][newT] = i+1;
         int k = affCount-i;
         for(int t = oldT; t < newT; t++) Ncount[a][t] = k;
         Ncount[a][newT] = k-1;
         oldT = newT+1;
      }  // end of loop i
      for(int t = oldT; t < TOTALCUT; t++){
         Pcount[a][t] = affCount;
         Ncount[a][t] = 0;
      }
   }  // end of loop a
   for(int t = 0; t < TOTALCUT; t++){
      specificity[t] = (Ncount[0][t] + Pcount[0][t])? 1.0 * Ncount[0][t] / (Ncount[0][t] + Pcount[0][t]): 1;
      sensitivity[t] = (Ncount[1][t] + Pcount[1][t])? 1.0 * Pcount[1][t] / (Ncount[1][t] + Pcount[1][t]): 1;
   }
   if(printFlag){
      double tempT[]={0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999};
      int indexcount = 8;
      IntArray index(0);
      for(int i = 0; i < indexcount; i++)
         index.Push((1-tempT[i])*TOTALCUT-1);
      printf("Cutoff Value");
      for(int i = 0; i < indexcount; i++)
         printf("\t%.3lf", printFlag==2? 1.0/TOTALCUT*(index[i]+1): 1-1.0/TOTALCUT*(index[i]+1));
      printf("\n");
      printf("TruePositives");
      for(int i = 0; i < indexcount; i++) printf("\t%d", Pcount[1][index[i]]);
      printf("\n");
      printf("FalsePositives");
      for(int i = 0; i < indexcount; i++) printf("\t%d", Pcount[0][index[i]]);
      printf("\n");
      printf("TrueNegatives");
      for(int i = 0; i < indexcount; i++) printf("\t%d", Ncount[0][index[i]]);
      printf("\n");
      printf("FalseNegatives");
      for(int i = 0; i < indexcount; i++) printf("\t%d", Ncount[1][index[i]]);
      printf("\n");
      printf("Sensitivity");
      for(int i = 0; i < indexcount; i++) printf("\t%.4lf", sensitivity[index[i]]);
      printf("\n");
      printf("Specificity");
      for(int i = 0; i < indexcount; i++) printf("\t%.4lf", specificity[index[i]]);
      printf("\n");
      if(prevalence != _NAN_){
         printf("Positive PV");
         for(int i = 0; i < indexcount; i++) {
            int t = index[i];
            double PPV = prevalence * sensitivity[t] / (prevalence * sensitivity[t] + (1-prevalence)*(1-specificity[t]));
            printf("\t%.4lf", PPV);
         }
         printf("\n");
         printf("Negative PV");
         for(int i = 0; i < indexcount; i++) {
            int t = index[i];
            double NPV = (1-prevalence)*specificity[t] / ((1-prevalence)*specificity[t] + prevalence * (1-sensitivity[t]));
            printf("\t%.4lf", NPV);
         }
         printf("\n");
      }
      printf("\n");
   }
   double AUC = (1-specificity[1])*sensitivity[0];
   for(int i = 1; i < TOTALCUT-1; i++)
      AUC += (specificity[i-1] - specificity[i+1])*sensitivity[i];
   AUC += (specificity[TOTALCUT-2] - specificity[TOTALCUT-1])*sensitivity[TOTALCUT-1];
   AUC *= 0.5;
   return AUC;
}

void Engine::PrintGeneticRiskScoreSNPMajor(const char* weightfile, bool noflipFlag)
{
   if(sortedSNP.Length()==0) error("No SNPs were extracted");
   printf("\nOptions in effect:\n");
   printf("\t--risk\n");
   printf("\t--model %s\n", (const char*)weightfile);
   if(noflipFlag)
      printf("\t--noflip\n");
   if(prevalence != _NAN_)
      printf("\t--prevalence %G\n", prevalence);
   if(prefix != "king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");

   if(ped.affectionCount==0 && prevalence != _NAN_)
      printf("--prevalence is not used in absence of disease status.\n");

   StringIntHash colHash;
   colHash.Clear();
   colHash.SetInteger("SNP", SNPCOL);
   colHash.SetInteger("SNPNAME", SNPCOL);
   colHash.SetInteger("MARKER", SNPCOL);
   colHash.SetInteger("EA", EACOL);
   colHash.SetInteger("EAL", EACOL);
   colHash.SetInteger("ALLELE", EACOL);
   colHash.SetInteger("AF", AFCOL);
   colHash.SetInteger("EAF", AFCOL);
   colHash.SetInteger("MAF", AFCOL);
   colHash.SetInteger("WT", WTCOL);
   colHash.SetInteger("WEIGHT", WTCOL);
   colHash.SetInteger("CHR", CHRCOL);
   colHash.SetInteger("CHROMOSOME", CHRCOL);
   colHash.SetInteger("POS", POSCOL);
   colHash.SetInteger("POSITION", POSCOL);
   colHash.SetInteger("AA", OACOL);
   colHash.SetInteger("OA", OACOL);
   colHash.SetInteger("AAL", OACOL);
   colHash.SetInteger("OAL", OACOL);
   colHash.SetInteger("ON", ONCOL);
   colHash.SetInteger("PROXY", ONCOL);
   String line;
   StringArray tokens;
   IFILE input;
   input=ifopen((const char*)weightfile, "rt");
   if(input==NULL) error("Cannot open %s to read", (const char*)weightfile);
   printf("Read in %s...\n", (const char*)weightfile);
   tokens.Clear();
   line.ReadLine(input);
   tokens.AddTokens(line);
   if(tokens.Length() < 4)
      error("A header needs to include: SNP EA AF WT\n  Or even better: SNP EA AF WT CHR POS OA ON");
   int col[8];
   for(int i = 0; i < 8; i++) col[i] = -1;
   for(int i = 0; i < tokens.Length(); i++){
      int k = colHash.Integer(tokens[i]);
      if(k > -1) col[k] = i;
   }
   if(col[SNPCOL]==-1 || col[EACOL]==-1 || col[AFCOL]==-1 || col[WTCOL]==-1)
      error("A header needs to include: SNP EA AF WT\n  Or even better: SNP EA AF WT CHR POS OA ON");
   IntArray snpidx(0), chrs(0), positions(0);
   StringArray snpnames(0), eals(0), aals(0);
   Vector eafs(0), weights(0);
   StringIntHash markerLookup;
   markerLookup.Clear();
   for(int i = 0; i < sortedSNP.Length(); i++)
      markerLookup.SetInteger(snpName[sortedSNP[i]], i);
   while(!ifeof(input)){
      tokens.Clear();
      line.ReadLine(input);
      tokens.AddTokens(line);
      if(tokens.Length() < 4) continue;
      int k = markerLookup.Integer(tokens[col[SNPCOL]]);
      if(k==-1) {
         if(col[ONCOL]!=-1)
            k = markerLookup.Integer(tokens[col[ONCOL]]);
         if(k==-1) continue;
         snpnames.Push(tokens[col[ONCOL]]);
      }else
         snpnames.Push(tokens[col[SNPCOL]]);
      snpidx.Push(k);
      eals.Push(tokens[col[EACOL]]);
      eafs.Push(double(tokens[col[AFCOL]]));
      weights.Push(double(tokens[col[WTCOL]]));
      if(col[CHRCOL]!=-1)
         chrs.Push(int(tokens[col[CHRCOL]]));
      if(col[POSCOL]!=-1)
         positions.Push(int(tokens[col[POSCOL]]));
      if(col[OACOL]!=-1)
         aals.Push(tokens[col[OACOL]]);
   }
   ifclose(input);

   printf("Compute genetic risk scores...\n");
   int snpcount=snpidx.Length();
   Vector GRS(idCount);
   Vector GRS_MissingVar(idCount);
   Vector GRS_MissingSNP(idCount);
   GRS.Zero();
   GRS_MissingVar.Zero();
   GRS_MissingSNP.Zero();
   double WT[4][256], missingVar[4][256];
   int missingSNP[4][256];
   String action;
   printf("\nSummary at %d SNPs in model vs. prediction data\n", snpcount);
   printf("%15s", "SNP_NAME");
   if(col[CHRCOL]!=-1) printf(" %5s", "CHR");
   printf(" %5s", "EA");
   if(col[OACOL]!=-1) printf(" %5s", "OA");
   printf(" %7s %7s %5s %5s %7s %10s\n", "AF", "WEIGHT", "A1", "A2", "AF1", "ACTION");
   bool freqFlag = true;
   if(idCount < 10) {
      freqFlag = false;
      if(!flipFlag)
         printf("Due to small sample size, manual strand alighment with --noflip option is highly recommended.\n");
   }
   int MISSCOUNT[256];
   for(int byte = 0; byte < 256; byte++)
      MISSCOUNT[byte] = oneCount[(byte & 0x55) & ~((byte & 0xAA)>>1)];
   int byteCount = (idCount+3)/4;
   int Acount = 0, misscount = 0;
   double allmean=0;
   double allVar = 0;
   int validSNPCount = 0;
   double temp;
   char FLIP[256];
   for(int i = 0; i < 256; i++)
      FLIP[i] = i;
   FLIP['a'] = FLIP['A'] = FLIP['t'] = FLIP['T'] = 'A';
   FLIP['c'] = FLIP['C'] = FLIP['g'] = FLIP['G'] = 'C';
   for(int m = 0; m < snpcount; m++){
      int s = snpidx[m];
      Acount = byteCount*8;   // N_A = N*2
      misscount = 0;
      for(int i  = 0; i < byteCount; i++){
         Acount -= oneCount[G_SNP[s][i]]; // misscount is substracted too
         misscount += MISSCOUNT[G_SNP[s][i]];
      }
      Acount -= misscount; // misscount*2 - misscount
      double freq = (misscount==byteCount*4)? _NAN_: Acount*0.5 / (byteCount*4 - misscount);
      action = "OK";
      if(FLIP[alleleLabel[0][sortedSNP[s]][0]] == FLIP[alleleLabel[1][sortedSNP[s]][0]]){ // ambiguous SNP
         if(noflipFlag){
            if(eals[m] == alleleLabel[1][sortedSNP[s]]) // Allele switched
               action = "SWITCHED";
         }else{
            if(!freqFlag || freq==_NAN_) action = "SKIPPED";   // too few samples
            else{
               if(eals[m] == alleleLabel[0][sortedSNP[s]]){ // Allele matched
                  if(freq - eafs[m] > 0.3 || eafs[m] - freq > 0.3) action = "SKIPPED";
               }else if(eals[m] == alleleLabel[1][sortedSNP[s]]){ // Allele switched
                  if(freq + eafs[m] < 0.3)   // 1 - freq - eafs[m] > 0.7
                     action = "FLIPPED";
                  else
                     action = "SKIPPED";
               }else // Allele not matched at all
                  action = "SKIPPED";
            }
         }
      }else if(alleleLabel[0][sortedSNP[s]]=="0"){ // one allele only
         if(noflipFlag && eals[m] == alleleLabel[1][sortedSNP[s]])
            action = "SWITCHED";
         else if(!noflipFlag || freq==_NAN_)
            action = "SKIPPED";
      }else{   // unambiguous SNP
         if(FLIP[eals[m][0]] == FLIP[alleleLabel[1][sortedSNP[s]][0]])
            action = "SWITCHED";
         else if(eals[m] != alleleLabel[0][sortedSNP[s]]){
            if(FLIP[eals[m][0]] == FLIP[alleleLabel[0][sortedSNP[s]][0]])
               action = "FLIPPED";
            else
               action = "SKIPPED";
         }
      }  // end of if
      printf("%15s", (const char*)snpnames[m]);
      if(col[CHRCOL]!=-1) printf(" %5d", chrs[m]);
      printf(" %5s",(const char*)eals[m]);
      if(col[OACOL]!=-1) printf(" %5s", (const char*)aals[m]);
      printf(" %7.3lf %7.3lf %5s %5s %7.3lf %10s\n",
         eafs[m], weights[m],
         (const char*)alleleLabel[0][sortedSNP[s]], (const char*)alleleLabel[1][sortedSNP[s]],
         freq==_NAN_?-9:freq, (const char*)action);
      if(action=="SKIPPED") continue;

      validSNPCount ++;
      if(action=="SWITCHED") weights[m] = -weights[m];
      temp = weights[m] * 2.0 * eafs[m];
      double snpvar = 2.0 * eafs[m] * (1-eafs[m])*weights[m]*weights[m];
      allVar += snpvar;
      for(int k = 0; k < 4; k++)
         for(int byte = 0; byte < 256; byte++){
            missingVar[k][byte] = 0;
            missingSNP[k][byte] = 0;
         }
      for(int byte = 0; byte < 256; byte++)
         for(int k = 0; k < 4; k++){
            int bit = (byte >> (k*2)) & 3;
            if(bit == 0) WT[k][byte] = weights[m] * 2;   // AA
            else if(bit == 2) WT[k][byte] = weights[m];  // Aa
            else if(bit == 1) {                          // missing
               WT[k][byte] = temp;
               missingVar[k][byte] = snpvar;
               missingSNP[k][byte] = 1;
            }else WT[k][byte] = 0;                        // aa
         }
      allmean += temp;
      for(int i = 0; i < idCount; i++){
         int g = G_SNP[s][i/4];
         GRS[i] += WT[i%4][g];
         GRS_MissingVar[i] += missingVar[i%4][g];
         GRS_MissingSNP[i] += missingSNP[i%4][g];
      }  // end of loop i for idCount
   }  // end of loop m for snpcount
   for(int i = 0; i < idCount; i++){
      GRS[i] -= allmean;
      GRS_MissingVar[i] /= allVar;
      GRS_MissingSNP[i] /= validSNPCount;
   }
   printf("\n");
   double sigma = sqrt(allVar);
   Vector risks(idCount);
   Vector zscores(idCount);
   Vector pvalues(idCount);
   for(int i = 0; i < idCount; i++){
      temp = exp(GRS[i]);
      risks[i] = temp / (temp+1);
      zscores[i] = GRS[i] / sigma;
      pvalues[i] = ndist(zscores[i], true);
   }
   String riskscorefile=prefix;
   riskscorefile.Add("grs.txt");
   FILE *fp=fopen((const char*)riskscorefile, "wt");
   if(fp==NULL) error("Cannot open %s to write", (const char*)riskscorefile);
   fprintf(fp, "%15s %15s %7s %7s %7s %7s %7s %9s",
      "FID", "IID", "InfoSNP", "InfoVar", "GRS", "Zscore", "Percent", "ScaledGRS");
   if(ped.affectionCount>0)
      fprintf(fp, " %7s", "Status");
   fprintf(fp, "\n");
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = 0; i < id[f].Length(); i++){
         int k = id[f][i];
         int j = geno[k];
         fprintf(fp, "%15s %15s %7.3lf %7.3lf %7.3lf %7.3lf %7.4lf %9.4lf",
         (const char*)ped[k].famid, (const char*)ped[k].pid,
         1 - GRS_MissingSNP[j], 1 - GRS_MissingVar[j],
         GRS[j], zscores[j], pvalues[j], risks[j]);
         if(ped.affectionCount>0){
            if(ped[k].affections[0]>0)
               fprintf(fp, " %7d", ped[k].affections[0]-1);
            else
               fprintf(fp, " %7s", "NA");
         }
         fprintf(fp, "\n");
   }
   fclose(fp);
   if(ped.affectionCount>0){
      IntArray diseaseStatus(idCount);
      for(int i = 0; i < idCount; i++)
         diseaseStatus[i] = ped[phenoid[i]].affections[0]-1;
      printf("Using estimated percentiles as cutoff values:\n");
      double AUC = ComputeAUC(pvalues, diseaseStatus, 2);
      printf("AUC (Area under the ROC curve) = %.4lf\n", AUC);
      printf("\nUsing scaled-GRS as cutoff values:\n");
      AUC = ComputeAUC(risks, diseaseStatus, 1);
      printf("AUC (Area under the ROC curve) = %.4lf\n", AUC);
      IntArray maleDiseaseStatus(0), femaleDiseaseStatus(0);
      Vector malerisks(0), femalerisks(0);
      for(int i = 0; i < idCount; i++){
         if(ped[phenoid[i]].sex==1) {
            maleDiseaseStatus.Push(ped[phenoid[i]].affections[0]-1);
            malerisks.Push(risks[i]);
         }else if(ped[phenoid[i]].sex==2){
            femaleDiseaseStatus.Push(ped[phenoid[i]].affections[0]-1);
            femalerisks.Push(risks[i]);
         }
      }
      if(maleDiseaseStatus.Length()){
         AUC = ComputeAUC(malerisks, maleDiseaseStatus, 0);
         printf("AUC among %d males is %.4lf\n", maleDiseaseStatus.Length(), AUC);
      }
      if(femaleDiseaseStatus.Length()){
         AUC = ComputeAUC(femalerisks, femaleDiseaseStatus, 0);
         printf("AUC among %d females is %.4lf\n", femaleDiseaseStatus.Length(), AUC);
      }
      printf("\n");
   }
   printf("Genetic risk scores are saved in file %s\n", (const char*)riskscorefile);
   for(int m = 0; m < sortedSNP.Length(); m++)
      delete[] G_SNP[m];
   delete[] G_SNP;
   printf("Genetic risk score prediction ends at %s", currentTime());
}







/*
void Engine::PrintGeneticRiskScore(const char* weightfile)
{
   printf("Autosome genotypes stored in %d words for each of %d individuals.\n",
      shortCount, idCount);

   printf("\nOptions in effect:\n");
   printf("\t--risk\n");
   printf("\t--model %s\n", (const char*)weightfile);
   if(prevalence!=_NAN_)
      printf("\t--prevalence %lf\n", prevalence);
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");
   String line;
   StringArray tokens;
   IFILE input;
   input=ifopen((const char*)weightfile, "rt");
   if(input==NULL) error("Cannot open %s to read", (const char*)weightfile);
   printf("Read in %s...\n", (const char*)weightfile);
   tokens.Clear();
   line.ReadLine(input);
   tokens.AddTokens(line);
   if(tokens.Length() < 7)
      error("Header: SNP CHR POS EAL AAL EAF WT");
   if(tokens[0]!="SNP" || (tokens[3]!="WEIGHT" && tokens[3]!="WT") )
      error("Header: SNP EAL AF WT CHR POS AAL");

   IntArray snpidx(0), chrs(0), positions(0);
   StringArray snpnames(0), eals(0), aals(0);
   Vector eafs(0), weights(0);
   StringIntHash markerLookup;
   markerLookup.Clear();
   for(int i = 0; i < markerCount; i++)
      markerLookup.SetInteger(snpName[i], i);
   while(!ifeof(input)){
      tokens.Clear();
      line.ReadLine(input);
      tokens.AddTokens(line);
      if(tokens.Length() < 7) continue;
      int k = markerLookup.Integer(tokens[SNPCOL]);
      if(k==-1){  // SNP not available
         printf("SNP %s from the risk model cannot be found in the data to be predicted.\n",
            (const char*)tokens[SNPCOL]);
         continue;
      }
      weights.Push(double(tokens[WTCOL]));
      snpidx.Push(k);
      snpnames.Push(tokens[SNPCOL]);
      chrs.Push(int(tokens[CHRCOL]));
      positions.Push(int(tokens[POSCOL]));
      eals.Push(tokens[EACOL]);
      aals.Push(tokens[OACOL]);
      eafs.Push(double(tokens[AFCOL]));
   }
   ifclose(input);

   printf("Compute genetic risk scores...\n");
   int snpcount=snpidx.Length();
   Vector GRS(idCount);
   GRS.Zero();
   double temp;
   double WT[2][2];

   String action;
   printf("\nSummary at %d SNPs in model vs. prediction data\n", snpcount);

   printf("SNP_NAME");
   if(CHRCOL!=-1) printf("\tCHR");
//   if(POSCOL!=-1) printf("\tPOSITION");
   printf("\tEA");
   if(OACOL!=-1) printf("\tOA");
   printf("\tEAF\tWEIGHT\tAL1\tAL2\tAF1\tACTION\n");
   for(int m = 0; m < snpcount; m++){
      int w = snpidx[m]/16;
      int b = snpidx[m]%16;

      double freq=0;
      int miss = 0;
      for(int i = 0; i < idCount; i++){
         freq += ((~GG[0][i][w] & GG[1][i][w] & (1<<b)) > 0) * 0.5; // Aa
         freq += ((GG[0][i][w] & GG[1][i][w] & (1<<b)) > 0);   // AA
         miss += ((~GG[0][i][w] & ~GG[1][i][w] & (1<<b)) > 0); // Missing
      }
      freq /= (idCount - miss);
      action = "OK";
      if((alleleLabel[0][snpidx[m]]=="A" && alleleLabel[1][snpidx[m]]=="T") ||
      (alleleLabel[0][snpidx[m]]=="T" && alleleLabel[1][snpidx[m]]=="A") ||
      (alleleLabel[0][snpidx[m]]=="C" && alleleLabel[1][snpidx[m]]=="G") ||
      (alleleLabel[0][snpidx[m]]=="G" && alleleLabel[1][snpidx[m]]=="A") ){
         if(eals[m] == alleleLabel[0][snpidx[m]]){ // Allele matched
            if(freq - eafs[m] > 0.5 || eafs[m] - freq > 0.5) action="SKIPPED";
         }else if(eals[m] == alleleLabel[0][snpidx[m]]){ // Allele not matched
            if(freq - eafs[m] > 0.5 || eafs[m] - freq > 0.5)
               action="FLIPPED";
            else
               action="SKIPPED";
         }
      }else{
         if(eals[m] == alleleLabel[1][snpidx[m]]) // Allele not matched
            action="FLIPPED";
      }

      printf("%s", (const char*)snpnames[m]);
      if(CHRCOL!=-1) printf("\t%d", chrs[m]);
//      if(POSCOL!=-1) printf("\t%d", positions[m]);
      printf("\t%s",(const char*)eals[m]);
      if(OACOL!=-1) printf("\t%s", (const char*)aals[m]);
      printf("\t%.3lf\t%.3lf\t%s\t%s\t%.3lf\t%s\n",
         eafs[m], weights[m],
         (const char*)alleleLabel[0][snpidx[m]], (const char*)alleleLabel[1][snpidx[m]],
         freq, (const char*)action);
      if(action=="SKIPPED") continue;
      if(action=="FLIPPED") weights[m] = -weights[m];
      temp = weights[m] * 2.0 * eafs[m];
      WT[0][0] = temp;  // 2pw for missing
      WT[0][1] = weights[m];                 // w for Aa
      WT[1][0] = 0;                          // 0 for aa
      WT[1][1] = weights[m] * 2.0;           // 2w for AA
      for(int i = 0; i < idCount; i++)
         GRS[i] += WT[(GG[0][i][w]>>b)&1][(GG[1][i][w]>>b)&1] - temp;
   }
   printf("\n");

   Vector risks(idCount);
   for(int i = 0; i < idCount; i++){
      temp = exp(GRS[i]);
      risks[i] = temp / (temp+1);
   }
   String riskscorefile=prefix;
   riskscorefile.Add("grs.txt");
   FILE *fp=fopen((const char*)riskscorefile, "wt");
   if(fp==NULL) error("Cannot open %s to write", (const char*)riskscorefile);
   fprintf(fp, "FID\tIID\tGRS");
   if(ped.affectionCount>0){
      fprintf(fp, "\tPr(%s)", (const char*)ped.affectionNames[0]);
      fprintf(fp, "\tObs_%s", (const char*)ped.affectionNames[0]);
   }else
      fprintf(fp, "\tPr(D)");
   fprintf(fp, "\n");
   for(int i = 0; i < idCount; i++){
      fprintf(fp, "%s\t%s\t%.3lf\t%.4lf",
         (const char*)ped[phenoid[i]].famid, (const char*)ped[phenoid[i]].pid, GRS[i], risks[i]);
      if(ped.affectionCount>0)
         fprintf(fp, "\t%d", ped[phenoid[i]].affections[0]-1);
      fprintf(fp, "\n");
   }
   fclose(fp);
   const int TOTALCUT=10000;
   if(ped.affectionCount>0){
      double cutoff[TOTALCUT], sensitivity[TOTALCUT], specificity[TOTALCUT];
      double PPV[TOTALCUT], NPV[TOTALCUT];
      int TP, TN, FP, FN;
      for(int t = 0; t < TOTALCUT; t++){
         cutoff[t] = 1-1.0/TOTALCUT*(t+1);
         TP = TN = FP = FN = 0;
         for(int i = 0; i < idCount; i++){
            TP += (risks[i] > cutoff[t] && ped[phenoid[i]].affections[0] == 2);
            TN += (risks[i] <= cutoff[t] && ped[phenoid[i]].affections[0] == 1);
            FP += (risks[i] > cutoff[t] && ped[phenoid[i]].affections[0] == 1);
            FN += (risks[i] <= cutoff[t] && ped[phenoid[i]].affections[0] == 2);
         }
         specificity[t] = TN+FP > 0? (double)TN / (TN+FP): 1;
         sensitivity[t] = TP+FN > 0? (double)TP / (TP+FN): 1;
         if(prevalence != _NAN_){
            PPV[t] = prevalence * sensitivity[t] / (prevalence * sensitivity[t] + (1-prevalence)*(1-specificity[t]));
            NPV[t] = (1-prevalence)*specificity[t] / ((1-prevalence)*specificity[t] + prevalence * (1-sensitivity[t]));
         }
      }
      printf("\nRisk Cutoff");
      for(int t = 8; t >= 0; t--) printf("\t%.1lf", cutoff[(t+1)*TOTALCUT/10-1]);
      printf("\n");
      printf("Sensitivity");
      for(int t = 8; t >= 0; t--) printf("\t%.4lf", sensitivity[(t+1)*TOTALCUT/10-1]);
      printf("\n");
      printf("Specificity");
      for(int t = 8; t >= 0; t--) printf("\t%.4lf", specificity[(t+1)*TOTALCUT/10-1]);
      printf("\n");
      if(prevalence != _NAN_){
         printf("Positive PV");
         for(int t = 8; t >= 0; t--) printf("\t%.4lf", PPV[(t+1)*TOTALCUT/10-1]);
         printf("\n");
         printf("Negative PV");
         for(int t = 8; t >= 0; t--) printf("\t%.4lf", NPV[(t+1)*TOTALCUT/10-1]);
         printf("\n");
      }
      printf("\n");
      double AUC = sensitivity[0] * (1-specificity[0])*0.5;
      for(int i = 1; i < TOTALCUT; i++)
         AUC += (specificity[i-1] - specificity[i]) * (sensitivity[i-1]+sensitivity[i])/2;
      printf("AUC (Area under the ROC curve) = %.4lf\n\n", AUC);
   }
   printf("Genetic risk scores are saved in file %s\n", (const char*)riskscorefile);
   printf("Genetic risk score prediction ends at %s", currentTime());

}
*/

         /*
      Vector pvalue(0);
      if(sigma!=_NAN_)
         for(int i = 0; i < indexcount; i++){
            double r = cutoff[index[i]];
            double z = log(r/(1-r))/sigma;
            pvalue.Push(ndist(z, true));
         } */
/*
   for(int i = 0; i < totalCount; i++){
      int aff = diseaseStatus[i];
      if((aff != 1 && aff != 0)|| risks[i] < -0.999) continue;
      int riskfold = int(risks[i]*TOTALCUT);
      if(printFlag==2){
         for(int t = TOTALCUT-1; t > riskfold; Pcount[aff][t]++, t--);  // risk < cutoff
         for(int t = riskfold; t >=0; Ncount[aff][t]++, t--);// risk >= cutoff
      }else{
         for(int t = TOTALCUT-1; t > TOTALCUT-1-riskfold; Pcount[aff][t]++, t--);  // cutoff <= risk
         for(int t = TOTALCUT-2-riskfold; t >=0; Ncount[aff][t]++, t--);// cutoff > risk
      }
   }
   */
//   double AUC = sensitivity[0] * (1-specificity[0])*0.5;
//   for(int i = 1; i < TOTALCUT; i++)
//      AUC += (specificity[i-1] - specificity[i]) * (sensitivity[i-1]+sensitivity[i])/2;


