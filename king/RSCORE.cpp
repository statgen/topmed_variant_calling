// RSCORE.cpp
// 4/6/2013 Wei-Min Chen

#include "analysis.h"
#include <math.h>
#include "MathStats.h"
#include "QuickIndex.h"
#include "MathCholesky.h"
#include "MathSVD.h"

#ifdef WITH_LAPACK
extern "C" void dgesdd_(char*, int*, int*, double*, int*, double*, double *, int*, double*, int*, double*, int*, int*, int*);
extern "C" void dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);
#endif

void Engine::RSCORE()
{
   if(detailFlag) countGenotype();
   if(geno.Length()==0) {
      individualInfo = true;
      if(shortFlag) BuildShortBinary();
      else BuildBinary();
   }
   if(FixedEff)
      PreScan_FixedEff();
   else
      PreScan_RSCORE();
   int N_Trait = traits.Length();
   printf("RSCORE genome scan starts at %s", currentTime());
   String lmmfile = prefix;
   lmmfile.Add("rscore.txt");
   FILE *fp = fopen(lmmfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)lmmfile);
   fprintf(fp, "SNP");
   if(N_Trait > 1)
      fprintf(fp, "\tTraitName");
   if(chromosomes.Length()) fprintf(fp, "\tChr");
   if(bp.Length()) fprintf(fp, "\tPos");
   fprintf(fp, "\tLabelA\tLabela\tFreqA\tN\tBeta\tSE\tChisq\tH2\tPvalue\n");
   ScoreScan(fp);
   fclose(fp);
   for(int t = 0; t < N_Trait; t++)
      delete missingpheno[t];
   delete []missingpheno;
   printf("RSCORE genome scan ends at %s", currentTime());
   printf("RSCORE genome scan results saved in file %s\n\n", (const char*)lmmfile);
   if(chisqFilter==_NAN_)
      PostVC();
}

void Engine::ScoreScan(FILE *fp)
{
   int N_ID = ID.Length();
   int N_Trait = traits.Length();
   double pvalue, h2, freq0[16], freq[16], statistic, score_local, var_snp, Tm;
   int missingCount0[16], missingCount[16], nonmissing, id, AA, Aa, missing, pos;
   if(chisqFilter==_NAN_){
      chisqs = new Vector[N_Trait];
      for(int t = 0; t < N_Trait; t++)
         chisqs[t].Dimension(0);
   }
   int *G2[16], *G1[16], *GM[16];
   for(int m = 0; m < 16; m++){
      G2[m] = new int[N_ID+1];
      G1[m] = new int[N_ID+1];
      GM[m] = new int[N_ID+1];
   }
   double *Sum_VR = new double[N_Trait];
   for(int t = 0; t < N_Trait; t++){
      Sum_VR[t] = 0;
      for(int i = 0; i < N_ID; i++)
         Sum_VR[t] += VR[t][i];
   }
   // parallel computing here
   for(int b = 0; b < shortCount; b++){ // loop over short integer
      for(int m = 0; m < 16; m++){
         freq0[m] = 0;
         missingCount0[m] = 0;
         G2[m][0] = G1[m][0] = GM[m][0] = 0;
      }
      for(int i = 0; i < N_ID; i++){
         id = geno[ID[i]];
         AA = GG[0][id][b] & GG[1][id][b];
         Aa = (~GG[0][id][b]) & GG[1][id][b];
         missing = (~GG[0][id][b]) & (~GG[1][id][b]) & 65535;
         for(int m = 0; m < 16; m++){
            if(missing & shortbase[m]){
               GM[m][0] ++;
               GM[m][GM[m][0]] = i;
               missingCount0[m] ++;
            }else if(AA & shortbase[m]){
               G2[m][0] ++;
               G2[m][G2[m][0]] = i;
               freq0[m] += 1.0;
            }else if(Aa & shortbase[m]){
               G1[m][0] ++;
               G1[m][G1[m][0]] = i;
               freq0[m] += 0.5;
            }
         }
      }
      for(int t = 0; t < N_Trait; t++){
         for(int m = 0; m < 16; m++){
            freq[m] = freq0[m];
            missingCount[m] = missingCount0[m];
         }
         for(int i = 0; i < missingpheno[t][0]; i++){
            id = geno[ID[missingpheno[t][i+1]]];
            AA = GG[0][id][b] & GG[1][id][b];
            Aa = (~GG[0][id][b]) & GG[1][id][b];
            nonmissing = GG[0][id][b] | GG[1][id][b];
            for(int m = 0; m < 16; m++)
               if(nonmissing & shortbase[m]){
                  missingCount[m] ++;
                  if(AA & shortbase[m])
                     freq[m] -= 1.0;
                  else if(Aa & shortbase[m])
                     freq[m] -= 0.5;
               }
         }      
         for(int m = 0; m < 16; m++){
            if(missingCount[m]==N_ID) continue;
            freq[m] /= (N_ID-missingCount[m]);
            var_snp = 2*freq[m]*(1-freq[m]);
            if(var_snp < 1E-10) continue;
            score_local = -Sum_VR[t];
            for(int i = 0; i < GM[m][0]; i++)
               score_local += VR[t][GM[m][i+1]];
            score_local *= freq[m];
            for(int i = 0; i < G2[m][0]; i++)
               score_local += VR[t][G2[m][i+1]];
            score_local += score_local;
            for(int i = 0; i < G1[m][0]; i++)
               score_local += VR[t][G1[m][i+1]];
            statistic = score_local * score_local / var_snp;
            if(chisqFilter==_NAN_){
               chisqs[t].Push(statistic);
            }else if(statistic < chisqFilter)
               continue;
            pvalue = ndist(sqrt(statistic))*2;
            h2 = score_local * score_local / varVR[t] / var_snp;
            h2 = h2/((1+lambda0[t])/tau0[t]+h2);
            pos = b*16 + m;
            if(snpName.Length())
               fprintf(fp, "%s", (const char*)snpName[pos]);
            else
               fprintf(fp, "SNP%d", pos+1);
            if(N_Trait>1)
               fprintf(fp, "\t%s", (const char*)ped.traitNames[traits[t]]);
            if(chromosomes.Length())
               fprintf(fp, "\t%d", chromosomes[pos]);
            if(bp.Length())
               fprintf(fp, "\t%.6lf", bp[pos]*0.000001);
            fprintf(fp, "\t%s\t%s\t%.3lf\t%d",
               (const char*)alleleLabel[0][pos], (const char*)alleleLabel[1][pos],
               freq[m], N_ID-missingCount[m]);
            if(statistic != _NAN_)
               fprintf(fp, "\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.3G\n",
                  score_local / sqrt(varVR[t]) / var_snp,
                  1/sqrt(varVR[t]*var_snp),
                  statistic, h2, pvalue);
            else
               fprintf(fp, "\t%s\t%s\n", "NA", "NA");
         }
      }

   }
   for(int t = 0; t < N_Trait; t++)
      delete VR[t];
   delete []VR;
   if(varVR) delete []varVR;

   if(Sum_VR) delete []Sum_VR;
   for(int m = 0; m < 16; m++){
      delete []G2[m];
      delete []G1[m];
      delete []GM[m];
   }
}

void Engine::PreScan_RSCORE()
{
   int N_Trait = traits.Length();
   if(N_Trait == 0) error("Quantitative traits are not found for variance component analysis.");
   int N_Covariate = covariates.Length();
   IntArray validFlag(ped.count);
   validFlag.Set(1);
   for(int i = 0; i < ped.count; i++){
      bool allmissing = true;
      for(int t = 0; t < N_Trait; t++)
         if(ped[i].isPhenotyped(traits[t])) {
            allmissing = false;
            break;
         }
      if(allmissing) validFlag[i] = 0;
      bool somemissing = false;
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
   ID0 = ID;
   int N_ID0 = ID.Length();
   Matrix Ys(N_ID0, N_Trait);
   Matrix X(N_ID0, 1+N_Covariate);
   for(int i = 0; i < N_ID0; i++)
      X[i][0] = 1.0;
   for(int t = 0; t < N_Trait; t++)
      for(int i = 0; i < N_ID0; i++)
         Ys[i][t] = ped[ID[i]].traits[traits[t]];
   for(int i = 0; i < N_ID0; i++)
      for(int j = 0; j < N_Covariate; j++)
         X[i][j+1] = ped[ID[i]].covariates[covariates[j]];
   printf("Genotypes stored in %d integers for each of %d individuals used in analysis.\n",
      shortCount, N_ID0);

   if(svdinfile!=""){
      if(N_Trait==1){
         EV.Dimension(N_ID0);
         UT = new double * [N_ID0];
         for(int i = 0; i < N_ID0; i++)
            UT[i] = new double [N_ID0];
       }
   if(D==NULL){
      D = new double * [N_ID0];
      for(int i = 0; i < N_ID0; i++)
         D[i] = new double [N_ID0];
      for(int i = 0; i < N_ID0; i++)
         for(int j = 0; j < N_ID0; j++)
            D[i][j] = 0;
   }
   StringArray tokens;
   String line;
   IFILE input = ifopen(svdinfile, "rt");
   if(input == NULL)
      error("Cannot open %s to read", (const char*)svdinfile);
   tokens.Clear();
   line.ReadLine(input);
   tokens.AddTokens(line);
   if(tokens.Length() < N_ID0)
      error("%d samples in SVD file %s; however there are %d samples in the data.\n",
            tokens.Length(), (const char*)svdinfile, N_ID0);
   String tempID;
   int N_tokens = tokens.Length();
   IntArray key(N_tokens);
   key.Set(-1);
   for(int i = 0; i < N_ID0; i++){
         tempID = ped[ID[i]].famid;
         tempID.Add("_");
         tempID.Add(ped[ID[i]].pid);
         if(tempID == tokens[i])
            key[i] = i;
         else{
            int k = tokens.Find(tempID);
            if(k < 0)
               error("Family %s ID %s cannot be found in SVD file",
               (const char*)ped[ID[i]].famid, (const char*)ped[ID[i]].pid);
            key[k] = i;
         }
   }
   tokens.Clear();
   line.ReadLine(input);
   tokens.AddTokens(line);
   if(tokens.Length() < N_tokens)
      error("There are only %d Eigenvalues (not %d) in the 2nd row of the SVD file",
         tokens.Length(), N_tokens);
   if(N_Trait==1)
   for(int i = 0; i < N_tokens; i++){
      if(key[i]==-1) continue;
      EV[key[i]] = tokens[i].AsDouble();
   }
   int row = 0;
   while(!ifeof(input)){
      if(row == N_tokens) break;
      tokens.Clear();
      line.ReadLine(input);
      tokens.AddTokens(line);
      if(tokens.Length() < N_tokens) continue;
      if(key[row]==-1) {
         row++;
         continue;
      }
      if(N_Trait==1)
      for(int i = 0; i < N_tokens; i++){
         if(key[i] == -1) continue;
         UT[key[row]][key[i]] = tokens[i].AsDouble();
      }
      row++;
   }
   if(row < N_tokens)
      error("There are only %d rows (not %d) in the Eigenvector data",
            row, N_tokens);
   while(!ifeof(input)){
      if(row == N_tokens*2) break;
      tokens.Clear();
      line.ReadLine(input);
      tokens.AddTokens(line);
      if(tokens.Length() < N_tokens) continue;
      if(key[row-N_tokens]==-1){
         row++;
         continue;
      }
      for(int i = 0; i < N_tokens; i++){
         if(key[i]==-1) continue;
         D[key[row-N_tokens]][key[i]] = tokens[i].AsDouble();
      }
      row++;
   }
   if(row < N_tokens*2)
      error("There are only %d rows (not %d) in the genetic similarity data",
            row-N_ID0, N_tokens);
   printf("Eigenvalues, eigenvectors, and genetic similarity are successfully loaded from file %s\n",
         (const char*)svdinfile);
   ifclose(input);
   
   }else
      ComputeSimilarity();
   if(normalization)
      printf("Inverse normal transformation is applied to phenotypes.\n");

   Vector Y;
   UX.Dimension(N_ID0, N_Covariate+1);
   UX.Zero();
   Vector tV;
   Cholesky chol;
   QuickIndex idx;
   double a, b, LB, T, x, temp, funcx;
   int numiter;
   int maxIter=10000;
   bool quiet=true;
   T = 0.001;
   lambda0.Dimension(N_Trait);
   lambda0.Zero();
   Vector score0(N_Covariate+1);
   Matrix Info0(N_Trait, N_Trait);
   printf("\nPolygenic parameter estimates\n");
   printf("%-15s %7s %7s %7s %9s",
      "TraitName", "N", "Herit", "Lambda", "Mu");
   for(int i = 0; i < N_Covariate; i++)
      printf(" %9s", (const char*)ped.covariateNames[covariates[i]]);
   printf("\n");

   Vector beta0(N_Covariate+1);
   Vector Resid0;
   double loglik0;
   tau0.Dimension(N_Trait);
   double R1R;
   VR = new double * [N_Trait];
   for(int t = 0; t < N_Trait; t++){
      VR[t] = new double [N_ID0];
      for(int i = 0; i < N_ID0; i++)
         VR[t][i] = 0;
   }
   varVR = new double [N_Trait];
   int N_ID;


   IntArray missingpheno0;
   missingpheno = new int *[N_Trait];
   for(int t = 0; t < N_Trait; t++){   // Precompute, by trait
      validpheno.Dimension(0);
      missingpheno0.Dimension(0);
      for(int i = 0; i < N_ID0; i++)
         if(Ys[i][t] != _NAN_)
            validpheno.Push(i);
         else
            missingpheno0.Push(i);
      missingpheno[t] = new int [missingpheno0.Length()+1];
      missingpheno[t][0] = missingpheno0.Length();
      for(int i = 0; i < missingpheno[t][0]; i++)
         missingpheno[t][i+1] = missingpheno0[i];
      ID.Dimension(0);
      for(int i = 0; i < validpheno.Length(); i++)
         ID.Push(ID0[validpheno[i]]);
      N_ID = ID.Length();
      Y.Dimension(N_ID);
      for(int i = 0; i < N_ID; i++)
         Y[i] = ped[ID[i]].traits[traits[t]];
      if(N_Trait>1 || svdinfile==""){
         EV.Dimension(N_ID);
         UT = new double * [N_ID];
         for(int i = 0; i < N_ID; i++)
            UT[i] = new double [N_ID];
         ComputeScoreSVD(quiet);
      }
      if(normalization){
         tV.Dimension(N_ID);
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

      UY.Dimension(N_ID);
      UY.Zero();
      for(int i = 0; i < N_ID; i++)
         for(int k = 0; k < N_ID; k++)
            UY[i] += UT[i][k] * Y[k];
      UX.Dimension(N_ID, N_Covariate+1);
      UX.Zero();
      for(int i = 0; i < N_ID; i++)
         for(int j = 0; j < N_Covariate+1; j++)
            for(int k = 0; k < N_ID; k++)
               UX[i][j] += UT[i][k] * X[validpheno[k]][j];

      currentT = t;
      loglik0 = 1E10;
      for(int i = 0; i < 100; i++){
         a = exp(i*0.1 - 5);
         b = exp(i*0.1 - 4.9);
         x = minimize(a, b, T, funcx, numiter, maxIter, quiet);
         if(funcx < loglik0) {
            loglik0 = funcx;
            lambda0[currentT] = 1/x;
         }
      }
      if(lambda0[currentT]<0.0067){  // border
         lambda0[currentT] = 0;
      }
      Info0.Dimension(N_Covariate+1, N_Covariate+1);
      Info0.Zero();
      for(int i = 0; i < N_Covariate+1; i++)
         for(int j = 0; j < N_Covariate+1; j++)
            for(int k = 0; k < N_ID; k++)
               Info0[i][j] += UX[k][i] * UX[k][j] / (lambda0[currentT] * EV[k] + 1);
      chol.Decompose(Info0);
      score0.Zero();
      for(int i = 0; i < N_Covariate+1; i++)
         for(int k = 0; k < N_ID; k++)
            score0[i] += UX[k][i] * UY[k] / (lambda0[currentT] * EV[k] + 1);
      chol.BackSubst(score0);
      beta0 = chol.x;
      R1R = 0;

      Resid0.Dimension(N_ID);
      for(int i = 0; i < N_ID; i++){
         Resid0[i] = UY[i];
         for(int j = 0; j < N_Covariate+1; j++)
            Resid0[i] -= UX[i][j] * beta0[j];
         R1R += Resid0[i] * Resid0[i] / (lambda0[currentT] * EV[i]+1);
      }
      tau0[currentT] = N_ID / R1R;

      printf("%-15s %7d %7.4lf %7.4lf",
         (const char*)ped.traitNames[traits[currentT]], N_ID,
         lambda0[currentT]/(1+lambda0[currentT]), lambda0[currentT]);
      for(int i = 0; i < N_Covariate+1; i++)
         printf(" %9.3lf", beta0[i]);
      printf("\n");

      varVR[t] = 0;
      for(int i = 0; i < N_ID0; i++)
         VR[t][i] = 0;
      for(int i = 0; i < N_ID; i++)
         for(int j = 0; j < N_ID; j++)
            VR[t][validpheno[i]] += UT[j][i] * Resid0[j] * tau0[t] / (lambda0[t]*EV[j]+1);
      for(int i = 0; i < N_ID; i++)
         varVR[t] += Resid0[i] * Resid0[i] * EV[i] * tau0[t] * tau0[t] / (lambda0[t]*EV[i]+1) / (lambda0[t]*EV[i]+1);
      for(int i = 0; i < N_ID; i++)
         VR[t][validpheno[i]] /= sqrt(varVR[t]);
      for(int i = 0; i < N_ID; i++)
         delete UT[i];
      delete []UT;
   }
   ID = ID0;
   for(int i = 0; i < N_ID0; i++)
      delete D[i];
   delete []D;
   chisqs = new Vector[N_Trait];
   for(int t = 0; t < N_Trait; t++)
      chisqs[t].Dimension(0);
}


void Engine::ComputeSimilarity()
{
   int IBS0Count, het1Count, notMissingCount;
   int id1, id2, count;
   int N_ID = ID.Length();

   printf("Calculating similarity matrix starts at %s", currentTime());
   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++){
      oneoneCount[i] = 0;
      for(int j = 0; j < 16; j++)
         if(i & shortbase[j]) oneoneCount[i]++;
   }                   

   if(D==NULL){
      D = new double * [N_ID];
      for(int i = 0; i < N_ID; i++)
         D[i] = new double [N_ID];
      for(int i = 0; i < N_ID; i++)
         for(int j = 0; j < N_ID; j++)
            D[i][j] = 0;
   }
   // parallel computing here
   for(int i = 0; i < N_ID; i++)
      for(int j = i+1; j < N_ID; j++){
         id1 = geno[ID[i]]; id2 = geno[ID[j]];
         if(id1 < 0 || id2 < 0) continue;
         IBS0Count = het1Count = notMissingCount = 0;
         for(int m = 0; m < shortCount; m++){
            IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
            notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
            het1Count += oneoneCount[(GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m]) |
                     (GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]) ];
         }
         if(notMissingCount)
            D[i][j] = D[j][i] = (het1Count + 4.0*IBS0Count) / notMissingCount;
      }

   Vector tempV(N_ID);
   tempV.Zero();
   for(int j = 0; j < N_ID; j++)
      for(int i = 0; i < N_ID; i++)
         tempV[j] += D[i][j];
   tempV.Multiply(1.0/N_ID);
   // (I-11'/N) * D
   for(int i = 0; i < N_ID; i++)
      for(int j = 0; j < N_ID; j++)
         D[i][j] -= tempV[j];
   tempV.Zero();
   for(int i = 0; i < N_ID; i++)
      for(int j = 0; j < N_ID; j++)
         tempV[i] += D[i][j];
   tempV.Multiply(1.0/N_ID);
   // D * (I-11'/N)
   for(int i = 0; i < N_ID; i++)
      for(int j = 0; j < N_ID; j++)
         D[i][j] -= tempV[i];
//   D.Multiply(-0.5);
   for(int i = 0; i < N_ID; i++)
      for(int j = 0; j < N_ID; j++)
         D[i][j] *= -0.5;


   for(int i = 0; i < N_ID; i++)
      D[i][i] = sqrt(D[i][i]); // save computation for below
   for(int i = 0; i < N_ID; i++)
      for(int j = i+1; j < N_ID; j++)
         D[i][j] = D[j][i] = D[i][j] / (D[i][i]*D[j][j]);
   for(int i = 0; i < N_ID; i++)
      D[i][i] = 1.0;

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
}

void Engine::ComputeScoreSVD(bool quiet)
{
   int N_ID = ID.Length();
   if(!quiet)
      printf("SVD starts at %s...", currentTime());
#ifdef WITH_LAPACK
   if(!quiet)
      printf("  LAPACK is used.\n");
   char JOBZ = 'A';
   int info;
   double *A = new double[N_ID*N_ID];
   int dimLA = N_ID;
   for(int i = 0; i < N_ID; i++)
     for(int j = 0; j < N_ID; j++)
       A[i*N_ID+j] = D[validpheno[j]][validpheno[i]];
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
   for(int i = 0; i < N_ID; i++) EV[i] = S[i];
   delete S;
   // VT[k*N_ID+j] stores jth eigenvector. k: id; j: marker
   for(int j = 0; j < N_ID; j++)
      for(int k = 0; k < N_ID; k++)
         UT[j][k] = VT[k*N_ID+j];
   delete []VT;
#else
   SVD svd;
   Matrix MD(N_ID, N_ID);
   for(int i = 0; i < N_ID; i++)
      for(int j = 0; j < N_ID; j++)
         MD[i][j] = D[validpheno[i]][validpheno[j]];
   svd.Decompose(MD);
   printf("done\n");
   if(svd.n == 0) return;
   for(int i = 0; i < N_ID; i++) EV[i] = svd.w[i];
   // svd.v[k][idx[N_ID-1-j]] stores jth eigenvector
   for(int j = 0; j < N_ID; j++)
      for(int k = 0; k < N_ID; k++)
         UT[j][k] = svd.v[k][j];
#endif
   if(!quiet){
      printf("SVD ends at %s", currentTime());
      printf("Eigenvalues range from %G to %G.\n", EV.Min(), EV.Max());
   }
   if(svdoutFlag){
      String svdoutfile = prefix;
      svdoutfile.Add(".svd");
      FILE *fp = fopen(svdoutfile, "wt");
      if(fp == NULL) error("Cannot open %s to write.", (const char*)svdoutfile);
      fprintf(fp, "%s_%s",
         (const char*)ped[ID[0]].famid, (const char*)ped[ID[0]].pid);
      for(int i = 1; i < N_ID; i++)
         fprintf(fp, "\t%s_%s", (const char*)ped[ID[i]].famid, (const char*)ped[ID[i]].pid);
      fprintf(fp, "\n");
      fprintf(fp, "%lf", EV[0]);
      for(int i = 1; i < N_ID; i++)
         if(EV[i] > 1E-6 || EV[i] < -1E-6)
            fprintf(fp, "\t%lf", EV[i]);
         else
            fprintf(fp, "\t%G", EV[i]);
      fprintf(fp, "\n");
      for(int i = 0; i < N_ID; i++){
         fprintf(fp, "%G", UT[i][0]);
         for(int j = 1; j < N_ID; j++)
            fprintf(fp, "\t%G", UT[i][j]);
         fprintf(fp, "\n");
      }
      for(int i = 0; i < N_ID; i++){
         fprintf(fp, "%G", D[i][0]);
         for(int j = 1; j < N_ID; j++)
            fprintf(fp, "\t%G", D[i][j]);
         fprintf(fp, "\n");
      }
      fclose(fp);
      if(!quiet)
      printf("Eigenvalues, eigenvectors and genetic similarity are saved in file %s\n",
         (const char*)svdoutfile);
      exit(0);
   }
}



