//////////////////////////////////////////////////////////////////////
// LMMSCORE.cpp
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

#include "analysis.h"
#include <math.h>
#include "MathStats.h"
#include "QuickIndex.h"
#include "MathCholesky.h"
#include "MathSVD.h"
#include "Random.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

#ifdef WITH_LAPACK
extern "C" void dgesdd_(char*, int*, int*, double*, int*, double*, double *, int*, double*, int*, double*, int*, int*, int*);
extern "C" void dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);
#endif

void Engine::GenomeScan64Bit()
{
   if(snpName.Length()==0 || chromosomes.Length()==0 || bp.Length()==0){
      warning("Incomplete SNP information.");
      return;
   }
   int N_ID = ID.Length();
   int N_Trait = traits.Length();
   double pvalue, h2, statistic, var_snp, Tm;
   int id, pos;
   if(chisqFilter==_NAN_){
      chisqs = new Vector[N_Trait];
      for(int t = 0; t < N_Trait; t++)
         chisqs[t].Dimension(0);
   }
   IntArray InvalidTrait(0);
   double *Sum_VR = new double[N_Trait];
   Matrix tVR(N_ID, N_Trait);
   for(int t = 0; t < N_Trait; t++){
      Sum_VR[t] = 0;
      for(int i = 0; i < N_ID; i++){
         if(VR[t][i] == _NAN_){
            InvalidTrait.Push(t);
            break;
         }
         Sum_VR[t] += VR[t][i];
         tVR[i][t] = VR[t][i];
      }
   }
   if(InvalidTrait.Length()){
      printf("The following traits have issues and are skipped in the anaysis:");
      for(int t = 0; t < InvalidTrait.Length(); t++)
         printf(" %s", (const char*)ped.traitNames[traits[t]]);
      printf("\n");
   }

int thread = 0;
FILE **fps;

#ifdef _OPENMP
   printf("Scanning autosomes with %d CPU cores...\n", defaultMaxCoreCount);
   fps = new FILE *[defaultMaxCoreCount];
   StringArray outfiles(defaultMaxCoreCount);
   for(int c = 1; c < defaultMaxCoreCount; c++){
      outfiles[c].Copy(prefix);
      outfiles[c] += (c+1);
      outfiles[c].Add("$$$.mtscore");
      fps[c] = fopen((const char*)outfiles[c], "wt");
   }
#else
   fps = new FILE *[1];
#endif
   double v0,v1, RR, RA;
   String outfile(prefix);
   outfile.Add("mtscore.txt");
   fps[0] = fopen((const char*)outfile, "wt");
   fprintf(fps[0], "SNP\tTrait\tChr\tPos\tLabelA\tLabela\tFreqA\tN\tBeta\tSE\tChisq\tH2\tPvalue\n");
   double maxR, lxy, lxx, vn;
   unsigned long long int word;
   unsigned char *pchar;
   unsigned char byte;
   char revbase[256];
   for(int i = 0; i < 8; i++)
      revbase[base[i]] = i;
   char rightmost[256];
   for(int i = 0; i < 256; i++)
      rightmost[i] = revbase[i&(-i)];
   const int CACHESIZE = 2;   // 128 SNPs at a time
   const int CACHESIZE_SNP = (CACHESIZE<<6);
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
   private(pos, var_snp, Tm, id, word, pchar, byte,\
   pvalue, statistic, h2, v0,v1,vn, RR, RA, maxR, lxy, lxx, thread)
{
#endif
   int *G2 = new int[CACHESIZE_SNP+1];
   int *G1 = new int[CACHESIZE_SNP+1];
   int *GM = new int[CACHESIZE_SNP+1];
   float *scores_local[CACHESIZE_SNP];
   for(int m = 0; m < CACHESIZE_SNP; m++)
      scores_local[m] = new float [N_Trait];
   float *vr = new float [N_Trait];
   float *vr2 = new float [N_Trait];
   double freq[CACHESIZE_SNP];
   int missingCount[CACHESIZE_SNP];
#ifdef _OPENMP
   thread = omp_get_thread_num();
   #pragma omp for
#endif
   for(int blockb = 0; blockb < longCount; blockb += CACHESIZE){
      int bMax = (blockb >= longCount-CACHESIZE) ? longCount: blockb+CACHESIZE;
      int mMax = (blockb >= longCount-CACHESIZE)? markerCount-(blockb<<6): ((bMax-blockb)<<6);
      for(int j = 0; j < CACHESIZE_SNP; j++){
         freq[j] = 0.0;
         missingCount[j] = 0;
      }
      for(int i = 0; i < N_ID; i++){
         id = geno[ID[i]];
         for(int b = blockb; b < bMax; b++){
            int bb = (b-blockb)<<6;
            word = LG[0][id][b] & LG[1][id][b];  // AA
            pchar = (unsigned char*)&word;
            for(int k = 0; k < 8; k++)
               for(byte = pchar[k]; byte; byte &= (byte-1))
                  freq[bb+(k<<3)+rightmost[byte]] += 1;
            word = (~LG[0][id][b]) & LG[1][id][b];  // Aa
            pchar = (unsigned char*)&word;
            for(int k = 0; k < 8; k++)
               for(byte = pchar[k]; byte; byte &= (byte-1))
                  freq[bb+(k<<3)+rightmost[byte]] += 0.5;
            word = (~LG[0][id][b]) & (~LG[1][id][b]);  // Miss
            pchar = (unsigned char*)&word;
            for(int k = 0; k < 8; k++)
               for(byte = pchar[k]; byte; byte &= (byte-1))
                  missingCount[bb+(k<<3)+rightmost[byte]] ++;
         }  // end of b
      }  // end of i
      for(int m = 0; m < mMax; m++){
         freq[m] /= (N_ID - missingCount[m]);
         for(int t = 0; t < N_Trait; t++)
            scores_local[m][t] = -Sum_VR[t];
      }
      for(int i = 0; i < N_ID; i++){   // most computationally intensive
         id = geno[ID[i]];
         for(int t = 0; t < N_Trait; t++){
            vr[t] = tVR[i][t];
            vr2[t] = tVR[i][t]*2;
         }
         GM[0] = G1[0] = G2[0] = 0;
         for(int b = blockb; b < bMax; b++){
            int bb = (b-blockb)<<6;
            word = LG[0][id][b] & LG[1][id][b];  // AA
            pchar = (unsigned char*)&word;
            for(int k = 0; k < 8; k++)
               for(byte = pchar[k]; byte; byte &= (byte-1))
                  G2[++G2[0]] = bb+(k<<3)+rightmost[byte];   // AA
            word = (~LG[0][id][b]) & LG[1][id][b];  // Aa
            pchar = (unsigned char*)&word;
            for(int k = 0; k < 8; k++)
               for(byte = pchar[k]; byte; byte &= (byte-1))
                  G1[++G1[0]] = bb+(k<<3)+rightmost[byte];   // Aa
            word = (~LG[0][id][b]) & (~LG[1][id][b]);  // Miss
            pchar = (unsigned char*)&word;
            for(int k = 0; k < 8; k++)
               for(byte = pchar[k]; byte; byte &= (byte-1))
                  GM[++GM[0]] = bb+(k<<3)+rightmost[byte];   // Miss
         }  // end of b
         for(int k = 0; k < GM[0]; k++){
            int m = GM[k+1];
            for(int t = 0; t < N_Trait; t++)
               scores_local[m][t] += vr[t]*freq[m];
         }
         for(int k = 0; k < G2[0]; k++){
            int m = G2[k+1];
            for(int t = 0; t < N_Trait; t++)
               scores_local[m][t] += vr2[t];
         }
         for(int k = 0; k < G1[0]; k++){
            int m = G1[k+1];
            for(int t = 0; t < N_Trait; t++) // trait loop
               scores_local[m][t] += vr[t];
         }  // end of marker loop
      }  // end of individual loop
      if(chisqFilter!=_NAN_){  // filter P-value, faster computation
         for(int m = 0; m < mMax; m++){
            pos = (blockb<<6) + m;
            var_snp = 2*freq[m]*(1-freq[m]);
            if(var_snp < 1E-10) continue;
            Tm = sqrt(chisqFilter*var_snp);
            maxR = -1;
            for(int t = 0; t < N_Trait; t++){// computationally intensive for large N_Trait
               if(Sum_VR[t]==_NAN_) continue;
               if((scores_local[m][t] > -Tm) && (scores_local[m][t] < Tm)) continue;
               statistic = scores_local[m][t] * scores_local[m][t] / var_snp;
               if(statistic==_NAN_) continue;
               if(statistic > maxR){
                  maxR = statistic;
               }
               pvalue = ndist(sqrt(statistic))*2;
               h2 = scores_local[m][t] * scores_local[m][t] / varVR[t] / var_snp;
               h2 = h2/((1+lambda0[t])/tau0[t]+h2);
               fprintf(fps[thread], "%s\t%s\t%d\t%d\t%s\t%s\t%.3lf\t%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.3G\n",
                  (const char*)snpName[pos],
                  (const char*)ped.traitNames[traits[t]],
                  chromosomes[pos], bp[pos],
                  (const char*)alleleLabel[0][pos],
                  (const char*)alleleLabel[1][pos],
                  freq[m], N_ID-missingCount[m],
                  scores_local[m][t] / sqrt(varVR[t]) / var_snp,
                  1/sqrt(varVR[t]*var_snp), statistic, h2, pvalue);
            }
         }  // end of m
      }
      else // print at all SNPs
         for(int m = 0; m < mMax; m++){
            pos = (blockb<<6) + m;
            var_snp = 2*freq[m]*(1-freq[m]);
            if(var_snp < 1E-10) continue;
            for(int t = 0; t < N_Trait; t++){
               statistic = scores_local[m][t] * scores_local[m][t] / var_snp;
               if(statistic==_NAN_) continue;
               pvalue = ndist(sqrt(statistic))*2;
               h2 = scores_local[m][t] * scores_local[m][t] / varVR[t] / var_snp;
               h2 = h2/((1+lambda0[t])/tau0[t]+h2);
               fprintf(fps[thread], "%s\t%s\t%d\t%d\t%s\t%s\t%.3lf\t%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.3G\n",
                  (const char*)snpName[pos],
                  (const char*)ped.traitNames[traits[t]],
                  chromosomes[pos],
                  bp[pos],
                  (const char*)alleleLabel[0][pos],
                  (const char*)alleleLabel[1][pos],
                  freq[m], N_ID-missingCount[m],
                  scores_local[m][t] / sqrt(varVR[t]) / var_snp,
                  1/sqrt(varVR[t]*var_snp),
                  statistic, h2, pvalue);
            }  // end of t
         }  // end of  m
   }  // end of blockb
   for(int m = 0; m < CACHESIZE_SNP; m++)
      delete []scores_local[m];
   delete []GM;
   delete []G1;
   delete []G2;
   delete []vr;
   delete []vr2;
   fclose(fps[thread]);
#ifdef _OPENMP
}  // extra bracket for omp
   fps[0] = fopen(outfile, "at");
   char buffer[60000];
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
   printf("Association scan results saved in file %s\n", (const char*)outfile);
   for(int t = 0; t < N_Trait; t++)
      delete VR[t];
   delete []VR;
   if(varVR) delete []varVR;
   if(Sum_VR) delete []Sum_VR;
}

void Engine::GenomeScan()
{
   if(snpName.Length()==0 || chromosomes.Length()==0 || bp.Length()==0){
      warning("Incomplete SNP information.");
      return;
   }
   int N_ID = ID.Length();
   int N_Trait = traits.Length();
   double pvalue, h2, statistic, var_snp, Tm;
   int id, pos;
   if(chisqFilter==_NAN_){
      chisqs = new Vector[N_Trait];
      for(int t = 0; t < N_Trait; t++)
         chisqs[t].Dimension(0);
   }
   IntArray InvalidTrait(0);
   double *Sum_VR = new double[N_Trait];
   Matrix tVR(N_ID, N_Trait);
   for(int t = 0; t < N_Trait; t++){
      Sum_VR[t] = 0;
      for(int i = 0; i < N_ID; i++){
         if(VR[t][i] == _NAN_){
            InvalidTrait.Push(t);
            break;
         }
         Sum_VR[t] += VR[t][i];
         tVR[i][t] = VR[t][i];
      }
   }
   if(InvalidTrait.Length()){
      printf("The following traits have issues and are skipped in the anaysis:");
      for(int t = 0; t < InvalidTrait.Length(); t++)
         printf(" %s", (const char*)ped.traitNames[traits[t]]);
      printf("\n");
   }

int thread = 0;
FILE **fps;
//FILE **ft;

#ifdef _OPENMP
   printf("Scanning autosomes with %d CPU cores...\n", defaultMaxCoreCount);
   fps = new FILE *[defaultMaxCoreCount];
   StringArray outfiles(defaultMaxCoreCount);
   for(int c = 1; c < defaultMaxCoreCount; c++){
      outfiles[c].Copy(prefix);
      outfiles[c] += (c+1);
      outfiles[c].Add("$$$.mtscore");
      fps[c] = fopen((const char*)outfiles[c], "wt");
   }
   /*
   StringArray tfiles(CoreInUse);
   if(chisqFilter!=_NAN_){
      ft = new FILE *[CoreInUse];
      for(int c = 1; c < CoreInUse; c++){
         tfiles[c].Copy(prefix);
         tfiles[c] += (c+1);
         tfiles[c].Add("$$$.tped");
         ft[c] = fopen((const char*)tfiles[c], "wt");
      }
   } */
#else
   fps = new FILE *[1];
//   if(chisqFilter!=_NAN_)
//      ft = new FILE *[1];
#endif
   double v0,v1, RR, RA;
   String outfile(prefix);
   outfile.Add("mtscore.txt");
   fps[0] = fopen((const char*)outfile, "wt");
   fprintf(fps[0], "SNP\tTrait\tChr\tPos\tLabelA\tLabela\tFreqA\tN\tBeta\tSE\tChisq\tH2\tPvalue\n");
/*
   String tfile(prefix);
   if(chisqFilter!=_NAN_){
      tfile.Add(".tped");
      ft[0] = fopen((const char*)tfile, "wt");
   }
*/
   int maxT;
   double maxR, lxy, lxx, vn;
   unsigned short int word;
   char revbase[65536];
   for(int i = 0; i < 16; i++)
      revbase[shortbase[i]] = i;
   char rightmost[65536];
   for(int i = 0; i < 65536; i++)
      rightmost[i] = revbase[i&(-i)];
   const int CACHESIZE = 4;   // 64 SNPs at a time
   const int CACHESIZE_SNP = CACHESIZE*16;
   double freq[CACHESIZE_SNP];
   int missingCount[CACHESIZE_SNP];
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
   private(freq, missingCount, pos, var_snp, Tm, id, word,\
   pvalue, statistic, h2, v0,v1,vn, RR, RA, maxT, maxR, lxy, lxx, thread)
{
#endif
   int *G2 = new int[CACHESIZE_SNP+1];
   int *G1 = new int[CACHESIZE_SNP+1];
   int *GM = new int[CACHESIZE_SNP+1];
   float *scores_local[CACHESIZE_SNP];
   for(int m = 0; m < CACHESIZE_SNP; m++)
      scores_local[m] = new float [N_Trait];
   float *vr = new float [N_Trait];
   float *vr2 = new float [N_Trait];
#ifdef _OPENMP
   thread = omp_get_thread_num();
   #pragma omp for
#endif
   for(int blockb = 0; blockb < shortCount; blockb += CACHESIZE){
      int bMax = (blockb > shortCount-CACHESIZE) ? shortCount: blockb+CACHESIZE;
      int mMax = (blockb > shortCount-CACHESIZE)? markerCount-blockb*16: (bMax-blockb)*16;
      for(int j = 0; j < CACHESIZE_SNP; j++){
         freq[j] = 0.0;
         missingCount[j] = 0;
      }
      for(int i = 0; i < N_ID; i++){
         id = geno[ID[i]];
         for(int b = blockb; b < bMax; b++){
            int bb = (b-blockb)*16;
            for(word = GG[0][id][b] & GG[1][id][b]; word; word &= (word-1))
               freq[bb+rightmost[word]] += 1;   // AA
            for(word = (~GG[0][id][b]) & GG[1][id][b]; word; word &= (word-1))
               freq[bb+rightmost[word]] += 0.5; // Aa
            for(word = (~GG[0][id][b]) & (~GG[1][id][b]) & 65535; word; word &= (word-1))
               missingCount[bb+rightmost[word]] ++;   // missing
         }  // end of b
      }  // end of i
      for(int m = 0; m < mMax; m++){
         freq[m] /= (N_ID - missingCount[m]);
         for(int t = 0; t < N_Trait; t++)
            scores_local[m][t] = -Sum_VR[t];
      }
      for(int i = 0; i < N_ID; i++){   // most computationally intensive
         id = geno[ID[i]];
         for(int t = 0; t < N_Trait; t++){
            vr[t] = tVR[i][t];
            vr2[t] = tVR[i][t]*2;
         }
         GM[0] = G1[0] = G2[0] = 0;
         for(int b = blockb; b < bMax; b++){
            int bb = (b-blockb)*16;
            for(word = GG[0][id][b] & GG[1][id][b]; word; word &= (word-1))
               G2[++G2[0]] = bb+rightmost[word];   // AA
            for(word = (~GG[0][id][b]) & GG[1][id][b]; word; word &= (word-1))
               G1[++G1[0]] = bb+rightmost[word];   // Aa
            for(word = (~GG[0][id][b]) & (~GG[1][id][b]) & 65535; word; word &= (word-1))
               GM[++GM[0]] = bb+rightmost[word];   // Missing
         }  // end of b
         for(int k = 0; k < GM[0]; k++){
            int m = GM[k+1];
            for(int t = 0; t < N_Trait; t++)
               scores_local[m][t] += vr[t]*freq[m];
         }
         for(int k = 0; k < G2[0]; k++){
            int m = G2[k+1];
            for(int t = 0; t < N_Trait; t++)
               scores_local[m][t] += vr2[t];
         }
         for(int k = 0; k < G1[0]; k++){
            int m = G1[k+1];
            for(int t = 0; t < N_Trait; t++) // trait loop
               scores_local[m][t] += vr[t];
         }  // end of marker loop
      }  // end of individual loop
      if(chisqFilter!=_NAN_){  // filter P-value, faster computation
         for(int m = 0; m < mMax; m++){
            pos = blockb*16 + m;
            var_snp = 2*freq[m]*(1-freq[m]);
            if(var_snp < 1E-10) continue;
            Tm = sqrt(chisqFilter*var_snp);
            maxT = -1;
            maxR = -1;
            for(int t = 0; t < N_Trait; t++){// computationally intensive for large N_Trait
               if(Sum_VR[t]==_NAN_) continue;
               if((scores_local[m][t] > -Tm) && (scores_local[m][t] < Tm)) continue;
               statistic = scores_local[m][t] * scores_local[m][t] / var_snp;
               if(statistic==_NAN_) continue;
               if(statistic > maxR){
                  maxR = statistic;
                  maxT = t;
               }
               pvalue = ndist(sqrt(statistic))*2;
               h2 = scores_local[m][t] * scores_local[m][t] / varVR[t] / var_snp;
               h2 = h2/((1+lambda0[t])/tau0[t]+h2);
               fprintf(fps[thread], "%s\t%s\t%d\t%d\t%s\t%s\t%.3lf\t%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.3G\n",
                  (const char*)snpName[pos],
                  (const char*)ped.traitNames[traits[t]],
                  chromosomes[pos], bp[pos],
                  (const char*)alleleLabel[0][pos],
                  (const char*)alleleLabel[1][pos],
                  freq[m], N_ID-missingCount[m],
                  scores_local[m][t] / sqrt(varVR[t]) / var_snp,
                  1/sqrt(varVR[t]*var_snp), statistic, h2, pvalue);
            }
         }  // end of m
      }
      else // print at all SNPs
         for(int m = 0; m < mMax; m++){
            pos = blockb*16 + m;
            var_snp = 2*freq[m]*(1-freq[m]);
            if(var_snp < 1E-10) continue;
            for(int t = 0; t < N_Trait; t++){
               statistic = scores_local[m][t] * scores_local[m][t] / var_snp;
               if(statistic==_NAN_) continue;
               pvalue = ndist(sqrt(statistic))*2;
               h2 = scores_local[m][t] * scores_local[m][t] / varVR[t] / var_snp;
               h2 = h2/((1+lambda0[t])/tau0[t]+h2);
               fprintf(fps[thread], "%s\t%s\t%d\t%d\t%s\t%s\t%.3lf\t%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.3G\n",
                  (const char*)snpName[pos],
                  (const char*)ped.traitNames[traits[t]],
                chromosomes[pos],
                  bp[pos],
                  (const char*)alleleLabel[0][pos],
                  (const char*)alleleLabel[1][pos],
                  freq[m], N_ID-missingCount[m],
                  scores_local[m][t] / sqrt(varVR[t]) / var_snp,
                  1/sqrt(varVR[t]*var_snp),
                  statistic, h2, pvalue);
            }  // end of t
         }  // end of  m
   }  // end of blockb
   for(int m = 0; m < CACHESIZE_SNP; m++)
      delete []scores_local[m];
   delete []GM;
   delete []G1;
   delete []G2;
   delete []vr;
   delete []vr2;
   fclose(fps[thread]);
//   if(chisqFilter!=_NAN_)
//      fclose(ft[thread]);
#ifdef _OPENMP
}  // extra bracket for omp
   fps[0] = fopen(outfile, "at");
   char buffer[60000];
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
//   if(chisqFilter!=_NAN_)
//      fclose(ft[0]);
   printf("Association scan results saved in file %s\n", (const char*)outfile);
   for(int t = 0; t < N_Trait; t++)
      delete VR[t];
   delete []VR;
   if(varVR) delete []varVR;
   if(Sum_VR) delete []Sum_VR;
}

void Engine::DosageScan(FILE *fp)
{
   int N_ID = ID.Length();
   int N_Trait = traits.Length();
   double pvalue, h2, freq, statistic, score_local, var_snp, Tm;
   double *dosage = new double[N_ID];
   double *Sum_VR = new double[N_Trait];
   StringArray tokens;
   String line;
   if(dfamfile==""){
      int m = dosagefile.Find(".");
      if(m == -1) error("Please use PLINK dosage format as input.");
      String filename = dosagefile.SubStr(0, m);
      dfamfile=filename;
      dfamfile.Add(".dfam");
   }
   IFILE input = ifopen(dfamfile, "rt");
   if(input == NULL){
      String secondtry(dfamfile);
      secondtry.Add(".gz");
      input = ifopen(secondtry, "rt");
      if(input == NULL)
         error("ID file for dosage data %s cannot be opened", (const char*)secondtry);
      dfamfile=secondtry;
   }
   printf("Read in ID file for dosage data %s...\n", (const char*)dfamfile);
   StringArray uniqueName(N_ID);
   for(int i = 0; i < N_ID; i++){
      uniqueName[i] = ped[ID[i]].famid;
      uniqueName[i].Add("_");
      uniqueName[i].Add(ped[ID[i]].pid);
   }
   String tempID;
   IntArray key(N_ID);
   key.Set(-1);
   int row = 0;
   int pos;
   int missingCount=0;
   while(!ifeof(input)){
      tokens.Clear();
      line.ReadLine(input);
      tokens.AddTokens(line);
      if(tokens.Length() < 2) continue;
      tempID = tokens[0];
      tempID.Add("_");
      tempID.Add(tokens[1]);
      pos = uniqueName.Find(tempID);
      if(pos!=-1)
         key[pos] = row;
      else
         missingCount++;
      row++;
   }
   ifclose(input);
   int N_famID = row;
   int **valid;
   valid = new int* [N_Trait];
   for(int t = 0; t < N_Trait; t++){
      valid[t] = new int [N_ID+1];
      valid[t][0] = 0;
   }
   IntArray SampleSize(N_Trait);
   if(FixedEff){
      SampleSize.Zero();
      for(int t = 0; t < N_Trait; t++)
         for(int i = 0; i < N_ID; i++)
            if(key[i]!=-1 && VR[t][i]!=0){
               SampleSize[t]++;
               valid[t][0] ++;
               valid[t][valid[t][0]] = i;
            }
   }else{ // for random effect, assuming non-missing phenotypes
      SampleSize.Set(N_ID - missingCount);
      for(int t = 0; t < N_Trait; t++)
         for(int i = 0; i < N_ID; i++)
            if(key[i]!=-1){
               valid[t][0] ++;
               valid[t][valid[t][0]] = i;
            }
   }
   for(int t = 0; t < N_Trait; t++){
      Sum_VR[t] = 0;
      for(int i = 0; i < valid[t][0]; i++)
         Sum_VR[t] += VR[t][valid[t][i+1]];
   }
   for(int i = 0; i < N_ID; i++)
      if(key[i]==-1)
         printf("Cannot find Fam %s ID %s in the ID file %s\n",
         (const char*)ped[ID[i]].pid, (const char*)ped[ID[i]].pid,
         (const char*)dfamfile);
   StringIntHash markerLookup;
   IntArray markerChr(0);
   IntArray markerPos(0);
   int chr;
   if(dmapfile!=""){
      input = ifopen(dmapfile, "rt");
      if(input == NULL)
         error("Cannot open %s to read", (const char*)dmapfile);
      printf("Read in map file for dosage data %s...\n", (const char*)dmapfile);
      markerCount = 0;
      while(!ifeof(input)){
         tokens.Clear();
         line.ReadLine(input);
         tokens.AddTokens(line);
         if(tokens.Length()>=4){
            markerLookup.SetInteger(tokens[1], markerCount);
            markerPos.Push(tokens[3].AsInteger());
            if(tokens[0]=="X")
               chr = SEXCHR;
            else if(tokens[0]=="Y")
               chr = SEXCHR+1;
            else if(tokens[0]=="XY")
               chr = SEXCHR+2;
            else if(tokens[0]=="MT")
               chr = SEXCHR+3;
            else
               chr = tokens[0].AsInteger();
            markerChr.Push(chr);
         }
         markerCount++;
      }
      ifclose(input);
   }
   input = ifopen(dosagefile, "rt");
   if(input == NULL)
      error("Cannot open %s to read", (const char*)dosagefile);
   tokens.Clear();
   line.ReadLine(input);
   tokens.AddTokens(line);
   StringArray dosagefiles(0);
   if(tokens.Length()==1){
      dosagefiles.Push(tokens[0]);
      while(!ifeof(input)){
         tokens.Clear();
         line.ReadLine(input);
         tokens.AddTokens(line);
         if(tokens.Length() && tokens[0]!="")
            dosagefiles.Push(tokens[0]);
      }
   }else{
      dosagefiles.Push(dosagefile);
      ifrewind(input);
   }
   ifclose(input);
   for(int f = 0; f < dosagefiles.Length(); f++){
      input = ifopen(dosagefiles[f], "rt");
      if(input == NULL)
         error("Cannot open %s to read", (const char*)dosagefiles[f]);
      printf("Read in dosage data %s and start genome scan...\n", (const char*)dosagefiles[f]);

   row = 0;
   while(!ifeof(input)){
      tokens.Clear();
      line.ReadLine(input);
      tokens.AddTokens(line);
      if(tokens.Length() < N_famID+3) continue;
      for(int i = 0; i < N_ID; i++)
         if(key[i]!=-1)
            dosage[i] = tokens[3+key[i]].AsDouble();
      row++;
      for(int t = 0; t < N_Trait; t++){// computationally intensive for large N_Trait
         freq = var_snp = score_local = 0;
         for(int i = 0; i < valid[t][0]; i++){
            pos = valid[t][i+1];
            score_local += dosage[pos] * VR[t][pos];
            var_snp += dosage[pos] * dosage[pos];
            freq += dosage[pos];
         }
         freq /= (valid[t][0]*2);
         var_snp /= valid[t][0];
         var_snp -= 4*freq*freq;
         if(var_snp < 1E-10 || freq < 1E-10 || freq > 1-1E-10) continue;
         score_local -= 2*freq*Sum_VR[t];
         statistic = score_local * score_local / var_snp;
         if((chisqFilter!=_NAN_) && statistic < chisqFilter) continue;
         pvalue = ndist(sqrt(statistic))*2;
         h2 = score_local * score_local / varVR[t] / var_snp;
         h2 = h2/((1+lambda0[t])/tau0[t]+h2);
         fprintf(fp, "%s", (const char*)tokens[0]);
         if(N_Trait>1)
            fprintf(fp, "\t%s", (const char*)ped.traitNames[traits[t]]);
         if(dmapfile!=""){
            int idx = markerLookup.Integer(tokens[0]);
            if(idx!=-1)
               fprintf(fp, "\t%d\t%.6lf", markerChr[idx], markerPos[idx]*0.000001);
            else
               fprintf(fp, "\tNA\tNA");
         }
         fprintf(fp, "\t%c\t%c\t%.4lf\t%.4lf\t%d",
            tokens[1][0], tokens[2][0], freq,
            var_snp/(2*freq*(1-freq)), SampleSize[t]);
         if(statistic != _NAN_)
            fprintf(fp, "\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.3G\n",
               score_local / sqrt(varVR[t]) / var_snp,
               1/sqrt(varVR[t]*var_snp),
               statistic, h2, pvalue);
         else
            fprintf(fp, "\t%s\t%s\n", "NA", "NA");
      }
   }
   ifclose(input);
   printf("%d SNPs processed\n", row);
   }
   if(Sum_VR) delete []Sum_VR;
   for(int t = 0; t < N_Trait; t++)
      delete VR[t];
   delete []VR;
   if(varVR) delete []varVR;
   delete []dosage;
   for(int t = 0; t < N_Trait; t++)
      delete valid[t];
   delete []valid;
}


void Engine::GenomeScanWithPermutation(FILE *fp)
{
   int N_ID = ID.Length();
   int N_Trait = traits.Length();
   Random rand;
   int **perm = new int * [permuteCount];
   for(int p = 0; p < permuteCount; p++)
      perm[p] = new int [N_ID];
   int pos, swap;
   for(int p = 0; p < permuteCount; p++){
      for(int j = 0; j < N_ID; j++) perm[p][j] = j;
      for(int j = N_ID; j > 1; j--){
         pos = rand.NextInt()%j;
         swap = perm[p][j-1];
         perm[p][j-1] = perm[p][pos];
         perm[p][pos] = swap;
      }
   }

   Matrix maxStat(N_Trait, permuteCount);
   maxStat.Zero();
   double pvalue, h2, freq[16], statistic, score_local, score_perm, var_snp, Tm, tempstat;
   int missingCount[16], id, AA, Aa, missing;
   chisqs = new Vector[N_Trait];
   IntArray *chiPos = new IntArray[N_Trait];
   for(int t = 0; t < N_Trait; t++){
      chisqs[t].Dimension(0);
      chiPos[t].Dimension(0);
   }                         
   int *G2[16], *G1[16], *GM[16];
   for(int m = 0; m < 16; m++){
      G2[m] = new int[N_ID+1];
      G1[m] = new int[N_ID+1];
      GM[m] = new int[N_ID+1];
   }
   double Sum_VR = 0;
   for(int i = 0; i < N_ID; i++)
      Sum_VR += VR[0][i];
   for(int b = 0; b < shortCount; b++){ // loop over short integer
      for(int m = 0; m < 16; m++){
         freq[m] = 0;
         missingCount[m] = 0;
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
               missingCount[m] ++;
            }else if(AA & shortbase[m]){
               G2[m][0] ++;
               G2[m][G2[m][0]] = i;
               freq[m] += 1.0;
            }else if(Aa & shortbase[m]){
               G1[m][0] ++;
               G1[m][G1[m][0]] = i;
               freq[m] += 0.5;
            }
         }
      }
      for(int m = 0; m < 16; m++)
         freq[m] /= (N_ID-missingCount[m]);
      for(int m = 0; m < 16; m++){
         pos = b*16 + m;
         if(pos >= markerCount) break;
         var_snp = 2*freq[m]*(1-freq[m]);
         if(var_snp < 1E-10) continue;
         if(chisqFilter!=_NAN_)
            Tm = sqrt(chisqFilter*var_snp);
         for(int t = 0; t < N_Trait; t++){
            score_local = -Sum_VR;
            for(int i = 0; i < GM[m][0]; i++)
               score_local += VR[t][GM[m][i+1]];
            score_local *= freq[m];
            for(int i = 0; i < G2[m][0]; i++)
               score_local += VR[t][G2[m][i+1]];
            score_local += score_local;
            for(int i = 0; i < G1[m][0]; i++)
               score_local += VR[t][G1[m][i+1]];
            for(int p = 0; p < permuteCount; p++){// computationally intensive for large N_Trait
               score_perm = -Sum_VR;
               for(int i = 0; i < GM[m][0]; i++)
                  score_perm += VR[t][perm[p][GM[m][i+1]]];
               score_perm *= freq[m];
               for(int i = 0; i < G2[m][0]; i++)
                  score_perm += VR[t][perm[p][G2[m][i+1]]];
               score_perm += score_perm;
               for(int i = 0; i < G1[m][0]; i++)
                  score_perm += VR[t][perm[p][G1[m][i+1]]];
               tempstat = score_perm * score_perm / var_snp;
               if(tempstat > maxStat[t][p])
                  maxStat[t][p] = tempstat;
            }
            if((chisqFilter!=_NAN_) && (score_local > -Tm) && (score_local < Tm))
               continue;
            statistic = score_local * score_local / var_snp;
            chisqs[t].Push(statistic);
            chiPos[t].Push(pos);
         } // end of t
      } // end of m
   } // end of b
   for(int m = 0; m < 16; m++){
      delete []G2[m];
      delete []G1[m];
      delete []GM[m];
   }
   for(int p = 0; p < permuteCount; p++)
      delete perm[p];
   delete []perm;

   Vector maxStats;
   if(N_Trait>1){
      maxStats.Dimension(permuteCount);
      for(int p = 0; p < permuteCount; p++){
         maxStats[p] = 0;
         for(int t = 0; t < N_Trait; t++)
            if(maxStat[t][p] > maxStats[p])
               maxStats[p] = maxStat[t][p];
      }
      maxStats.Sort();
      printf("P value cutoff for genome-wide & pheno-wide significance is %.3G\n",
         ndist(sqrt(maxStats[int(0.95*permuteCount)]))*2);
   }
   for(int t = 0; t < N_Trait; t++){
      maxStat[t].Sort();
      printf("P value cutoff for genome-wide significance of %s is %.3G\n",
         (const char*)ped.traitNames[traits[t]],
         ndist(sqrt(maxStat[t][int(0.95*permuteCount)]))*2);
      for(int m = 0; m < chiPos[t].Length(); m++){
         int pos = chiPos[t][m];
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
         fprintf(fp, "\t%s\t%s",
         (const char*)alleleLabel[0][pos], (const char*)alleleLabel[1][pos]);
         fprintf(fp, "\t%.4lf\t%.3G", chisqs[t][m], ndist(sqrt(chisqs[t][m]))*2);
         pvalue = EmpP(chisqs[t][m], maxStat[t]);
         fprintf(fp, "\t%.3G", pvalue);
         if(N_Trait>1){
            pvalue = EmpP(chisqs[t][m], maxStats);
            fprintf(fp, "\t%.3G", pvalue);
         }
         fprintf(fp, "\n");
      }
   }
}

double Engine::EmpP(double stat, Vector & dist)
{
   int left, right, middle, rank;
   double pvalue;
   int permuteCount = dist.Length();
   if(stat > dist[permuteCount-1]) pvalue = 0.5 / permuteCount;
   else if(stat < dist[0]) pvalue = 1-0.5 / permuteCount;
   else{
         left = 0;
         right = permuteCount-1;
         middle = (left+right)/2;
         for(; left < right-1; middle=(left+right)/2){
            if(dist[middle] > stat) right = middle;
            else if(dist[middle] < stat) left = middle;
            else {
               rank = middle;
               break;
            }
         }
         if(left==right || left < right-1){
            if(left==right) rank = left;
            pvalue = 1-(rank + 0.5)/permuteCount;
         }else
            pvalue = 1-right*1.0/permuteCount;
   }
   return pvalue;
}


void Engine::PreScan_FixedEff()
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
   int N_ID = ID.Length();

   IntArray *missingpheno = new IntArray[N_Trait];
   for(int t = 0; t < N_Trait; t++)
      missingpheno[t].Dimension(0);
   IntArray *validpheno = new IntArray[N_Trait];
   for(int t = 0; t < N_Trait; t++)
      validpheno[t].Dimension(0);

   printf("Genotypes stored in %d integers for each of %d individuals used in analysis.\n",
      shortCount, N_ID);
   EV.Dimension(N_ID);
   UT = new double * [N_ID];
   for(int i = 0; i < N_ID; i++)
      UT[i] = new double [N_ID];
   if(normalization)
      printf("Inverse normal transformation is applied to phenotypes.\n");

   if(svdinfile!="")
      ReadSVD();
   lambda0.Dimension(N_Trait);
   lambda0.Zero();
   tau0.Dimension(N_Trait);

   Matrix Y(N_ID, N_Trait);
   int N_X = 1+N_Covariate;
   if(svdinfile!="")
      N_X += FixedEff;
   Matrix X(N_ID, N_X);
   for(int i = 0; i < N_ID; i++)
      X[i][0] = 1.0;
   for(int t = 0; t < N_Trait; t++)
      for(int i = 0; i < N_ID; i++)
         Y[i][t] = ped[ID[i]].traits[traits[t]];
   for(int i = 0; i < N_ID; i++)
      for(int j = 0; j < N_Covariate; j++)
         X[i][j+1] = ped[ID[i]].covariates[covariates[j]];
   if(svdinfile!="")
   for(int i = 0; i < N_ID; i++)
      for(int j = 0; j < FixedEff; j++)
         X[i][N_Covariate+1+j] = UT[j][i];
   for(int i = 0; i < N_ID; i++)
      delete UT[i];
   delete []UT;
   for(int t = 0; t < N_Trait; t++)
      for(int i = 0; i < N_ID; i++)
         if(Y[i][t] == _NAN_)
            missingpheno[t].Push(i);
         else
            validpheno[t].Push(i);
   Vector tV(N_ID);
   Cholesky chol;
   QuickIndex idx;
   if(normalization)
      for(int t = 0; t < N_Trait; t++){
         tV.Dimension(validpheno[t].Length());
         for(int i = 0; i < validpheno[t].Length(); i++)
            tV[i] = Y[validpheno[t][i]][t];
         idx.Index(tV);
         for(int i = 0; i < validpheno[t].Length();){
            int start = i, end = i + 1;
            while(end < validpheno[t].Length() && tV[idx[end]] == tV[idx[start]])
               end++;
            end --;
            double q = ninv((start + (end-start)/2.0 + 0.5)/validpheno[t].Length());
            for(int j = start; j <= end; j++)
               Y[validpheno[t][idx[j]]][t] = q;
            i = end + 1;
         }
      }

   Matrix *Info0 = new Matrix[N_Trait];
   Matrix beta0(N_Trait, N_X);
   for(int t = 0; t < N_Trait; t++){
      int N = validpheno[t].Length();
      Info0[t].Dimension(N_X, N_X);
      Info0[t].Zero();
      for(int i = 0; i < N_X; i++)
         for(int j = 0; j < N_X; j++)
            for(int k = 0; k < N; k++)
               Info0[t][i][j] += X[validpheno[t][k]][i] * X[validpheno[t][k]][j];
      if(chol.TryDecompose(Info0[t])==0){
         beta0[t].Zero();
         continue;
      }
      chol.Decompose(Info0[t]);
      tV.Dimension(N_X);
      tV.Zero();
      for(int i = 0; i < N_X; i++)
         for(int k = 0; k < N; k++)
            tV[i] += X[validpheno[t][k]][i] * Y[validpheno[t][k]][t];
      chol.BackSubst(tV);
      beta0[t] = chol.x;
   }
   if(ped.traitCount < 100){
      printf("Regression coefficients\n");
      printf("%-15s%10s", "Trait", "Mu");
      for(int i = 0; i < N_Covariate; i++)
         printf(" %9s", (const char*)ped.covariateNames[covariates[i]]);
      if(svdinfile!="")
         for(int i = 0; i < FixedEff; i++)
            printf("       PC%d", i+1);
      printf("\n");
      for(int t = 0; t < N_Trait; t++){
         printf("%-15s", (const char*)ped.traitNames[traits[t]]);
         for(int i = 0; i < N_Covariate+1; i++)
            printf(" %9.3lf", beta0[t][i]);
         if(svdinfile!="")
         for(int i = N_Covariate+1; i < N_X; i++)
            printf(" %9.3lf", beta0[t][i]);
         printf("\n");
      }
   }

   delete []Info0;
   VR = new double * [N_Trait];
   for(int t = 0; t < N_Trait; t++){
      VR[t] = new double [N_ID];
      for(int i = 0; i < N_ID; i++)
         VR[t][i] = 0;
   }
   varVR = new double [N_Trait];
   for(int t = 0; t < N_Trait; t++){
      for(int i = 0; i < validpheno[t].Length(); i++){
         VR[t][validpheno[t][i]] = Y[validpheno[t][i]][t];
         for(int j = 0; j < N_X; j++)
            VR[t][validpheno[t][i]] -= beta0[t][j]*X[validpheno[t][i]][j];
      }
      for(int i = 0; i < missingpheno[t].Length(); i++)
         VR[t][missingpheno[t][i]] = 0;
      varVR[t] = 0;
      for(int i = 0; i < N_ID; i++)
         varVR[t] += VR[t][i] * VR[t][i];

      if(varVR[t] < 1E-200){
         for(int i = 0; i < N_ID; i++)
            VR[t][i] = _NAN_;
         tau0[t] = _NAN_;
      }else{
         tau0[t] = validpheno[t].Length()/varVR[t];
         for(int i = 0; i < N_ID; i++)
            VR[t][i] *= tau0[t];
         varVR[t] *= tau0[t] * tau0[t];
         for(int i = 0; i < N_ID; i++)
            VR[t][i] /= sqrt(varVR[t]);
      }
   }
   if(missingpheno) delete []missingpheno;
   if(validpheno) delete []validpheno;
}


void Engine::PreScan_MTSCORE()
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
   int N_ID = ID.Length();

   IntArray *missingpheno = new IntArray[N_Trait];
   for(int t = 0; t < N_Trait; t++)
      missingpheno[t].Dimension(0);
   IntArray *validpheno = new IntArray[N_Trait];
   for(int t = 0; t < N_Trait; t++)
      validpheno[t].Dimension(0);

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
         if(Y[i][t] == _NAN_){
            missingpheno[t].Push(i);
            Y[i][t] = 0;
         }else
            validpheno[t].Push(i);
   }
   Vector tV(N_ID);
   Cholesky chol;
   QuickIndex idx;
   if(normalization)
      for(int t = 0; t < N_Trait; t++){
         tV.Dimension(validpheno[t].Length());
         for(int i = 0; i < validpheno[t].Length(); i++)
            tV[i] = Y[validpheno[t][i]][t];
         idx.Index(tV);
         for(int i = 0; i < validpheno[t].Length();){
            int start = i, end = i + 1;
            while(end < validpheno[t].Length() && tV[idx[end]] == tV[idx[start]])
               end++;
            end --;
            double q = ninv((start + (end-start)/2.0 + 0.5)/validpheno[t].Length());
            for(int j = start; j <= end; j++)
               Y[validpheno[t][idx[j]]][t] = q;
            i = end + 1;
         }
      }

   Matrix *Info0 = new Matrix[N_Trait];
   Matrix beta0(N_Trait, N_Covariate+1);
   if(N_Covariate==0) // chol is not needed
      for(int t = 0; t < N_Trait; t++){
         int N = validpheno[t].Length();
         Info0[t].Dimension(1, 1);
         Info0[t][0] = 0;
         beta0[t][0] = 0;
         for(int k = 0; k < N; k++)
            beta0[t][0] += Y[validpheno[t][k]][t];
         beta0[t][0] /= N;
      }
   else
   for(int t = 0; t < N_Trait; t++){
      int N = validpheno[t].Length();
      Info0[t].Dimension(N_Covariate+1, N_Covariate+1);
      for(int i = 0; i < N_Covariate+1; i++)
         for(int j = 0; j < N_Covariate+1; j++)
            for(int k = 0; k < N; k++)
               Info0[t][i][j] += X[validpheno[t][k]][i] * X[validpheno[t][k]][j];
      if(chol.TryDecompose(Info0[t])==0){
         beta0[t].Zero();
         continue;
      }
      chol.Decompose(Info0[t]);
      tV.Dimension(N_Covariate+1);
      tV.Zero();
      for(int i = 0; i < N_Covariate+1; i++)
         for(int k = 0; k < N; k++)
            tV[i] += X[validpheno[t][k]][i] * Y[validpheno[t][k]][t];
      chol.BackSubst(tV);
      beta0[t] = chol.x;
   }

   Matrix Resid0(N_ID, N_Trait);
   for(int t = 0; t < N_Trait; t++)
      for(int i = 0; i < N_ID; i++){
         Resid0[i][t] = Y[i][t];
         for(int j = 0; j < N_Covariate+1; j++)
            Resid0[i][t] -= beta0[t][j] * X[i][j];
      }

   printf("Genotypes stored in %d integers for each of %d individuals used in analysis.\n",
      shortCount, N_ID);
   EV.Dimension(N_ID);
   UT = new double * [N_ID];
   for(int i = 0; i < N_ID; i++)
      UT[i] = new double [N_ID];
   if(normalization)
      printf("Inverse normal transformation is applied to phenotypes.\n");

   ComputeMTSVD();

   lambda0.Dimension(N_Trait);
   tau0.Dimension(N_Trait);
   if(ped.traitCount < 100){
      printf("\nPolygenic parameter estimates\n");
      printf("%-15s %7s %7s %7s %7s %7s %9s",
         "TraitName", "N", "N_Iter", "Herit", "Lambda", "Tau", "Mu");
      for(int i = 0; i < N_Covariate; i++)
         printf(" %9s", (const char*)ped.covariateNames[covariates[i]]);
      printf("\n");
   }
   
   UX.Dimension(N_ID, N_Covariate+1);
   UX.Zero();
   for(int i = 0; i < N_ID; i++)
      for(int j = 0; j < N_Covariate+1; j++)
         for(int k = 0; k < N_ID; k++)
            UX[i][j] += UT[i][k] * X[k][j];

   Matrix UYs(N_ID, N_Trait);
   UYs.Zero();
   for(int i = 0; i < N_ID; i++)
      for(int t = 0; t < N_Trait; t++)
         for(int k = 0; k < validpheno[t].Length(); k++)
            UYs.data[i]->data[t] += UT[i][validpheno[t][k]] * Y.data[validpheno[t][k]]->data[t];

   Vector resid, uy(N_ID);
   double lxy, lxx, temp, diff, diffsq, similarity, meandiffsq, meansimilarity;
   double beta, mu;
   int N;
   Vector mV;
   for(int t = 0; t < N_Trait; t++){
      N = validpheno[t].Length();
      lxx = lxy = 0;
      resid.Dimension(N);
      tV.Dimension(N_Covariate+1);
      for(int i = 0; i < N; i++)
         resid[i] = Resid0[validpheno[t][i]][t];
      int iter=0, totalpair;
      if(noiterFlag){
         meandiffsq = meansimilarity = lxx = lxy = 0;
         for(int i = 0; i < N; i++)
            for(int j = i+1; j < N; j++){
               diff = resid.data[i] - resid.data[j];
               diffsq = diff * diff;
               similarity = D[validpheno[t][i]][validpheno[t][j]];
               meandiffsq += diffsq;
               meansimilarity += similarity;
               lxy += similarity * diffsq;
               lxx += similarity * similarity;
            }
         totalpair = N*(N-1)/2;
         beta = (lxy - meansimilarity*meandiffsq/totalpair) / (lxx - meansimilarity*meansimilarity/totalpair);
         meandiffsq /= totalpair;
         meansimilarity /= totalpair;
         if(beta > 0) {// trait not genetic
            beta = 0;
            mu = meandiffsq;
            lambda0[t] = 0;
            tau0[t] = 2/mu;
         }else{
            mu = meandiffsq - beta*meansimilarity;
            lambda0[t] = -beta/(mu+beta);
            tau0[t] = 2/(mu+beta);
         }
      }else{       // iterations for best H2 and Beta estimates
         double delta_local = 100;
         mV.Dimension(missingpheno[t].Length());
         for(iter = 0; (iter < 20) && ((delta_local < -0.0001) || (delta_local > 0.0001)); iter++){
            meandiffsq = meansimilarity = lxx = lxy = 0;
            for(int i = 0; i < N; i++)
               for(int j = i+1; j < N; j++){
                  diff = resid[i] - resid[j];
                  diffsq = diff * diff;
                  similarity = D[validpheno[t][i]][validpheno[t][j]];
                  meandiffsq += diffsq;
                  meansimilarity += similarity;
                  lxy += similarity * diffsq;
                  lxx += similarity * similarity;
               }
            totalpair = N*(N-1)/2;
            beta = (lxy - meansimilarity*meandiffsq/totalpair) / (lxx - meansimilarity*meansimilarity/totalpair);
            meandiffsq /= totalpair;
            meansimilarity /= totalpair;
            if(beta > 0) {// trait not genetic
               beta = 0;
               mu = meandiffsq;
               delta_local = 0;
               lambda0[t] = 0;
               tau0[t] = 2/mu;
            }else{
               mu = meandiffsq - beta*meansimilarity;
               delta_local = -beta/(mu+beta)-lambda0[t];
               lambda0[t] = -beta/(mu+beta);
               tau0[t] = 2/(mu+beta);
            }
            for(int i = 0; i < N_ID; i++)
               uy[i] = UYs[i][t];
            mV.Zero();
            for(int i = 0; i < missingpheno[t].Length(); i++)
               for(int j = 0; j < N_Covariate+1; j++)
                  mV[i] += beta0[t][j] * X[missingpheno[t][i]][j];
            for(int i = 0; i < N_ID; i++)
               for(int k = 0; k < missingpheno[t].Length(); k++)
                  uy[i] += UT[i][missingpheno[t][k]] * mV[k];
            Info0[t].Zero();
            if(N_Covariate==0){
               for(int k = 0; k < N_ID; k++)
                  Info0[t].data[0]->data[0] += UX.data[k]->data[0] * UX.data[k]->data[0] / (lambda0[t] * EV[k] + 1);
               beta0[t][0] = 0;
               for(int k = 0; k < N_ID; k++)
                  beta0[t][0] += UX.data[k]->data[0] * uy.data[k] / (lambda0[t] * EV[k] + 1);
               beta0[t][0] /= Info0[t][0][0];
            }else{
               for(int i = 0; i < N_Covariate+1; i++)
                  for(int j = 0; j < N_Covariate+1; j++)
                     for(int k = 0; k < N_ID; k++)
                        Info0[t].data[i]->data[j] += UX.data[k]->data[i] * UX.data[k]->data[j] / (lambda0[t] * EV[k] + 1);
               if(chol.TryDecompose(Info0[t])==0) break;
               chol.Decompose(Info0[t]);
               tV.Zero();
               for(int i = 0; i < N_Covariate+1; i++)
                  for(int k = 0; k < N_ID; k++)
                     tV.data[i] += UX.data[k]->data[i] * uy.data[k] / (lambda0[t] * EV[k] + 1);
               chol.BackSubst(tV);
               beta0[t] = chol.x;
            }
            for(int i = 0; i < N; i++){
               resid[i] = Y[validpheno[t][i]][t];
               for(int j = 0; j < N_Covariate+1; j++)
                  resid[i] -= beta0.data[t]->data[j] * X.data[validpheno[t][i]]->data[j];
            }
         }// iteration ends
      }// if ends
      for(int i = 0; i < N; i++)
         Resid0[validpheno[t][i]][t] = resid[i];
      for(int i = 0; i < missingpheno[t].Length(); i++)
         Resid0[missingpheno[t][i]][t] = 0;

      if(ped.traitCount < 100){
         printf("%-15s %7d %7d %7.4lf %7.4lf %7.3lf",
         (const char*)ped.traitNames[traits[t]], N, iter,
         lambda0[t]/(1+lambda0[t]),lambda0[t], tau0[t]);
         for(int i = 0; i < N_Covariate+1; i++)
            printf(" %9.3lf", beta0[t][i]);
         printf("\n");
      }
   }


   if(Info0) delete []Info0;
   if(missingpheno) delete []missingpheno;
   if(validpheno) delete []validpheno;
   for(int i = 0; i < N_ID; i++)
      delete D[i];
   delete []D;
   UYs.Zero(); // reusing matrix UYs for U'R
   for(int i = 0; i < N_ID; i++)
      for(int t = 0; t < N_Trait; t++)
         for(int k = 0; k < N_ID; k++)
            UYs.data[i]->data[t] += UT[i][k] * Resid0.data[k]->data[t];
   VR = new double * [N_Trait];
   for(int t = 0; t < N_Trait; t++){
      VR[t] = new double [N_ID];
      for(int i = 0; i < N_ID; i++)
         VR[t][i] = 0;
   }
   varVR = new double [N_Trait];
   for(int t = 0; t < N_Trait; t++){
      varVR[t] = 0;
      for(int i = 0; i < N_ID; i++)
         for(int j = 0; j < N_ID; j++)
            VR[t][i] += UT[j][i] * UYs.data[j]->data[t] / (lambda0[t]*EV[j]+1) * tau0[t];
      for(int i = 0; i < N_ID; i++)
         varVR[t] += UYs[i][t] * UYs[i][t] * EV[i] * tau0[t] * tau0[t] / (lambda0[t]*EV[i]+1) / (lambda0[t]*EV[i]+1);
      for(int i = 0; i < N_ID; i++)
         VR[t][i] /= sqrt(varVR[t]);
   }
   for(int i = 0; i < N_ID; i++)
      delete UT[i];
   delete []UT;
}

void Engine::LMMSCORE()
{
   printf("\nOptions in effect:\n");
   printf("\t--mtscore\n");
   if(normalization)
      printf("\t--InvNorm\n");
   if(CoreCount)
      printf("\t--cpus %d\n", CoreCount);
   if(Bit64Flag)
      printf("\t--sysbit 64\n");
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");

   if(FixedEff)
      PreScan_FixedEff();
   else
      PreScan_MTSCORE(); // everything before the scan

   printf("Score-MT genome scan starts at %s", currentTime());
   if(Bit64==64)
      GenomeScan64Bit();
   else
      GenomeScan();   // genome scan
   printf("Score-MT genome scan ends at %s", currentTime());
}


void Engine::ComputeMTSVD()
{
   int N_ID = ID.Length();
   D = new double * [N_ID];
   for(int i = 0; i < N_ID; i++)
      D[i] = new double [N_ID];
   for(int i = 0; i < N_ID; i++)
      for(int j = 0; j < N_ID; j++)
         D[i][j] = 0;

   if(svdinfile!=""){
   StringArray tokens;
   String line;
   IFILE input = ifopen(svdinfile, "rt");
   if(input == NULL)
      error("Cannot open %s to read", (const char*)svdinfile);
   tokens.Clear();
   line.ReadLine(input);
   tokens.AddTokens(line);
   if(tokens.Length() < N_ID)
      error("%d samples in SVD file %s; however there are %d samples in the data.\n",
            tokens.Length(), (const char*)svdinfile, N_ID);
   String tempID;
   int N_tokens = tokens.Length();
   IntArray key(N_tokens);
   key.Set(-1);
   for(int i = 0; i < N_ID; i++){
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
            row-N_ID, N_tokens);
   printf("Eigenvalues, eigenvectors, and genetic similarity are successfully loaded from file %s\n",
         (const char*)svdinfile);
   ifclose(input);
   }else{
   int IBS0Count, het1Count, notMissingCount;
   int id1, id2, count;
   printf("Calculating similarity matrix starts at %s", currentTime());
   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++){
      oneoneCount[i] = 0;
      for(int j = 0; j < 16; j++)
         if(i & shortbase[j]) oneoneCount[i]++;
   }

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
   printf("                                Ends at %s", currentTime());

   printf("SVD starts at %s...", currentTime());
#ifdef WITH_LAPACK
   printf("  LAPACK is used.\n");
   char JOBZ = 'A';
   int info;
   double *A = new double[N_ID*N_ID];
   int dimLA = N_ID;
   for(int i = 0; i < N_ID; i++)
     for(int j = 0; j < N_ID; j++)
       A[i*N_ID+j] = D[j][i];
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
         MD[i][j] = D[i][j];
   svd.Decompose(MD);
   printf("done\n");
   if(svd.n == 0) return;
   for(int i = 0; i < N_ID; i++) EV[i] = svd.w[i];
   // svd.v[k][idx[N_ID-1-j]] stores jth eigenvector
   for(int j = 0; j < N_ID; j++)
      for(int k = 0; k < N_ID; k++)
         UT[j][k] = svd.v[k][j];
#endif
   printf("SVD ends at %s", currentTime());
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
      printf("Eigenvalues, eigenvectors and genetic similarity are saved in file %s\n",
         (const char*)svdoutfile);
      exit(0);
   }
   }   // SVD ends
}


