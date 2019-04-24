//////////////////////////////////////////////////////////////////////
// assoc.cpp
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
#include "MathStats.h"
#include "Kinship.h"
#include "QuickIndex.h"
#include <math.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

void Engine::CountTDT64Bit(IntArray &trioList, IntArray &tcounts, IntArray & ntcounts)
{
   int N_Trio = trioList.Length() / 3;
   tcounts.Dimension(longCount<<6);
   ntcounts.Dimension(longCount<<6);
   tcounts.Zero();
   ntcounts.Zero();
   if(N_Trio==0) return;
   const int BLOCKSIZE=127;
   const int CACHESIZE=16; // 2^(6+5-2+7+2) = 2^17
   const int CACHESIZE_BITPAR = CACHESIZE<<3;
   int trio[3];
   unsigned long long int nomiss, informative, informativePO, word, word2, T[CACHESIZE_BITPAR], NT[CACHESIZE_BITPAR];
   unsigned char *pchar;

   int thread = 0;
   int **localcounts[2];
   for(int k = 0; k < 2; k++){
      localcounts[k] = new int *[defaultMaxCoreCount];
      for(int i = 0; i < defaultMaxCoreCount; i++)
         localcounts[k][i] = NULL;
   }
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(nomiss, informative, informativePO, word, word2, trio, pchar, T, NT, thread)
{
   thread = omp_get_thread_num();
   for(int k = 0; k < 2; k++){
      localcounts[k][thread] = new int [longCount<<6];   // allocate for effective thread
      for(int m = 0; m < markerCount; m++)
         localcounts[k][thread][m] = 0;
   }
   #pragma omp for
#else
   for(int k = 0; k < 2; k++){
      localcounts[k][thread] = new int [longCount<<6];   // allocate for effective thread
      for(int m = 0; m < markerCount; m++)
         localcounts[k][thread][m] = 0;
   }
#endif
   for(int bi = 0; bi < N_Trio; bi+=BLOCKSIZE){ // blocked trio list
      int iMax = (bi > N_Trio-BLOCKSIZE) ? N_Trio: bi+BLOCKSIZE;
      for(int blockb = 0; blockb < longCount; blockb += CACHESIZE){  // blocked SNP words
         int bMax = (blockb > longCount-CACHESIZE) ? longCount: blockb+CACHESIZE;
         for(int j = 0; j < CACHESIZE_BITPAR; j++)
            T[j] = NT[j] = 0;
         for(int i = bi; i < iMax; i++){
            trio[0] = trioList[i*3+1]; // parent 1
            trio[1] = trioList[i*3+2]; // parent 2
            trio[2] = trioList[i*3]; // trio[2]=child
            for(int b = blockb; b < bMax; b++){
               nomiss = (LG[0][trio[0]][b] | LG[1][trio[0]][b]) &
                  (LG[0][trio[1]][b] | LG[1][trio[1]][b]) &
                  (LG[0][trio[2]][b] | LG[1][trio[2]][b]);
               informative = nomiss & ~(LG[0][trio[0]][b] & LG[0][trio[1]][b]);
               int bbb = ((b-blockb)<<3);
               word = informative & // Aa x Aa -> Aa
                  ~(LG[0][trio[0]][b] | LG[0][trio[1]][b] | LG[0][trio[2]][b]);
               word2 = (word&0x0101010101010101);
               T[bbb] += word2;
               NT[bbb] += word2;
               for(int k = 1; k < 8; k++){
                  word2 = ((word>>k)&0x0101010101010101);
                  T[bbb|k] += word2;
                  NT[bbb|k] += word2;
               }
               for(int p = 0; p < 2; p++){   // loop over two parents
                  informativePO = informative & (LG[0][trio[p]][b] ^ LG[0][trio[2]][b]);
                  word = informativePO & ((LG[0][trio[2]][b] & LG[1][trio[2]][b]) | (~LG[1][trio[p]][b]));
                  T[bbb] += (word&0x0101010101010101);
                  for(int k = 1; k < 8; k++)
                     T[bbb|k] += ((word>>k)&0x0101010101010101);
                  word = informativePO & ((LG[0][trio[p]][b] & LG[1][trio[p]][b]) | (~LG[1][trio[2]][b]));
                  NT[bbb] += (word&0x0101010101010101);
                  for(int k = 1; k < 8; k++)
                     NT[bbb|k] += ((word>>k)&0x0101010101010101);
               }  // end of parent
            }  // end of SNP words
         }  // end of trio
         for(int b = blockb; b < bMax; b++){
            int bbb = ((b-blockb)<<3);
            for(int k = 0; k < 8; k++){   // shift
               int pos = ((b<<6) | k);
               pchar = (unsigned char *)&T[bbb|k];
               for(int j = 0; j < 8; j++) // bit parallel
                  localcounts[0][thread][pos | (j<<3)] += pchar[j];
               pchar = (unsigned char*)&NT[bbb|k];
               for(int j = 0; j < 8; j++) // bit parallel
                  localcounts[1][thread][pos | (j<<3)] += pchar[j];
            }  // end of 8 shifts
         }  // end of b
      }  // end of blockb
   }  // end of blocked triolist
#ifdef _OPENMP
}
#endif
   for(int i = 0; i < defaultMaxCoreCount; i++)
      if(localcounts[0][i]){  // thread effectively used 
         for(int m = 0; m < markerCount; m++){
            tcounts[m] += localcounts[0][i][m];
            ntcounts[m] += localcounts[1][i][m];
         }
         delete []localcounts[0][i];
         delete []localcounts[1][i];
      }
   delete []localcounts[0];
   delete []localcounts[1];
}


//   int TCount[CACHESIZE_SNP], NTCount[CACHESIZE_SNP];
//         int mMax = (blockb > longCount-CACHESIZE)? markerCount-(blockb<<6): ((bMax-blockb)<<6);
/*
         for(int b = blockb; b < bMax; b++){
            int bbb = ((b-blockb)<<3);
            for(int k = 0; k < 8; k++){   // shift
               int pos = (bbb<<3) | k;
               pchar = (unsigned char *)&T[bbb|k];
               for(int j = 0; j < 8; j++) // bit parallel
                  TCount[pos | (j<<3)] += pchar[j];
               pchar = (unsigned char*)&NT[bbb|k];
               for(int j = 0; j < 8; j++) // bit parallel
                  NTCount[pos | (j<<3)] += pchar[j];
            }  // end of 8 shifts
         }  // end of b
*/
/*
      for(int m = 0; m < mMax; m++){
         tcounts[(blockb<<6)+m] = TCount[m];
         ntcounts[(blockb<<6)+m] = NTCount[m];
      }
*/

/*
      for(int i = 0; i < N_Trio; i++){
         trio[0] = trioList[i*3+1]; // parent 1
         trio[1] = trioList[i*3+2]; // parent 2
         trio[2] = trioList[i*3]; // trio[2]=child
         for(int b = blockb; b < bMax; b++){
            int bb = (b-blockb)<<6;
            nomiss = (LG[0][trio[0]][b] | LG[1][trio[0]][b]) &
               (LG[0][trio[1]][b] | LG[1][trio[1]][b]) &
               (LG[0][trio[2]][b] | LG[1][trio[2]][b]);
            informative = nomiss & ~(LG[0][trio[0]][b] & LG[0][trio[1]][b]);
            word = informative & // Aa x Aa -> Aa
               ~(LG[0][trio[0]][b] | LG[0][trio[1]][b] | LG[0][trio[2]][b]);
            pchar = (unsigned char*)&word;
            for(int k = 0; k < 8; k++)
               for(byte = pchar[k]; byte; byte &= (byte-1)){
                  pos = bb+(k<<3)+rightmost[byte];
                  TCount[pos] ++;
                  NTCount[pos] ++;
               }
            for(int p = 0; p < 2; p++){   // loop over two parents
               informativePO = informative & (LG[0][trio[p]][b] ^ LG[0][trio[2]][b]);
               word = informativePO & ((LG[0][trio[2]][b] & LG[1][trio[2]][b]) | (~LG[1][trio[p]][b]));
               pchar = (unsigned char*)&word;   // O=AA or P=aa
               for(int k = 0; k < 8; k++)
                  for(byte = pchar[k]; byte; byte &= (byte-1))
                     TCount[bb+(k<<3)+rightmost[byte]] ++;
               word = informativePO & ((LG[0][trio[p]][b] & LG[1][trio[p]][b]) | (~LG[1][trio[2]][b]));
               pchar = (unsigned char*)&word;   // P=AA or O=aa
               for(int k = 0; k < 8; k++)
                  for(byte = pchar[k]; byte; byte &= (byte-1))
                     NTCount[bb+(k<<3)+rightmost[byte]] ++;
            }  // end of parent
         }  // end of pair
      }  // end of trio
      int mMax = (blockb > longCount-CACHESIZE)? markerCount-(blockb<<6): ((bMax-blockb)<<6);
      for(int m = 0; m < mMax; m++){
         ncounts[(blockb<<6)+m] += TCount[m];
         ntcounts[(blockb<<6)+m] += NTCount[m];
      }
*/

void Engine::CountTDT(IntArray &trioList, IntArray &ncounts, IntArray & ntcounts)
{
   int N_Trio = trioList.Length() / 3;
   ncounts.Dimension(shortCount*16);
   ncounts.Zero();
   ntcounts.Dimension(shortCount*16);
   ntcounts.Zero();
   if(N_Trio==0) return;

   char revbase[65536];
   for(int i = 0; i < 16; i++)
      revbase[shortbase[i]] = i;
   char rightmost[65536];
   for(int i = 0; i < 65536; i++)
      rightmost[i] = revbase[i&(-i)];

   const int CACHESIZE=256;
   const int CACHESIZE_SNP = CACHESIZE*16;
   unsigned int TCount[CACHESIZE_SNP], NTCount[CACHESIZE_SNP];
   int trio[3], pos;
   unsigned short int nomiss, informative, word, informativePO;

#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
      private(nomiss, word, informative, informativePO, TCount, NTCount, trio, pos)
#endif
   for(int blockb = 0; blockb < shortCount; blockb += CACHESIZE){
      int bMax = (blockb > shortCount-CACHESIZE) ? shortCount: blockb+CACHESIZE;
      for(int j = 0; j < CACHESIZE_SNP; j++)
         TCount[j] = NTCount[j] = 0;
      for(int i = 0; i < N_Trio; i++){
         trio[0] = trioList[i*3+1]; // parent 1
         trio[1] = trioList[i*3+2]; // parent 2
         trio[2] = trioList[i*3]; // trio[2]=child
         for(int b = blockb; b < bMax; b++){
            int bb = (b-blockb)*16;
            nomiss = (GG[0][trio[0]][b] | GG[1][trio[0]][b]) &
               (GG[0][trio[1]][b] | GG[1][trio[1]][b]) &
               (GG[0][trio[2]][b] | GG[1][trio[2]][b]);
            informative = nomiss & ~(GG[0][trio[0]][b] & GG[0][trio[1]][b]);
            for(word = informative & // Aa x Aa -> Aa
               ~(GG[0][trio[0]][b] | GG[0][trio[1]][b] | GG[0][trio[2]][b]);
               word; word &= (word-1)){
               pos = bb+rightmost[word];
               TCount[pos] ++;
               NTCount[pos] ++;
            }
            for(int p = 0; p < 2; p++){   // loop over two parents
               informativePO = informative & (GG[0][trio[p]][b] ^ GG[0][trio[2]][b]);
               for(word = informativePO & ((GG[0][trio[2]][b] & GG[1][trio[2]][b]) | (~GG[1][trio[p]][b]));
                  word; word &= (word-1)){   // O=AA or P=aa
                  pos = bb+rightmost[word];
                  TCount[pos] ++;
               }
               for(word = informativePO & ((GG[0][trio[p]][b] & GG[1][trio[p]][b]) | (~GG[1][trio[2]][b]));
                  word; word &= (word-1)){   // P=AA or O=aa
                  pos = bb+rightmost[word];
                  NTCount[pos] ++;
               }
            }  // end of parent
         }  // end of pair
      }  // end of trio
      int mMax = (blockb > shortCount-CACHESIZE)? markerCount-blockb*16: (bMax-blockb)*16;
      for(int m = 0; m < mMax; m++){
         ncounts[blockb*16+m] += TCount[m];
         ntcounts[blockb*16+m] += NTCount[m];
      }
   }  // end of blockb
}

void Engine::CountTDTinX(IntArray &trioList, IntArray &ncounts, IntArray & ntcounts)
{
   int N_Trio = trioList.Length() / 3;
   ncounts.Dimension(xshortCount*16);
   ncounts.Zero();
   ntcounts.Dimension(xshortCount*16);
   ntcounts.Zero();
   if(N_Trio==0) return;

   char revbase[65536];
   for(int i = 0; i < 16; i++)
      revbase[shortbase[i]] = i;
   char rightmost[65536];
   for(int i = 0; i < 65536; i++)
      rightmost[i] = revbase[i&(-i)];

   const int CACHESIZE=256;
   const int CACHESIZE_SNP = CACHESIZE*16;
   unsigned int TCount[CACHESIZE_SNP], NTCount[CACHESIZE_SNP];//, informativeCount[4096];
   int trio[3], pos;
   unsigned short int nomiss, informative, word, informativePO;

#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
      private(nomiss, word, informative, informativePO, TCount, NTCount, trio, pos)
#endif
   for(int blockb = 0; blockb < xshortCount; blockb += CACHESIZE){
      int bMax = (blockb > xshortCount-CACHESIZE) ? xshortCount: blockb+CACHESIZE;
      for(int j = 0; j < CACHESIZE_SNP; j++)
         TCount[j] = NTCount[j] = 0;
      for(int i = 0; i < N_Trio; i++){
         trio[0] = trioList[i*3+1]; // parent 1
         trio[1] = trioList[i*3+2]; // parent 2
         trio[2] = trioList[i*3]; // trio[2]=child
         for(int b = blockb; b < bMax; b++){
            int bb = (b-blockb)*16;
            nomiss = (XG[0][trio[0]][b] | XG[1][trio[0]][b]) &
               (XG[0][trio[1]][b] | XG[1][trio[1]][b]) &
               (XG[0][trio[2]][b] | XG[1][trio[2]][b]);
            informative = nomiss & ~(XG[0][trio[0]][b] & XG[0][trio[1]][b]);
            for(word = informative & // Aa x Aa -> Aa
               ~(XG[0][trio[0]][b] | XG[0][trio[1]][b] | XG[0][trio[2]][b]);
               word; word &= (word-1)){
               pos = bb+rightmost[word];
               TCount[pos] ++;
               NTCount[pos] ++;
            }
            for(int p = 0; p < 2; p++){   // loop over two parents
               informativePO = informative & (XG[0][trio[p]][b] ^ XG[0][trio[2]][b]);
               for(word = informativePO & ((XG[0][trio[2]][b] & XG[1][trio[2]][b]) | (~XG[1][trio[p]][b]));
                  word; word &= (word-1)){   // O=AA or P=aa
                  pos = bb+rightmost[word];
                  TCount[pos] ++;
               }
               for(word = informativePO & ((XG[0][trio[p]][b] & XG[1][trio[p]][b]) | (~XG[1][trio[2]][b]));
                  word; word &= (word-1)){   // P=AA or O=aa
                  pos = bb+rightmost[word];
                  NTCount[pos] ++;
               }
            }  // end of parent
         }  // end of pair
      }  // end of trio
      int mMax = (blockb > xshortCount-CACHESIZE)? xmarkerCount-blockb*16: (bMax-blockb)*16;
      for(int m = 0; m < mMax; m++){
         ncounts[blockb*16+m] += TCount[m];
         ntcounts[blockb*16+m] += NTCount[m];
      }
   }  // end of blockb
}

void Engine::TDT()
{
   printf("\nOptions in effect:\n");
   printf("\t--tdt\n");
   if(CoreCount)
      printf("\t--cpus %d\n", CoreCount);
   if(Bit64Flag)
      printf("\t--sysbit 64\n");
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");

   printf("Genome-wide TDT scan starts at %s", currentTime());
   MakeTrioForTDT();
   int N_Trio = Ltrio.Length() / 3;
   if(N_Trio==0){
      warning("TDT analysis requires parent-affected-offspring trios.\n");
      return;
   }

   // uninformative families
   IntArray infoID(ped.count);
   infoID.Zero();
   IntArray ifam(ped.familyCount);
   ifam.Zero();
   for(int i = 0; i < Ltrio.Length(); i++)
         infoID[i] = 1;
   int uninformativeFamilyCount = 0;
   for(int f = 0; f < ped.familyCount; f++){
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(infoID[i]) ifam[f]=1;
      if(!ifam[f]) uninformativeFamilyCount++;
   }
   if(uninformativeFamilyCount){
      String uninfofile = prefix;
      uninfofile.Add("TDTuninfo.txt");
      printf("%d parent-affected-child trios in %d pedigrees will be used in TDT analysis.\n",
         N_Trio, ped.familyCount - uninformativeFamilyCount);
      FILE *ufp = fopen((const char*)uninfofile, "wt");
      for(int f = 0; f < ped.familyCount; f++)
         if(!ifam[f])
            for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
               fprintf(ufp, "%s %s\n",
               (const char*)ped[i].famid, (const char*)ped[i].pid);
      fclose(ufp);
      printf("IDs of %d TDT uninformative pedigrees saved in file %s.\n",
         uninformativeFamilyCount, (const char*)uninfofile);
   }

   printf("Scanning autosomes with %d CPU cores...\n", defaultMaxCoreCount);
   IntArray trioList(0);
   for(int i = 0; i < N_Trio; i++){
      trioList.Push(geno[Ltrio[i*3]]);
      trioList.Push(geno[Ltrio[i*3+1]]);
      trioList.Push(geno[Ltrio[i*3+2]]);
   }
   IntArray tcounts, ntcounts;
   if(Bit64 == 64)
      CountTDT64Bit(trioList, tcounts, ntcounts);
   else
      CountTDT(trioList, tcounts, ntcounts);
   String pedfile = prefix;
   pedfile.Add("tdt.txt");
   FILE *fp = fopen(pedfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
   fprintf(fp, "SNP");
   if(chromosomes.Length()) fprintf(fp, "\tChr");
   if(bp.Length()) fprintf(fp, "\tPos");
   fprintf(fp, "\tA1\tA2\tT\tNT\tOR\tChisq\tP\n");
   double statistic, p;
   for(int m = 0; m < markerCount; m++){
      p = tcounts[m] + ntcounts[m];
      if(p>0.5){
         statistic = tcounts[m] - ntcounts[m];
         statistic = statistic*statistic/p;
         p = chidist(statistic, 1);
      }else{
         statistic = 0.0;
         p = 1.0;
      }
      fprintf(fp, "%s\t%d\t%d\t%s\t%s\t%d\t%d\t%.3lf\t%.3lf\t%.3G\n",
         (const char*)snpName[m], chromosomes[m], bp[m],
         (const char*)alleleLabel[0][m], (const char*)alleleLabel[1][m],
         tcounts[m], ntcounts[m],
         ntcounts[m]? tcounts[m]*1.0 / ntcounts[m]: -9,
         statistic, p);
   }  // end of m
   if(xmarkerCount)
      printf("Scanning X-chromosome with %d CPU cores...\n", defaultMaxCoreCount);
   CountTDTinX(trioList, tcounts, ntcounts);
   for(int m = 0; m < xmarkerCount; m++){
      p = tcounts[m] + ntcounts[m];
      if(p>0.5){
         statistic = tcounts[m] - ntcounts[m];
         statistic = statistic*statistic/p;
         p = chidist(statistic, 1);
      }else{
         statistic = 0.0;
         p = 1.0;
      }
      if(xsnpName.Length())
         fprintf(fp, "%s", (const char*)xsnpName[m]);
      else
         fprintf(fp, "XSNP%d", m+1);
      fprintf(fp, "\tX");
      if(xbp.Length())
         fprintf(fp, "\t%d", bp[m]);
      fprintf(fp, "\t%s\t%s\t%d\t%d\t%.3lf\t%.3lf\t%.3G\n",
         (const char*)xalleleLabel[0][m], (const char*)xalleleLabel[1][m],
         tcounts[m], ntcounts[m],
         ntcounts[m]? tcounts[m]*1.0 / ntcounts[m]: -9,
         statistic, p);
   }  // end of m
   fclose(fp);
   printf("TDT scan results for %s saved in file %s\n",
      (const char*)ped.affectionNames[0], (const char*)pedfile);
   printf("Genome-wide TDT scan ends at %s\n", currentTime());
}

void Engine::ExtractUnrelated()
{
   IntArray unrelatedCount, ID_unrelated;
   QuickIndex idx;

   SemifamilyKinship();
   unrelatedList.Dimension(0);
   for(int f = 0; f < ped.familyCount; f++){
      if(id[f].Length() < 2) {
         if(id[f].Length() == 1)
            unrelatedList.Push(id[f][0]);
         continue;
      }
      unrelatedCount.Dimension(id[f].Length());
      unrelatedCount.Zero();
      for(int i = 0; i < id[f].Length(); i++)
         for(int j = i+1; j < id[f].Length(); j++)
            if( (semifamilyFlag && pedKin[f][i][j] < 0.022) || (pedKin[f][i][j] < 0.001) ){
               unrelatedCount[i] ++;
               unrelatedCount[j] ++;
            }
      idx.Index(unrelatedCount);
      ID_unrelated.Dimension(0);
      ID_unrelated.Push(idx[unrelatedCount.Length()-1]);
      for(int i = unrelatedCount.Length()-2; (i >= 0) && (unrelatedCount[idx[i]] > 0); i--){
         bool isUnrelated = true;
         for(int j = 0; j < ID_unrelated.Length(); j++)
            if( (semifamilyFlag && pedKin[f][idx[i]][ID_unrelated[j]] >= 0.022)
               || (!semifamilyFlag && pedKin[f][idx[i]][ID_unrelated[j]] >= 0.001) ) {
               isUnrelated = false;
               break;
            }
         if(isUnrelated) ID_unrelated.Push(idx[i]);
      }
      for(int i = 0; i < ID_unrelated.Length(); i++)
         unrelatedList.Push(id[f][ID_unrelated[i]]);
   }
   if(pedKin) {
      delete []pedKin;
      pedKin = NULL;
   }
}

void Engine::HetSharing()
{
   if(geno.Length()==0) BuildShortBinary();

   String pedfile = prefix;
   pedfile.Add("het.txt");
   FILE *fp = fopen(pedfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
   int pos, k;

   int hom1, hom2, rareAllele, homCount[16];
   int HetHetCount[16], HetCount[16], nonmissingCount[16];
   IntArray hetRatio[16];
   int HetHet, het1, het2, het, nonmissing, b;
   int id1, id2, no;

   fprintf(fp, "SNP\tChr\tPosition\tN\tHetHet\tHet");
   no = 0;
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = 0; i < id[f].Length(); i++)
         for(int j = i+1; j < id[f].Length(); j++){
            fprintf(fp, "\t%s_%d%d", (const char*)ped.families[f]->famid,
               i+1, j+1);
            no++;
         }

   fprintf(fp, "\n");
   for(int k = 0; k < 16; k++)
      hetRatio[k].Dimension(no);

   for(int m = 0; m < markerCount; m+=16){
      b = m/16;
      for(int k = 0; k < 16; k++){
         HetHetCount[k] = HetCount[k] = nonmissingCount[k] = 0;
         hetRatio[k].Zero();
      }
      no = 0;
      for(int f = 0; f < ped.familyCount; f++){
         for(int i = 0; i < id[f].Length(); i++)
            for(int j = i+1; j < id[f].Length(); j++){
               id1 = geno[id[f][i]]; id2 = geno[id[f][j]];
               HetHet = (~GG[0][id1][b]) & GG[1][id1][b] & (~GG[0][id2][b]) & GG[1][id2][b];
               het1 = (GG[0][id2][b] | GG[1][id2][b]) & (~GG[0][id1][b]) & GG[1][id1][b];
               het2 = (GG[0][id1][b] | GG[1][id1][b]) & (~GG[0][id2][b]) & GG[1][id2][b];
               nonmissing = (GG[0][id1][b] | GG[1][id1][b]) & (GG[0][id2][b] | GG[1][id2][b]);

               het = het1 | het2;
               for(int k = 0; k < 16; k++){
                  if(het & shortbase[k]) // at least one is heterozygote
                     HetCount[k]++;
                  else
                     hetRatio[k][no] = -1;
                  if(HetHet & shortbase[k]){ // both are heterozygotes
                     HetHetCount[k]++;
                     hetRatio[k][no] = 1;
                  }
                  if(nonmissing & shortbase[k])
                     nonmissingCount[k] ++;
               }
               no ++;
            }
      }
      for(int k = 0; k < 16; k++){
         if(m+k >= markerCount) continue;
         if(snpName.Length())
            fprintf(fp, "%s", (const char*)snpName[m+k]);
         else
            fprintf(fp, "SNP%d", m+k+1);
         if(chromosomes.Length())
            fprintf(fp, " %d", chromosomes[m+k]);
         if(bp.Length())
            fprintf(fp, " %d", bp[m+k]);

         fprintf(fp, "\t%d\t%d\t%d", nonmissingCount[k], HetHetCount[k], HetCount[k]);
         no = 0;
         for(int f = 0; f < ped.familyCount; f++)
            for(int i = 0; i < id[f].Length(); i++)
               for(int j = i+1; j < id[f].Length(); j++){
                  fprintf(fp, "\t%d", hetRatio[k][no]);
                  no++;
               }
         fprintf(fp, "\n");
      }
   }

   fclose(fp);
   printf("Heterozygosity-sharing statistics by SNPs saved in file %s\n\n", (const char*)pedfile);
}


void Engine::MakeTrioForTDT()
{
   Ltrio.Dimension(0);
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(geno[i]!=-1 && ped[i].affections[0]==2
            && ped[i].father && geno[ped[i].father->serial]!=-1
            && ped[i].mother && geno[ped[i].mother->serial]!=-1){
            Ltrio.Push(i);
            Ltrio.Push(ped[i].father->serial);
            Ltrio.Push(ped[i].mother->serial);
         }
   if(Ltrio.Length())
      printf("There are %d parent-affected-offspring trios according to the pedigree.\n",
         Ltrio.Length()/3);
   if(!semifamilyFlag) { // not semi-pedigree
      if(Ltrio.Length()==0)
         printf("There are no parent-affected-offspring trios in the data.\n");
      return;
   }
   int LtrioCount = Ltrio.Length()/3;
   double kinship;
   Kinship kin;
   int id1, id2;
   int HetHetCount, IBS0Count, het1Count, notMissingCount;

   Lpo.Dimension(0);
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
         if(geno[i]==-1 || ped[i].affections[0] != 2) continue;
         if(ped[i].father && geno[ped[i].father->serial]!=-1){
            Lpo.Push(i); Lpo.Push(ped[i].father->serial);
         }
         if(ped[i].mother && geno[ped[i].mother->serial]!=-1){
            Lpo.Push(i); Lpo.Push(ped[i].mother->serial);
         }
      }
   IntArray Lpo1F;

   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
      printf("1st-degree relatives are treated as parent-offspring if IBS0 < %.4lf\n",
         errorrateCutoff);
      IntArray connections[5000];
      for(int f = 0; f < ped.familyCount; f++){
         if(id[f].Length() < 2) continue;
         kin.Setup(*ped.families[f]);
         Lpo1F.Dimension(0);

         for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
            if(geno[i]==-1 || ped[i].affections[0] != 2) continue;
            if(ped[i].father && geno[ped[i].father->serial]!=-1
               && ped[i].mother && geno[ped[i].mother->serial]!=-1);
            else if(ped[i].father && geno[ped[i].father->serial]!=-1){
               Lpo1F.Push(i); Lpo1F.Push(ped[i].father->serial);
            }else if(ped[i].mother && geno[ped[i].mother->serial]!=-1){
               Lpo1F.Push(i); Lpo1F.Push(ped[i].mother->serial);
            }
         }

         for(int i = 0; i < id[f].Length(); i++){
            id1 = geno[id[f][i]];
            for(int j = i+1; j < id[f].Length(); j++){
               id2 = geno[id[f][j]];
               if(kin(ped[id[f][i]], ped[id[f][j]]) > 0.005) continue;
               HetHetCount = IBS0Count = het1Count = 0;
               for(int m = 0; m < shortCount; m++){
                  HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
                  IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
               }
               if(HetHetCount < IBS0Count*2) continue;   // unrelated
               for(int m = 0; m < shortCount; m++)
                  het1Count += (oneoneCount[(GG[0][id2][m] | GG[1][id2][m]) & (~GG[0][id1][m]) & GG[1][id1][m]]
                     + oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]]);
               kinship = het1Count? (HetHetCount - IBS0Count*2)*1.0/het1Count: 0;
               if(kinship > 0.15 && kinship < 0.38){ // PO
                  notMissingCount = 0;
                  for(int m = 0; m < shortCount; m++)
                     notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
                  if(IBS0Count < notMissingCount * errorrateCutoff){
                     Lpo.Push(id[f][i]); Lpo.Push(id[f][j]);
                     Lpo1F.Push(id[f][i]); Lpo1F.Push(id[f][j]);
                  }
               }
            }
         }
         for(int i = 0; i < ped.families[f]->count; i++) connections[i].Dimension(0);
         for(int i = 0; i < Lpo1F.Length()/2; i++){
            connections[Lpo1F[i*2]-ped.families[f]->first].Push(Lpo1F[i*2+1]);
            connections[Lpo1F[i*2+1]-ped.families[f]->first].Push(Lpo1F[i*2]);
         }
         for(int i = 0; i < ped.families[f]->count; i++)
            if(ped[ped.families[f]->first+i].affections[0]==2 && connections[i].Length() >= 2){
               for(int j = 0; j < connections[i].Length(); j++){
                  id1 = geno[connections[i][j]];
                  for(int k = j+1; k < connections[i].Length(); k++){
                     id2 = geno[connections[i][k]];
                     HetHetCount = IBS0Count = het1Count = 0;
                     for(int m = 0; m < shortCount; m++){
                        HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
                        het1Count += (oneoneCount[(GG[0][id2][m] | GG[1][id2][m]) & (~GG[0][id1][m]) & GG[1][id1][m]]
                           + oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]]);
                        IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
                     }
                     if(het1Count){
                        kinship = (HetHetCount - IBS0Count*2)*1.0/het1Count;
                        if(kinship < 0.0625){   // not 2nd-degree, likely unrelated
                           Ltrio.Push(ped.families[f]->first+i);
                           Ltrio.Push(connections[i][j]);
                           Ltrio.Push(connections[i][k]);
                           j = k = connections[i].Length();
                           break;
                        }
                     }
                  }
               }
            }
      }
   if(Ltrio.Length() > LtrioCount*3)
      printf("There are additional %d parent-affected-offspring trios by inference.\n",
         Ltrio.Length()/3 - LtrioCount);
   if(Ltrio.Length()==0)
      printf("There are no parent-affected-offspring trios in the data.\n");
}


