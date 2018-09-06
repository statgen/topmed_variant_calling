//////////////////////////////////////////////////////////////////////
// assoc.cpp
// (c) 2010-2017 Wei-Min Chen
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
// May 25, 2017

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
   if(positions.Length()) fprintf(fp, "\tPos");
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
/*
      if(snpName.Length())
         fprintf(fp, "%s", (const char*)snpName[m]);
      else
         fprintf(fp, "SNP%d", m+1);
      if(chromosomes.Length())
         fprintf(fp, "\t%d", chromosomes[m]);
      if(positions.Length())
         fprintf(fp, "\t%d", int(positions[m]*1000000+0.5));
      fprintf(fp, "\t%s\t%s\t%d\t%d\t%.3lf\t%.3lf\t%.3G\n",
         (const char*)alleleLabel[0][m], (const char*)alleleLabel[1][m],
         tcounts[m], ntcounts[m],
         ntcounts[m]? tcounts[m]*1.0 / ntcounts[m]: -9,
         statistic, p);
*/
      fprintf(fp, "%s\t%d\t%d\t%s\t%s\t%d\t%d\t%.3lf\t%.3lf\t%.3G\n",
         (const char*)snpName[m], chromosomes[m], int(positions[m]*1000000+0.5),
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
      if(xpositions.Length())
         fprintf(fp, "\t%d", int(xpositions[m]*1000000+0.5));
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
         if(positions.Length())
            fprintf(fp, " %d", int(positions[m+k]*1000000+0.5));

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

/*
void Engine::MakeTrioForTDT_64bit()
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
}
*/

/*
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
      private(missingCount, freq, AA, Aa, missing, word,\
      informative, TCount, NTCount, informativeCount, \
      id1, id2, id3, p, statistic, buffer)
#endif
   for(int blockb = 0; blockb < shortCount; blockb += CACHESIZE){
      int bMax = (blockb > shortCount-CACHESIZE) ? shortCount: blockb+CACHESIZE;
      for(int i = 0; i < 4096; i++){
         missingCount[i] = 0;
         freq[i] = 0.0;
      }
      for(int u = 0; u < unrelatedList.Length(); u++)
         for(int b = blockb; b < bMax; b++){
            int bb = (b-blockb)*16;
            int i = geno[unrelatedList[u]];
            AA = GG[0][i][b] & GG[1][i][b];
            Aa = (~GG[0][i][b]) & GG[1][i][b];
            missing = (~GG[0][i][b]) & (~GG[1][i][b]) & 65535;
            for(word = AA; word; word &= (word-1))
               freq[bb+rightmost[word]] += 1.0;
            for(word = Aa; word; word &= (word-1))
               freq[bb+rightmost[word]] += 0.5;
            for(word = missing; word; word &= (word-1))
               missingCount[bb+rightmost[word]] ++;
         }
      for(int j = 0; j < 4096; j++)
         if(missingCount[j] < unrelatedList.Length())
            freq[j] /= (unrelatedList.Length() - missingCount[j]);
      for(int j = 0; j < 4096; j++)
         TCount[j] = NTCount[j] = informativeCount[j] = 0;
      for(int i = 0; i < N_Trio; i++){
         id1 = geno[Ltrio[i*3+1]];
         id2 = geno[Ltrio[i*3+2]];
         id3 = geno[Ltrio[i*3]];   // id3 is the child
         for(int b = blockb; b < bMax; b++){
            int bb = (b-blockb)*16;
            informative = (GG[0][id3][b] | GG[1][id3][b]) & // offspring not missing
            ( (~GG[0][id1][b] & GG[1][id1][b] & (GG[0][id2][b] | GG[1][id2][b])) |
              (~GG[0][id2][b] & GG[1][id2][b] & (GG[0][id1][b] | GG[1][id1][b])) );
            for(word = informative; word; word &= (word-1))
               informativeCount[bb+rightmost[word]] ++;
            for(word = ~GG[0][id1][b] & GG[1][id1][b] &
               ~GG[0][id2][b] & GG[1][id2][b] &
               ~GG[0][id3][b] & GG[1][id3][b]; // Aa x Aa -> Aa
               word; word &= (word-1)){
               Aa = bb+rightmost[word];
               TCount[Aa] ++;
               NTCount[Aa] ++;
            }
            for(word = informative &   // Aa->AA or aa->Aa
              ( (~GG[0][id1][b] & GG[1][id1][b] & GG[0][id3][b] & GG[1][id3][b]) |
                (GG[0][id1][b] & ~GG[1][id1][b] & ~GG[0][id3][b] & GG[1][id3][b]) );
               word; word &= (word-1))
                  TCount[bb+rightmost[word]] ++;
            for(word = informative &   // AA->Aa or Aa->aa
              ( (~GG[0][id3][b] & GG[1][id3][b] & GG[0][id1][b] & GG[1][id1][b]) |
                (GG[0][id3][b] & ~GG[1][id3][b] & ~GG[0][id1][b] & GG[1][id1][b]) );
               word; word &= (word-1))
                  NTCount[bb+rightmost[word]] ++;

            for(word = informative &   // Aa->AA or aa->Aa
              ( (~GG[0][id2][b] & GG[1][id2][b] & GG[0][id3][b] & GG[1][id3][b]) |
                (GG[0][id2][b] & ~GG[1][id2][b] & ~GG[0][id3][b] & GG[1][id3][b]) );
               word; word &= (word-1))
                  TCount[bb+rightmost[word]] ++;
            for(word = informative &   // AA->Aa or Aa->aa
              ( (~GG[0][id3][b] & GG[1][id3][b] & GG[0][id2][b] & GG[1][id2][b]) |
                (GG[0][id3][b] & ~GG[1][id3][b] & ~GG[0][id2][b] & GG[1][id2][b]) );
               word; word &= (word-1))
                  NTCount[bb+rightmost[word]] ++;
         }  // end of pair
      }  // end of trio
      for(int b = blockb; b < bMax; b++)
         for(int j = 0; j < 16; j++){
            int m = b*16+j;
            if(m >= markerCount) continue;
            int  jj = (b-blockb)*16+j;
            p = TCount[jj]+NTCount[jj];
            if(p > 0.5){
               statistic = TCount[jj] - NTCount[jj];
               statistic = statistic*statistic/p;
            }else
               statistic = 0;
            p = chidist(statistic, 1);
            if(snpName.Length())
               sprintf(buffer, "%s", (const char*)snpName[m]);
            else
               sprintf(buffer, "SNP%d", m+1);
            buffers[b].Add(buffer);
            if(chromosomes.Length()){
               sprintf(buffer, "\t%d", chromosomes[m]);
               buffers[b].Add(buffer);
            }
            if(positions.Length()){
               sprintf(buffer, "\t%d", int(positions[m]*1000000+0.5));
               buffers[b].Add(buffer);
            }
            if(freq[jj] < 0.5){
               sprintf(buffer, "\t%s\t%s",
                  (const char*)alleleLabel[0][m], (const char*)alleleLabel[1][m]);
            }else{
               sprintf(buffer, "\t%s\t%s",
                  (const char*)alleleLabel[1][m], (const char*)alleleLabel[0][m]);
            }
            buffers[b].Add(buffer);
            if(freq[jj] < 0.5){
               sprintf(buffer, "\t%.3lf\t%d\t%d\t%d\t%.3lf\t%.3lf\t%.3G\n",
                  freq[jj], informativeCount[jj], TCount[jj], NTCount[jj],
                  NTCount[jj]? TCount[jj]*1.0/NTCount[jj]: -9, statistic, p);
            }else{
               sprintf(buffer, "\t%.3lf\t%d\t%d\t%d\t%.3lf\t%.3lf\t%.3G\n",
                  1-freq[jj], informativeCount[jj], NTCount[jj], TCount[jj],
                  TCount[jj]? NTCount[jj]*1.0/TCount[jj]: -9, statistic, p);
            }
            buffers[b].Add(buffer);
         }
   }
   String pedfile = prefix;
   pedfile.Add("tdt.txt");
   FILE *fp = fopen(pedfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
   fprintf(fp, "SNP");
   if(chromosomes.Length()) fprintf(fp, "\tChr");
   if(positions.Length()) fprintf(fp, "\tPos");
   fprintf(fp, "\tA1\tA2\tFreq1\tN_Info\tT\tNT\tOR\tChisq\tP\n");
   for(int b = 0; b < buffers.Length(); b++)
      buffers[b].Write(fp);
*/
/*

   for(int m = 0; m < xmarkerCount; m+=16){
      int b = m/16;
      for(int i = 0; i < 16; i++){
         missingCount[i] = 0;
         xfreq[i] = 0.0;
      }
      for(int u = 0; u < unrelatedList.Length(); u++){
         int i = geno[unrelatedList[u]];
         AA = XG[0][i][b] & XG[1][i][b];
         Aa = (~XG[0][i][b]) & XG[1][i][b];
         missing = (~XG[0][i][b]) & (~XG[1][i][b]);
         for(int j = 0; j < 16; j++){
            if(AA & shortbase[j])
               xfreq[j] ++;
            else if(Aa & shortbase[j])
               xfreq[j] += 0.5;
            else if(missing & shortbase[j])
               missingCount[j] ++;
         }
      }
      for(int j = 0; j < 16; j++)
         if(missingCount[j] < unrelatedList.Length())
            xfreq[j] /= (unrelatedList.Length() - missingCount[j]);
      for(int j = 0; j < 16; j++)
         TCount[j] = NTCount[j] = informativeCount[j] = 0;
      for(int i = 0; i < N_Trio; i++){
         id1 = geno[Ltrio[i*3+1]];
         id2 = geno[Ltrio[i*3+2]];
         id3 = geno[Ltrio[i*3]];   // id3 is the child
         informative = (XG[0][id3][b] | XG[1][id3][b]) & // offspring not missing
            ( (~XG[0][id1][b] & XG[1][id1][b] & (XG[0][id2][b] | XG[1][id2][b])) |
              (~XG[0][id2][b] & XG[1][id2][b] & (XG[0][id1][b] | XG[1][id1][b])) );
         if(!informative) continue; // no informative trios at any SNPs
         for(int j = 0; j < 16; j++)
            if(informative & shortbase[j])
               informativeCount[j]++;
         for(int k = 0; k < 2; k++){
            id1 = geno[Ltrio[i*3+1+k]];
            T = informative &   // Aa->AA or aa->Aa
              ( (~XG[0][id1][b] & XG[1][id1][b] & XG[0][id3][b] & XG[1][id3][b]) |
                (XG[0][id1][b] & ~XG[1][id1][b] & ~XG[0][id3][b] & XG[1][id3][b]) );
            NT = informative &   // AA->Aa or Aa->aa
              ( (~XG[0][id3][b] & XG[1][id3][b] & XG[0][id1][b] & XG[1][id1][b]) |
                (XG[0][id3][b] & ~XG[1][id3][b] & ~XG[0][id1][b] & XG[1][id1][b]) );
            if(T)
               for(int j = 0; j < 16; j++)
                  if(T & shortbase[j])
                     TCount[j]++;
            if(NT)
               for(int j = 0; j < 16; j++)
                  if(NT & shortbase[j])
                     NTCount[j]++;
         }
      }
      for(int j = 0; j < 16; j++){
         if(m+j >= xmarkerCount) continue;
         p = TCount[j]+NTCount[j];
         if(p > 0.5){
            statistic = TCount[j] - NTCount[j];
            statistic = statistic*statistic/p;
         }else
            statistic = 0;
         p = chidist(statistic, 1);
         if(xsnpName.Length())
            fprintf(fp, "%s", (const char*)xsnpName[m+j]);
         else
            fprintf(fp, "SNPX%d", m+j+1);
         fprintf(fp, "\tX");
         if(xpositions.Length())
            fprintf(fp, "\t%.6lf", xpositions[m+j]);
         if(xfreq[j] < 0.5){
            fprintf(fp, "\t%s\t%s\t%.3lf",
               (const char*)xalleleLabel[0][m+j], (const char*)xalleleLabel[1][m+j], xfreq[j]);
            fprintf(fp, "\t%d\t%d\t%d\t%.3lf\t%.3lf\t%.3G\n",
               informativeCount[j], TCount[j], NTCount[j],
               TCount[j]? NTCount[j]*1.0/TCount[j]: -9, statistic, p);
         }else{
            fprintf(fp, "\t%s\t%s\t%.3lf",
               (const char*)xalleleLabel[1][m+j], (const char*)xalleleLabel[0][m+j], 1-xfreq[j]);
            fprintf(fp, "\t%d\t%d\t%d\t%.3lf\t%.3lf\t%.3G\n",
               informativeCount[j], NTCount[j], TCount[j],
               TCount[j]? NTCount[j]*1.0/TCount[j]: -9, statistic, p);
         }
      }
   }
   fclose(fp);
   printf("\n");
   */

//   if(geno.Length()==0) BuildShortBinary();

/*
void Engine::TDT_64bit()
{
   unsigned long long int informative, T, NT, HHH;
   int missingCount[64], id1, id2, id3;
   int TCount[64], NTCount[64], informativeCount[64];
   double freq[64], xfreq[64], statistic, p, temp;
   if(geno.Length()==0) error("Genotype data not available");//BuildShortBinary();

   printf("\nOptions in effect:\n");
   printf("\t--tdt\n");
   if(CoreCount)
      printf("\t--cpus %d\n", CoreCount);
   if(Bit64)
      printf("\t--bit64\n");
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");

   printf("Genome-wide TDT scan starts at %s", currentTime());
   MakeTrioForTDT_64bit();
   int N_Trio = Ltrio.Length() / 3;
   if(N_Trio==0){
      warning("TDT analysis requires parent-affected-offspring trios.\n");
      return;
   }

   char buffer[1024];
   StringArray buffers;
   buffers.Dimension(longCount);
   for(int i = 0; i < buffers.Length(); i++)
      buffers[i].Clear();

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
   unsigned long long int longbase[64];
   for(int i = 0; i < Bit64; i++)
      longbase[i] = ((unsigned long long int)1 << i);

   int *trioID[3];
   for(int j = 0; j < 3; j++)
      trioID[j] = new int[N_Trio];
   for(int i = 0; i < N_Trio; i++){
      trioID[0][i] = geno[Ltrio[i*3+1]];
      trioID[1][i] = geno[Ltrio[i*3+2]];
      trioID[2][i] = geno[Ltrio[i*3]]; // child
   }
#ifdef _OPENMP
   printf("Scanning autosomes with %d CPU cores...",
      CoreCount?CoreCount:defaultEfficientCoreCount);
   #pragma omp parallel num_threads(CoreCount?CoreCount:defaultEfficientCoreCount) \
   private(missingCount, freq, TCount, NTCount, informativeCount, p, statistic, buffer)
{
#endif
   unsigned long long int *IsAA, *IsAa, *Isaa, *IsMissing;
   IsAA = new unsigned long long int [idCount];
   IsAa = new unsigned long long int [idCount];
   Isaa = new unsigned long long int [idCount];
   IsMissing = new unsigned long long int [idCount];
   unsigned long long int word;
   unsigned long long int CS[16];
   int tempCount[64];
#ifdef _OPENMP
   #pragma omp for
#endif
   for(int b = 0; b < longCount; b++){
      for(int i = 0; i < idCount; i++){
         IsAA[i] = LG[0][i][b] & LG[1][i][b];
         IsAa[i] = (~LG[0][i][b]) & LG[1][i][b];
         Isaa[i] = LG[0][i][b] & (~LG[1][i][b]);
         IsMissing[i] = (~LG[0][i][b]) & (~LG[1][i][b]);
      }

      // Missing Count
      for(int i = 0; i < 16; i++) CS[i] = 0;
      for(int i = 0; i < idCount; i++){
         word = IsMissing[i];
         for(int bit = 0; bit < 16; bit++)
            CS[bit] += ((word>>bit)&0x0001000100010001);
      }
      for(int bit = 0; bit < 16; bit++)
         for(int repeat = 0; repeat < 4; repeat++)
            missingCount[repeat*16+bit] = (CS[bit]>>(16*repeat))&0xFFFF;

      // Aa Count
      for(int i = 0; i < 16; i++) CS[i] = 0;
      for(int i = 0; i < idCount; i++){
         word = IsAa[i];
         for(int bit = 0; bit < 16; bit++)
            CS[bit] += ((word>>bit)&0x0001000100010001);
      }
      for(int bit = 0; bit < 16; bit++)
         for(int repeat = 0; repeat < 4; repeat++)
            tempCount[repeat*16+bit] = (CS[bit]>>(16*repeat))&0xFFFF;

      // AA Count
      for(int i = 0; i < 16; i++) CS[i] = 0;
      for(int i = 0; i < idCount; i++){
         word = IsAA[i];
         for(int bit = 0; bit < 16; bit++)
            CS[bit] += ((word>>bit)&0x0001000100010001);
      }
      for(int bit = 0; bit < 16; bit++)
         for(int repeat = 0; repeat < 4; repeat++)
            tempCount[repeat*16+bit] += ((CS[bit]>>(16*repeat))&0xFFFF)<<1;

      for(int i = 0; i < 64; i++)
         if(missingCount[i] < idCount)
            freq[i] = tempCount[i] * 0.5 / (idCount - missingCount[i]);

      // Informative Count
      for(int i = 0; i < 16; i++) CS[i] = 0;
      for(int i = 0; i < N_Trio; i++){
         word = (~IsMissing[trioID[2][i]]) & // child not missing
            ( (IsAa[trioID[0][i]] & (~IsMissing[trioID[1][i]])) | ((~IsMissing[trioID[0][i]]) & IsAa[trioID[1][i]]) );
         for(int bit = 0; bit < 16; bit++)
            CS[bit] += ((word>>bit)&0x0001000100010001);
      }
      for(int bit = 0; bit < 16; bit++)
         for(int repeat = 0; repeat < 4; repeat++)
            informativeCount[repeat*16+bit] = (CS[bit]>>(16*repeat))&0xFFFF;

      // T Count
      for(int i = 0; i < 16; i++) CS[i] = 0;
      for(int i = 0; i < N_Trio; i++){
         word = (
            (IsAa[trioID[0][i]] & IsAA[trioID[1][i]] & IsAA[trioID[2][i]]) | // AaxAA->AA
            (IsAa[trioID[0][i]] & Isaa[trioID[1][i]] & IsAa[trioID[2][i]]) | // Aaxaa->Aa
            (IsAA[trioID[0][i]] & IsAa[trioID[1][i]] & IsAA[trioID[2][i]]) | // AAxAa->AA
            (Isaa[trioID[0][i]] & IsAa[trioID[1][i]] & IsAa[trioID[2][i]])  // aaxAa->Aa
            );
         for(int bit = 0; bit < 16; bit++)
            CS[bit] += ((word>>bit)&0x0001000100010001);
      }
      for(int bit = 0; bit < 16; bit++)
         for(int repeat = 0; repeat < 4; repeat++)
            TCount[repeat*16+bit] = (CS[bit]>>(16*repeat))&0xFFFF;

      // NT Count
      for(int i = 0; i < 16; i++) CS[i] = 0;
      for(int i = 0; i < N_Trio; i++){
         word = (
            (IsAa[trioID[0][i]] & IsAA[trioID[1][i]] & IsAa[trioID[2][i]]) | // AaxAA->Aa
            (IsAa[trioID[0][i]] & Isaa[trioID[1][i]] & Isaa[trioID[2][i]]) | // Aaxaa->aa
            (IsAA[trioID[0][i]] & IsAa[trioID[1][i]] & IsAa[trioID[2][i]]) | // AAxAa->Aa
            (Isaa[trioID[0][i]] & IsAa[trioID[1][i]] & Isaa[trioID[2][i]])  // aaxAa->aa
            );
         for(int bit = 0; bit < 16; bit++)
            CS[bit] += ((word>>bit)&0x0001000100010001);
      }
      for(int bit = 0; bit < 16; bit++)
         for(int repeat = 0; repeat < 4; repeat++)
            NTCount[repeat*16+bit] = (CS[bit]>>(16*repeat))&0xFFFF;

      // HHH Count
      for(int i = 0; i < 16; i++) CS[i] = 0;
      for(int i = 0; i < N_Trio; i++){
         word = IsAa[trioID[0][i]] & IsAa[trioID[1][i]] & IsAa[trioID[2][i]];  // AaxAa->Aa
         for(int bit = 0; bit < 16; bit++)
            CS[bit] += ((word>>bit)&0x0001000100010001);
      }
      for(int bit = 0; bit < 16; bit++)
         for(int repeat = 0; repeat < 4; repeat++){
            word = (CS[bit]>>(16*repeat))&0xFFFF;
            TCount[repeat*16+bit] += (int)word;
            NTCount[repeat*16+bit] += (int)word;
         }

      // double T count
      for(int i = 0; i < 16; i++) CS[i] = 0;
      for(int i = 0; i < N_Trio; i++){
         word = IsAa[trioID[0][i]] & IsAa[trioID[1][i]] & IsAA[trioID[2][i]];  // AaxAa->AA
         for(int bit = 0; bit < 16; bit++)
            CS[bit] += ((word>>bit)&0x0001000100010001);
      }
      for(int bit = 0; bit < 16; bit++)
         for(int repeat = 0; repeat < 4; repeat++){
            word = (CS[bit]>>(16*repeat))&0xFFFF;
            TCount[repeat*16+bit] += (int)(word<<1);
         }

      // double NT count
      for(int i = 0; i < 16; i++) CS[i] = 0;
      for(int i = 0; i < N_Trio; i++){
         word = IsAa[trioID[0][i]] & IsAa[trioID[1][i]] & Isaa[trioID[2][i]];  // AaxAa->aa
         for(int bit = 0; bit < 16; bit++)
            CS[bit] += ((word>>bit)&0x0001000100010001);
      }
      for(int bit = 0; bit < 16; bit++)
         for(int repeat = 0; repeat < 4; repeat++){
            word = (CS[bit]>>(16*repeat))&0xFFFF;
            NTCount[repeat*16+bit] += (int)(word<<1);
         }

      for(int j = 0; j < Bit64; j++){
         if( (b==longCount-1) && (b*Bit64+j >= markerCount) ) continue;
         p = TCount[j]+NTCount[j];
         if(p > 0.5)
            statistic = (TCount[j] - NTCount[j])*(TCount[j] - NTCount[j])/p;
         else
            statistic = 0;
         p = chidist(statistic, 1);

         int k = bigdataIdx[markerCount-1-b*Bit64-j];
         if(chromosomes.Length()){
            sprintf(buffer, "%d", chromosomes[k]);
            buffers[b].Add(buffer);
         }
         if(snpName.Length())
            sprintf(buffer, "\t%20s", (const char*)snpName[k]);
         else
            sprintf(buffer, "\tSNP%d", k);
         buffers[b].Add(buffer);
         if(positions.Length()){
            sprintf(buffer, "\t%10d", int(positions[k]*1000000+0.5));
            buffers[b].Add(buffer);
         }
         if(freq[j] < 0.5){
            sprintf(buffer, "\t%s\t%s",
               (const char*)alleleLabel[0][k], (const char*)alleleLabel[1][k]);
         }else{
            sprintf(buffer, "\t%s\t%s",
               (const char*)alleleLabel[1][k], (const char*)alleleLabel[0][k]);
         }
         buffers[b].Add(buffer);

         if(freq[j] < 0.5){
            sprintf(buffer, "\t%.3lf\t%d\t%d\t%d\t%.3lf\t%.3lf\t%.3G\n",
               freq[j], informativeCount[j], TCount[j], NTCount[j],
               NTCount[j]? TCount[j]*1.0/NTCount[j]: -9, statistic, p);
         }else{
            sprintf(buffer, "\t%.3lf\t%d\t%d\t%d\t%.3lf\t%.3lf\t%.3G\n",
               1-freq[j], informativeCount[j], NTCount[j], TCount[j],
               TCount[j]? NTCount[j]*1.0/TCount[j]: -9, statistic, p);
         }
         buffers[b].Add(buffer);
      }
   }  // end of omp for loop
   delete []IsAA;
   delete []IsAa;
   delete []Isaa;
   delete []IsMissing;
#ifdef _OPENMP
}  // extra bracket for omp
#endif
   for(int i = 0; i < 3; i++)
      delete []trioID[i];

   String pedfile = prefix;
   pedfile.Add("tdt.txt");
   FILE *fp = fopen(pedfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
   if(chromosomes.Length()) fprintf(fp, "Chr");
   fprintf(fp, "\t%20s", "SNP");
   if(positions.Length()) fprintf(fp, "\t%10s", "Pos");
   fprintf(fp, "\tA1\tA2\tFreq1\tN_Info\tT\tNT\tOR\tChisq\tP\n");
   for(int b = 0; b < buffers.Length(); b++)
      buffers[b].Write(fp);

   if(xmarkerCount)
      printf("\nScanning X-chromosome...");
   int AA, Aa, missing;
   for(int m = 0; m < xmarkerCount; m+=16){
      int b = m/16;
      for(int i = 0; i < 16; i++){
         missingCount[i] = 0;
         xfreq[i] = 0.0;
      }
      for(int u = 0; u < unrelatedList.Length(); u++){
         int i = geno[unrelatedList[u]];
         AA = XG[0][i][b] & XG[1][i][b];
         Aa = (~XG[0][i][b]) & XG[1][i][b];
         missing = (~XG[0][i][b]) & (~XG[1][i][b]);
         for(int j = 0; j < 16; j++){
            if(AA & shortbase[j])
               xfreq[j] ++;
            else if(Aa & shortbase[j])
               xfreq[j] += 0.5;
            else if(missing & shortbase[j])
               missingCount[j] ++;
         }
      }
      for(int j = 0; j < 16; j++)
         if(missingCount[j] < unrelatedList.Length())
            xfreq[j] /= (unrelatedList.Length() - missingCount[j]);
      for(int j = 0; j < 16; j++)
         TCount[j] = NTCount[j] = informativeCount[j] = 0;
      for(int i = 0; i < N_Trio; i++){
         id1 = geno[Ltrio[i*3+1]];
         id2 = geno[Ltrio[i*3+2]];
         id3 = geno[Ltrio[i*3]];   // id3 is the child
         informative = (XG[0][id3][b] | XG[1][id3][b]) & // offspring not missing
            ( (~XG[0][id1][b] & XG[1][id1][b] & (XG[0][id2][b] | XG[1][id2][b])) |
              (~XG[0][id2][b] & XG[1][id2][b] & (XG[0][id1][b] | XG[1][id1][b])) );
         if(!informative) continue; // no informative trios at any SNPs
         for(int j = 0; j < 16; j++)
            if(informative & shortbase[j])
               informativeCount[j]++;
         for(int k = 0; k < 2; k++){
            id1 = geno[Ltrio[i*3+1+k]];
            T = informative &   // Aa->AA or aa->Aa
              ( (~XG[0][id1][b] & XG[1][id1][b] & XG[0][id3][b] & XG[1][id3][b]) |
                (XG[0][id1][b] & ~XG[1][id1][b] & ~XG[0][id3][b] & XG[1][id3][b]) );
            NT = informative &   // AA->Aa or Aa->aa
              ( (~XG[0][id3][b] & XG[1][id3][b] & XG[0][id1][b] & XG[1][id1][b]) |
                (XG[0][id3][b] & ~XG[1][id3][b] & ~XG[0][id1][b] & XG[1][id1][b]) );
            if(T)
               for(int j = 0; j < 16; j++)
                  if(T & shortbase[j])
                     TCount[j]++;
            if(NT)
               for(int j = 0; j < 16; j++)
                  if(NT & shortbase[j])
                     NTCount[j]++;
         }
      }
      for(int j = 0; j < 16; j++){
         if(m+j >= xmarkerCount) continue;
         p = TCount[j]+NTCount[j];
         if(p > 0.5){
            temp = TCount[j] - NTCount[j];
            statistic = temp*temp/p;
         }else
            statistic = 0;
         p = chidist(statistic, 1);
         fprintf(fp, "X");
         if(xsnpName.Length())
            fprintf(fp, "\t%20s", (const char*)xsnpName[m+j]);
         else
            fprintf(fp, "SNPX%d", m+j+1);
         if(xpositions.Length())
            fprintf(fp, "\t%10d", int(xpositions[m+j]*1000000+0.5));
         if(xfreq[j] < 0.5){
            fprintf(fp, "\t%s\t%s\t%.3lf",
               (const char*)xalleleLabel[0][m+j], (const char*)xalleleLabel[1][m+j], xfreq[j]);
            fprintf(fp, "\t%d\t%d\t%d\t%.3lf\t%.3lf\t%.3G\n",
               informativeCount[j], TCount[j], NTCount[j],
               TCount[j]? NTCount[j]*1.0/TCount[j]: -9, statistic, p);
         }else{
            fprintf(fp, "\t%s\t%s\t%.3lf",
               (const char*)xalleleLabel[1][m+j], (const char*)xalleleLabel[0][m+j], 1-xfreq[j]);
            fprintf(fp, "\t%d\t%d\t%d\t%.3lf\t%.3lf\t%.3G\n",
               informativeCount[j], NTCount[j], TCount[j],
               TCount[j]? NTCount[j]*1.0/TCount[j]: -9, statistic, p);
         }
      }
   }
   printf("\n");
   fclose(fp);

   printf("TDT scan results for disease %s saved in file %s\n",
      (const char*)ped.affectionNames[0], (const char*)pedfile);
   printf("Genome-wide TDT scan ends at %s\n", currentTime());
}
*/
/*
void Engine::TDT()
{
   int missingCount[16], AA, Aa, missing, id1, id2, id3, informative, T, NT, HHH;
   int TCount[16], NTCount[16], informativeCount[16];
   double freq[16], xfreq[16], statistic, p, temp;
   if(geno.Length()==0) BuildShortBinary();

   printf("\nOptions in effect:\n");
   printf("\t--tdt\n");
   if(CoreCount)
      printf("\t--cpus %d\n", CoreCount);
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

   char buffer[1024];
   StringArray buffers;
   buffers.Dimension(shortCount);
   for(int i = 0; i < buffers.Length(); i++)
      buffers[i].Clear();

   int pos, k;

   ExtractUnrelated();
   printf("A subset of %d unrelated individuals are used to calculate allele frequencies\n",
      unrelatedList.Length());
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

   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   printf("Scanning autosomes with %d CPU cores...", CoreCount?CoreCount:defaultEfficientCoreCount);
   unsigned short int word;
#ifdef _OPENMP
   #pragma omp parallel for num_threads(CoreCount?CoreCount:defaultEfficientCoreCount) \
      private(missingCount, freq, AA, Aa, missing, word, \
      informative, TCount, NTCount, informativeCount, \
      id1, id2, id3, p, temp, statistic, buffer)
#endif
   for(int b = 0; b < shortCount; b++){
      for(int i = 0; i < 16; i++){
         missingCount[i] = 0;
         freq[i] = 0.0;
      }
      for(int u = 0; u < unrelatedList.Length(); u++){
         int i = geno[unrelatedList[u]];
         AA = GG[0][i][b] & GG[1][i][b];
         Aa = (~GG[0][i][b]) & GG[1][i][b];
         missing = (~GG[0][i][b]) & (~GG[1][i][b]) & 65535;
         for(word = AA; word; word &= (word-1))
            freq[oneoneCount[(word&(-word))-1]] += 1.0;
         for(word = Aa; word; word &= (word-1))
            freq[oneoneCount[(word&(-word))-1]] += 0.5;
         for(word = missing; word; word &= (word-1))
            missingCount[oneoneCount[(word&(-word))-1]] ++;
      }
      for(int j = 0; j < 16; j++)
         if(missingCount[j] < unrelatedList.Length())
            freq[j] /= (unrelatedList.Length() - missingCount[j]);
      for(int j = 0; j < 16; j++)
         TCount[j] = NTCount[j] = informativeCount[j] = 0;
      for(int i = 0; i < N_Trio; i++){
         id1 = geno[Ltrio[i*3+1]];
         id2 = geno[Ltrio[i*3+2]];
         id3 = geno[Ltrio[i*3]];   // id3 is the child
         informative = (GG[0][id3][b] | GG[1][id3][b]) & // offspring not missing
            ( (~GG[0][id1][b] & GG[1][id1][b] & (GG[0][id2][b] | GG[1][id2][b])) |
              (~GG[0][id2][b] & GG[1][id2][b] & (GG[0][id1][b] | GG[1][id1][b])) );
         if(!informative) continue; // no informative trios at any SNPs
         for(word = informative; word; word &= (word-1))
            informativeCount[oneoneCount[(word&(-word))-1]] ++;
         for(word = ~GG[0][id1][b] & GG[1][id1][b] &
               ~GG[0][id2][b] & GG[1][id2][b] &
               ~GG[0][id3][b] & GG[1][id3][b]; // Aa x Aa -> Aa
            word; word &= (word-1)){
            Aa = oneoneCount[(word&(-word))-1];
            TCount[Aa] ++;
            NTCount[Aa] ++;
         }
         for(int k = 0; k < 2; k++){
            id1 = geno[Ltrio[i*3+1+k]];
            for(word = informative &   // Aa->AA or aa->Aa
              ( (~GG[0][id1][b] & GG[1][id1][b] & GG[0][id3][b] & GG[1][id3][b]) |
                (GG[0][id1][b] & ~GG[1][id1][b] & ~GG[0][id3][b] & GG[1][id3][b]) );
               word; word &= (word-1))
               TCount[oneoneCount[(word&(-word))-1]] ++;
            for(word = informative &   // AA->Aa or Aa->aa
              ( (~GG[0][id3][b] & GG[1][id3][b] & GG[0][id1][b] & GG[1][id1][b]) |
                (GG[0][id3][b] & ~GG[1][id3][b] & ~GG[0][id1][b] & GG[1][id1][b]) );
               word; word &= (word-1))
               NTCount[oneoneCount[(word&(-word))-1]] ++;
         }
      }
      for(int j = 0; j < 16; j++){
         if( (b==shortCount-1) && (b*16+j >= markerCount) ) continue;
         p = TCount[j]+NTCount[j];
         if(p > 0.5){
            temp = TCount[j] - NTCount[j];
            statistic = temp*temp/p;
         }else
            statistic = 0;
         p = chidist(statistic, 1);
            if(snpName.Length())
               sprintf(buffer, "%s", (const char*)snpName[b*16+j]);
            else
               sprintf(buffer, "SNP%d", b*16+j+1);
            buffers[b].Add(buffer);
            if(chromosomes.Length()){
               sprintf(buffer, "\t%d", chromosomes[b*16+j]);
               buffers[b].Add(buffer);
            }
            if(positions.Length()){
               sprintf(buffer, "\t%d", int(positions[b*16+j]*1000000+0.5));
               buffers[b].Add(buffer);
            }
            if(freq[j] < 0.5){
               sprintf(buffer, "\t%s\t%s",
                  (const char*)alleleLabel[0][b*16+j], (const char*)alleleLabel[1][b*16+j]);
            }else{
               sprintf(buffer, "\t%s\t%s",
                  (const char*)alleleLabel[1][b*16+j], (const char*)alleleLabel[0][b*16+j]);
            }
            buffers[b].Add(buffer);
         if(freq[j] < 0.5){
            sprintf(buffer, "\t%.3lf\t%d\t%d\t%d\t%.3lf\t%.3lf\t%.3G\n",
               freq[j], informativeCount[j], TCount[j], NTCount[j],
               NTCount[j]? TCount[j]*1.0/NTCount[j]: -9, statistic, p);
         }else{
            sprintf(buffer, "\t%.3lf\t%d\t%d\t%d\t%.3lf\t%.3lf\t%.3G\n",
               1-freq[j], informativeCount[j], NTCount[j], TCount[j],
               TCount[j]? NTCount[j]*1.0/TCount[j]: -9, statistic, p);
         }
         buffers[b].Add(buffer);
      }
   }
   String pedfile = prefix;
   pedfile.Add("tdt.txt");
   FILE *fp = fopen(pedfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
   fprintf(fp, "SNP");
   if(chromosomes.Length()) fprintf(fp, "\tChr");
   if(positions.Length()) fprintf(fp, "\tPos");
   fprintf(fp, "\tA1\tA2\tFreq1\tN_Info\tT\tNT\tOR\tChisq\tP\n");
   for(int b = 0; b < buffers.Length(); b++)
      buffers[b].Write(fp);

   if(xmarkerCount)
      printf("\nScanning X-chromosome...");
   for(int m = 0; m < xmarkerCount; m+=16){
      int b = m/16;
      for(int i = 0; i < 16; i++){
         missingCount[i] = 0;
         xfreq[i] = 0.0;
      }
      for(int u = 0; u < unrelatedList.Length(); u++){
         int i = geno[unrelatedList[u]];
         AA = XG[0][i][b] & XG[1][i][b];
         Aa = (~XG[0][i][b]) & XG[1][i][b];
         missing = (~XG[0][i][b]) & (~XG[1][i][b]);
         for(int j = 0; j < 16; j++){
            if(AA & shortbase[j])
               xfreq[j] ++;
            else if(Aa & shortbase[j])
               xfreq[j] += 0.5;
            else if(missing & shortbase[j])
               missingCount[j] ++;
         }
      }
      for(int j = 0; j < 16; j++)
         if(missingCount[j] < unrelatedList.Length())
            xfreq[j] /= (unrelatedList.Length() - missingCount[j]);
      for(int j = 0; j < 16; j++)
         TCount[j] = NTCount[j] = informativeCount[j] = 0;
      for(int i = 0; i < N_Trio; i++){
         id1 = geno[Ltrio[i*3+1]];
         id2 = geno[Ltrio[i*3+2]];
         id3 = geno[Ltrio[i*3]];   // id3 is the child
         informative = (XG[0][id3][b] | XG[1][id3][b]) & // offspring not missing
            ( (~XG[0][id1][b] & XG[1][id1][b] & (XG[0][id2][b] | XG[1][id2][b])) |
              (~XG[0][id2][b] & XG[1][id2][b] & (XG[0][id1][b] | XG[1][id1][b])) );
         if(!informative) continue; // no informative trios at any SNPs
         for(int j = 0; j < 16; j++)
            if(informative & shortbase[j])
               informativeCount[j]++;
         for(int k = 0; k < 2; k++){
            id1 = geno[Ltrio[i*3+1+k]];
            T = informative &   // Aa->AA or aa->Aa
              ( (~XG[0][id1][b] & XG[1][id1][b] & XG[0][id3][b] & XG[1][id3][b]) |
                (XG[0][id1][b] & ~XG[1][id1][b] & ~XG[0][id3][b] & XG[1][id3][b]) );
            NT = informative &   // AA->Aa or Aa->aa
              ( (~XG[0][id3][b] & XG[1][id3][b] & XG[0][id1][b] & XG[1][id1][b]) |
                (XG[0][id3][b] & ~XG[1][id3][b] & ~XG[0][id1][b] & XG[1][id1][b]) );
            if(T)
               for(int j = 0; j < 16; j++)
                  if(T & shortbase[j])
                     TCount[j]++;
            if(NT)
               for(int j = 0; j < 16; j++)
                  if(NT & shortbase[j])
                     NTCount[j]++;
         }
      }
      for(int j = 0; j < 16; j++){
         if(m+j >= xmarkerCount) continue;
         p = TCount[j]+NTCount[j];
         if(p > 0.5){
            temp = TCount[j] - NTCount[j];
            statistic = temp*temp/p;
         }else
            statistic = 0;
         p = chidist(statistic, 1);
         if(xsnpName.Length())
            fprintf(fp, "%s", (const char*)xsnpName[m+j]);
         else
            fprintf(fp, "SNPX%d", m+j+1);
         fprintf(fp, "\tX");
         if(xpositions.Length())
            fprintf(fp, "\t%.6lf", xpositions[m+j]);
         if(xfreq[j] < 0.5){
            fprintf(fp, "\t%s\t%s\t%.3lf",
               (const char*)xalleleLabel[0][m+j], (const char*)xalleleLabel[1][m+j], xfreq[j]);
            fprintf(fp, "\t%d\t%d\t%d\t%.3lf\t%.3lf\t%.3G\n",
               informativeCount[j], TCount[j], NTCount[j],
               TCount[j]? NTCount[j]*1.0/TCount[j]: -9, statistic, p);
         }else{
            fprintf(fp, "\t%s\t%s\t%.3lf",
               (const char*)xalleleLabel[1][m+j], (const char*)xalleleLabel[0][m+j], 1-xfreq[j]);
            fprintf(fp, "\t%d\t%d\t%d\t%.3lf\t%.3lf\t%.3G\n",
               informativeCount[j], NTCount[j], TCount[j],
               TCount[j]? NTCount[j]*1.0/TCount[j]: -9, statistic, p);
         }
      }
   }
   fclose(fp);
   printf("\n");                    
   printf("TDT scan results for disease %s saved in file %s\n",
      (const char*)ped.affectionNames[0], (const char*)pedfile);
   printf("Genome-wide TDT scan ends at %s\n", currentTime());
}

void Engine::TDT_64bit()
{
   unsigned long long int informative, T, NT, HHH;
   int missingCount[64], id1, id2, id3;
   int TCount[64], NTCount[64], informativeCount[64];
   double freq[64], xfreq[64], statistic, p, temp;
   if(geno.Length()==0) error("Genotype data not available");//BuildShortBinary();

   printf("\nOptions in effect:\n");
   printf("\t--tdt\n");
   if(CoreCount)
      printf("\t--cpus %d\n", CoreCount);
   if(Bit64)
      printf("\t--bit64\n");
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");

   printf("Genome-wide TDT scan starts at %s", currentTime());
   MakeTrioForTDT_64bit();
   int N_Trio = Ltrio.Length() / 3;
   if(N_Trio==0){
      warning("TDT analysis requires parent-affected-offspring trios.\n");
      return;
   }

   char buffer[1024];
   StringArray buffers;
   buffers.Dimension(longCount);
   for(int i = 0; i < buffers.Length(); i++)
      buffers[i].Clear();

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
   unsigned long long int longbase[64];
   for(int i = 0; i < Bit64; i++)
      longbase[i] = ((unsigned long long int)1 << i);

   int *trioID[3];
   for(int j = 0; j < 3; j++)
      trioID[j] = new int[N_Trio];
   for(int i = 0; i < N_Trio; i++){
      trioID[0][i] = geno[Ltrio[i*3+1]];
      trioID[1][i] = geno[Ltrio[i*3+2]];
      trioID[2][i] = geno[Ltrio[i*3]]; // child
   }
#ifdef _OPENMP
   printf("Scanning autosomes with %d CPU cores...",
      CoreCount?CoreCount:defaultEfficientCoreCount);
   #pragma omp parallel num_threads(CoreCount?CoreCount:defaultEfficientCoreCount) \
   private(missingCount, freq, TCount, NTCount, informativeCount, p, statistic, buffer)
{
#endif
   unsigned long long int *IsAA, *IsAa, *Isaa, *IsMissing;
   IsAA = new unsigned long long int [idCount];
   IsAa = new unsigned long long int [idCount];
   Isaa = new unsigned long long int [idCount];
   IsMissing = new unsigned long long int [idCount];
   unsigned long long int word;
   unsigned long long int CS[16];
   int tempCount[64];
#ifdef _OPENMP
   #pragma omp for
#endif
   for(int b = 0; b < longCount; b++){
      for(int i = 0; i < idCount; i++){
         IsAA[i] = LG[0][i][b] & LG[1][i][b];
         IsAa[i] = (~LG[0][i][b]) & LG[1][i][b];
         Isaa[i] = LG[0][i][b] & (~LG[1][i][b]);
         IsMissing[i] = (~LG[0][i][b]) & (~LG[1][i][b]);
      }

      // Missing Count
      for(int i = 0; i < 16; i++) CS[i] = 0;
      for(int i = 0; i < idCount; i++){
         word = IsMissing[i];
         for(int bit = 0; bit < 16; bit++)
            CS[bit] += ((word>>bit)&0x0001000100010001);
      }
      for(int bit = 0; bit < 16; bit++)
         for(int repeat = 0; repeat < 4; repeat++)
            missingCount[repeat*16+bit] = (CS[bit]>>(16*repeat))&0xFFFF;

      // Aa Count
      for(int i = 0; i < 16; i++) CS[i] = 0;
      for(int i = 0; i < idCount; i++){
         word = IsAa[i];
         for(int bit = 0; bit < 16; bit++)
            CS[bit] += ((word>>bit)&0x0001000100010001);
      }
      for(int bit = 0; bit < 16; bit++)
         for(int repeat = 0; repeat < 4; repeat++)
            tempCount[repeat*16+bit] = (CS[bit]>>(16*repeat))&0xFFFF;

      // AA Count
      for(int i = 0; i < 16; i++) CS[i] = 0;
      for(int i = 0; i < idCount; i++){
         word = IsAA[i];
         for(int bit = 0; bit < 16; bit++)
            CS[bit] += ((word>>bit)&0x0001000100010001);
      }
      for(int bit = 0; bit < 16; bit++)
         for(int repeat = 0; repeat < 4; repeat++)
            tempCount[repeat*16+bit] += ((CS[bit]>>(16*repeat))&0xFFFF)<<1;

      for(int i = 0; i < 64; i++)
         if(missingCount[i] < idCount)
            freq[i] = tempCount[i] * 0.5 / (idCount - missingCount[i]);

      // Informative Count
      for(int i = 0; i < 16; i++) CS[i] = 0;
      for(int i = 0; i < N_Trio; i++){
         word = (~IsMissing[trioID[2][i]]) & // child not missing
            ( (IsAa[trioID[0][i]] & (~IsMissing[trioID[1][i]])) | ((~IsMissing[trioID[0][i]]) & IsAa[trioID[1][i]]) );
         for(int bit = 0; bit < 16; bit++)
            CS[bit] += ((word>>bit)&0x0001000100010001);
      }
      for(int bit = 0; bit < 16; bit++)
         for(int repeat = 0; repeat < 4; repeat++)
            informativeCount[repeat*16+bit] = (CS[bit]>>(16*repeat))&0xFFFF;

      // T Count
      for(int i = 0; i < 16; i++) CS[i] = 0;
      for(int i = 0; i < N_Trio; i++){
         word = (
            (IsAa[trioID[0][i]] & IsAA[trioID[1][i]] & IsAA[trioID[2][i]]) | // AaxAA->AA
            (IsAa[trioID[0][i]] & Isaa[trioID[1][i]] & IsAa[trioID[2][i]]) | // Aaxaa->Aa
            (IsAA[trioID[0][i]] & IsAa[trioID[1][i]] & IsAA[trioID[2][i]]) | // AAxAa->AA
            (Isaa[trioID[0][i]] & IsAa[trioID[1][i]] & IsAa[trioID[2][i]])  // aaxAa->Aa
            );
         for(int bit = 0; bit < 16; bit++)
            CS[bit] += ((word>>bit)&0x0001000100010001);
      }
      for(int bit = 0; bit < 16; bit++)
         for(int repeat = 0; repeat < 4; repeat++)
            TCount[repeat*16+bit] = (CS[bit]>>(16*repeat))&0xFFFF;

      // NT Count
      for(int i = 0; i < 16; i++) CS[i] = 0;
      for(int i = 0; i < N_Trio; i++){
         word = (
            (IsAa[trioID[0][i]] & IsAA[trioID[1][i]] & IsAa[trioID[2][i]]) | // AaxAA->Aa
            (IsAa[trioID[0][i]] & Isaa[trioID[1][i]] & Isaa[trioID[2][i]]) | // Aaxaa->aa
            (IsAA[trioID[0][i]] & IsAa[trioID[1][i]] & IsAa[trioID[2][i]]) | // AAxAa->Aa
            (Isaa[trioID[0][i]] & IsAa[trioID[1][i]] & Isaa[trioID[2][i]])  // aaxAa->aa
            );
         for(int bit = 0; bit < 16; bit++)
            CS[bit] += ((word>>bit)&0x0001000100010001);
      }
      for(int bit = 0; bit < 16; bit++)
         for(int repeat = 0; repeat < 4; repeat++)
            NTCount[repeat*16+bit] = (CS[bit]>>(16*repeat))&0xFFFF;

      // HHH Count
      for(int i = 0; i < 16; i++) CS[i] = 0;
      for(int i = 0; i < N_Trio; i++){
         word = IsAa[trioID[0][i]] & IsAa[trioID[1][i]] & IsAa[trioID[2][i]];  // AaxAa->Aa
         for(int bit = 0; bit < 16; bit++)
            CS[bit] += ((word>>bit)&0x0001000100010001);
      }
      for(int bit = 0; bit < 16; bit++)
         for(int repeat = 0; repeat < 4; repeat++){
            word = (CS[bit]>>(16*repeat))&0xFFFF;
            TCount[repeat*16+bit] += (int)word;
            NTCount[repeat*16+bit] += (int)word;
         }

      // double T count
      for(int i = 0; i < 16; i++) CS[i] = 0;
      for(int i = 0; i < N_Trio; i++){
         word = IsAa[trioID[0][i]] & IsAa[trioID[1][i]] & IsAA[trioID[2][i]];  // AaxAa->AA
         for(int bit = 0; bit < 16; bit++)
            CS[bit] += ((word>>bit)&0x0001000100010001);
      }
      for(int bit = 0; bit < 16; bit++)
         for(int repeat = 0; repeat < 4; repeat++){
            word = (CS[bit]>>(16*repeat))&0xFFFF;
            TCount[repeat*16+bit] += (int)(word<<1);
         }

      // double NT count
      for(int i = 0; i < 16; i++) CS[i] = 0;
      for(int i = 0; i < N_Trio; i++){
         word = IsAa[trioID[0][i]] & IsAa[trioID[1][i]] & Isaa[trioID[2][i]];  // AaxAa->aa
         for(int bit = 0; bit < 16; bit++)
            CS[bit] += ((word>>bit)&0x0001000100010001);
      }
      for(int bit = 0; bit < 16; bit++)
         for(int repeat = 0; repeat < 4; repeat++){
            word = (CS[bit]>>(16*repeat))&0xFFFF;
            NTCount[repeat*16+bit] += (int)(word<<1);
         }

      for(int j = 0; j < Bit64; j++){
         if( (b==longCount-1) && (b*Bit64+j >= markerCount) ) continue;
         p = TCount[j]+NTCount[j];
         if(p > 0.5)
            statistic = (TCount[j] - NTCount[j])*(TCount[j] - NTCount[j])/p;
         else
            statistic = 0;
         p = chidist(statistic, 1);

         int k = bigdataIdx[markerCount-1-b*Bit64-j];
         if(chromosomes.Length()){
            sprintf(buffer, "%d", chromosomes[k]);
            buffers[b].Add(buffer);
         }
         if(snpName.Length())
            sprintf(buffer, "\t%20s", (const char*)snpName[k]);
         else
            sprintf(buffer, "\tSNP%d", k);
         buffers[b].Add(buffer);
         if(positions.Length()){
            sprintf(buffer, "\t%10d", int(positions[k]*1000000+0.5));
            buffers[b].Add(buffer);
         }
         if(freq[j] < 0.5){
            sprintf(buffer, "\t%s\t%s",
               (const char*)alleleLabel[0][k], (const char*)alleleLabel[1][k]);
         }else{
            sprintf(buffer, "\t%s\t%s",
               (const char*)alleleLabel[1][k], (const char*)alleleLabel[0][k]);
         }
         buffers[b].Add(buffer);

         if(freq[j] < 0.5){
            sprintf(buffer, "\t%.3lf\t%d\t%d\t%d\t%.3lf\t%.3lf\t%.3G\n",
               freq[j], informativeCount[j], TCount[j], NTCount[j],
               NTCount[j]? TCount[j]*1.0/NTCount[j]: -9, statistic, p);
         }else{
            sprintf(buffer, "\t%.3lf\t%d\t%d\t%d\t%.3lf\t%.3lf\t%.3G\n",
               1-freq[j], informativeCount[j], NTCount[j], TCount[j],
               TCount[j]? NTCount[j]*1.0/TCount[j]: -9, statistic, p);
         }
         buffers[b].Add(buffer);
      }
   }  // end of omp for loop
   delete []IsAA;
   delete []IsAa;
   delete []Isaa;
   delete []IsMissing;
#ifdef _OPENMP
}  // extra bracket for omp
#endif
   for(int i = 0; i < 3; i++)
      delete []trioID[i];

   String pedfile = prefix;
   pedfile.Add("tdt.txt");
   FILE *fp = fopen(pedfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
   if(chromosomes.Length()) fprintf(fp, "Chr");
   fprintf(fp, "\t%20s", "SNP");
   if(positions.Length()) fprintf(fp, "\t%10s", "Pos");
   fprintf(fp, "\tA1\tA2\tFreq1\tN_Info\tT\tNT\tOR\tChisq\tP\n");
   for(int b = 0; b < buffers.Length(); b++)
      buffers[b].Write(fp);

   if(xmarkerCount)
      printf("\nScanning X-chromosome...");
   int AA, Aa, missing;
   for(int m = 0; m < xmarkerCount; m+=16){
      int b = m/16;
      for(int i = 0; i < 16; i++){
         missingCount[i] = 0;
         xfreq[i] = 0.0;
      }
      for(int u = 0; u < unrelatedList.Length(); u++){
         int i = geno[unrelatedList[u]];
         AA = XG[0][i][b] & XG[1][i][b];
         Aa = (~XG[0][i][b]) & XG[1][i][b];
         missing = (~XG[0][i][b]) & (~XG[1][i][b]);
         for(int j = 0; j < 16; j++){
            if(AA & shortbase[j])
               xfreq[j] ++;
            else if(Aa & shortbase[j])
               xfreq[j] += 0.5;
            else if(missing & shortbase[j])
               missingCount[j] ++;
         }
      }
      for(int j = 0; j < 16; j++)
         if(missingCount[j] < unrelatedList.Length())
            xfreq[j] /= (unrelatedList.Length() - missingCount[j]);
      for(int j = 0; j < 16; j++)
         TCount[j] = NTCount[j] = informativeCount[j] = 0;
      for(int i = 0; i < N_Trio; i++){
         id1 = geno[Ltrio[i*3+1]];
         id2 = geno[Ltrio[i*3+2]];
         id3 = geno[Ltrio[i*3]];   // id3 is the child
         informative = (XG[0][id3][b] | XG[1][id3][b]) & // offspring not missing
            ( (~XG[0][id1][b] & XG[1][id1][b] & (XG[0][id2][b] | XG[1][id2][b])) |
              (~XG[0][id2][b] & XG[1][id2][b] & (XG[0][id1][b] | XG[1][id1][b])) );
         if(!informative) continue; // no informative trios at any SNPs
         for(int j = 0; j < 16; j++)
            if(informative & shortbase[j])
               informativeCount[j]++;
         for(int k = 0; k < 2; k++){
            id1 = geno[Ltrio[i*3+1+k]];
            T = informative &   // Aa->AA or aa->Aa
              ( (~XG[0][id1][b] & XG[1][id1][b] & XG[0][id3][b] & XG[1][id3][b]) |
                (XG[0][id1][b] & ~XG[1][id1][b] & ~XG[0][id3][b] & XG[1][id3][b]) );
            NT = informative &   // AA->Aa or Aa->aa
              ( (~XG[0][id3][b] & XG[1][id3][b] & XG[0][id1][b] & XG[1][id1][b]) |
                (XG[0][id3][b] & ~XG[1][id3][b] & ~XG[0][id1][b] & XG[1][id1][b]) );
            if(T)
               for(int j = 0; j < 16; j++)
                  if(T & shortbase[j])
                     TCount[j]++;
            if(NT)
               for(int j = 0; j < 16; j++)
                  if(NT & shortbase[j])
                     NTCount[j]++;
         }
      }
      for(int j = 0; j < 16; j++){
         if(m+j >= xmarkerCount) continue;
         p = TCount[j]+NTCount[j];
         if(p > 0.5){
            temp = TCount[j] - NTCount[j];
            statistic = temp*temp/p;
         }else
            statistic = 0;
         p = chidist(statistic, 1);
         fprintf(fp, "X");
         if(xsnpName.Length())
            fprintf(fp, "\t%20s", (const char*)xsnpName[m+j]);
         else
            fprintf(fp, "SNPX%d", m+j+1);
         if(xpositions.Length())
            fprintf(fp, "\t%10d", int(xpositions[m+j]*1000000+0.5));
         if(xfreq[j] < 0.5){
            fprintf(fp, "\t%s\t%s\t%.3lf",
               (const char*)xalleleLabel[0][m+j], (const char*)xalleleLabel[1][m+j], xfreq[j]);
            fprintf(fp, "\t%d\t%d\t%d\t%.3lf\t%.3lf\t%.3G\n",
               informativeCount[j], TCount[j], NTCount[j],
               TCount[j]? NTCount[j]*1.0/TCount[j]: -9, statistic, p);
         }else{
            fprintf(fp, "\t%s\t%s\t%.3lf",
               (const char*)xalleleLabel[1][m+j], (const char*)xalleleLabel[0][m+j], 1-xfreq[j]);
            fprintf(fp, "\t%d\t%d\t%d\t%.3lf\t%.3lf\t%.3G\n",
               informativeCount[j], NTCount[j], TCount[j],
               TCount[j]? NTCount[j]*1.0/TCount[j]: -9, statistic, p);
         }
      }
   }
   printf("\n");
   fclose(fp);

   printf("TDT scan results for disease %s saved in file %s\n",
      (const char*)ped.affectionNames[0], (const char*)pedfile);
   printf("Genome-wide TDT scan ends at %s\n", currentTime());
}

*/

/*
         }else{
            if(chromosomes.Length()){
               sprintf(buffer, "%d", chromosomes[b*Bit64+j]);
               buffers[b].Add(buffer);
            }
            if(snpName.Length())
               sprintf(buffer, "%20s", (const char*)snpName[b*Bit64+j]);
            else
               sprintf(buffer, "SNP%d", b*Bit64+j+1);
            buffers[b].Add(buffer);
            if(positions.Length()){
               sprintf(buffer, "\t%10d", int(positions[b*Bit64+j]*1000000+0.5));
               buffers[b].Add(buffer);
            }
            if(freq[j] < 0.5){
               sprintf(buffer, "\t%s\t%s",
                  (const char*)alleleLabel[0][b*Bit64+j], (const char*)alleleLabel[1][b*Bit64+j]);
            }else{
               sprintf(buffer, "\t%s\t%s",
                  (const char*)alleleLabel[1][b*Bit64+j], (const char*)alleleLabel[0][b*Bit64+j]);
            }
            buffers[b].Add(buffer);
         }
   */




