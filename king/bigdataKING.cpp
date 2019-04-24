//////////////////////////////////////////////////////////////////////
// bigdataKING.cpp
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
// Dec 12, 2018

#include <math.h>
#include "analysis.h"
#include "Kinship.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

void Engine::ComputeInnerProduct64Bit(IntArray & subset, IntArray * counts)
{
   unsigned long long int word, word0, word1, word2;
   int subsetCount = subset.Length();
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) private(word, word0)
#endif
   for(int i = 0; i < subsetCount; i++){
      int id = subset[i];
      word0 = 0;
      for(int m = 0; m < longCount; m++){
         word = LG[0][id][m];          // Hom
         word = word - ((word>>1)&0x5555555555555555);
         word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
         word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
         word = (word+(word>>8)) & 0x00FF00FF00FF00FF;
         word0 += (word+(word>>16)) & 0x0000FFFF0000FFFF;  // Hom
      }
      counts[i][i] = (word0+(word0>>32)) & 0xFFFFFFFF;
   }
   const int BLOCKSIZE=8;
   const int CACHESIZE=512;   // 2^(9+4+4)=128K total cache
   int localD[BLOCKSIZE][BLOCKSIZE];
   IntArray loopIndex[2];
   loopIndex[0].Dimension(0); loopIndex[1].Dimension(0);
   for(int i = 0; i < subsetCount; i += BLOCKSIZE)
      for(int j = i; j < subsetCount; j += BLOCKSIZE){
         loopIndex[0].Push(i);
         loopIndex[1].Push(j);
      }
   int loopIndexLength = loopIndex[0].Length();
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
      private(word0, word1, word2, word, localD)
#endif
   for(int k = 0; k < loopIndexLength; k++){
      int i = loopIndex[0][k];
      int iMax = i<subsetCount-BLOCKSIZE? i+BLOCKSIZE: subsetCount;
      int j = loopIndex[1][k];
      int jMax = j<subsetCount-BLOCKSIZE? j+BLOCKSIZE: subsetCount;
      int jMin = j;
      for(int ii = 0; ii < BLOCKSIZE; ii++)
         for(int jj = 0; jj < BLOCKSIZE; jj++)
            localD[ii][jj] = 0;
      for(int mm = 0; mm < longCount; mm += CACHESIZE){
         int mMax = (mm > longCount-CACHESIZE) ? longCount: mm+CACHESIZE;
         for(int i1 = i; i1 < iMax; i1++){
            char ii = i1 - i;
            int id1 = subset[i1];
            for(int i2 = j; i2 < jMax; i2++){
               char jj = i2 - j;
               int id2 = subset[i2];
               word1 = word2 = 0;
               for(int m = mm; m < mMax; m++){
                  word0 = LG[0][id1][m] & LG[0][id2][m];          // HomHom
                  word = word0 - ((word0>>1)&0x5555555555555555);
                  word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
                  word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
                  word1 += (word+(word>>8)) & 0x00FF00FF00FF00FF;
                  word0 &= (LG[1][id1][m] ^ LG[1][id2][m]);       // IBS0
                  word = word0 - ((word0>>1)&0x5555555555555555);  // IBS0
                  word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
                  word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
                  word2 += (word+(word>>8)) & 0x00FF00FF00FF00FF;
               }
               word1 = (word1+(word1>>16)) & 0x0000FFFF0000FFFF;  // HomHom
               localD[ii][jj] += (word1+(word1>>32)) & 0xFFFFFFFF;
               word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;  // IBS0
               localD[ii][jj] -= (((word2+(word2>>32)) & 0xFFFFFFFF)<<1);
            }  // end of id2
         }  // end of id1
      }// end of mm
      for(int id1 = i; id1 < iMax; id1++){
         char ii = id1 - i;
         if(i == j) jMin = id1 + 1;
         for(int id2 = jMin; id2 < jMax; id2++){
            char jj = id2 - j;
            counts[id2][id1] = counts[id1][id2] = localD[ii][jj];
         }  // end of id2
      }  // end of id1
   }  // end of pairs
}

void Engine::ComputeDistanceMatrix64Bit(Matrix & Dist)
{
   int notMissingCount, id1, id2, sum;
   unsigned long long int word, word1, word2;
   int *missingInOnePersonCount = new int[idCount];
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
      private(sum, word)
#endif
   for(int i = 0; i < idCount; i++){
      sum = 0;
      for(int m = 0; m < longCount; m++)   // not all non-missing
         for(word = ~(LG[0][i][m] | LG[1][i][m]); word; word &= (word-1), sum++);
      missingInOnePersonCount[i] = sum;
   }  // parallel among individuals ends
   const int cutoffMissingCount = markerCount-MINSNPCOUNT;  // almost all missing
   IntArray outindex(idCount);
   outindex.Set(-1);
   sum = 0;
   for(int i = 0; i < idCount; i++)
      if(missingInOnePersonCount[i] < cutoffMissingCount)
         outindex[i] = sum++;
   Dist.Dimension(sum, sum);
   Dist.Zero();
#ifdef _OPENMP
   printf("%d CPU cores are used.\n", defaultMaxCoreCount);
#endif
   const int BLOCKSIZE=8;
   const int CACHESIZE=512;   // 2^(9+4+4)=128K total cache
   int localD[2][BLOCKSIZE][BLOCKSIZE];
   IntArray loopIndex[2];
   loopIndex[0].Dimension(0); loopIndex[1].Dimension(0);
   for(int i = 0; i < idCount; i += BLOCKSIZE)
      for(int j = i; j < idCount; j += BLOCKSIZE){
         loopIndex[0].Push(i);
         loopIndex[1].Push(j);
      }
   int loopIndexLength = loopIndex[0].Length();
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
      private(notMissingCount, id1, id2, word1, word2, word, localD)
#endif
   for(int k = 0; k < loopIndexLength; k++){
      int i = loopIndex[0][k];
      int iMax = i<idCount-BLOCKSIZE? i+BLOCKSIZE: idCount;
      int j = loopIndex[1][k];
      int jMax = j<idCount-BLOCKSIZE? j+BLOCKSIZE: idCount;
      int jMin = j;
      for(int c = 0; c < 2; c++)
         for(int ii = 0; ii < BLOCKSIZE; ii++)
            for(int jj = 0; jj < BLOCKSIZE; jj++)
               localD[c][ii][jj] = 0;
      for(int mm = 0; mm < longCount; mm += CACHESIZE){
         int mMax = (mm > longCount-CACHESIZE) ? longCount: mm+CACHESIZE;
         for(id1 = i; id1 < iMax; id1++){
            char ii = id1 - i;
            for(id2 = j; id2 < jMax; id2++){
               char jj = id2 - j;
               word1 = word2 = notMissingCount = 0;
               for(int m = mm; m < mMax; m++){
                  for(word = ~(LG[0][id1][m] | LG[0][id2][m] | LG[1][id1][m] | LG[1][id2][m]);
                     word; word &= (word-1), notMissingCount++);  // MissMiss
                  word = (LG[0][id1][m]^LG[0][id2][m]) & (LG[0][id1][m] | LG[1][id1][m] ) & (LG[0][id2][m]| LG[1][id2][m]); // HetHom
                  word = word - ((word>>1)&0x5555555555555555);
                  word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
                  word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
                  word1 += (word+(word>>8)) & 0x00FF00FF00FF00FF;
                  word = LG[0][id1][m] & LG[0][id2][m] & (LG[1][id1][m] ^ LG[1][id2][m]);
                  word = word - ((word>>1)&0x5555555555555555);   // IBS0
                  word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
                  word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
                  word2 += (word+(word>>8)) & 0x00FF00FF00FF00FF;
               }
               word1 = (word1+(word1>>16)) & 0x0000FFFF0000FFFF;  // HetHom
               localD[0][ii][jj] += (word1+(word1>>32)) & 0xFFFFFFFF;
               word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;  // IBS0
               localD[0][ii][jj] += ((word2+(word2>>32)) & 0xFFFFFFFF)<<2;
               localD[1][ii][jj] += notMissingCount;  // MissMiss
            }  // end of id2
         }  // end of id1
      }// end of mm
      for(id1 = i; id1 < iMax; id1++){
         if(missingInOnePersonCount[id1]>=cutoffMissingCount) continue;
         char ii = id1 - i;
         if(i == j) jMin = id1 + 1;
         for(id2 = jMin; id2 < jMax; id2++){
            char jj = id2 - j;
            if(missingInOnePersonCount[id2]>=cutoffMissingCount) continue;
            localD[1][ii][jj] += markerCount - missingInOnePersonCount[id1] - missingInOnePersonCount[id2];
            Dist[outindex[id2]][outindex[id1]] = Dist[outindex[id1]][outindex[id2]] = (localD[0][ii][jj]*1.0/localD[1][ii][jj]);
         }  // end of id2
      }  // end of id1
   }  // end of pairs
   delete []missingInOnePersonCount;
}

long long int Engine::ScreenDuplicates64Bit(IntArray rpList[])
{
   if(longCount < 8 || !bigdataFlag) return -1;  // screening is not needed
   int long_prescan = 8; // 512 SNPs
   int prescanMarkerCount = long_prescan==longCount? markerCount: (long_prescan<<6);
   int **missingInOnePersonCount = new int * [idCount];
   for(int i = 0; i < idCount; i++)
      missingInOnePersonCount[i] = new int [2];
   int *hetInOnePersonCount = new int[idCount];
   for(int i = 0; i < idCount; i++)
      missingInOnePersonCount[i][0] = missingInOnePersonCount[i][1] = hetInOnePersonCount[i] = 0;
   unsigned long long int *masks = new unsigned long long int[long_prescan];
   for(int i = 0; i < long_prescan;  i++)
      masks[i] = 0xFFFFFFFFFFFFFFFF;
   if(long_prescan == longCount)
      masks[long_prescan-1] = markerCount%64? ((1 << (markerCount % 64))-1): 0xFFFFFFFFFFFFFFFF;
   unsigned long long int word, word1, word2;
   int m1, m2, m3, sum;
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
      private(sum, m1, word, word2)
#endif
   for(int i = 0; i < idCount; i++){
      for(sum=0, word = ~(SLG[0][i][0] | SLG[1][i][0]); word; word &= (word-1), sum++);
      missingInOnePersonCount[i][0] = sum;
      for(m1 = 1; m1 < long_prescan; m1++)
         for(word = ~(SLG[0][i][m1] | SLG[1][i][m1]) & masks[m1];
            word; word &= (word-1), sum++);
      missingInOnePersonCount[i][1] = sum;
      for(word2 = m1 = 0; m1 < long_prescan; m1++){
         word = (~SLG[0][i][m1]) & SLG[1][i][m1];  // Het
         word = word - ((word>>1)&0x5555555555555555);
         word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
         word2 += (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
      }
      word2 = (word2+(word2>>8)) & 0x00FF00FF00FF00FF;
      word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;
      hetInOnePersonCount[i] = (word2+(word2>>32)) & 0xFFFFFFFF;
   }  // parallel among individuals ends
   delete []masks;
   int id1, id2, HetHomCount, HomHomCount, notMissingHetCount;
   for(int i = 0; i < defaultMaxCoreCount; i++)
      rpList[i].Dimension(0);
   short int LOOPBLOCKINGSIZE=256;  // use 64KB cache
   char SHIFTSIZE=8; // 2^8 * 2 * 2^9 / 2^2 = 2^16 = 64kB
   if(idCount > (1<<23)){  // sample size over 8 million
      int count = (idCount >> 23);
      int F = 0;
      for(; count; F++) count >>= 1;
      LOOPBLOCKINGSIZE <<= F;
      SHIFTSIZE += F;
   }
   unsigned int blockCount = (idCount-1)/LOOPBLOCKINGSIZE+1;
   unsigned int loopIndexLength = (blockCount & 1)? blockCount * ((blockCount+1)/2):
      (blockCount/2) * (blockCount+1);
   unsigned int *loopIndex = new unsigned int [loopIndexLength];
   unsigned int index = 0;
   for(unsigned int i = 0; i < idCount; i += LOOPBLOCKINGSIZE){
      unsigned int iShifted = (i<<(16-SHIFTSIZE));
      for(int j = i; j < idCount; j += LOOPBLOCKINGSIZE)
         loopIndex[index++] = iShifted | (j>>SHIFTSIZE);
   }
   int firstCount = loopIndexLength / defaultMaxCoreCount;
   int onepercent = firstCount? (firstCount+99) / 100: 1;
   int thread = 0;
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(HetHomCount, notMissingHetCount, HomHomCount, \
      id1, id2, m1, m2, m3, thread, word, word1, word2)
{
   thread = omp_get_thread_num();
   #pragma omp for
#endif
   for(int k = 0; k < loopIndexLength; k++){
      if(k < firstCount && (k % onepercent) == 0) {
         printf("%d%%\r", k/onepercent);
         fflush(stdout);
      }
      int i = (loopIndex[k]>>16)<<SHIFTSIZE;
      int iMax = (i > idCount - LOOPBLOCKINGSIZE? idCount: i + LOOPBLOCKINGSIZE);
      int j = (loopIndex[k]&0xFFFF)<<SHIFTSIZE;
      int jMax = (j > idCount - LOOPBLOCKINGSIZE? idCount: j + LOOPBLOCKINGSIZE);
      for(id1 = i; id1 < iMax; id1++){
         for(id2 = j; id2 < jMax; id2++){ // Computationally intensive here
            // Stage 1: all pairs at few SNPs
            // Stage 1.1: Quick and dirty, almost all computation goes here
            notMissingHetCount = hetInOnePersonCount[id1] + hetInOnePersonCount[id2];  // Upper Bound
            if(notMissingHetCount>100){   // informative enough to make inference
               double cutoff = notMissingHetCount * 0.1428571; // not too many HetHom
               m2 = missingInOnePersonCount[id1][0]+missingInOnePersonCount[id2][0];
               word = SLG[0][id1][0]^SLG[0][id2][0];
               word = word - ((word>>1)&0x5555555555555555);
               word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
               word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
               word = (word+(word>>8)) & 0x00FF00FF00FF00FF;
               word = (word+(word>>16)) & 0x0000FFFF0000FFFF;
               HetHomCount = ((word+(word>>32)) & 0xFFFFFFFF) - m2;
               if(HetHomCount*6 >= cutoff) continue; // C_het <= 0.75

               for(word2 = 0, m1 = 1; m1 < 4; m1++){ // Lower bound of HetHom Count
                  word = SLG[0][id1][m1]^SLG[0][id2][m1];
                  word = word - ((word>>1)&0x5555555555555555);
                  word2 += (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
               }
               word2 = (word2 & 0x0F0F0F0F0F0F0F0F) + ((word2>>4) & 0x0F0F0F0F0F0F0F0F);
               word2 = (word2+(word2>>8)) & 0x00FF00FF00FF00FF;
               word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;
               HetHomCount += ((word2+(word2>>32)) & 0xFFFFFFFF) + m2 -
                  missingInOnePersonCount[id1][1] - missingInOnePersonCount[id2][1];
               if(HetHomCount >= cutoff*0.5) continue; // C_het <= 0.75
               for(word2 = 0, m1 = 4; m1 < long_prescan; m1++){ // HetHom Count
                  word = SLG[0][id1][m1]^SLG[0][id2][m1];
                  word = word - ((word>>1)&0x5555555555555555);
                  word2 += (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
               }// Rare chance of hitting LB, when all 4 words are 1111
               word2 = (word2 & 0x0F0F0F0F0F0F0F0F) + ((word2>>4) & 0x0F0F0F0F0F0F0F0F);
               word2 = (word2+(word2>>8)) & 0x00FF00FF00FF00FF;
               word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;
               HetHomCount += (word2+(word2>>32)) & 0xFFFFFFFF;
               if(HetHomCount >= cutoff) continue; // C_het <= 0.75
               else if(HetHomCount < cutoff*0.125){  // C_het >= 0.97
                  if(j > i || id2 > id1){
                     rpList[thread].Push(id1);
                     rpList[thread].Push(id2);
                     continue;
                  }
               }                  

               // Stage 1.2: Standard heterozygote concordance
               for(m1 = m2 = m3 = 0; m1 < long_prescan; m1++){
                  word1 = ~(SLG[0][id1][m1] | SLG[1][id1][m1]);
                  word2 = ~(SLG[0][id2][m1] | SLG[1][id2][m1]);
                  for(word = word1 & word2; word; word &= (word-1), m2++);
                  for(word = (~SLG[0][id1][m1]) & SLG[1][id1][m1] & word2; word; word &= (word-1), m3++);
                  for(word = word1 & (~SLG[0][id2][m1]) & SLG[1][id2][m1]; word; word &= (word-1), m3++);
               }
               HetHomCount += (m2<<1) + m3;
               notMissingHetCount -= m3;
               if(HetHomCount >= notMissingHetCount * 0.1428571) continue; // C_het <= 0.75
               HomHomCount = prescanMarkerCount - (notMissingHetCount+HetHomCount)/2
               - missingInOnePersonCount[id1][1] - missingInOnePersonCount[id2][1] + m2;
            }else{
               for(word2 = m1 = 0; m1 < long_prescan; m1++){
                  word = SLG[0][id1][m1] & SLG[0][id2][m1]; // HomHom
                  word = word - ((word>>1)&0x5555555555555555);
                  word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
                  word2 += (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
               }
               word2 = (word2+(word2>>8)) & 0x00FF00FF00FF00FF;
               word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;
               HomHomCount = (word2+(word2>>32)) & 0xFFFFFFFF;
            }

            // Stage 1.3: Standard homozygote concordance
            if(HomHomCount > 100){
               for(m3 = m1 = 0; m1 < long_prescan; m1++)   // IBS0
                  for(word = SLG[0][id1][m1] & SLG[0][id2][m1] & (SLG[1][id1][m1] ^ SLG[1][id2][m1]);
                     word; word &= (word-1), m3++);
               if(m3 > HomHomCount * 0.1) continue;   // C_hom < 0.9
            }
            if(i==j && id2 <= id1) continue;
            rpList[thread].Push(id1);
            rpList[thread].Push(id2);
         }  // end of id2 (from j to j+99)
      }  // end of id1 (from i to i+99)
   }  // end of (i, j) block
#ifdef _OPENMP
}
#endif
   delete []loopIndex;
   for(int i = 0; i < idCount; i++)
      delete []missingInOnePersonCount[i];
   delete []missingInOnePersonCount;
   delete []hetInOnePersonCount;
   long long int count = 0;
   for(int i = 0; i < defaultMaxCoreCount; i++)
      count += rpList[i].Length();
   return (count>>1);
}

int Engine::ScreenDuplicates(IntArray & rpList)
{
   if(shortCount < 32 || !bigdataFlag) return -1;  // screening is not needed
   int short_prescan = 32; // 512 SNPs
   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   int prescanMarkerCount = short_prescan==shortCount? markerCount: short_prescan*16;
   int **missingInOnePersonCount = new int * [idCount];
   for(int i = 0; i < idCount; i++)
      missingInOnePersonCount[i] = new int [2];
   int *hetInOnePersonCount;
   hetInOnePersonCount = new int[idCount];
   for(int i = 0; i < idCount; i++)
      missingInOnePersonCount[i][0] = missingInOnePersonCount[i][1] = hetInOnePersonCount[i] = 0;
   unsigned short int *masks = new unsigned short int[short_prescan];
   for(int i = 0; i < short_prescan;  i++)
      masks[i] = 0xFFFF;
   if(short_prescan == shortCount)
      masks[short_prescan-1] = markerCount%16? ((1 << (markerCount % 16))-1): 0xFFFF;
   unsigned short int **missingWordIndex = new unsigned short int *[idCount];
   for(int i = 0; i < idCount; i++){
      missingWordIndex[i] = new unsigned short int [2];
      missingWordIndex[i][0] = missingWordIndex[i][1] = 0;
   }
   unsigned short int word, word2;
   int m1, m2, m3, sum;
   for(int i = 0; i < idCount; i++){
      for(sum = short_prescan*2, m1 = 0; m1 < short_prescan/8; m1++)
         sum -= oneoneCount[SG[0][i][m1] | SG[1][i][m1]];
      missingInOnePersonCount[i][0] = sum;
      for(sum = m1 = 0; m1 < short_prescan; m1++){
         word = (~SG[0][i][m1]) & (~SG[1][i][m1]) & masks[m1];
         if(word){   // some missing
            missingWordIndex[i][m1/16] |= (1<<(m1%16));
            sum += oneoneCount[word];
         }
      }
      missingInOnePersonCount[i][1] = sum;
      for(sum = m1 = 0; m1 < short_prescan; m1++)
         sum += oneoneCount[(~SG[0][i][m1]) & SG[1][i][m1]];
      hetInOnePersonCount[i] = sum;
   }
   delete []masks;
   char revbase[65536];
   for(int i = 0; i < 16; i++)
      revbase[shortbase[i]] = i;
   char rightmost[65536];
   for(int i = 0; i < 65536; i++)
      rightmost[i] = revbase[i&(-i)];
   int id1, id2, HetHomCount, DiffHomCount, HomHomCount, notMissingHetCount;
   const int cutoffMissingCount = short_prescan*16-100;//MissingCount < 412
   IntArray *ibuffer = new IntArray[defaultMaxCoreCount];
   for(int i = 0; i < defaultMaxCoreCount; i++)
      ibuffer[i].Dimension(0);
   const char LOOPBLOCKINGSIZE=64;
   const char SHIFTSIZE=6;
   unsigned int loopIndexLength = (unsigned int)((idCount-1)/LOOPBLOCKINGSIZE+1)*(((idCount-1)/LOOPBLOCKINGSIZE+1)+1)/2;
   unsigned int *loopIndex = new unsigned int [loopIndexLength];
   unsigned int index = 0;
   for(unsigned int i = 0; i < idCount; i += LOOPBLOCKINGSIZE){
      unsigned int iShifted = (i<<(16-SHIFTSIZE));
      for(int j = i; j < idCount; j += LOOPBLOCKINGSIZE)
         loopIndex[index++] = iShifted | (j>>SHIFTSIZE);
   }
   int thread = 0;
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
   private(HetHomCount, notMissingHetCount, HomHomCount, DiffHomCount, \
   id1, id2, m1, m2, m3, thread, word, word2)
{
   thread = omp_get_thread_num();
   #pragma omp for
#endif
   for(int k = 0; k < loopIndexLength; k++){
      int i = (loopIndex[k]>>16)<<SHIFTSIZE;
      int iMax = (i > idCount - LOOPBLOCKINGSIZE? idCount: i + LOOPBLOCKINGSIZE);
      int j = (loopIndex[k]&0xFFFF)<<SHIFTSIZE;
      int jMax = (j > idCount - LOOPBLOCKINGSIZE? idCount: j + LOOPBLOCKINGSIZE);
      int jMin = j;
      for(id1 = i; id1 < iMax; id1++){
         if(missingInOnePersonCount[id1][1]>=cutoffMissingCount) continue;
         if(i == j) jMin = id1 + 1; // diagonal blocks only
         for(id2 = jMin; id2 < jMax; id2++){
            if(missingInOnePersonCount[id2][1]>=cutoffMissingCount) continue;
            // Stage 1: all pairs at few SNPs
            // Stage 1.1: Quick and dirty, almost all computation goes here
            notMissingHetCount = hetInOnePersonCount[id1] + hetInOnePersonCount[id2];  // Upper Bound
            if(notMissingHetCount>100){   // informative enough to make inference
               double cutoff = notMissingHetCount * 0.1428571; // not too many HetHom
               m2 = missingInOnePersonCount[id1][0]+missingInOnePersonCount[id2][0];
               for(HetHomCount = -m2, m1 = 0; m1 < short_prescan/8; m1++) // Lower bound of HetHom Count
                  HetHomCount += oneoneCount[SG[0][id1][m1]^SG[0][id2][m1]];
               if(HetHomCount*6 >= cutoff) continue; // C_het <= 0.75
               for(HetHomCount += m2-missingInOnePersonCount[id1][1]-missingInOnePersonCount[id2][1],
               m1 = short_prescan/8; m1 < short_prescan/2; m1++) // Lower bound of HetHom Count
                  HetHomCount += oneoneCount[SG[0][id1][m1]^SG[0][id2][m1]];
               if(HetHomCount*2 >= cutoff) continue; // C_het <= 0.75
               for(m1 = short_prescan/2; (HetHomCount < cutoff) && (m1 < short_prescan); m1++) // HetHom Count
                  HetHomCount += oneoneCount[SG[0][id1][m1]^SG[0][id2][m1]];
               if(HetHomCount >= cutoff) continue; // C_het <= 0.75

               // Stage 1.2: Standard heterozygote concordance
               for(HomHomCount = m3 = m1 = 0; m1 < 2; m1++){
                  for(word = missingWordIndex[id1][m1] & missingWordIndex[id2][m1];
                     word; word &= (word-1)){   // MissMiss
                     m2 = m1*16+rightmost[word];   // m2 is the word with nonmissingness
                     for(word2 = ~(SG[0][id1][m2] | SG[0][id2][m2] | SG[1][id1][m2] | SG[1][id2][m2])&0xFFFF;
                        word2; word2 &= (word2-1))
                        HomHomCount ++;
                  }
                  for(word = missingWordIndex[id2][m1];  // HetMiss
                     word; word &= (word-1)){
                     m2 = m1*16+rightmost[word];
                     for(word2 = ~(SG[0][id1][m2] | SG[0][id2][m2] | SG[1][id2][m2]) & SG[1][id1][m2];
                         word2; word2 &= (word2-1))
                        m3 ++;
                  }
                  for(word = missingWordIndex[id1][m1];  // MissHet
                     word; word &= (word-1)){
                     m2 = m1*16+rightmost[word];
                     for(word2 = ~(SG[0][id1][m2] | SG[0][id2][m2] | SG[1][id1][m2]) & SG[1][id2][m2];
                         word2; word2 &= (word2-1))
                        m3 ++;
                  }
               }
               HetHomCount += HomHomCount*2 + m3;
               notMissingHetCount -= m3;
               if(HetHomCount >= notMissingHetCount * 0.1428571) continue; // C_het <= 0.75
               HomHomCount = prescanMarkerCount - (notMissingHetCount+HetHomCount)/2
               - missingInOnePersonCount[id1][1] - missingInOnePersonCount[id2][1] + HomHomCount;
            }else
               for(HomHomCount = m1 = 0; m1 < short_prescan; m1++)
                  HomHomCount += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1]];
            // Stage 1.3: Standard homozygote concordance
            if(HomHomCount > 100){
               for(DiffHomCount = m1 = 0; m1 < short_prescan; m1++)
                  DiffHomCount += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1] & (SG[1][id1][m1] ^ SG[1][id2][m1])];
               if(DiffHomCount > HomHomCount * 0.1) continue;   // C_hom < 0.9
            }
            ibuffer[thread].Push(id1);
            ibuffer[thread].Push(id2);
         }  // end of id2 (from j to j+99)
      }  // end of id1 (from i to i+99)
   }  // end of (i, j) block
#ifdef _OPENMP
}
#endif
   delete []loopIndex;
   for(int i = 0; i < idCount; i++){
      delete []missingWordIndex[i];
      delete []missingInOnePersonCount[i];
   }
   delete []missingWordIndex;
   delete []missingInOnePersonCount;
   delete []hetInOnePersonCount;
   rpList.Dimension(0);
   for(int i = 0; i < defaultMaxCoreCount; i++)
      rpList.Append(ibuffer[i]);
   delete []ibuffer;
   return rpList.Length()/2;
}

void Engine::ComputeDistanceMatrix(Matrix & Dist)
{
   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   int IBS0Count, notMissingCount, HetHomCount, id1, id2;
   int m1, m2, m3;
   IntArray tArray;
   unsigned int **missingWordInOnePerson = new unsigned int *[idCount];
   int *missingWordInOnePersonCount = new int[idCount];
   int *missingInOnePersonCount = new int[idCount];
   for(int i = 0; i < idCount; i++)
      missingWordInOnePersonCount[i] = missingInOnePersonCount[i] = 0;
   unsigned short int *masks = new unsigned short int[shortCount];
   for(int i = 0; i < shortCount;  i++)
      masks[i] = 0xFFFF;
   masks[shortCount-1] = markerCount%16? ((1 << (markerCount % 16))-1): 0xFFFF;
   int sum, sum2;
   for(int i = 0; i < idCount; i++){
      tArray.Dimension(0);
      sum = sum2 = 0;
      for(int m = 0; m < shortCount; m++)
         if(((~GG[0][i][m]) & (~GG[1][i][m]) & masks[m])!=0){
            sum ++;
            sum2 += oneoneCount[(~GG[0][i][m]) & (~GG[1][i][m]) & masks[m]];
            tArray.Push(m);
          }
      missingWordInOnePersonCount[i] = sum;
      missingInOnePersonCount[i] = sum2;
      if(tArray.Length()){
         missingWordInOnePerson[i] = new unsigned int [tArray.Length()];
         for(int m = 0; m < tArray.Length(); m++)
            missingWordInOnePerson[i][m] = tArray[m];
      }else
         missingWordInOnePerson[i] = NULL;
   }
   delete []masks;
   IntArray outindex(idCount);
   outindex.Set(-1);
   int N = 0;
   const int cutoffMissingCount = markerCount-MINSNPCOUNT;  // almost all missing
   for(int i = 0; i < idCount; i++)
      if(missingInOnePersonCount[i] < cutoffMissingCount)
         outindex[i] = N++;
   Dist.Dimension(N,N);
   Dist.Zero();

//   const double distancePar = (allflags&(1<<IBSFLAG))? 2.0: 4.0;

#ifdef _OPENMP
   printf("%d CPU cores are used.\n", defaultMaxCoreCount);
#endif

const int BLOCKSIZE=16;
const int CACHESIZE=512;
IntArray loopIndex[2];
loopIndex[0].Dimension(0); loopIndex[1].Dimension(0);
for(int i = 0; i < idCount; i += BLOCKSIZE)
   for(int j = i; j < idCount; j += BLOCKSIZE){
      loopIndex[0].Push(i);
      loopIndex[1].Push(j);
   }
int loopIndexLength = loopIndex[0].Length();
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(HetHomCount, IBS0Count, notMissingCount, id1, id2, m1, m2, m3)
{
#endif
   int missingCount[2][BLOCKSIZE];
   int mmMin[2][BLOCKSIZE], mmMax[2][BLOCKSIZE];
   double localD[2][BLOCKSIZE][BLOCKSIZE];
#ifdef _OPENMP
   #pragma omp for
#endif
   for(int k = 0; k < loopIndexLength; k++){
      int i = loopIndex[0][k];
      int iMax = i<idCount-BLOCKSIZE? i+BLOCKSIZE: idCount;
      int j = loopIndex[1][k];
      int jMax = j<idCount-BLOCKSIZE? j+BLOCKSIZE: idCount;
      int jMin = j;
      for(int r = 0; r < 2; r++)
         for(int c = 0; c < BLOCKSIZE; c++)
            mmMin[r][c] = 0;
      for(int c = 0; c < BLOCKSIZE; c++){
         missingCount[0][c] = missingWordInOnePersonCount[i+c];
         missingCount[1][c] = missingWordInOnePersonCount[j+c];
      }
      for(int c = 0; c < 2; c++)
         for(int ii = 0; ii < BLOCKSIZE; ii++)
            for(int jj = 0; jj < BLOCKSIZE; jj++)
               localD[c][ii][jj] = 0.0;
      for(int mm = 0; mm < shortCount; mm += CACHESIZE){
         int mMax = (mm > shortCount-CACHESIZE) ? shortCount: mm+CACHESIZE;
         for(id1 = i; id1 < iMax; id1++){
            mmMin[0][id1-i] = (mm == 0? 0: mmMax[0][id1-i]);
            for(m1 = mmMin[0][id1-i]; m1 < missingCount[0][id1-i] && missingWordInOnePerson[id1][m1] < mMax; m1++);
            mmMax[0][id1-i] = m1;
         }
         for(id2 = j; id2 < jMax; id2++){
            mmMin[1][id2-j] = (mm == 0? 0: mmMax[1][id2-j]);
            for(m1 = mmMin[1][id2-j]; m1 < missingCount[1][id2-j] && missingWordInOnePerson[id2][m1] < mMax; m1++);
            mmMax[1][id2-j] = m1;
         }
         for(id1 = i; id1 < iMax; id1++){
            if(missingInOnePersonCount[id1] >= cutoffMissingCount) continue;
            char ii = id1 - i;
            if(i == j) jMin = id1 + 1;
            for(id2 = jMin; id2 < jMax; id2++){
               if(missingInOnePersonCount[id2] >= cutoffMissingCount) continue;
               char jj = id2 - j;
               for(IBS0Count=HetHomCount=0, m1 = mm; m1 < mMax; m1++){
                  HetHomCount += oneoneCount[GG[0][id1][m1]^GG[0][id2][m1]];
                  IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
               } // Lower bound of HetHom Count, minus HetMiss and MissMiss
               for(notMissingCount = 0, m1 = mmMin[0][ii], m2 = mmMin[1][jj]; m1 < mmMax[0][ii]; m1++){
                  m3 = missingWordInOnePerson[id1][m1]; // m for id1
                  for(; m2 < mmMax[1][jj] && missingWordInOnePerson[id2][m2] < m3; m2++);
                  if(m2 == mmMax[1][jj]) break;
                  if(m3 == missingWordInOnePerson[id2][m2])   // id1 and id2 are both missing
                     notMissingCount += oneoneCount[~(GG[0][id1][m3] | GG[1][id1][m3] | GG[0][id2][m3] | GG[1][id2][m3]) & 0xFFFF];
               }  // MissMissCount
               for(m1 = mmMin[1][jj]; m1 < mmMax[1][jj]; m1++){
                  m2 = missingWordInOnePerson[id2][m1];
                  HetHomCount += oneoneCount[~(GG[0][id1][m2] | GG[0][id2][m2] | GG[1][id2][m2]) & GG[1][id1][m2]];
               }  // Add HetMiss
               for(m1 = mmMin[0][ii]; m1 < mmMax[0][ii]; m1++){
                  m2 = missingWordInOnePerson[id1][m1];
                  HetHomCount += oneoneCount[~(GG[0][id1][m2] | GG[0][id2][m2] | GG[1][id1][m2]) & GG[1][id2][m2]];
               }  // Add MissHet
               localD[0][ii][jj] += HetHomCount + (notMissingCount<<1) + (IBS0Count<<2);
               localD[1][ii][jj] += notMissingCount;
            }  // end of id2
         }  // end of id1
      }// end of mm
      for(id1 = i; id1 < iMax; id1++){
         if(missingInOnePersonCount[id1]>=cutoffMissingCount) continue;
         char ii = id1 - i;
         if(i == j) jMin = id1 + 1;
         for(id2 = jMin; id2 < jMax; id2++){
            char jj = id2 - j;
            if(missingInOnePersonCount[id2]>=cutoffMissingCount) continue;
            localD[0][ii][jj] -= (missingInOnePersonCount[id1] + missingInOnePersonCount[id2]);
            localD[1][ii][jj] += markerCount - missingInOnePersonCount[id1] - missingInOnePersonCount[id2];
            localD[0][ii][jj] /= localD[1][ii][jj];
            Dist[outindex[id2]][outindex[id1]] = Dist[outindex[id1]][outindex[id2]] = localD[0][ii][jj];
         }
      }
   }
#ifdef _OPENMP
}  // extra bracket for omp
#endif
   for(int i = 0; i < idCount; i++)
      if(missingWordInOnePerson[i]) delete []missingWordInOnePerson[i];
   if(missingWordInOnePersonCount) delete []missingWordInOnePersonCount;
   if(missingInOnePersonCount) delete []missingInOnePersonCount;
}

void Engine::ComputeBigDataDuplicate()
{
   if(shortCount == 0) {
      printf("No genotype data.");
      return;
   }
   printf("Autosome genotypes stored in %d", shortCount);
   printf(" words for each of %d individuals.\n", idCount);

   printf("\nOptions in effect:\n");
   printf("\t--duplicate\n");
   if(bigdataFlag){
      if(faster)
         printf("\t--faster %d\n", faster);
      if(slower)
         printf("\t--slower %d\n", slower);
   }
   if(CoreCount)
      printf("\t--cpus %d\n", CoreCount);
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);

   printf("\n");
#ifdef _OPENMP
   printf("%d CPU cores are used to sort autosomes...\n", defaultMaxCoreCount);
#endif
   if(SG[0]==NULL)
   //*****************Make sorted genotype data****************************
   ConvertGGtoSG(GG, markerCount, SG, (32 < shortCount)? 512: markerCount);
   //*****************Make sorted genotype data****************************

   printf("Computing pairwise genotype concordance starts at %s", currentTime());
#ifdef _OPENMP
   printf("  %d CPU cores are used.\n", defaultMaxCoreCount);
#endif
   IntArray allpairs0;

   //*******************Screen Duplicates*******************************
   int rawdupCount = ScreenDuplicates(allpairs0);
   //*******************Screen Duplicates*******************************

   if(rawdupCount==0){
      printf("                                          ends at %s", currentTime());
      printf("No duplicates are found with heterozygote concordance rate > 80%%.\n\n");
      return;
   }
   IntArray allpairs;
   if(rawdupCount > 0){
      printf("        Stage 1 (with %d SNPs) screening ends at %s", 512, currentTime());
      Vector pM(rawdupCount);
      for(int i = 0; i < rawdupCount; i++){
         pM[i] = allpairs0[i*2];
         pM[i] *= idCount;
         pM[i] += allpairs0[i*2+1];
      }
      QuickIndex index;
      index.Index(pM);
      allpairs.Dimension(rawdupCount*2);
      for(int i = 0; i < rawdupCount; i++){
         int k = index[i];
         allpairs[i*2] = allpairs0[k*2];
         allpairs[i*2+1] = allpairs0[k*2+1];
      }
   }else{ // rawdupCount == -1
      allpairs.Dimension(0);
      for(int i = 0; i < idCount; i++)
         for(int j = i+1; j < idCount; j++){
            allpairs.Push(i);
            allpairs.Push(j);
         }
      rawdupCount = allpairs.Length() / 2;
   }
   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   char buffer[1024];
   StringArray buffers;
   buffers.Dimension(defaultMaxCoreCount);
   for(int i = 0; i < buffers.Length(); i++)
      buffers[i].Clear();
   int rpCount = 0;
   int thread = 0;
   int id1, id2, HetHetCount, notMissingHetCount, notMissingCount, HomHomCount, DiffHomCount, SameHomCount, m1;
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
         private(HetHetCount, notMissingHetCount, notMissingCount, \
         HomHomCount, DiffHomCount, SameHomCount, id1, id2, buffer, m1, thread)\
         reduction(+:rpCount)
{
   thread = omp_get_thread_num();
   #pragma omp for
#endif
      for(int i = 0; i < rawdupCount; i++){
         id1 = allpairs[i*2];
         id2 = allpairs[i*2+1];
         for(HetHetCount=notMissingHetCount=m1=0; m1 < shortCount; m1++){
            HetHetCount += oneoneCount[~(GG[0][id1][m1] | GG[0][id2][m1]) & (GG[1][id1][m1]) & GG[1][id2][m1]];
            notMissingHetCount += oneoneCount[(GG[0][id1][m1] | GG[1][id1][m1]) &
            (GG[0][id2][m1] | GG[1][id2][m1]) & ~(GG[0][id1][m1] & GG[0][id2][m1])];
         }
         if(HetHetCount <= notMissingHetCount*0.8) continue; // not duplicates, heterozygote concordance rate <= 80%
         for(HomHomCount=DiffHomCount=m1=0; m1 < shortCount; m1++){
            HomHomCount += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1]];
            DiffHomCount += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
         }
         notMissingCount = notMissingHetCount + HomHomCount;
         SameHomCount = HomHomCount - DiffHomCount;
         sprintf(buffer, "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.5lf\t%.5lf\t%.5lf\n",
               (const char*)ped[phenoid[id1]].famid, (const char*)ped[phenoid[id1]].pid,
               (const char*)ped[phenoid[id2]].famid, (const char*)ped[phenoid[id2]].pid,
               notMissingCount, DiffHomCount, notMissingHetCount - HetHetCount,
               HetHetCount+SameHomCount,
               (HetHetCount+SameHomCount)*1.0/notMissingCount,
               1 - DiffHomCount * 1.0 / HomHomCount,
               HetHetCount * 1.0 / notMissingHetCount);
         buffers[thread].Add(buffer);
         rpCount ++;
      }
#ifdef _OPENMP
}
#endif
   printf("        Stage 2 (with all SNPs) inference ends at %s", currentTime());
   if(rpCount){
      String outfile;
      outfile.Copy(prefix);
      outfile.Add(".con");
      FILE *fp = fopen(outfile, "wt");
      fprintf(fp, "FID1\tID1\tFID2\tID2\tN\tN_IBS0\tN_IBS1\tN_IBS2\tConcord\tHomConc\tHetConc\n");
      for(int b = 0; b < buffers.Length(); b++)
         buffers[b].Write(fp);
      fclose(fp);
      printf("%d pairs of duplicates with heterozygote concordance rate > 80%% are saved in file %s\n\n",
         rpCount, (const char*)outfile);
   }else
      printf("No duplicates are found with heterozygote concordance rate > 80%%.\n\n");
   if(rpCount < rawdupCount){
      long long int pMax = idCount*(idCount-1)/2;
      if(rawdupCount < pMax)  // with screening
         printf("  %d additional pairs from screening stage cannot be not confirmed in the final stage\n\n", rawdupCount - rpCount);
   }
}

void Engine::ComputeBigDataDuplicate64Bit()
{
   if(longCount == 0) {
      printf("No genotype data.");
      return;
   }
   printf("Autosome genotypes stored in %d", longCount);
   printf(" words for each of %d individuals.\n", idCount);

   printf("\nOptions in effect:\n");
   printf("\t--duplicate\n");
   if(CoreCount)
      printf("\t--cpus %d\n", CoreCount);
   if(Bit64Flag)
      printf("\t--sysbit 64\n");
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);

   printf("\n");
   printf("Sorting autosomes...\n");
   if(SLG[0]==NULL)
      //*****************Make sorted genotype data****************************
      ConvertLGtoSLG(LG, markerCount, SLG, (longCount > 8)? 512: markerCount);
      //*****************Make sorted genotype data****************************

   printf("Computing pairwise genotype concordance starts at %s", currentTime());
#ifdef _OPENMP
   printf("  %d CPU cores are used...\n", defaultMaxCoreCount);
#endif
   IntArray *allpairs = new IntArray[defaultMaxCoreCount];

   //*******************Screen Duplicates*******************************
   long long int rawdupCount = ScreenDuplicates64Bit(allpairs);
   //*******************Screen Duplicates*******************************

   if(rawdupCount==0){
      printf("                                          ends at %s", currentTime());
      printf("No duplicates are found with heterozygote concordance rate > 80%%.\n\n");
      return;
   }
   if(rawdupCount > 0)
      printf("        Stage 1 (with %d SNPs) screening ends at %s", 512, currentTime());
   else{ // rawdupCount == -1 due to few SNPs
      for(int t = 0; t < defaultMaxCoreCount; t++)
         allpairs[t].Dimension(0);
      rawdupCount = (idCount & 1)? (idCount-1)/2*idCount: idCount/2*(idCount-1);
      long long int blockSize = (rawdupCount-1)/defaultMaxCoreCount+1;
      long long int index = 0;
      for(int i = 0; i < idCount; i++)
         for(int j = i+1; j < idCount; j++){
            int thread = index / blockSize;
            allpairs[thread].Push(i);
            allpairs[thread].Push(j);
            index ++;
         }
   }
   String outfile;
   outfile.Copy(prefix);
   outfile.Add(".con");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "FID1\tID1\tFID2\tID2\tN\tN_IBS0\tN_IBS1\tN_IBS2\tConcord\tHomConc\tHetConc\n");
   int fpCount = 0;  // false pair Count
   IntArray HetHetCounts, DiffHomCounts, HomHomCounts, notMissingHetCounts;
   for(int t = 0; t < defaultMaxCoreCount; t++){
      DuplicateInSubset64Bit(allpairs[t], HetHetCounts, DiffHomCounts, HomHomCounts, notMissingHetCounts);
      int pairCount = allpairs[t].Length()/2;
      for(int p = 0; p < pairCount; p++){
         int HetHetCount = HetHetCounts[p];
         int DiffHomCount = DiffHomCounts[p];
         int HomHomCount = HomHomCounts[p];
         int notMissingHetCount = notMissingHetCounts[p];
         if(HetHetCount <= notMissingHetCount*0.8) {
            fpCount ++; // false positive count adds
            continue; // not duplicates, heterozygote concordance rate <= 80%
         }
         int notMissingCount = notMissingHetCount + HomHomCount;
         int SameHomCount = HomHomCount - DiffHomCount;
         int id1 = phenoid[allpairs[t][p*2]];
         int id2 = phenoid[allpairs[t][p*2+1]];
         fprintf(fp, "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.5lf\t%.5lf\t%.5lf\n",
            (const char*)ped[id1].famid, (const char*)ped[id1].pid,
            (const char*)ped[id2].famid, (const char*)ped[id2].pid,
            notMissingCount, DiffHomCount,   // N_SNP, IBS0Count
            notMissingHetCount - HetHetCount,   // IBS1Count
            HetHetCount+SameHomCount,  // IBS2Count
            (HetHetCount+SameHomCount)*1.0/notMissingCount, // Concordance
            SameHomCount * 1.0 / HomHomCount,  // HomConc
            HetHetCount * 1.0 / notMissingHetCount);  // HetConc
         }
   }
   fclose(fp);
   printf("        Stage 2 (with all SNPs) inference ends at %s", currentTime());
   if(rawdupCount > fpCount)
      printf("%lli pairs of duplicates with heterozygote concordance rate > 80%% are saved in file %s\n\n",
         rawdupCount - fpCount, (const char*)outfile);
   else
      printf("No duplicates are found with heterozygote concordance rate > 80%%.\n\n");
   if(fpCount && rawdupCount>0)
      printf("  %d additional pairs from screening stage not confirmed in the final stage\n\n",
         fpCount);
   delete []allpairs;
}


/*
   char buffer[1024];
   StringArray buffers(defaultMaxCoreCount);
   for(int i = 0; i < defaultMaxCoreCount; i++)
      buffers[i].Clear();
   int fpCount = 0;  // false pair Count
   int HetHetCount, DiffHomCount, m1;
   unsigned long long int word, word1, word2;
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
   private(HetHetCount, DiffHomCount, buffer, m1, word, word1, word2) reduction(+:fpCount)
#endif
   for(int t = 0; t < defaultMaxCoreCount; t++){
      int pairCount = allpairs[t].Length()/2;
      for(int i = 0; i < pairCount; i++){
         int id1 = allpairs[t][i*2];
         int id2 = allpairs[t][i*2+1];
         for(word1=word2=HetHetCount=DiffHomCount=m1=0; m1 < longCount; m1++){
            word = (LG[0][id1][m1] | LG[1][id1][m1]) &   // HetNomiss
            (LG[0][id2][m1] | LG[1][id2][m1]) & ~(LG[0][id1][m1] & LG[0][id2][m1]);
            word = word - ((word>>1)&0x5555555555555555);
            word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
            word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
            word = (word+(word>>8)) & 0x00FF00FF00FF00FF;
            word1 += (word+(word>>16)) & 0x0000FFFF0000FFFF;
            word = LG[0][id1][m1] & LG[0][id2][m1];   // HomHom
            word = word - ((word>>1)&0x5555555555555555);
            word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
            word = (word + (word>>4)) & 0x0F0F0F0F0F0F0F0F;
            word = (word + (word>>8)) & 0x00FF00FF00FF00FF;
            word2 += (word + (word>>16)) & 0x0000FFFF0000FFFF;
            for(word = (~LG[0][id1][m1] & LG[1][id1][m1] & LG[0][id2][m1]) |  // HetHom
               (~LG[0][id2][m1] & LG[1][id2][m1] & LG[0][id1][m1]);  // HomHet
               word; word &= (word-1), HetHetCount++);
            for(word = LG[0][id1][m1] & LG[0][id2][m1] & (LG[1][id1][m1] ^ LG[1][id2][m1]);
               word; word &= (word-1), DiffHomCount++);  // IBS0
         }
         int notMissingHetCount = (word1+(word1>>32)) & 0xFFFFFFFF;
         int HomHomCount = (word2 + (word2>>32)) & 0xFFFFFFFF;
         HetHetCount = notMissingHetCount - HetHetCount;
         if(HetHetCount <= notMissingHetCount*0.8) {
            fpCount ++; // false positive count adds
            continue; // not duplicates, heterozygote concordance rate <= 80%
         }
         int notMissingCount = notMissingHetCount + HomHomCount;
         int SameHomCount = HomHomCount - DiffHomCount;
         sprintf(buffer, "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.5lf\t%.5lf\t%.5lf\n",
            (const char*)ped[phenoid[id1]].famid, (const char*)ped[phenoid[id1]].pid,
            (const char*)ped[phenoid[id2]].famid, (const char*)ped[phenoid[id2]].pid,
            notMissingCount, DiffHomCount,   // N_SNP, IBS0Count
            notMissingHetCount - HetHetCount,   // IBS1Count
            HetHetCount+SameHomCount,  // IBS2Count
            (HetHetCount+SameHomCount)*1.0/notMissingCount, // Concordance
            SameHomCount * 1.0 / HomHomCount,  // HomConc
            HetHetCount * 1.0 / notMissingHetCount);  // HetConc
         buffers[t].Add(buffer);
      }
   }  // openMP with thread loop ends
   delete []allpairs;
*/


