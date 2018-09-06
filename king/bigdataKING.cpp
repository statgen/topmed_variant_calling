//////////////////////////////////////////////////////////////////////
// bigdataKING.cpp
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
// March 22, 2018

#include <math.h>
#include "analysis.h"
#include "Kinship.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

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

int Engine::ScreenDuplicates64Bit(IntArray & rpList)
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
   IntArray *ibuffer = new IntArray[defaultMaxCoreCount];
   for(int i = 0; i < defaultMaxCoreCount; i++)
      ibuffer[i].Dimension(0);
   const short int LOOPBLOCKINGSIZE=256;  // use 64KB cache
   const char SHIFTSIZE=8; // 2^8 * 2 * 2^9 / 2^2 = 2^16 = 64kB
   unsigned int loopIndexLength = (unsigned int)((idCount-1)/LOOPBLOCKINGSIZE+1)*(((idCount-1)/LOOPBLOCKINGSIZE+1)+1)/2;
   unsigned short int **loopIndex = new unsigned short int * [loopIndexLength];
   for(int k = 0; k < loopIndexLength; k++)
      loopIndex[k] = new unsigned short int [2];
   unsigned int index = 0;
   for(int i = 0; i < idCount; i += LOOPBLOCKINGSIZE)
      for(int j = i; j < idCount; j += LOOPBLOCKINGSIZE){
         loopIndex[index][0] = (i>>SHIFTSIZE);
         loopIndex[index][1] = (j>>SHIFTSIZE);
         index ++;
      }
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
      int i = int(loopIndex[k][0])<<SHIFTSIZE;
      int iMax = (i > idCount - LOOPBLOCKINGSIZE? idCount: i + LOOPBLOCKINGSIZE);
      int j = int(loopIndex[k][1])<<SHIFTSIZE;
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
                     ibuffer[thread].Push(id1);
                     ibuffer[thread].Push(id2);
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
            // if(i == j) jMin = id1 + 1; // diagonal blocks only
            if(i==j && id2 <= id1) continue;
            ibuffer[thread].Push(id1);
            ibuffer[thread].Push(id2);
         }  // end of id2 (from j to j+99)
      }  // end of id1 (from i to i+99)
   }  // end of (i, j) block
#ifdef _OPENMP
}
#endif
   for(int k = 0; k < loopIndexLength; k++)
      delete []loopIndex[k];
   delete []loopIndex;
   for(int i = 0; i < idCount; i++)
      delete []missingInOnePersonCount[i];
   delete []missingInOnePersonCount;
   delete []hetInOnePersonCount;
   rpList.Dimension(0);
   for(int i = 0; i < defaultMaxCoreCount; i++)
      rpList.Append(ibuffer[i]);
   delete []ibuffer;
   return rpList.Length()/2;
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
   unsigned short int **loopIndex = new unsigned short int * [loopIndexLength];
   for(int k = 0; k < loopIndexLength; k++)
      loopIndex[k] = new unsigned short int [2];
   unsigned int index = 0;
   for(int i = 0; i < idCount; i += LOOPBLOCKINGSIZE)
      for(int j = i; j < idCount; j += LOOPBLOCKINGSIZE){
         loopIndex[index][0] = (i>>SHIFTSIZE);
         loopIndex[index][1] = (j>>SHIFTSIZE);
         index ++;
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
      int i = int(loopIndex[k][0])<<SHIFTSIZE;
      int iMax = (i > idCount - LOOPBLOCKINGSIZE? idCount: i + LOOPBLOCKINGSIZE);
      int j = int(loopIndex[k][1])<<SHIFTSIZE;
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
   for(int k = 0; k < loopIndexLength; k++)
      delete []loopIndex[k];
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
   printf("  %d CPU cores are used.\n", defaultMaxCoreCount);
#endif
   IntArray allpairs0;

   //*******************Screen Duplicates*******************************
   int rawdupCount = ScreenDuplicates64Bit(allpairs0);
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

   char buffer[1024];
   StringArray buffers;
   buffers.Dimension(defaultMaxCoreCount);
   for(int i = 0; i < buffers.Length(); i++)
      buffers[i].Clear();
   int rpCount = 0;
   int thread = 0;
   int id1, id2, HetHetCount, notMissingHetCount, notMissingCount, HomHomCount, DiffHomCount, SameHomCount, m1;
   unsigned long long int word, word1, word2;
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
   private(HetHetCount, notMissingHetCount, notMissingCount, HomHomCount, \
   DiffHomCount, SameHomCount, id1, id2, buffer, m1, thread, word, word1, word2)\
   reduction(+:rpCount)
{
   thread = omp_get_thread_num();
   #pragma omp for
#endif
      for(int i = 0; i < rawdupCount; i++){
         id1 = allpairs[i*2];
         id2 = allpairs[i*2+1];
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
         notMissingHetCount = (word1+(word1>>32)) & 0xFFFFFFFF;
         HomHomCount = (word2 + (word2>>32)) & 0xFFFFFFFF;
         HetHetCount = notMissingHetCount - HetHetCount;
         if(HetHetCount <= notMissingHetCount*0.8) continue; // not duplicates, heterozygote concordance rate <= 80%
         notMissingCount = notMissingHetCount + HomHomCount;
         SameHomCount = HomHomCount - DiffHomCount;
         sprintf(buffer, "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.5lf\t%.5lf\t%.5lf\n",
               (const char*)ped[phenoid[id1]].famid, (const char*)ped[phenoid[id1]].pid,
               (const char*)ped[phenoid[id2]].famid, (const char*)ped[phenoid[id2]].pid,
               notMissingCount, DiffHomCount,   // N_SNP, IBS0Count
               notMissingHetCount - HetHetCount,   // IBS1Count
               HetHetCount+SameHomCount,  // IBS2Count
               (HetHetCount+SameHomCount)*1.0/notMissingCount, // Concordance
               SameHomCount * 1.0 / HomHomCount,  // HomConc
               HetHetCount * 1.0 / notMissingHetCount);  // HetConc
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
// to be updated, for --noibdseg only
void Engine::ComputeBigDataKinship()
{
   if(shortCount==0) error("No genotype data");
   printf("Autosome genotypes stored in %d", shortCount);
   printf(" words for each of %d individuals.\n", idCount);

   printf("\nOptions in effect:\n");
   if(errorrateCutoff!=_NAN_)
      printf("\t--errorrate %G\n", errorrateCutoff);
   printf("\t--related\n");
   if(relativedegree>1)
      printf("\t--degree %d\n", relativedegree);
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

   int prescanCount[10];
   for(int i = 0; i < 10; i++)
      prescanCount[i] = 256 * (1<<(i*3));
   int short_prescan = relativedegree < 10? prescanCount[relativedegree-1]: prescanCount[9];
   if(short_prescan > shortCount || (!bigdataFlag)) short_prescan=shortCount;

   if(bigdataFlag){
      if(faster > 1){
         short_prescan /= faster;
         if(short_prescan < 10) short_prescan = 10;
      }
      if(slower > 1){
         short_prescan *= slower;
         if(short_prescan > shortCount) short_prescan = shortCount;
      }
   }
   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   String outfile;
   outfile.Copy(prefix);
   outfile.Add(".kin");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "FID\tID1\tID2\tN_SNP\tZ0\tPhi\tHetHet\tIBS0\tHetConc\tHomConc\tKinship\tError\n");
   double kinship, inflation, smaller, larger;
   int id1, id2;
   Kinship kin;
   int het1, het2, homhom, unequal;
   int HetHetCount, IBS0Count, het1Count, het2Count, notMissingCount;
   double phi, pi0, errorFlag;
   int beforeCount[6], afterCount[6];
   Vector IBS0PO(0), IBS0FS(0), IBS0L1(0);
   int degree; double ibs0;
   for(int i = 0; i < 6; i++) beforeCount[i] = afterCount[i] = 0;
   for(int f = 0; f < ped.familyCount; f++){
      if(id[f].Length()<2) continue;
      kin.Setup(*ped.families[f]);
      for(int i = 0; i < id[f].Length(); i++)
         for(int j = i+1; j < id[f].Length(); j++){
            id1 = geno[id[f][i]]; id2 = geno[id[f][j]];
            HetHetCount = IBS0Count = het1Count = notMissingCount = 0;
            for(int m = 0; m < shortCount; m++){
               HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
               IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
               notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
               het1Count += oneoneCount[(GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m]) |
                  (GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m])]; // HomHet
            }
            if(het1Count+HetHetCount==0) continue;
            kinship = (HetHetCount - IBS0Count*2.0) / (HetHetCount*2+het1Count);
            phi = kin(ped[id[f][i]], ped[id[f][j]]);
            pi0 = 0.0;
            if(phi < 0.2)
               pi0 = 1-4*phi;
            else if(phi < 0.3 && ped[id[f][i]].isSib(ped[id[f][j]]))
               pi0 = 0.25;
            errorFlag = 0;
            inflation = (phi > 0)? kinship / phi: -1;
            if(phi > 0.03){  // up to 4th-degree relative
               if(inflation > 2 || inflation < 0.5)
                  errorFlag = 1;
               else if(inflation > 1.4142 || inflation < 0.70711)
                  errorFlag = 0.5;
            }else if(phi < 0.005){  // unrelated pair
               if(kinship > 0.0442) errorFlag = 1;
               else if(kinship > 0.0221) errorFlag = 0.5;
            }else{   // distant relatives
               if(kinship < -0.0221) errorFlag = 1;
               else errorFlag = 0.5;
            }
            ibs0 = IBS0Count*1.0/notMissingCount;
            degree = 4;
            if(phi > 0.0442)
               degree = int(-log(phi)/log(2.0) - 0.5);
            if(degree < 4){
               beforeCount[degree] ++;
               if(degree == 1)
                  if(pi0==0) beforeCount[5] ++;
            }else
               beforeCount[4] ++;
            degree = 4;
            if(kinship > 0.0442)
               degree = int(-log(kinship)/log(2.0) - 0.5);
            if(degree < 4){
               afterCount[degree] ++;
               if(degree == 1 && errorrateCutoff != _NAN_ && ibs0 < errorrateCutoff)
                  afterCount[5]++;
            }else
               afterCount[4] ++;

            if(degree==1 && phi==0.25){
               if(pi0==0)
                  IBS0PO.Push(ibs0);
               else
                  IBS0FS.Push(ibs0);
            }
            if(degree==1)
               IBS0L1.Push(ibs0);
            fprintf(fp, "%s\t%s\t%s\t%d\t%.3lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%G\n",
               (const char*)ped[id[f][i]].famid, (const char*)ped[id[f][i]].pid,
               (const char*)ped[id[f][j]].pid, notMissingCount, pi0, phi,
               HetHetCount*1.0/notMissingCount, ibs0,
               HetHetCount * 1.0 / (HetHetCount+het1Count), // CHet
               1 - IBS0Count * 1.0 / (notMissingCount - het1Count - HetHetCount),   // CHom
               kinship, errorFlag);
         }
   }
   fclose(fp);

   bool pedigreeFlag = false;
   for(int i = 0; i < 6; i++)
      if(beforeCount[i]) pedigreeFlag = true;
   if(pedigreeFlag){
      if(errorrateCutoff==_NAN_){
         if(IBS0PO.Length() && IBS0FS.Length())
            errorrateCutoff = (IBS0PO.Sum()/IBS0PO.Length() + IBS0FS.Sum()/IBS0FS.Length())*0.5;
         else if(IBS0FS.Length())
            errorrateCutoff = IBS0FS.Sum()/IBS0FS.Length()*0.5;
         else if(IBS0PO.Length())
            errorrateCutoff = 0.005;
         for(int i = 0; i < IBS0L1.Length(); i++)
            if(IBS0L1[i] < errorrateCutoff)
               afterCount[5]++;
      }
      printf("Within-family kinship data saved in file %s\n", (const char*)outfile);
      printRelationship(beforeCount, afterCount);
   }else
      printf("Each family consists of one individual.\n");

   if(ped.familyCount < 2) {
      if(ped.familyCount==1 && ped.families[0]->famid=="0")
         warning("All individuals with family ID 0 are considered as relatives.\n");
      printf("There is only one family.\n");
      return;
   }
   printf("Relationship inference across families starts at %s", currentTime());
#ifdef _OPENMP
   printf("%d CPU cores are used.\n", defaultMaxCoreCount);
#endif
   IntArray allpairs0, allpairs1(0), allpairs;

   // most computationally intensive here: SCREENING RELATIVES
   long int rawrelativeCount = ScreenCloseRelatives(allpairs0);
   //****************SCREENING RELATIVES end************************

   for(int i = 0; i < allpairs0.Length()/2; i++){
      id1 = allpairs0[i*2];
      id2 = allpairs0[i*2+1];
      if(ped[phenoid[id1]].famid != ped[phenoid[id2]].famid){
         allpairs1.Push(id1);
         allpairs1.Push(id2);
      }
   }
   unsigned int midrelativeCount = allpairs1.Length() / 2;

   if(rawrelativeCount==0 || midrelativeCount==0){
      printf("                                           ends at %s", currentTime());
      printf("No close relatives are inferred.\n\n");
      return;
   }

   unsigned short int **OG[2];   // original G
//   IntArray chrSeg;
//   double totalLength;
//   String segmessage;
   bool IBSvalidFlag = PreSegment(/*chrSeg, totalLength, segmessage*/);
   int segCount = chrSeg.Length()/2;

   if(IBSvalidFlag){
      IntArray allid(idCount);
      allid.Zero();
      for(int i = 0; i < midrelativeCount; i++){
         allid[allpairs1[i*2]] = 1;
         allid[allpairs1[i*2+1]] = 1;
      }
      for(int j = 0; j < 2; j++){
         OG[j] = new unsigned short int * [idCount];
         for(int i = 0; i < idCount; i++){
            if(allid[i]){
               OG[j][i] = new unsigned short int [shortCount];
               for(int m = 0; m < shortCount; m++)
                  OG[j][i][m] = 0;
            }else
               OG[j][i] = NULL;
         }
      }

   unsigned short int word[2];
   int m, oldm, oldbyte, oldbit;
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(word, m, oldm, oldbyte, oldbit)
{
#endif
   unsigned short int *aGG[2];
   for(int k = 0; k < 2; k++)
      aGG[k] = new unsigned short int [shortCount];
#ifdef _OPENMP
   #pragma omp for
#endif
   for(int i = 0; i < idCount; i++){
      if(OG[0][i]==NULL) continue;
      for(int k = 0; k < 2; k++)
         for(int w = 0; w < shortCount; w++)
            aGG[k][w] = 0;
      for(int w = 0; w < shortCount; w++){
         for(int k = 0; k < 2; k++)
            word[k] = GG[k][i][w];
         for(int j = 0; j < 16; j++){
            m = w*16+j;
            if(m>=markerCount) break;
            oldm = bigdataIdx[markerCount-1-m];
            oldbyte = oldm/16;
            oldbit = oldm%16;
            for(int k = 0; k < 2; k++)
               aGG[k][oldbyte] |= ((word[k] >> j)&1)<<oldbit;
         }
      }
      for(int k = 0; k < 2; k++)
         for(int w = 0; w < shortCount; w++)
            OG[k][i][w] = aGG[k][w];
   }
#ifdef _OPENMP
   for(int k = 0; k < 2; k++)
      delete []aGG[k];
}
#endif
   }else
      printf("%s\n", (const char*)segmessage);
   //   printf("Chromosome positions need to be sorted for the chromosome segment analysis to work\n");

   int relativeCount = 0;
   int stop1, stop2;
   stop1 = short_prescan;
   if(stop1 > shortCount) stop1 = shortCount;
   stop2 = stop1 + (short_prescan<<3);
   if(stop2 > shortCount) stop2 = shortCount;
   double lowerbound = pow(2.0, -(relativedegree+1.5));
   double kincutoff1 = 0.08838835; // 1-degree
   if(short_prescan==shortCount)
      kincutoff1 = lowerbound;
   else if(relativedegree == 2)
      kincutoff1 = 0.04419417;
   else if(relativedegree == 3)
      kincutoff1 = 0.015625;
   else if(relativedegree > 3)
      kincutoff1 = 0;
   double kincutoff2 = sqrt(lowerbound * kincutoff1);
   if(rawrelativeCount > 0){
      printf("Relationship inference speed-up by screening top informative SNPs:\n");
      printf("  Stage 1 (with %d SNPs): %u relative pairs detected (with Kinship > %.4lf)\n",
         stop1*16>markerCount? markerCount: stop1*16,
         (unsigned int)rawrelativeCount, kincutoff1);
      printf("  Stage 2 (with %d SNPs): %u out of %u remain related (with Kinship > %.4lf)\n",
         stop2*16>markerCount? markerCount: stop2*16,
         midrelativeCount, (unsigned int)rawrelativeCount, kincutoff2);
      printf("                               Screening ends at %s", currentTime());

      Vector pM(midrelativeCount);
      for(int i = 0; i < midrelativeCount; i++){
         pM[i] = allpairs1[i*2];
         pM[i] *= idCount;
         pM[i] += allpairs1[i*2+1];
      }
      QuickIndex index;
      index.Index(pM);
      allpairs.Dimension(midrelativeCount*2);
      for(int i = 0; i < midrelativeCount; i++){
         int k = index[i];
         allpairs[i*2] = allpairs1[k*2];
         allpairs[i*2+1] = allpairs1[k*2+1];
      }
   }else{ // not screened
      allpairs.Dimension(0);
      for(int i = 0; i < idCount; i++)
         for(int j = i+1; j < idCount; j++){
            allpairs.Push(i);
            allpairs.Push(j);
         }
      midrelativeCount = allpairs.Length() / 2;
   }
   char buffer[1024];
   StringArray buffers;
   buffers.Dimension(defaultMaxCoreCount);
   for(int i = 0; i < buffers.Length(); i++)
      buffers[i].Clear();
//   int blockSize = (midrelativeCount-1) / defaultMaxCoreCount + 1;
   double ibd2prop, maxLength;
   double CHet, CHom;
   int thread = 0;
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(HetHetCount, IBS0Count, het1Count, het2Count, notMissingCount,\
      id1, id2, kinship, smaller, buffer, ibd2prop, maxLength, CHet, CHom, thread) \
      reduction(+:relativeCount)
{
#endif
   IntArray startPos, stopPos;
#ifdef _OPENMP
   thread = omp_get_thread_num();
   #pragma omp for
#endif
   for(int i = 0; i < midrelativeCount; i++){
      id1 = allpairs[i*2];
      id2 = allpairs[i*2+1];
      HetHetCount = IBS0Count = het1Count = het2Count = notMissingCount = 0;
      for(int m = 0; m < stop1; m++){
         HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
         IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
         het1Count += oneoneCount[(~GG[0][id1][m]) & GG[0][id2][m] & GG[1][id1][m]]; // Het1Hom2
         het2Count += oneoneCount[GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]]; // Hom1Het2
         notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
      }
      CHom = 1-IBS0Count*1.0/(notMissingCount-het1Count-het2Count-HetHetCount);
      for(int m = stop1; m < shortCount; m++){
         HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
         IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
         het1Count += oneoneCount[(~GG[0][id1][m]) & GG[0][id2][m] & GG[1][id1][m]]; // Het1Hom2
         het2Count += oneoneCount[GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]]; // Hom1Het2
         notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
      }
      smaller = HetHetCount + (het1Count < het2Count? het1Count: het2Count);
      kinship = 0.5 - ((het1Count+het2Count)*0.25+IBS0Count)/smaller;
      CHet = HetHetCount*1.0/(het1Count+het2Count+HetHetCount);

      if(IBSvalidFlag && CHet<0.8 && kinship >= 0.125){ // IBD2 segment
         int localHetHetCount, localHetCount, totalHetHetCount, totalHetCount;
         int start, stop, blockstart, blockstop;
         totalHetHetCount = totalHetCount = 0;
         startPos.Dimension(0);
         stopPos.Dimension(0);
         for(int seg = 0; seg < segCount; seg++){
            localHetHetCount = localHetCount = 0;
            start = -1;
            blockstart = chrSeg[seg*2];
            for(int m = chrSeg[seg*2]; m <= chrSeg[seg*2+1]; m++){
               localHetHetCount += oneoneCount[~OG[0][id1][m] & OG[1][id1][m] & ~OG[0][id2][m] & OG[1][id2][m]];
               localHetCount += oneoneCount[(~OG[0][id1][m] & OG[1][id1][m] & (OG[0][id2][m] | OG[1][id2][m])) | // HetNomiss
                  (~OG[0][id2][m] & OG[1][id2][m] & (OG[0][id1][m] | OG[1][id1][m])) ];   // NomissHet
               if(localHetCount >= 100){ // start a new block
                  if(localHetHetCount >= 95){  // IBD2 segment
                     stop = m;
                     if(start == -1){  // a new IBD2 segment
                        start = blockstart;
                        for(start--; start>=chrSeg[seg*2] &&
                        ((~OG[0][id1][start] & OG[1][id1][start] & OG[0][id2][start])
                        | (~OG[0][id2][start] & OG[1][id2][start] & OG[0][id1][start]))==0; start--);
                        if(start>=chrSeg[seg*2] && // one error only
                        oneoneCount[(~OG[0][id1][start] & OG[1][id1][start] & OG[0][id2][start])
                        | (~OG[0][id2][start] & OG[1][id2][start] & OG[0][id1][start])]==1)
                           for(start--; start>=chrSeg[seg*2] &&
                           ((~OG[0][id1][start] & OG[1][id1][start] & OG[0][id2][start])
                           | (~OG[0][id2][start] & OG[1][id2][start] & OG[0][id1][start]))==0; start--);
                        start++;
                     }
                  }else if(start != -1){ // complete an old IBD2 segment
                     startPos.Push(start);
                     for(stop++; ((~OG[0][id1][stop] & OG[1][id1][stop] & OG[0][id2][stop])
                     | (~OG[0][id2][stop] & OG[1][id2][stop] & OG[0][id1][stop]))==0; stop++);
                     if(oneoneCount[(~OG[0][id1][stop] & OG[1][id1][stop] & OG[0][id2][stop])
                     | (~OG[0][id2][stop] & OG[1][id2][stop] & OG[0][id1][stop])]==1)  // one error only
                        for(stop++; ((~OG[0][id1][stop] & OG[1][id1][stop] & OG[0][id2][stop])
                        | (~OG[0][id2][stop] & OG[1][id2][stop] & OG[0][id1][stop]))==0; stop++);
                     stop--;
                     stopPos.Push(stop);
                     start = -1;
                  }
                  totalHetHetCount += localHetHetCount;
                  totalHetCount += localHetCount;
                  localHetHetCount = localHetCount = 0;
                  blockstop = m;
                  blockstart = m+1;
               }else if(m==chrSeg[seg*2+1] && start != -1){   // complete an old IBD2 segment
                  startPos.Push(start);
                  if(localHetHetCount > localHetCount-2) // add this block to the IBD2 segment
                     stopPos.Push(m);
                  else{  // without this block
                     for(stop++; ((~OG[0][id1][stop] & OG[1][id1][stop] & OG[0][id2][stop])
                     | (~OG[0][id2][stop] & OG[1][id2][stop] & OG[0][id1][stop]))==0; stop++);
                     stop--;
                     stopPos.Push(stop);
                  }
                  totalHetHetCount += localHetHetCount;
                  totalHetCount += localHetCount;
               }
            }  // end of word m
         }// end of seg

         int newsegCount = startPos.Length();
         ibd2prop = maxLength = 0.0;
         if(newsegCount){
            double length;
            for(int seg = 0; seg < newsegCount; seg++){
               if(seg+1 < newsegCount && // two IBD2 Segments can be merged
                  chromosomes[startPos[seg+1]*16] == chromosomes[stopPos[seg]*16+15] &&
                  positions[startPos[seg+1]*16]-positions[stopPos[seg]*16+15]<1.0){
                  startPos[seg+1] = startPos[seg];
                  continue;
               }
               if(stopPos[seg] == shortCount-1)
                  length = positions[markerCount-1] - positions[startPos[seg]*16];
               else
                  length = positions[stopPos[seg]*16+15] - positions[startPos[seg]*16];
               ibd2prop += length;
               if(length > maxLength)
                  maxLength = length;
            }
            if(maxLength < 1) ibd2prop=0.0;
            else
               ibd2prop /= totalLength;
         }
      }  // end of IBD2 segment

      if(relativedegree == 1){
         if(kinship < 0.125 || (kinship < lowerbound && CHom < 0.97 && ibd2prop < 0.05) ) continue;
      }else{
         if(kinship < lowerbound) continue;
      }
      
      relativeCount ++;
      sprintf(buffer, "%s\t%s\t%s\t%s\t%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf",
         (const char*)ped[phenoid[id1]].famid, (const char*)ped[phenoid[id1]].pid,
         (const char*)ped[phenoid[id2]].famid, (const char*)ped[phenoid[id2]].pid,
         notMissingCount, HetHetCount*1.0/notMissingCount, IBS0Count*1.0/notMissingCount,
         CHet, // CHet
         CHom, // CHom
         kinship);
      buffers[thread].Add(buffer);
      if(IBSvalidFlag){
         if(CHet<0.8){
            if(kinship < 0.125)
               sprintf(buffer, "\t-9\t-9\t%s", "REL");
            else if(CHom>0.99)
               sprintf(buffer, "\t%.4lf\t%.3lf\t%s", ibd2prop, maxLength, "PO");
            else if(CHom>0.96 && ibd2prop <= 0.05) // cannot be FS for small IBD2
               sprintf(buffer, "\t%.4lf\t%.3lf\t%s", ibd2prop, maxLength, "PO");
            else if(ibd2prop>0.05)  // FS must be Pr(IBD2)>5%
               sprintf(buffer, "\t%.4lf\t%.3lf\t%s", ibd2prop, maxLength, "FS");
            else
               sprintf(buffer, "\t%.4lf\t%.3lf\t%s", ibd2prop, maxLength, "2nd");
         }else // Duplicate
            sprintf(buffer, "\t-9\t-9\t%s", "DUP/MZ");
         buffers[thread].Add(buffer);
      }
      sprintf(buffer, "\n");
      buffers[thread].Add(buffer);
   }  // end  of loop over possible relative pairs
#ifdef _OPENMP
}
#endif

   printf("  Final Stage (with %d SNPs): %d pairs of relatives up to %d-degree are identified\n",
      markerCount, relativeCount, relativedegree);
   printf("                               Inference ends at %s", currentTime());

   outfile.Copy(prefix);
   outfile.Add(".kin0");
   fp = fopen(outfile, "wt");
   fprintf(fp, "FID1\tID1\tFID2\tID2\tN_SNP\tHetHet\tIBS0\tHetConc\tHomConc\tKinship");
   if(IBSvalidFlag) fprintf(fp, "\tPr_IBD2\tMaxIBD2\tRelType");
   fprintf(fp,"\n");
   for(int b = 0; b < buffers.Length(); b++)
      buffers[b].Write(fp);
   fclose(fp);
   printf("\nBetween-family relatives (kinship >= %.5lf) saved in file %s\n",
      lowerbound, (const char*)outfile);
   if(relativedegree == 1){
      printf("\nNote only duplicates and 1st-degree relative pairs are included here.\n");
      printf("  Specifying '--degree 2' for a higher degree relationship inference.\n\n");
   }
}

long int Engine::ScreenCloseRelatives(IntArray & rpList)
{
   if(shortCount < 128 || !bigdataFlag) return -1;  // screening is not needed
   int prescanCount[10];
   for(int i = 0; i < 10; i++)
      prescanCount[i] = 256 * (1<<(i*3));
   int short_prescan = relativedegree < 10? prescanCount[relativedegree-1]: prescanCount[9];
   if(short_prescan > shortCount) short_prescan=shortCount;
   if(faster > 1){
      short_prescan /= faster;
      if(short_prescan < 10) short_prescan = 10;
   }
   if(slower > 1){
      short_prescan *= slower;
      if(short_prescan > shortCount) short_prescan = shortCount;
   }
   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   int start1, stop1, start2, stop2;
   start1 = 0;
   stop1 = short_prescan;
   if(stop1 > shortCount) stop1 = shortCount;
   start2 = stop1;
   stop2 = stop1 + (short_prescan<<3);
   if(stop2 > shortCount) stop2 = shortCount;
   double lowerbound = pow(2.0, -(relativedegree+1.5));
   double kincutoff1 = 0.08838835; // 1-degree
   if(short_prescan==shortCount)
      kincutoff1 = lowerbound;
   else if(relativedegree == 2)
      kincutoff1 = 0.04419417;
   else if(relativedegree == 3)
      kincutoff1 = 0.015625;
   else if(relativedegree > 3)
      kincutoff1 = 0;
   double kincutoff2 = sqrt(lowerbound * kincutoff1);
   IntArray tArray;
   int m1, m2, m3;
   unsigned int **missingWordInOnePerson = new unsigned int *[idCount];
   int **missingWordInOnePersonCount = new int * [idCount];
   int **missingInOnePersonCount = new int * [idCount];
   int **hetInOnePersonCount = new int * [idCount];
   for(int i = 0; i < idCount; i++){
      missingWordInOnePersonCount[i] = new int[2];
      missingInOnePersonCount[i] = new int [2];
      hetInOnePersonCount[i] = new int [2];
   }
   for(int i = 0; i < idCount; i++)
      for(int j = 0; j < 2; j++)
         missingWordInOnePersonCount[i][j] = missingInOnePersonCount[i][j] = hetInOnePersonCount[i][j] = 0;
   unsigned short int *masks = new unsigned short int[stop2];
   for(int i = 0; i < stop2;  i++)
      masks[i] = 0xFFFF;
   if(stop2 == shortCount)
      masks[stop2-1] = markerCount%16? ((1 << (markerCount % 16))-1): 0xFFFF;
   for(int i = 0; i < idCount; i++){
      tArray.Dimension(0);
      for(int m = 0; m < stop1; m++)
         if(((~GG[0][i][m]) & (~GG[1][i][m]) & masks[m])!=0){
            missingWordInOnePersonCount[i][0] ++;
            missingInOnePersonCount[i][0] += oneoneCount[(~GG[0][i][m]) & (~GG[1][i][m]) & masks[m]];
            tArray.Push(m);
          }
      for(int m = start2; m < stop2; m++)
         if(((~GG[0][i][m]) & (~GG[1][i][m]) & masks[m])!=0){
            missingWordInOnePersonCount[i][1] ++;
            missingInOnePersonCount[i][1] += oneoneCount[(~GG[0][i][m]) & (~GG[1][i][m]) & masks[m]];
            tArray.Push(m);
          }
      if(tArray.Length()){
         missingWordInOnePerson[i] = new unsigned int [tArray.Length()];
         for(int m = 0; m < tArray.Length(); m++)
            missingWordInOnePerson[i][m] = tArray[m];
      }else
         missingWordInOnePerson[i] = NULL;
      for(int m = 0; m < stop1; m++)
         hetInOnePersonCount[i][0] += oneoneCount[(~GG[0][i][m]) & GG[1][i][m]];
      for(int m = start2; m < stop2; m++)
         hetInOnePersonCount[i][1] += oneoneCount[(~GG[0][i][m]) & GG[1][i][m]];
   }
   delete []masks;
   const int cutoffMissingCount = stop1*16-MINSNPCOUNT;
   IntArray *ibuffer = new IntArray[defaultMaxCoreCount];
   for(int i = 0; i < defaultMaxCoreCount; i++)
      ibuffer[i].Dimension(0);
   double threshold, smaller, kinship;
   int id1, id2, IBS0Count, het1Count, het2Count, HetHomCount, HomHomCount, HetHetCount;
   const int LOOPBLOCKINGSIZE=32;
   IntArray loopIndex[2];
   loopIndex[0].Dimension(0); loopIndex[1].Dimension(0);
   for(int i = 0; i < idCount; i += LOOPBLOCKINGSIZE)
      for(int j = i; j < idCount; j += LOOPBLOCKINGSIZE){
         loopIndex[0].Push(i);
         loopIndex[1].Push(j);
      }
   int loopIndexLength = loopIndex[0].Length();
   long int *rawrelativeCount = new long int [defaultMaxCoreCount];
   for(int t = 0; t < defaultMaxCoreCount; t++)
      rawrelativeCount[t] = 0;

   int thread = 0;
   if(relativedegree == 1)
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(HomHomCount, IBS0Count, HetHetCount, HetHomCount, het1Count, het2Count,\
      threshold, id1, id2, kinship, smaller, m1, m2, m3, thread)
{
   thread = omp_get_thread_num();
   #pragma omp for
#endif
   for(int k = 0; k < loopIndexLength; k++){
      int i = loopIndex[0][k];
      int iMax = (i > idCount - LOOPBLOCKINGSIZE? idCount: i + LOOPBLOCKINGSIZE);
      int j = loopIndex[1][k];
      int jMax = (j > idCount - LOOPBLOCKINGSIZE? idCount: j + LOOPBLOCKINGSIZE);
      int jMin = j;
      for(id1 = i; id1 < iMax; id1++){
         if(missingInOnePersonCount[id1][0]>=cutoffMissingCount) continue;
         if(i == j) jMin = id1 + 1; // diagonal blocks only
         for(id2 = jMin; id2 < jMax; id2++){
            if(missingInOnePersonCount[id2][0]>=cutoffMissingCount) continue;
            // Stage 1: all pairs
            // Stop rule 1: C_Hom < 65%
            // Stop rule 2: kinship < 0.0625
            // Stop rule 3: kinship < 0.0884 && C_Hom < 70%

            // 1st quarter
            for(HomHomCount = IBS0Count = m1 = 0, m2 = stop1/4; m1 < m2; m1++){
               HomHomCount += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1]];
               IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
            }
            if(IBS0Count >= HomHomCount * 0.5) continue; // C_Hom < 50%
            // 2nd quarter
            for(m1 = stop1/4, m2 = stop1/2; m1 < m2; m1++){
               HomHomCount += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1]];
               IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
            }
            if(IBS0Count >= HomHomCount * 0.4) continue; // C_Hom < 60%
            // 3rd quarter
            for(m1 = stop1/2, m2 = stop1*3/4; m1 < m2; m1++){
               HomHomCount += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1]];
               IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
            }
            if(IBS0Count >= HomHomCount * 0.375) continue; // C_Hom < 62.5%
            // 4th quarter
            for(m1 = stop1*3/4; m1 < stop1; m1++){
               HomHomCount += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1]];
               IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
            }
            if(IBS0Count >= HomHomCount * 0.35) continue; // C_Hom < 65%

            m3 = (hetInOnePersonCount[id1][0] < hetInOnePersonCount[id2][0])?
                  hetInOnePersonCount[id1][0]: hetInOnePersonCount[id2][0];
            HetHomCount = stop1 * 32 - HomHomCount * 2 - hetInOnePersonCount[id1][0] - hetInOnePersonCount[id2][0]
               - missingInOnePersonCount[id1][0] - missingInOnePersonCount[id2][0];
            threshold = (0.5-kincutoff1*0.707107)*m3 - HetHomCount*0.25;   // kinship < 0.0625
            double cutoffDist = (0.5-kincutoff1)*m3 - HetHomCount*0.25;  // C_Hom < 70% && kinship < 0.0884
            if(HomHomCount * 0.3 > cutoffDist) cutoffDist = HomHomCount * 0.3;
            if(cutoffDist < threshold) threshold = cutoffDist;
            if(IBS0Count >= threshold) continue;

            int m1Max = missingWordInOnePersonCount[id1][0];
            int m2Max = missingWordInOnePersonCount[id2][0];
            for(HetHomCount = m1 = m2 = 0; m1 < m1Max; m1++){
               m3 = missingWordInOnePerson[id1][m1]; // m for id1
               for(; m2 < m2Max && missingWordInOnePerson[id2][m2] < m3; m2++);
               if( m2 == m2Max ) break;
               if(m3 == missingWordInOnePerson[id2][m2])   // id1 and id2 are both missing
                  HetHomCount += oneoneCount[~(GG[0][id1][m3] | GG[0][id2][m3] | GG[1][id1][m3] | GG[1][id2][m3])&0xFFFF];
            }  // Store MissMissCount in HetHomCount
            for(m1 = m3 = 0; m1 < m2Max; m1++){
               m2 = missingWordInOnePerson[id2][m1];
               m3 += oneoneCount[~(GG[0][id1][m2] | GG[0][id2][m2] | GG[1][id2][m2]) & GG[1][id1][m2]];
            }
            het1Count = hetInOnePersonCount[id1][0] - m3;
            for(m1 = m3 = 0; m1 < m1Max; m1++){
               m2 = missingWordInOnePerson[id1][m1];
               m3 += oneoneCount[~(GG[0][id1][m2] | GG[0][id2][m2] | GG[1][id1][m2]) & GG[1][id2][m2]];
            }
            het2Count = hetInOnePersonCount[id2][0] - m3;
            HetHomCount = stop1 * 32 - HomHomCount * 2 - het1Count - het2Count +
               (-missingInOnePersonCount[id1][0] - missingInOnePersonCount[id2][0] + HetHomCount) * 2;
            m3 = (het1Count < het2Count? het1Count: het2Count);
            threshold = (0.5-kincutoff1*0.707107)*m3 - HetHomCount*0.25;  // kinship < 0.0625
            cutoffDist = (0.5-kincutoff1)*m3 - HetHomCount*0.25;  // C_Hom < 70% && kinship < 0.0884
            if(HomHomCount * 0.7 > cutoffDist) cutoffDist = HomHomCount * 0.7;
            if(cutoffDist < threshold) threshold = cutoffDist;
            if(IBS0Count >= threshold) continue;
            rawrelativeCount[thread] ++;

            // Stage 2: few pairs
            for(HetHomCount -= (missingInOnePersonCount[id1][1]+missingInOnePersonCount[id2][1]),  // different non-missing genotypes
               m1 = start2; m1 < stop2; m1++) // different-genotpe Count
                  HetHomCount += oneoneCount[GG[0][id1][m1]^GG[0][id2][m1]];
            het1Count += hetInOnePersonCount[id1][1];
            het2Count += hetInOnePersonCount[id2][1];
            cutoffDist = (2-kincutoff2*4) * (het1Count < het2Count? het1Count: het2Count)-IBS0Count*4;
            if(HetHomCount >= cutoffDist) continue; // too big distance already
            cutoffDist = (cutoffDist+IBS0Count*4-HetHomCount)*0.25;
            for(m1 = start2; m1 < (start2+stop2)/2; m1++)
               IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
            if(IBS0Count >= cutoffDist * 1.1) continue;
            for(m1 = (start2+stop2)/2; (m1 < stop2) && (IBS0Count < cutoffDist); m1++)
               IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
            if(IBS0Count >= cutoffDist) continue;   // too big distance
            // now add the MissMissCount back
            int mm1Max = m1Max+missingWordInOnePersonCount[id1][1];
            int mm2Max = m2Max+missingWordInOnePersonCount[id2][1];
            for(m1 = m1Max, m2 = m2Max; m1 < mm1Max; m1++){
               m3 = missingWordInOnePerson[id1][m1]; // m for id1
               for(; m2 < mm2Max && missingWordInOnePerson[id2][m2] < m3; m2++);
               if( m2 == mm2Max ) break;
               if(m3 == missingWordInOnePerson[id2][m2])   // id1 and id2 are both missing
                  HetHomCount += 2 * oneoneCount[~(GG[0][id1][m3] | GG[0][id2][m3] | GG[1][id1][m3] | GG[1][id2][m3]) & 0xFFFF];
            }
            for(m1 = m2Max, m3 = 0; m1 < mm2Max; m1++){
               m2 = missingWordInOnePerson[id2][m1];
               m3 += oneoneCount[~(GG[0][id2][m2] | GG[0][id1][m2] | GG[1][id2][m2]) & GG[1][id1][m2]];
            }
            het1Count -= m3;
            HetHomCount += m3;
            for(m1 = m1Max, m3=0; m1 < mm1Max; m1++){
               m2 = missingWordInOnePerson[id1][m1];
               m3 += oneoneCount[~(GG[0][id1][m2] | GG[0][id2][m2] | GG[1][id1][m2]) & GG[1][id2][m2]];
            }
            het2Count -= m3;
            HetHomCount += m3;
            HetHetCount = (het1Count + het2Count - HetHomCount)/2;
            het1Count -= HetHetCount;
            het2Count -= HetHetCount;
            smaller = HetHetCount + (het1Count < het2Count? het1Count: het2Count);
            kinship = 0.5 - ((het1Count+het2Count)*0.25+IBS0Count)/smaller;
            if(kinship < kincutoff2) continue;
            ibuffer[thread].Push(id1);
            ibuffer[thread].Push(id2);
         }  // end of id2
      }  // end of id1
   }  // end of loop blocking
#ifdef _OPENMP
}
#endif
   else // more distant than first-degree
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(HomHomCount, IBS0Count, HetHetCount, HetHomCount, het1Count, het2Count,\
      threshold, id1, id2, kinship, smaller, m1, m2, m3, thread)
{
   thread = omp_get_thread_num();
   #pragma omp for
#endif
   for(int k = 0; k < loopIndexLength; k++){
      int i = loopIndex[0][k];
      int iMax = (i > idCount - LOOPBLOCKINGSIZE? idCount: i + LOOPBLOCKINGSIZE);
      int j = loopIndex[1][k];
      int jMax = (j > idCount - LOOPBLOCKINGSIZE? idCount: j + LOOPBLOCKINGSIZE);
      int jMin = j;
      for(id1 = i; id1 < iMax; id1++){
         if(missingInOnePersonCount[id1][0]>=cutoffMissingCount) continue;
         if(i == j) jMin = id1 + 1; // diagonal blocks only
         for(id2 = jMin; id2 < jMax; id2++){
            if(missingInOnePersonCount[id2][0]>=cutoffMissingCount) continue;
            // Stage 1: all pairs
            if(relativedegree == 2){
               // 1st quarter
               for(HomHomCount = IBS0Count = m1 = 0, m2 = stop1/16; m1 < m2; m1++){
                  HomHomCount += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1]];
                  IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
               }
               if(IBS0Count >= HomHomCount * 0.5) continue; // C_Hom < 50%
               // 2nd quarter
               for(m1 = stop1/16, m2 = stop1/8; m1 < m2; m1++){
                  HomHomCount += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1]];
                  IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
               }
               if(IBS0Count >= HomHomCount * 0.45) continue; // C_Hom < 55%
               // 3rd quarter
               for(m1 = stop1/8, m2 = stop1/4; m1 < m2; m1++){
                  HomHomCount += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1]];
                  IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
               }
               if(IBS0Count >= HomHomCount * 0.42) continue; // C_Hom < 58%

               // 4th quarter
               for(m1 = stop1/4; m1 < stop1; m1++){
                  HomHomCount += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1]];
                  IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
               }
               if(IBS0Count >= HomHomCount * 0.375) continue; // C_Hom < 62.5%
            }else{   // more distant
               // 1st quarter
               for(HomHomCount = IBS0Count = m1 = 0, m2 = stop1/16; m1 < m2; m1++){
                  HomHomCount += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1]];
                  IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
               }
               if(IBS0Count >= HomHomCount * 0.55) continue; // C_Hom < 45%
               // 2nd quarter
               for(m1 = stop1/16, m2 = stop1/8; m1 < m2; m1++){
                  HomHomCount += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1]];
                  IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
               }
               if(IBS0Count >= HomHomCount * 0.5) continue; // C_Hom < 50%
               // 3rd quarter
               for(m1 = stop1/8, m2 = stop1/4; m1 < m2; m1++){
                  HomHomCount += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1]];
                  IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
               }
               if(IBS0Count >= HomHomCount * 0.45) continue; // C_Hom < 55%

               // 4th quarter
               for(m1 = stop1/4; m1 < stop1; m1++){
                  HomHomCount += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1]];
                  IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
               }
               if(IBS0Count >= HomHomCount * 0.45) continue; // C_Hom < 55%
            }

            m3 = (hetInOnePersonCount[id1][0] < hetInOnePersonCount[id2][0])?
                  hetInOnePersonCount[id1][0]: hetInOnePersonCount[id2][0];
            HetHomCount = stop1 * 32 - HomHomCount * 2 - hetInOnePersonCount[id1][0] - hetInOnePersonCount[id2][0]
               - missingInOnePersonCount[id1][0] - missingInOnePersonCount[id2][0];
            threshold = (0.5-kincutoff1)*m3 - HetHomCount*0.25;   // kinship < 0.0625
            if(IBS0Count >= threshold) continue;

            int m1Max = missingWordInOnePersonCount[id1][0];
            int m2Max = missingWordInOnePersonCount[id2][0];
            for(HetHomCount = m1 = m2 = 0; m1 < m1Max; m1++){
               m3 = missingWordInOnePerson[id1][m1]; // m for id1
               for(; m2 < m2Max && missingWordInOnePerson[id2][m2] < m3; m2++);
               if( m2 == m2Max ) break;
               if(m3 == missingWordInOnePerson[id2][m2])   // id1 and id2 are both missing
                  HetHomCount += oneoneCount[~(GG[0][id1][m3] | GG[0][id2][m3] | GG[1][id1][m3] | GG[1][id2][m3])&0xFFFF];
            }  // Store MissMissCount in HetHomCount
            for(m1 = m3 = 0; m1 < m2Max; m1++){
               m2 = missingWordInOnePerson[id2][m1];
               m3 += oneoneCount[~(GG[0][id1][m2] | GG[0][id2][m2] | GG[1][id2][m2]) & GG[1][id1][m2]];
            }
            het1Count = hetInOnePersonCount[id1][0] - m3;
            for(m1 = m3 = 0; m1 < m1Max; m1++){
               m2 = missingWordInOnePerson[id1][m1];
               m3 += oneoneCount[~(GG[0][id1][m2] | GG[0][id2][m2] | GG[1][id1][m2]) & GG[1][id2][m2]];
            }
            het2Count = hetInOnePersonCount[id2][0] - m3;
            HetHomCount = stop1 * 32 - HomHomCount * 2 - het1Count - het2Count +
               (-missingInOnePersonCount[id1][0] - missingInOnePersonCount[id2][0] + HetHomCount) * 2;
            m3 = (het1Count < het2Count? het1Count: het2Count);
            threshold = (0.5-kincutoff1)*m3 - HetHomCount*0.25;  // kinship < 0.0625
            if(IBS0Count >= threshold) continue;
            rawrelativeCount[thread] ++;

            // Stage 2: few pairs
            for(HetHomCount -= (missingInOnePersonCount[id1][1]+missingInOnePersonCount[id2][1]),  // different non-missing genotypes
               m1 = start2; m1 < stop2; m1++) // different-genotpe Count
                  HetHomCount += oneoneCount[GG[0][id1][m1]^GG[0][id2][m1]];
            het1Count += hetInOnePersonCount[id1][1];
            het2Count += hetInOnePersonCount[id2][1];
            double cutoffDist = (2-kincutoff2*4) * (het1Count < het2Count? het1Count: het2Count)-IBS0Count*4;
            if(HetHomCount >= cutoffDist) continue; // too big distance already
            cutoffDist = (cutoffDist+IBS0Count*4-HetHomCount)*0.25;
            for(m1 = start2; m1 < (start2+stop2)/2; m1++)
               IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
            if(IBS0Count >= cutoffDist * 1.1) continue;
            for(m1 = (start2+stop2)/2; (m1 < stop2) && (IBS0Count < cutoffDist); m1++)
               IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
            if(IBS0Count >= cutoffDist) continue;   // too big distance
            // now add the MissMissCount back
            int mm1Max = m1Max+missingWordInOnePersonCount[id1][1];
            int mm2Max = m2Max+missingWordInOnePersonCount[id2][1];
            for(m1 = m1Max, m2 = m2Max; m1 < mm1Max; m1++){
               m3 = missingWordInOnePerson[id1][m1]; // m for id1
               for(; m2 < mm2Max && missingWordInOnePerson[id2][m2] < m3; m2++);
               if( m2 == mm2Max ) break;
               if(m3 == missingWordInOnePerson[id2][m2])   // id1 and id2 are both missing
                  HetHomCount += 2 * oneoneCount[~(GG[0][id1][m3] | GG[0][id2][m3] | GG[1][id1][m3] | GG[1][id2][m3]) & 0xFFFF];
            }
            for(m1 = m2Max, m3 = 0; m1 < mm2Max; m1++){
               m2 = missingWordInOnePerson[id2][m1];
               m3 += oneoneCount[~(GG[0][id2][m2] | GG[0][id1][m2] | GG[1][id2][m2]) & GG[1][id1][m2]];
            }
            het1Count -= m3;
            HetHomCount += m3;
            for(m1 = m1Max, m3=0; m1 < mm1Max; m1++){
               m2 = missingWordInOnePerson[id1][m1];
               m3 += oneoneCount[~(GG[0][id1][m2] | GG[0][id2][m2] | GG[1][id1][m2]) & GG[1][id2][m2]];
            }
            het2Count -= m3;
            HetHomCount += m3;
            HetHetCount = (het1Count + het2Count - HetHomCount)/2;
            het1Count -= HetHetCount;
            het2Count -= HetHetCount;
            smaller = HetHetCount + (het1Count < het2Count? het1Count: het2Count);
            kinship = 0.5 - ((het1Count+het2Count)*0.25+IBS0Count)/smaller;
            if(kinship < kincutoff2) continue;
            ibuffer[thread].Push(id1);
            ibuffer[thread].Push(id2);
         }  // end of id2
      }  // end of id1
   }  // end of loop blocking
#ifdef _OPENMP
}
#endif
   for(int i = 0; i < idCount; i++){
      if(missingWordInOnePerson[i]) delete []missingWordInOnePerson[i];
      delete []missingWordInOnePersonCount[i];
      delete []missingInOnePersonCount[i];
      delete []hetInOnePersonCount[i];
   }
   delete []missingWordInOnePersonCount;
   delete []missingInOnePersonCount;
   delete []hetInOnePersonCount;
   rpList.Dimension(0);
   long int rawcount = 0;
   for(int i = 0; i < defaultMaxCoreCount; i++){
      rawcount += rawrelativeCount[i];
      rpList.Append(ibuffer[i]);
   }
   delete []ibuffer;
   return rawcount;
}









/*
void Engine::ComputeBigDataKinship()
{
   if(shortCount==0) error("No genotype data");
   printf("Autosome genotypes stored in %d", shortCount);
   printf(" words for each of %d individuals.\n", idCount);

   printf("\nOptions in effect:\n");
   if(errorrateCutoff!=_NAN_)
      printf("\t--errorrate %G\n", errorrateCutoff);
   printf("\t--related\n");
   if(relativedegree>1)
      printf("\t--degree %d\n", relativedegree);
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

   int prescanCount[10];
   for(int i = 0; i < 10; i++)
      prescanCount[i] = 256 * (1<<(i*3));
   //      prescanCount[i] = 128 * (1<<(i*3));
   int short_prescan = relativedegree < 10? prescanCount[relativedegree-1]: prescanCount[9];
   if(short_prescan > shortCount || (!bigdataFlag)) short_prescan=shortCount;

   if(bigdataFlag){
      if(faster > 1){
         short_prescan /= faster;
         if(short_prescan < 10) short_prescan = 10;
      }
      if(slower > 1){
         short_prescan *= slower;
         if(short_prescan > shortCount) short_prescan = shortCount;
      }
   }

   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   String outfile;
   outfile.Copy(prefix);
   outfile.Add(".kin");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "FID\tID1\tID2\tN_SNP\tZ0\tPhi\tHetHet\tIBS0\tHetConc\tHomConc\tKinship\tError\n");
   double kinship, inflation, smaller, larger;
   int id1, id2;
   Kinship kin;
   int het1, het2, homhom, unequal;
   int HomHomCount, HetHetCount, IBS0Count, het1Count, het2Count, notMissingCount, HetHomCount;
   double phi, pi0, errorFlag;
   int beforeCount[6], afterCount[6];
   Vector IBS0PO(0), IBS0FS(0), IBS0L1(0);
   int degree; double ibs0;
   for(int i = 0; i < 6; i++) beforeCount[i] = afterCount[i] = 0;
   for(int f = 0; f < ped.familyCount; f++){
      if(id[f].Length()<2) continue;
      kin.Setup(*ped.families[f]);
      for(int i = 0; i < id[f].Length(); i++)
         for(int j = i+1; j < id[f].Length(); j++){
            id1 = geno[id[f][i]]; id2 = geno[id[f][j]];
            HetHetCount = IBS0Count = het1Count = notMissingCount = 0;
            for(int m = 0; m < shortCount; m++){
               HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
               IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
               notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
               het1Count += oneoneCount[(GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m]) |
                  (GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m])]; // HomHet
            }
            if(het1Count+HetHetCount==0) continue;
            kinship = (HetHetCount - IBS0Count*2.0) / (HetHetCount*2+het1Count);

            phi = kin(ped[id[f][i]], ped[id[f][j]]);
            pi0 = 0.0;
            if(phi < 0.2)
               pi0 = 1-4*phi;
            else if(phi < 0.3 && ped[id[f][i]].isSib(ped[id[f][j]]))
               pi0 = 0.25;
            errorFlag = 0;
            inflation = (phi > 0)? kinship / phi: -1;
            if(phi > 0.03){  // up to 4th-degree relative
               if(inflation > 2 || inflation < 0.5)
                  errorFlag = 1;
               else if(inflation > 1.4142 || inflation < 0.70711)
                  errorFlag = 0.5;
            }else if(phi < 0.005){  // unrelated pair
               if(kinship > 0.0442) errorFlag = 1;
               else if(kinship > 0.0221) errorFlag = 0.5;
            }else{   // distant relatives
               if(kinship < -0.0221) errorFlag = 1;
               else errorFlag = 0.5;
            }
            ibs0 = IBS0Count*1.0/notMissingCount;
            degree = 4;
            if(phi > 0.0442)
               degree = int(-log(phi)/log(2.0) - 0.5);
            if(degree < 4){
               beforeCount[degree] ++;
               if(degree == 1)
                  if(pi0==0) beforeCount[5] ++;
            }else
               beforeCount[4] ++;
            degree = 4;
            if(kinship > 0.0442)
               degree = int(-log(kinship)/log(2.0) - 0.5);
            if(degree < 4){
               afterCount[degree] ++;
               if(degree == 1 && errorrateCutoff != _NAN_ && ibs0 < errorrateCutoff)
                  afterCount[5]++;
            }else
               afterCount[4] ++;

            if(degree==1 && phi==0.25){
               if(pi0==0)
                  IBS0PO.Push(ibs0);
               else
                  IBS0FS.Push(ibs0);
            }
            if(degree==1)
               IBS0L1.Push(ibs0);
            fprintf(fp, "%s\t%s\t%s\t%d\t%.3lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%G\n",
               (const char*)ped[id[f][i]].famid, (const char*)ped[id[f][i]].pid,
               (const char*)ped[id[f][j]].pid, notMissingCount, pi0, phi,
               HetHetCount * 1.0 / notMissingCount, ibs0,
               HetHetCount * 1.0 / (HetHetCount+het1Count), // CHet
               1 - IBS0Count * 1.0 / (notMissingCount - het1Count - HetHetCount),   // CHom
               kinship, errorFlag);
         }
   }
   fclose(fp);

   bool pedigreeFlag = false;
   for(int i = 0; i < 6; i++)
      if(beforeCount[i]) pedigreeFlag = true;
   if(pedigreeFlag){
      if(errorrateCutoff==_NAN_){
         if(IBS0PO.Length() && IBS0FS.Length())
            errorrateCutoff = (IBS0PO.Sum()/IBS0PO.Length() + IBS0FS.Sum()/IBS0FS.Length())*0.5;
         else if(IBS0FS.Length())
            errorrateCutoff = IBS0FS.Sum()/IBS0FS.Length()*0.5;
         else if(IBS0PO.Length())
            errorrateCutoff = 0.005;
         for(int i = 0; i < IBS0L1.Length(); i++)
            if(IBS0L1[i] < errorrateCutoff)
               afterCount[5]++;
      }

      printf("Within-family kinship data saved in file %s\n", (const char*)outfile);
      printRelationship(beforeCount, afterCount);
   }else
      printf("Each family consists of one individual.\n");

   if(ped.familyCount < 2) {
      if(ped.familyCount==1 && ped.families[0]->famid=="0")
         warning("All individuals with family ID 0 are considered as relatives.\n");
      printf("There is only one family.\n");
      return;
   }
   printf("Relationship inference across families starts at %s", currentTime());
   int start1, stop1, start2, stop2, start3, stop3;
   start1 = 0;
   stop1 = short_prescan;
   if(stop1 > shortCount) stop1 = shortCount;
   start2 = stop1;
   stop2 = stop1 + (short_prescan<<3);
   if(stop2 > shortCount) stop2 = shortCount;
   start3 = stop2;
   stop3 = shortCount;

   long int relativeCount = 0;
   long long int rawrelativeCount = 0;
   long int midrelativeCount = 0;
   double lowerbound = pow(2.0, -(relativedegree+1.5));
   double kincutoff1 = 0.08838835; // 1-degree
   if(short_prescan==shortCount)
      kincutoff1 = lowerbound;
   else if(relativedegree == 2)
      kincutoff1 = 0.04419417;
   else if(relativedegree == 3)
      kincutoff1 = 0.015625;
   else if(relativedegree > 3)
      kincutoff1 = 0;
   double kincutoff2 = sqrt(lowerbound * kincutoff1);
   IntArray tArray;
   int m1, m2, m3;
   unsigned int **missingWordInOnePerson = new unsigned int *[idCount];
   int *missingWordInOnePersonCount[2];
   int *missingInOnePersonCount[2];
   int **hetInOnePersonCount = new int * [idCount];
   for(int i = 0; i < 2; i++){
      missingWordInOnePersonCount[i] = new int[idCount];
      missingInOnePersonCount[i] = new int[idCount];
   }
   for(int i = 0; i < idCount; i++)
      hetInOnePersonCount[i] = new int [5];
   for(int i = 0; i < idCount; i++)
      for(int j = 0; j < 2; j++)
         missingWordInOnePersonCount[j][i] = missingInOnePersonCount[j][i] = 0;
   for(int i = 0; i < idCount; i++)
      for(int j = 0; j < 5; j++)
         hetInOnePersonCount[i][j] = 0;
   unsigned short int *masks = new unsigned short int[stop2];
   for(int i = 0; i < stop2;  i++)
      masks[i] = 0xFFFF;
   if(stop2 == shortCount)
      masks[stop2-1] = markerCount%16? ((1 << (markerCount % 16))-1): 0xFFFF;
   for(int i = 0; i < idCount; i++){
      tArray.Dimension(0);
      for(int m = 0; m < stop1; m++)
         if(((~GG[0][i][m]) & (~GG[1][i][m]) & masks[m])!=0){
            missingWordInOnePersonCount[0][i] ++;
            missingInOnePersonCount[0][i] += oneoneCount[(~GG[0][i][m]) & (~GG[1][i][m]) & masks[m]];
            tArray.Push(m);
          }
      for(int m = start2; m < stop2; m++)
         if(((~GG[0][i][m]) & (~GG[1][i][m]) & masks[m])!=0){
            missingWordInOnePersonCount[1][i] ++;
            missingInOnePersonCount[1][i] += oneoneCount[(~GG[0][i][m]) & (~GG[1][i][m]) & masks[m]];
            tArray.Push(m);
          }
      if(tArray.Length()){
         missingWordInOnePerson[i] = new unsigned int [tArray.Length()];
         for(int m = 0; m < tArray.Length(); m++)
            missingWordInOnePerson[i][m] = tArray[m];
      }else
         missingWordInOnePerson[i] = NULL;
   }
   if(relativedegree == 1)
      for(int i = 0; i < idCount; i++){
         for(int mm = 0; mm < 4; mm++){
            for(int m = mm*(stop1/4); m < (mm+1)*(stop1/4); m++)
               hetInOnePersonCount[i][mm] += oneoneCount[(~GG[0][i][m]) & GG[1][i][m]];
            if(mm>0) hetInOnePersonCount[i][mm] += hetInOnePersonCount[i][mm-1];
         }
         for(int m = start2; m < stop2; m++)
            hetInOnePersonCount[i][4] += oneoneCount[(~GG[0][i][m]) & GG[1][i][m]];
   }else{ // second-degree or more distant
      int fr[5]={0, 1, 2, 4, 16};
      for(int i = 0; i < idCount; i++){
         for(int mm = 0; mm < 4; mm++){
            for(int m = fr[mm]*(stop1/16); m < fr[mm+1]*(stop1/16); m++)
               hetInOnePersonCount[i][mm] += oneoneCount[(~GG[0][i][m]) & GG[1][i][m]];
            if(mm>0) hetInOnePersonCount[i][mm] += hetInOnePersonCount[i][mm-1];
         }
         for(int m = start2; m < stop2; m++)
            hetInOnePersonCount[i][4] += oneoneCount[(~GG[0][i][m]) & GG[1][i][m]];
      }
   }

   delete []masks;
   const int cutoffMissingCount = stop1*16-MINSNPCOUNT;

   char buffer[1024];
   StringArray buffers;
   buffers.Dimension(defaultMaxCoreCount);
   for(int i = 0; i < buffers.Length(); i++)
      buffers[i].Clear();

#ifdef _OPENMP
   printf("%d CPU cores are used.\n", defaultMaxCoreCount);
#endif
   double threshold;
   const int LOOPBLOCKINGSIZE=100;
   IntArray loopIndex[2];
   loopIndex[0].Dimension(0); loopIndex[1].Dimension(0);
   for(int i = 0; i < idCount; i += LOOPBLOCKINGSIZE)
      for(int j = i; j < idCount; j += LOOPBLOCKINGSIZE){
         loopIndex[0].Push(i);
         loopIndex[1].Push(j);
      }
   int loopIndexLength = loopIndex[0].Length();
   int blockSize = (loopIndexLength-1) / defaultMaxCoreCount + 1;
   if(relativedegree == 1)
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
      private(HomHomCount, HetHetCount, IBS0Count, het1Count, het2Count, notMissingCount, threshold,\
      HetHomCount, id1, id2, kinship, smaller, buffer, m1, m2, m3) \
      reduction(+:rawrelativeCount, midrelativeCount, relativeCount)
#endif
   for(int k = 0; k < loopIndexLength; k++){
      int i = loopIndex[0][k];
      int iMax = (i > idCount - LOOPBLOCKINGSIZE? idCount: i + LOOPBLOCKINGSIZE);
      int j = loopIndex[1][k];
      int jMax = (j > idCount - LOOPBLOCKINGSIZE? idCount: j + LOOPBLOCKINGSIZE);
      int jMin = j;
      for(id1 = i; id1 < iMax; id1++){
         if(missingInOnePersonCount[0][id1]>=cutoffMissingCount) continue;
         if(i == j) jMin = id1 + 1; // diagonal blocks only
         for(id2 = jMin; id2 < jMax; id2++){
            if(missingInOnePersonCount[0][id2]>=cutoffMissingCount) continue;
            if(ped[phenoid[id1]].famid == ped[phenoid[id2]].famid) continue;
            // Stage 1: all pairs
            // Stop rule 1: kinship < 0.0625
            // Stop rule 2: C_Hom < 70%
            // Stop rule 3: kinship < 0.0884 && C_Hom < 75%

            // 1st quarter
            for(HomHomCount = IBS0Count = m1 = 0, m2 = stop1/4; m1 < m2; m1++){
               HomHomCount += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1]];
               IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
            }
            if(IBS0Count >= HomHomCount * 0.5) continue; // C_Hom < 50%
            // 2nd quarter
            for(m1 = stop1/4, m2 = stop1/2; m1 < m2; m1++){
               HomHomCount += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1]];
               IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
            }
            if(IBS0Count >= HomHomCount * 0.4) continue; // C_Hom < 60%
            // 3rd quarter
            for(m1 = stop1/2, m2 = stop1*3/4; m1 < m2; m1++){
               HomHomCount += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1]];
               IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
            }
            if(IBS0Count >= HomHomCount * 0.375) continue; // C_Hom < 62.5%
            m3 = (hetInOnePersonCount[id1][2] < hetInOnePersonCount[id2][2]?
                  hetInOnePersonCount[id1][2]: hetInOnePersonCount[id2][2]);
            HetHomCount = stop1 * 24 - HomHomCount * 2 - hetInOnePersonCount[id1][2] - hetInOnePersonCount[id2][2]
               - (missingInOnePersonCount[0][id1] + missingInOnePersonCount[0][id2])*2;
            threshold = (0.5-kincutoff1*0.5)*m3 - HetHomCount*0.25;
            if(IBS0Count >= threshold) continue;   // kinship < 0.0442

            // 4th quarter
            for(m1 = stop1*3/4; m1 < stop1; m1++){
               HomHomCount += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1]];
               IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
            }
            if(IBS0Count >= HomHomCount * 0.35) continue; // C_Hom < 65%
            m3 = (hetInOnePersonCount[id1][3] < hetInOnePersonCount[id2][3])?
                  hetInOnePersonCount[id1][3]: hetInOnePersonCount[id2][3];
            HetHomCount = stop1 * 32 - HomHomCount * 2 - hetInOnePersonCount[id1][3] - hetInOnePersonCount[id2][3]
               - missingInOnePersonCount[0][id1] - missingInOnePersonCount[0][id2];
            threshold = (0.5-kincutoff1*0.707107)*m3 - HetHomCount*0.25;   // kinship < 0.0625
            double cutoffDist = (0.5-kincutoff1)*m3 - HetHomCount*0.25;  // C_Hom < 70% && kinship < 0.0884
            if(HomHomCount * 0.3 > cutoffDist) cutoffDist = HomHomCount * 0.3;
            if(cutoffDist < threshold) threshold = cutoffDist;
            if(IBS0Count >= threshold) continue;

            int m1Max = missingWordInOnePersonCount[0][id1];
            int m2Max = missingWordInOnePersonCount[0][id2];
            for(HetHomCount = m1 = m2 = 0; m1 < m1Max; m1++){
               m3 = missingWordInOnePerson[id1][m1]; // m for id1
               for(; m2 < m2Max && missingWordInOnePerson[id2][m2] < m3; m2++);
               if( m2 == m2Max ) break;
               if(m3 == missingWordInOnePerson[id2][m2])   // id1 and id2 are both missing
                  HetHomCount += oneoneCount[~(GG[0][id1][m3] | GG[0][id2][m3] | GG[1][id1][m3] | GG[1][id2][m3])&0xFFFF];
            }  // Store MissMissCount in HetHomCount
            for(m1 = m3 = 0; m1 < m2Max; m1++){
               m2 = missingWordInOnePerson[id2][m1];
               m3 += oneoneCount[~(GG[0][id1][m2] | GG[0][id2][m2] | GG[1][id2][m2]) & GG[1][id1][m2]];
            }
            het1Count = hetInOnePersonCount[id1][3] - m3;
            for(m1 = m3 = 0; m1 < m1Max; m1++){
               m2 = missingWordInOnePerson[id1][m1];
               m3 += oneoneCount[~(GG[0][id1][m2] | GG[0][id2][m2] | GG[1][id1][m2]) & GG[1][id2][m2]];
            }
            het2Count = hetInOnePersonCount[id2][3] - m3;
            HetHomCount = stop1 * 32 - HomHomCount * 2 - het1Count - het2Count +
               (-missingInOnePersonCount[0][id1] - missingInOnePersonCount[0][id2] + HetHomCount) * 2;
            m3 = (het1Count < het2Count? het1Count: het2Count);
            threshold = (0.5-kincutoff1*0.707107)*m3 - HetHomCount*0.25;  // kinship < 0.0625
            cutoffDist = (0.5-kincutoff1)*m3 - HetHomCount*0.25;  // C_Hom < 75% && kinship < 0.0884
            if(HomHomCount * 0.7 > cutoffDist) cutoffDist = HomHomCount * 0.7;
            if(cutoffDist < threshold) threshold = cutoffDist;
            if(IBS0Count >= threshold) continue;
            rawrelativeCount ++;

            // Stage 2: few pairs
            for(HetHomCount -= (missingInOnePersonCount[1][id1]+missingInOnePersonCount[1][id2]),  // different non-missing genotypes
               m1 = start2; m1 < stop2; m1++) // different-genotpe Count
                  HetHomCount += oneoneCount[GG[0][id1][m1]^GG[0][id2][m1]];
            het1Count += hetInOnePersonCount[id1][4];
            het2Count += hetInOnePersonCount[id2][4];
            cutoffDist = (2-kincutoff2*4) * (het1Count < het2Count? het1Count: het2Count)-IBS0Count*4;
            if(HetHomCount >= cutoffDist) continue; // too big distance already
            cutoffDist = (cutoffDist+IBS0Count*4-HetHomCount)*0.25;
            for(m1 = start2; m1 < (start2+stop2)/2; m1++)
               IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
            if(IBS0Count >= cutoffDist * 1.1) continue;
            for(m1 = (start2+stop2)/2; (m1 < stop2) && (IBS0Count < cutoffDist); m1++)
               IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
            if(IBS0Count >= cutoffDist) continue;   // too big distance
            // now add the MissMissCount back
            int mm1Max = m1Max+missingWordInOnePersonCount[1][id1];
            int mm2Max = m2Max+missingWordInOnePersonCount[1][id2];
            for(m1 = m1Max, m2 = m2Max; m1 < mm1Max; m1++){
               m3 = missingWordInOnePerson[id1][m1]; // m for id1
               for(; m2 < mm2Max && missingWordInOnePerson[id2][m2] < m3; m2++);
               if( m2 == mm2Max ) break;
               if(m3 == missingWordInOnePerson[id2][m2])   // id1 and id2 are both missing
                  HetHomCount += 2 * oneoneCount[~(GG[0][id1][m3] | GG[0][id2][m3] | GG[1][id1][m3] | GG[1][id2][m3]) & 0xFFFF];
            }
            for(m1 = m2Max, m3 = 0; m1 < mm2Max; m1++){
               m2 = missingWordInOnePerson[id2][m1];
               m3 += oneoneCount[~(GG[0][id2][m2] | GG[0][id1][m2] | GG[1][id2][m2]) & GG[1][id1][m2]];
            }
            het1Count -= m3;
            HetHomCount += m3;
            for(m1 = m1Max, m3=0; m1 < mm1Max; m1++){
               m2 = missingWordInOnePerson[id1][m1];
               m3 += oneoneCount[~(GG[0][id1][m2] | GG[0][id2][m2] | GG[1][id1][m2]) & GG[1][id2][m2]];
            }
            het2Count -= m3;
            HetHomCount += m3;
            HetHetCount = (het1Count + het2Count - HetHomCount)/2;
            het1Count -= HetHetCount;
            het2Count -= HetHetCount;
            smaller = HetHetCount + (het1Count < het2Count? het1Count: het2Count);
            kinship = 0.5 - ((het1Count+het2Count)*0.25+IBS0Count)/smaller;
            if(kinship < kincutoff2) continue;

            midrelativeCount ++;
            // Stage 3: exact
            for(int m = start3; m < stop3; m++){
               HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
               IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
               het1Count += oneoneCount[GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m]]; // Het1Hom2
               het2Count += oneoneCount[GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]]; // Hom1Het2
            }
            smaller = HetHetCount + (het1Count < het2Count? het1Count: het2Count);
            kinship = 0.5 - ((het1Count+het2Count)*0.25+IBS0Count)/smaller;
            if(kinship < lowerbound) continue;
            relativeCount ++;
            notMissingCount = 0;
            for(int m = 0; m < shortCount; m++)
               notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
            sprintf(buffer, "%s\t%s\t%s\t%s\t%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\n",
               (const char*)ped[phenoid[id1]].famid, (const char*)ped[phenoid[id1]].pid,
               (const char*)ped[phenoid[id2]].famid, (const char*)ped[phenoid[id2]].pid,
               notMissingCount, HetHetCount*1.0/notMissingCount, IBS0Count*1.0/notMissingCount,
               HetHetCount*1.0/(het1Count+het2Count+HetHetCount), // CHet
               1-IBS0Count*1.0/(notMissingCount-het1Count-het2Count-HetHetCount), // CHom
               kinship);
            buffers[k/blockSize].Add(buffer);
         }
      }
   }
else // more distant than first-degree
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
      private(HomHomCount, HetHetCount, IBS0Count, het1Count, het2Count, notMissingCount, threshold,\
      HetHomCount, id1, id2, kinship, smaller, buffer, m1, m2, m3) \
      reduction(+:rawrelativeCount, midrelativeCount, relativeCount)
#endif
   for(int k = 0; k < loopIndexLength; k++){
      int i = loopIndex[0][k];
      int iMax = (i > idCount - LOOPBLOCKINGSIZE? idCount: i + LOOPBLOCKINGSIZE);
      int j = loopIndex[1][k];
      int jMax = (j > idCount - LOOPBLOCKINGSIZE? idCount: j + LOOPBLOCKINGSIZE);
      int jMin = j;
      for(id1 = i; id1 < iMax; id1++){
         if(missingInOnePersonCount[0][id1]>=cutoffMissingCount) continue;
         if(i == j) jMin = id1 + 1; // diagonal blocks only
         for(id2 = jMin; id2 < jMax; id2++){
            if(missingInOnePersonCount[0][id2]>=cutoffMissingCount) continue;
            if(ped[phenoid[id1]].famid == ped[phenoid[id2]].famid) continue;
            // Stage 1: all pairs
            if(relativedegree == 2){
               // 1st quarter
               for(HomHomCount = IBS0Count = m1 = 0, m2 = stop1/16; m1 < m2; m1++){
                  HomHomCount += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1]];
                  IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
               }
               if(IBS0Count >= HomHomCount * 0.5) continue; // C_Hom < 50%
               // 2nd quarter
               for(m1 = stop1/16, m2 = stop1/8; m1 < m2; m1++){
                  HomHomCount += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1]];
                  IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
               }
               if(IBS0Count >= HomHomCount * 0.45) continue; // C_Hom < 55%
               // 3rd quarter
               for(m1 = stop1/8, m2 = stop1/4; m1 < m2; m1++){
                  HomHomCount += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1]];
                  IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
               }
               if(IBS0Count >= HomHomCount * 0.42) continue; // C_Hom < 58%

               // 4th quarter
               for(m1 = stop1/4; m1 < stop1; m1++){
                  HomHomCount += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1]];
                  IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
               }
               if(IBS0Count >= HomHomCount * 0.375) continue; // C_Hom < 62.5%
            }else{   // more distant
               // 1st quarter
               for(HomHomCount = IBS0Count = m1 = 0, m2 = stop1/16; m1 < m2; m1++){
                  HomHomCount += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1]];
                  IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
               }
               if(IBS0Count >= HomHomCount * 0.55) continue; // C_Hom < 45%
               // 2nd quarter
               for(m1 = stop1/16, m2 = stop1/8; m1 < m2; m1++){
                  HomHomCount += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1]];
                  IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
               }
               if(IBS0Count >= HomHomCount * 0.5) continue; // C_Hom < 50%
               // 3rd quarter
               for(m1 = stop1/8, m2 = stop1/4; m1 < m2; m1++){
                  HomHomCount += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1]];
                  IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
               }
               if(IBS0Count >= HomHomCount * 0.45) continue; // C_Hom < 55%

               // 4th quarter
               for(m1 = stop1/4; m1 < stop1; m1++){
                  HomHomCount += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1]];
                  IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
               }
               if(IBS0Count >= HomHomCount * 0.45) continue; // C_Hom < 55%
            }

            m3 = (hetInOnePersonCount[id1][3] < hetInOnePersonCount[id2][3])?
                  hetInOnePersonCount[id1][3]: hetInOnePersonCount[id2][3];
            HetHomCount = stop1 * 32 - HomHomCount * 2 - hetInOnePersonCount[id1][3] - hetInOnePersonCount[id2][3]
               - missingInOnePersonCount[0][id1] - missingInOnePersonCount[0][id2];
            threshold = (0.5-kincutoff1)*m3 - HetHomCount*0.25;   // kinship < 0.0625
            if(IBS0Count >= threshold) continue;

            int m1Max = missingWordInOnePersonCount[0][id1];
            int m2Max = missingWordInOnePersonCount[0][id2];
            for(HetHomCount = m1 = m2 = 0; m1 < m1Max; m1++){
               m3 = missingWordInOnePerson[id1][m1]; // m for id1
               for(; m2 < m2Max && missingWordInOnePerson[id2][m2] < m3; m2++);
               if( m2 == m2Max ) break;
               if(m3 == missingWordInOnePerson[id2][m2])   // id1 and id2 are both missing
                  HetHomCount += oneoneCount[~(GG[0][id1][m3] | GG[0][id2][m3] | GG[1][id1][m3] | GG[1][id2][m3])&0xFFFF];
            }  // Store MissMissCount in HetHomCount
            for(m1 = m3 = 0; m1 < m2Max; m1++){
               m2 = missingWordInOnePerson[id2][m1];
               m3 += oneoneCount[~(GG[0][id1][m2] | GG[0][id2][m2] | GG[1][id2][m2]) & GG[1][id1][m2]];
            }
            het1Count = hetInOnePersonCount[id1][3] - m3;
            for(m1 = m3 = 0; m1 < m1Max; m1++){
               m2 = missingWordInOnePerson[id1][m1];
               m3 += oneoneCount[~(GG[0][id1][m2] | GG[0][id2][m2] | GG[1][id1][m2]) & GG[1][id2][m2]];
            }
            het2Count = hetInOnePersonCount[id2][3] - m3;
            HetHomCount = stop1 * 32 - HomHomCount * 2 - het1Count - het2Count +
               (-missingInOnePersonCount[0][id1] - missingInOnePersonCount[0][id2] + HetHomCount) * 2;
            m3 = (het1Count < het2Count? het1Count: het2Count);
            threshold = (0.5-kincutoff1)*m3 - HetHomCount*0.25;  // kinship < 0.0625
            if(IBS0Count >= threshold) continue;

            rawrelativeCount ++;

            // Stage 2: few pairs
            for(HetHomCount -= (missingInOnePersonCount[1][id1]+missingInOnePersonCount[1][id2]),  // different non-missing genotypes
               m1 = start2; m1 < stop2; m1++) // different-genotpe Count
                  HetHomCount += oneoneCount[GG[0][id1][m1]^GG[0][id2][m1]];
            het1Count += hetInOnePersonCount[id1][4];
            het2Count += hetInOnePersonCount[id2][4];
            double cutoffDist = (2-kincutoff2*4) * (het1Count < het2Count? het1Count: het2Count)-IBS0Count*4;
            if(HetHomCount >= cutoffDist) continue; // too big distance already
            cutoffDist = (cutoffDist+IBS0Count*4-HetHomCount)*0.25;
            for(m1 = start2; m1 < (start2+stop2)/2; m1++)
               IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
            if(IBS0Count >= cutoffDist * 1.1) continue;
            for(m1 = (start2+stop2)/2; (m1 < stop2) && (IBS0Count < cutoffDist); m1++)
               IBS0Count += oneoneCount[GG[0][id1][m1] & GG[0][id2][m1] & (GG[1][id1][m1] ^ GG[1][id2][m1])];
            if(IBS0Count >= cutoffDist) continue;   // too big distance
            // now add the MissMissCount back
            int mm1Max = m1Max+missingWordInOnePersonCount[1][id1];
            int mm2Max = m2Max+missingWordInOnePersonCount[1][id2];
            for(m1 = m1Max, m2 = m2Max; m1 < mm1Max; m1++){
               m3 = missingWordInOnePerson[id1][m1]; // m for id1
               for(; m2 < mm2Max && missingWordInOnePerson[id2][m2] < m3; m2++);
               if( m2 == mm2Max ) break;
               if(m3 == missingWordInOnePerson[id2][m2])   // id1 and id2 are both missing
                  HetHomCount += 2 * oneoneCount[~(GG[0][id1][m3] | GG[0][id2][m3] | GG[1][id1][m3] | GG[1][id2][m3]) & 0xFFFF];
            }
            for(m1 = m2Max, m3 = 0; m1 < mm2Max; m1++){
               m2 = missingWordInOnePerson[id2][m1];
               m3 += oneoneCount[~(GG[0][id2][m2] | GG[0][id1][m2] | GG[1][id2][m2]) & GG[1][id1][m2]];
            }
            het1Count -= m3;
            HetHomCount += m3;
            for(m1 = m1Max, m3=0; m1 < mm1Max; m1++){
               m2 = missingWordInOnePerson[id1][m1];
               m3 += oneoneCount[~(GG[0][id1][m2] | GG[0][id2][m2] | GG[1][id1][m2]) & GG[1][id2][m2]];
            }
            het2Count -= m3;
            HetHomCount += m3;
            HetHetCount = (het1Count + het2Count - HetHomCount)/2;
            het1Count -= HetHetCount;
            het2Count -= HetHetCount;
            smaller = HetHetCount + (het1Count < het2Count? het1Count: het2Count);
            kinship = 0.5 - ((het1Count+het2Count)*0.25+IBS0Count)/smaller;
            if(kinship < kincutoff2) continue;
            midrelativeCount ++;
            // Stage 3: exact
            for(int m = start3; m < stop3; m++){
               HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
               IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
               het1Count += oneoneCount[GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m]]; // Het1Hom2
               het2Count += oneoneCount[GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]]; // Hom1Het2
            }
            smaller = HetHetCount + (het1Count < het2Count? het1Count: het2Count);
            kinship = 0.5 - ((het1Count+het2Count)*0.25+IBS0Count)/smaller;
            if(kinship < lowerbound) continue;
            relativeCount ++;
            notMissingCount = 0;
            for(int m = 0; m < shortCount; m++)
               notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
            sprintf(buffer, "%s\t%s\t%s\t%s\t%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\n",
               (const char*)ped[phenoid[id1]].famid, (const char*)ped[phenoid[id1]].pid,
               (const char*)ped[phenoid[id2]].famid, (const char*)ped[phenoid[id2]].pid,
               notMissingCount, HetHetCount*1.0/notMissingCount, IBS0Count*1.0/notMissingCount,
               HetHetCount*1.0/(het1Count+het2Count+HetHetCount), // CHet
               1-IBS0Count*1.0/(notMissingCount-het1Count-het2Count-HetHetCount), // CHom
               kinship);
            buffers[k/blockSize].Add(buffer);
         }
      }
   }

   for(int i = 0; i < idCount; i++)
      if(missingWordInOnePerson[i]) delete []missingWordInOnePerson[i];
   for(int i = 0; i < 2; i++){
      if(missingWordInOnePersonCount[i]) delete []missingWordInOnePersonCount[i];
      if(missingInOnePersonCount[i]) delete []missingInOnePersonCount[i];
   }
   for(int i = 0; i < idCount; i++)
      delete []hetInOnePersonCount[i];
   delete []hetInOnePersonCount;


   outfile.Copy(prefix);
   outfile.Add(".kin0");
   fp = fopen(outfile, "wt");
   fprintf(fp, "FID1\tID1\tFID2\tID2\tN_SNP\tHetHet\tIBS0\tHetConc\tHomConc\tKinship\n");
   for(int b = 0; b < buffers.Length(); b++)
      buffers[b].Write(fp);
   fclose(fp);
   printf("                                         ends at %s", currentTime());
   if(relativedegree){
      printf("Relationship inference speed-up by screening top informative SNPs:\n");
      printf("  Stage 1 (with %d SNPs): %lld relative pairs detected (with Kinship > %.5lf)\n",
         stop1*16>markerCount? markerCount: stop1*16,
         rawrelativeCount, kincutoff1);
      printf("  Stage 2 (with %d SNPs): %ld out of %lld remain related (with Kinship > %.5lf)\n",
         stop2*16>markerCount? markerCount: stop2*16,
         midrelativeCount, rawrelativeCount, kincutoff2);
      printf("  Final Stage (with %d SNPs): %ld pairs of relatives up to %d-degree are identified\n",
         markerCount, relativeCount, relativedegree);
      printf("Between-family relatives (kinship >= %.5lf) saved in file %s\n",
         lowerbound, (const char*)outfile);
      if(relativedegree == 1){
         printf("\nNote only duplicates and 1st-degree relative pairs are included here.\n");
         printf("  Specifying '--degree 2' for a higher degree relationship inference.\n\n");
      }
   }
}
  */
  

/*
               for(m1 = 4; m1 < long_prescan; m1++) // Lower bound of HetHom Count
                  for(word = SLG[0][id1][m1]^SLG[0][id2][m1]; word; word &= (word-1), HetHomCount++);
               if(HetHomCount >= cutoff) continue; // C_het <= 0.75
*/

/*
               for(word2 = 0, m1 = 4; m1 < long_prescan; m1++){ // HetHom Count
                  word = SLG[0][id1][m1]^SLG[0][id2][m1];
                  word = word - ((word>>1)&0x5555555555555555);
                  word2 += (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
               }// Rare chance of hitting LB, when all 4 words are 1111
               word2 = (word2+(word2>>4)) & 0x0F0F0F0F0F0F0F0F;
               word2 = (word2+(word2>>8)) & 0x00FF00FF00FF00FF;
               word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;
               HetHomCount += (word2+(word2>>32)) & 0xFFFFFFFF;
               if(HetHomCount >= cutoff) continue; // C_het <= 0.75
*/
                  /*
               for(m1 = 1; m1 < 4; m1++) // Lower bound of HetHom Count
                  for(word = SLG[0][id1][m1]^SLG[0][id2][m1]; word; word &= (word-1), HetHomCount++);
               if(HetHomCount*2 >= cutoff) continue; // C_het <= 0.75
               for(m1 = 4; m1 < long_prescan; m1++) // Lower bound of HetHom Count
                  for(word = SLG[0][id1][m1]^SLG[0][id2][m1]; word; word &= (word-1), HetHomCount++);
               if(HetHomCount >= cutoff) continue; // C_het <= 0.75
                    */

/*
               word = SLG[0][id1][0]^SLG[0][id2][0];
               word = word - ((word>>1)&0x5555555555555555);
               word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
               word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
               word = (word+(word>>8)) & 0x00FF00FF00FF00FF;
               word = (word+(word>>16)) & 0x0000FFFF0000FFFF;
               HetHomCount = (word+(word>>32)) & 0xFFFFFFFF - m2;
               if(HetHomCount*6 >= cutoff) continue; // C_het <= 0.75
               for(m1 = 1, word2 = 0; m1 < 4; m1++){ // Lower bound of HetHom Count
                  word = SLG[0][id1][m1]^SLG[0][id2][m1];
                  word = word - ((word>>1)&0x5555555555555555);
                  word2 += (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
               }
               word2 = (word2+(word2>>4)) & 0x0F0F0F0F0F0F0F0F;
               word2 = (word2+(word2>>8)) & 0x00FF00FF00FF00FF;
               word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;
               HetHomCount += (word2+(word2>>32)) & 0xFFFFFFFF + m2 -
                  missingInOnePersonCount[id1][1] - missingInOnePersonCount[id2][1];
               if(HetHomCount*2 >= cutoff) continue; // C_het <= 0.75
               for(word2 = 0, m1 = 4; m1 < long_prescan; m1++){ // HetHom Count
                  word = SLG[0][id1][m1]^SLG[0][id2][m1];
                  word = word - ((word>>1)&0x5555555555555555);
                  word2 += (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
               }// Rare chance of hitting LB, when all 4 words are 1111
               word2 = (word2+(word2>>4)) & 0x0F0F0F0F0F0F0F0F;
               word2 = (word2+(word2>>8)) & 0x00FF00FF00FF00FF;
               word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;
               HetHomCount += (word2+(word2>>32)) & 0xFFFFFFFF;
               if(HetHomCount >= cutoff) continue; // C_het <= 0.75
*/

/*
               for(word = SLG[0][id1][0]^SLG[0][id2][0], HetHomCount = -m2*6;
                  word && HetHomCount < cutoff; word &= (word-1), HetHomCount += 6);
               if(HetHomCount >= cutoff) continue;
               HetHomCount = HetHomCount/6 + m2 -
                  missingInOnePersonCount[id1][1] - missingInOnePersonCount[id2][1];
*/
/*
               for(word2 = m1 = 0; m1 < long_prescan; m1++){   // IBS0
                  word = SLG[0][id1][m1] & SLG[0][id2][m1] & (SLG[1][id1][m1] ^ SLG[1][id2][m1]);
                  word = word - ((word>>1)&0x5555555555555555);
                  word2 += (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
               }  // Rare chance of hitting LB
               word2 = (word2+(word2>>4)) & 0x0F0F0F0F0F0F0F0F;
               word2 = (word2+(word2>>8)) & 0x00FF00FF00FF00FF;
               word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;
               m3 = (word2+(word2>>32)) & 0xFFFFFFFF;   // IBS0Count
*/

//   IntArray tArray;
//   unsigned int **missingWordInOnePerson = new unsigned int *[idCount];
//   int *missingWordInOnePersonCount;
//   missingWordInOnePersonCount = new int[idCount];

//      if(tArray.Length()){
//         missingWordInOnePerson[i] = new unsigned int [tArray.Length()];
//         for(int m = 0; m < tArray.Length(); m++)
//            missingWordInOnePerson[i][m] = tArray[m];
//      }else
//         missingWordInOnePerson[i] = NULL;

/*
   if(faster > 1){
      if(faster > 4) faster = 4;
      short_prescan /= faster;
   }
   if(slower > 1){
      short_prescan *= slower;
      if(short_prescan > shortCount) short_prescan = shortCount;
   }
*/

/*
               int m1Max = missingWordInOnePersonCount[id1];
               int m2Max = missingWordInOnePersonCount[id2];
               for(HomHomCount = m1 = m2 = 0; m1 < m1Max; m1++){
                  m3 = missingWordInOnePerson[id1][m1]; // m for id1
                  for(; m2<m2Max && missingWordInOnePerson[id2][m2]<m3; m2++);
                  if( m2 == m2Max ) break;
                  if(m3 == missingWordInOnePerson[id2][m2])   // id1 and id2 are both missing
                     HomHomCount += oneoneCount[~(SG[0][id1][m3] | SG[0][id2][m3] | SG[1][id1][m3] | SG[1][id2][m3]) & 0xFFFF];
               }
               if(HetHomCount+HomHomCount*2 >= cutoff) continue; // C_het <= 0.75
               for(m3 = m1 = 0; m1 < m2Max; m1++){
                  m2 = missingWordInOnePerson[id2][m1];
                  m3 += oneoneCount[~(SG[0][id2][m2] | SG[0][id1][m2] | SG[1][id2][m2]) & SG[1][id1][m2]];
               }
               for(m1 = 0; m1 < m1Max; m1++){
                  m2 = missingWordInOnePerson[id1][m1];
                  m3 += oneoneCount[~(SG[0][id1][m2] | SG[0][id2][m2] | SG[1][id1][m2]) & SG[1][id2][m2]];
               }
*/

/*
   delete []missingInOnePersonCount;
   delete []hetInOnePersonCount;
   delete []missingWordInOnePersonCount;
   for(int i = 0; i < idCount; i++)
      if(missingWordInOnePerson[i]) delete []missingWordInOnePerson[i];
*/

/*
         for(word2=m1=0; m1 < longCount; m1++){
            word = LG[0][id1][m1] & LG[0][id2][m1];   // HomHom
            word = word - ((word>>1)&0x5555555555555555);
            word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
            word = (word + (word>>4)) & 0x0F0F0F0F0F0F0F0F;
            word = (word + (word>>8)) & 0x00FF00FF00FF00FF;
            word2 += (word + (word>>16)) & 0x0000FFFF0000FFFF;
         }
         HomHomCount = (word2 + (word2>>32)) & 0xFFFFFFFF;
*/

/*
         for(word2=m1=0; m1 < longCount; m1++){ // HetHet
            word = ~(LG[0][id1][m1] | LG[0][id2][m1]) & LG[1][id1][m1] & LG[1][id2][m1];
            word = word - ((word>>1)&0x5555555555555555);
            word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
            word = (word + (word>>4)) & 0x0F0F0F0F0F0F0F0F;
            word = (word + (word>>8)) & 0x00FF00FF00FF00FF;
            word2 += (word + (word>>16)) & 0x0000FFFF0000FFFF;
         }
         HetHetCount = (word2 + (word2>>32)) & 0xFFFFFFFF;
*/


/*
         for(word2=m1=0; m1 < longCount; m1++){ // IBS0
            word = LG[0][id1][m1] & LG[0][id2][m1] & (LG[1][id1][m1] ^ LG[1][id2][m1]);
            word = word - ((word>>1)&0x5555555555555555);
            word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
            word = (word + (word>>4)) & 0x0F0F0F0F0F0F0F0F;
            word = (word + (word>>8)) & 0x00FF00FF00FF00FF;
            word2 += (word + (word>>16)) & 0x0000FFFF0000FFFF;
         }
         DiffHomCount = (word2 + (word2>>32)) & 0xFFFFFFFF;
*/

/*
         HetHetCount = notMissingHetCount = HomHomCount = DiffHomCount = 0;
         for(m1 = 0; m1 < longCount; m1++){
            word = ~(LG[0][id1][m1] | LG[0][id2][m1]) & LG[1][id1][m1] & LG[1][id2][m1];  // HetHet
            HetHetCount += oneoneCount[word & 0xFFFF] +
                           oneoneCount[(word>>16) & 0xFFFF] +
                           oneoneCount[(word>>32) & 0xFFFF] +
                           oneoneCount[(word>>48) & 0xFFFF];
            word = (LG[0][id1][m1] | LG[1][id1][m1]) &   // HetNomiss
            (LG[0][id2][m1] | LG[1][id2][m1]) & ~(LG[0][id1][m1] & LG[0][id2][m1]);
            notMissingHetCount += oneoneCount[word & 0xFFFF] +
                           oneoneCount[(word>>16) & 0xFFFF] +
                           oneoneCount[(word>>32) & 0xFFFF] +
                           oneoneCount[(word>>48) & 0xFFFF];
         }
         if(HetHetCount <= notMissingHetCount*0.8) continue; // not duplicates, heterozygote concordance rate <= 80%
         for(m1 = 0; m1 < longCount; m1++){
            word = LG[0][id1][m1] & LG[0][id2][m1];   // HomHom
            HomHomCount += oneoneCount[word & 0xFFFF] +
                           oneoneCount[(word>>16) & 0xFFFF] +
                           oneoneCount[(word>>32) & 0xFFFF] +
                           oneoneCount[(word>>48) & 0xFFFF];
            word = LG[0][id1][m1] & LG[0][id2][m1] & (LG[1][id1][m1] ^ LG[1][id2][m1]);   // IBS0
            DiffHomCount += oneoneCount[word & 0xFFFF] +
                           oneoneCount[(word>>16) & 0xFFFF] +
                           oneoneCount[(word>>32) & 0xFFFF] +
                           oneoneCount[(word>>48) & 0xFFFF];
         }
*/


/*
   int m1, m2, m3;
//   IntArray tArray;
   int *missingInOnePersonCount = new int[idCount];
   for(int i = 0; i < idCount; i++)
      missingInOnePersonCount[i] = 0;
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
   */

//#include "KingCore.h"
//#include "KinshipX.h"
//#include "MathStats.h"
//#include "MathSVD.h"
//#include "QuickIndex.h"
//#include "MathCholesky.h"


