//////////////////////////////////////////////////////////////////////
// integrated.cpp
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
// Dec 20, 2018

#include <math.h>
#include "analysis.h"
#include "Kinship.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

void Engine::ScreenCloseRelativesInSubset64Bit(IntArray rpList[])
{  // --related
   if(longCount < 65 || !bigdataFlag) return;  // screening is not needed
   int long_prescan = 64;
   int stop1 = long_prescan;
//   int stop2 = stop1 + (long_prescan<<3);
   int stop2 = (long_prescan<<3);   // 4096*8 = 32768 SNPs
   double kincutoff2 = 0.125;
   double kincutoff1 = kincutoff2 * 0.7071068;
   if(relativedegree==2){
      stop1 <<= 3;
//      stop2 <<= 3;
      kincutoff1 *= 0.5;   // 0.0442
      kincutoff2 *= 0.5;   // 0.0625
   }
   if(stop1 > longCount) stop1 = longCount;
   if(stop2 > longCount) stop2 = longCount;
   int m1, m2, m3;
   unsigned long long int word, word1, word2, ibs0;
   unsigned short int **missingInOnePersonCount = new unsigned short int * [idCount];
   unsigned short int **hetInOnePersonCount = new unsigned short int * [idCount];
   if(relativedegree == 1)
      for(int i = 0; i < idCount; i++){
         missingInOnePersonCount[i] = new unsigned short int [2];
         hetInOnePersonCount[i] = new unsigned short int [2];
      }
   else  // degree 2
      for(int i = 0; i < idCount; i++){
         missingInOnePersonCount[i] = new unsigned short int [1];
         hetInOnePersonCount[i] = new unsigned short int [1];
      }
   unsigned long long int mask = markerCount%64? (((unsigned long long int)1 << (markerCount % 64))-1): 0xFFFFFFFFFFFFFFFF;
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
      private(m1, m3, word, word2)
#endif
   for(int i = 0; i < idCount; i++){
      for(m3 = m1 = 0; m1 < stop1; m1++)   // not all non-missing
         for(word = ~(SLG[0][i][m1] | SLG[1][i][m1]); word; word &= (word-1), m3++);
      missingInOnePersonCount[i][0] = m3;
      if(relativedegree == 1){
         for(m3 = 0, m1 = stop1; m1 < stop2-1; m1++)   // not all non-missing
            for(word = ~(SLG[0][i][m1] | SLG[1][i][m1]); word; word &= (word-1), m3++);
         for(word = ~(SLG[0][i][stop2-1] | SLG[1][i][stop2-1]) & mask; word; word &= (word-1), m3++);
         missingInOnePersonCount[i][1] = m3;
      }
      for(word2 = m1 = 0; m1 < stop1; m1++){
         word = (~SLG[0][i][m1]) & SLG[1][i][m1];  // Het
         word = word - ((word>>1)&0x5555555555555555);
         word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
         word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
         word2 += (word+(word>>8)) & 0x00FF00FF00FF00FF;
      }
      word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;
      hetInOnePersonCount[i][0] = (word2+(word2>>32)) & 0xFFFFFFFF;
      if(relativedegree == 1){
         for(word2 = 0, m1 = stop1; m1 < stop2; m1++){
            word = (~SLG[0][i][m1]) & SLG[1][i][m1];  // Het
            word = word - ((word>>1)&0x5555555555555555);
            word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
            word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
            word = (word+(word>>8)) & 0x00FF00FF00FF00FF;
            word2 += (word+(word>>16)) & 0x0000FFFF0000FFFF;
         }
         hetInOnePersonCount[i][1] = (word2+(word2>>32)) & 0xFFFFFFFF;
      }
   }  // parallel among individuals ends
   for(int i = 0; i < defaultMaxCoreCount; i++)
      rpList[i].Dimension(0);
   double threshold;
   int id1, id2, IBS0Count, het1Count, het2Count, HetHomCount, HomHomCount;
   int LOOPBLOCKINGSIZE=(relativedegree==1?32:16);
   char SHIFTSIZE=(relativedegree==1?5:4);
   if(idCount > (1<<20)){  // sample size over a million
      int count = (idCount >> 20);
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
if(relativedegree==1){
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(HomHomCount, IBS0Count, HetHomCount, het1Count, het2Count,\
      threshold, id1, id2, m1, m2, m3, thread, word, word1, word2, ibs0)
{
   thread = omp_get_thread_num();
   #pragma omp for
#endif
   for(unsigned long int k = 0; k < loopIndexLength; k++){
      if(k < firstCount && (k % onepercent) == 0) {
         printf("%d%%\r", k/onepercent);
         fflush(stdout);
      }
      int i = (loopIndex[k]>>16)<<SHIFTSIZE;
      int iMax = (i > idCount - LOOPBLOCKINGSIZE? idCount: i + LOOPBLOCKINGSIZE);
      int j = (loopIndex[k]&0xFFFF)<<SHIFTSIZE;
      int jMax = (j > idCount - LOOPBLOCKINGSIZE? idCount: j + LOOPBLOCKINGSIZE);
      for(id1 = i; id1 < iMax; id1++){
         for(id2 = j; id2 < jMax; id2++){
            // Stage 1: all pairs
            // Stop rule 1: C_Hom < 65%
            // Stop rule 2: kinship < 0.0625
            // Stop rule 3: kinship < 0.0884 && C_Hom < 70%
            // 1st quarter
            word1 = word2 = 0;
            for(int m = 0; m < 16; m++){
               word = SLG[0][id1][m] & SLG[0][id2][m]; // HomHom
               ibs0 = word & (SLG[1][id1][m] ^ SLG[1][id2][m]);   // IBS0
               word = word - ((word>>1)&0x5555555555555555);
               word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
               word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
               word1 += (word+(word>>8)) & 0x00FF00FF00FF00FF;
               ibs0 = ibs0 - ((ibs0>>1)&0x5555555555555555);   // IBS0
               ibs0 = (ibs0&0x3333333333333333) + ((ibs0>>2)&0x3333333333333333);
               ibs0 = (ibs0+(ibs0>>4)) & 0x0F0F0F0F0F0F0F0F;
               word2 += (ibs0+(ibs0>>8)) & 0x00FF00FF00FF00FF;
            }
            word1 = (word1+(word1>>16)) & 0x0000FFFF0000FFFF;
            HomHomCount = (word1+(word1>>32)) & 0xFFFFFFFF;
            word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;
            IBS0Count = (word2+(word2>>32)) & 0xFFFFFFFF;
            if(IBS0Count >= HomHomCount * 0.5) continue; // C_Hom < 50%
            // 2nd quarter
            word1 = word2 = 0;
            for(int m = 16; m < 32; m++){
               word = SLG[0][id1][m] & SLG[0][id2][m]; // HomHom
               ibs0 = word & (SLG[1][id1][m] ^ SLG[1][id2][m]);   // IBS0
               word = word - ((word>>1)&0x5555555555555555);
               word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
               word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
               word1 += (word+(word>>8)) & 0x00FF00FF00FF00FF;
               ibs0 = ibs0 - ((ibs0>>1)&0x5555555555555555);   // IBS0
               ibs0 = (ibs0&0x3333333333333333) + ((ibs0>>2)&0x3333333333333333);
               ibs0 = (ibs0+(ibs0>>4)) & 0x0F0F0F0F0F0F0F0F;
               word2 += (ibs0+(ibs0>>8)) & 0x00FF00FF00FF00FF;
            }
            word1 = (word1+(word1>>16)) & 0x0000FFFF0000FFFF;
            HomHomCount += (word1+(word1>>32)) & 0xFFFFFFFF;
            word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;
            IBS0Count += (word2+(word2>>32)) & 0xFFFFFFFF;
            if(IBS0Count >= HomHomCount * 0.4) continue; // C_Hom < 60%
            // 3rd quarter
            word1 = word2 = 0;
            for(int m = 32; m < 48; m++){
               word = SLG[0][id1][m] & SLG[0][id2][m]; // HomHom
               ibs0 = word & (SLG[1][id1][m] ^ SLG[1][id2][m]);   // IBS0
               word = word - ((word>>1)&0x5555555555555555);
               word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
               word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
               word1 += (word+(word>>8)) & 0x00FF00FF00FF00FF;
               ibs0 = ibs0 - ((ibs0>>1)&0x5555555555555555);   // IBS0
               ibs0 = (ibs0&0x3333333333333333) + ((ibs0>>2)&0x3333333333333333);
               ibs0 = (ibs0+(ibs0>>4)) & 0x0F0F0F0F0F0F0F0F;
               word2 += (ibs0+(ibs0>>8)) & 0x00FF00FF00FF00FF;
            }
            word1 = (word1+(word1>>16)) & 0x0000FFFF0000FFFF;
            HomHomCount += (word1+(word1>>32)) & 0xFFFFFFFF;
            word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;
            IBS0Count += (word2+(word2>>32)) & 0xFFFFFFFF;
            if(IBS0Count >= HomHomCount * 0.375) continue; // C_Hom < 62.5%

            // 4th quarter
            word1 = word2 = 0;
            for(int m = 48; m < stop1; m++){
               word = SLG[0][id1][m] & SLG[0][id2][m]; // HomHom
               ibs0 = word & (SLG[1][id1][m] ^ SLG[1][id2][m]);   // IBS0
               word = word - ((word>>1)&0x5555555555555555);
               word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
               word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
               word1 += (word+(word>>8)) & 0x00FF00FF00FF00FF;
               ibs0 = ibs0 - ((ibs0>>1)&0x5555555555555555);   // IBS0
               ibs0 = (ibs0&0x3333333333333333) + ((ibs0>>2)&0x3333333333333333);
               ibs0 = (ibs0+(ibs0>>4)) & 0x0F0F0F0F0F0F0F0F;
               word2 += (ibs0+(ibs0>>8)) & 0x00FF00FF00FF00FF;
            }
            word1 = (word1+(word1>>16)) & 0x0000FFFF0000FFFF;
            HomHomCount += (word1+(word1>>32)) & 0xFFFFFFFF;
            word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;
            IBS0Count += (word2+(word2>>32)) & 0xFFFFFFFF;
            if(IBS0Count >= HomHomCount * 0.35) continue; // C_Hom < 65%

            if(i==j && id2 <= id1) continue;

            if(IBS0Count < HomHomCount * 0.01){   // C_Hom > 99%
               rpList[thread].Push(id1);
               rpList[thread].Push(id2);
               continue;
            }

            het1Count = hetInOnePersonCount[id1][0];
            het2Count = hetInOnePersonCount[id2][0];
            HetHomCount = (stop1 << 7) - het1Count - het2Count -
               ((HomHomCount + missingInOnePersonCount[id1][0] + missingInOnePersonCount[id2][0])<<1);
//            m3 = (het1Count < het2Count? het1Count: het2Count);
//            threshold = (0.5-0.125)*m3 - HetHomCount*0.25;   // kinship < 0.0625
            if((IBS0Count >= HomHomCount * 0.2) && (het1Count+het2Count < HetHomCount*1.857143))
                 continue; // pass if CHet < 0.3 && CHom < 0.8
            m1 = m2 = m3 = 0;
            for(int m = 0; m < stop1; m++){
               word1 = ~(SLG[0][id1][m] | SLG[1][id1][m]);
               word2 = ~(SLG[0][id2][m] | SLG[1][id2][m]);
               for(word = word1 & word2; word; word &= (word-1), m3++);
               for(word = (~SLG[0][id1][m]) & SLG[1][id1][m] & word2; word; word &= (word-1), m1++);
               for(word = word1 & (~SLG[0][id2][m]) & SLG[1][id2][m]; word; word &= (word-1), m2++);
            }
            het1Count -= m1;
            het2Count -= m2;
            HetHomCount += m1 + m2 + (m3<<1);
//            m3 = (het1Count < het2Count? het1Count: het2Count);
//            threshold = (0.5-0.125)*m3 - HetHomCount*0.25;  // kinship < 0.0625
            if((IBS0Count >= HomHomCount * 0.2) && (het1Count+het2Count < HetHomCount*1.857143))
               continue;   // pass if CHet < 0.3 AND CHom < 0.8

            // Stage 2: few pairs
            word1 = word2 = 0;
            for(int m = stop1; m < stop2; m++){
               word = SLG[0][id1][m] & SLG[0][id2][m]; // HomHom
               ibs0 = word & (SLG[1][id1][m] ^ SLG[1][id2][m]);   // IBS0
               word = word - ((word>>1)&0x5555555555555555);
               word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
               word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
               word1 += (word+(word>>8)) & 0x00FF00FF00FF00FF;
               ibs0 = ibs0 - ((ibs0>>1)&0x5555555555555555);   // IBS0
               ibs0 = (ibs0&0x3333333333333333) + ((ibs0>>2)&0x3333333333333333);
               ibs0 = (ibs0+(ibs0>>4)) & 0x0F0F0F0F0F0F0F0F;
               word2 += (ibs0+(ibs0>>8)) & 0x00FF00FF00FF00FF;
            }
            word1 = (word1+(word1>>16)) & 0x0000FFFF0000FFFF;
// No, this is not a typo
            HomHomCount = (word1+(word1>>32)) & 0xFFFFFFFF;
            word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;
            IBS0Count += (word2+(word2>>32)) & 0xFFFFFFFF;
            het1Count += hetInOnePersonCount[id1][1];
            het2Count += hetInOnePersonCount[id2][1];
            HetHomCount += ((stop2-stop1) << 7) - hetInOnePersonCount[id1][1] - hetInOnePersonCount[id2][1]
               - ((HomHomCount + missingInOnePersonCount[id1][1] + missingInOnePersonCount[id2][1])<<1);
            m3 = (het1Count < het2Count)? het1Count: het2Count;
            threshold = (0.5-kincutoff2)*m3 - HetHomCount*0.25;   // kinship < 0.125
            if(IBS0Count >= threshold) continue;
            m1 = m2 = m3 = 0;
            for(int m = stop1; m < stop2; m++){
               word1 = ~(SLG[0][id1][m] | SLG[1][id1][m]);
               word2 = ~(SLG[0][id2][m] | SLG[1][id2][m]);
               for(word = word1 & word2; word; word &= (word-1), m3++);
               for(word = (~SLG[0][id1][m]) & SLG[1][id1][m] & word2; word; word &= (word-1), m1++);
               for(word = word1 & (~SLG[0][id2][m]) & SLG[1][id2][m]; word; word &= (word-1), m2++);
            }
            het1Count -= m1;
            het2Count -= m2;
            HetHomCount += m1 + m2 + (m3<<1);
            m3 = (het1Count < het2Count? het1Count: het2Count);
            threshold = (0.5-kincutoff2)*m3 - HetHomCount*0.25;  // kinship < 0.125
            if(IBS0Count >= threshold) continue;
            rpList[thread].Push(id1);
            rpList[thread].Push(id2);
         }  // end of id2
      }  // end of id1
   }  // end of loop blocking
#ifdef _OPENMP
}
#endif
}else{   // 2nd-degree
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(HomHomCount, IBS0Count, HetHomCount, het1Count, het2Count,\
      threshold, id1, id2, m1, m2, m3, thread, word, word1, word2, ibs0)
{
   thread = omp_get_thread_num();
   #pragma omp for
#endif
   for(unsigned long int k = 0; k < loopIndexLength; k++){
      if(k < firstCount && (k % onepercent) == 0) {
         printf("%d%%\r", k/onepercent);
         fflush(stdout);
      }
      int i = (loopIndex[k]>>16)<<SHIFTSIZE;
      int iMax = (i > idCount - LOOPBLOCKINGSIZE? idCount: i + LOOPBLOCKINGSIZE);
      int j = (loopIndex[k]&0xFFFF)<<SHIFTSIZE;
      int jMax = (j > idCount - LOOPBLOCKINGSIZE? idCount: j + LOOPBLOCKINGSIZE);
      for(id1 = i; id1 < iMax; id1++){
         for(id2 = j; id2 < jMax; id2++){
            // Stage 1 only: all pairs
            // Stop rule 1: C_Hom < 62.5%
            // Stop rule 2: kinship < 0.0625
            // 1st quarter
            word1 = word2 = 0;
            m2 = (stop1>>4);  // 2098/64=32 words
            for(int m = 0; m < m2; m++){
               word = SLG[0][id1][m] & SLG[0][id2][m]; // HomHom
               ibs0 = word & (SLG[1][id1][m] ^ SLG[1][id2][m]);   // IBS0
               word = word - ((word>>1)&0x5555555555555555);
               word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
               word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
               word1 += (word+(word>>8)) & 0x00FF00FF00FF00FF;
               ibs0 = ibs0 - ((ibs0>>1)&0x5555555555555555);   // IBS0
               ibs0 = (ibs0&0x3333333333333333) + ((ibs0>>2)&0x3333333333333333);
               ibs0 = (ibs0+(ibs0>>4)) & 0x0F0F0F0F0F0F0F0F;
               word2 += (ibs0+(ibs0>>8)) & 0x00FF00FF00FF00FF;
            }
            word1 = (word1+(word1>>16)) & 0x0000FFFF0000FFFF;
            HomHomCount = (word1+(word1>>32)) & 0xFFFFFFFF;
            word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;
            IBS0Count = (word2+(word2>>32)) & 0xFFFFFFFF;
            if(IBS0Count >= HomHomCount * 0.5) continue; // C_Hom < 50%
            // 2nd quarter
            word1 = word2 = 0;
            m2 = (stop1>>3);
            for(int m = (stop1>>4); m < m2; m++){
               word = SLG[0][id1][m] & SLG[0][id2][m]; // HomHom
               ibs0 = word & (SLG[1][id1][m] ^ SLG[1][id2][m]);   // IBS0
               word = word - ((word>>1)&0x5555555555555555);
               word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
               word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
               word1 += (word+(word>>8)) & 0x00FF00FF00FF00FF;
               ibs0 = ibs0 - ((ibs0>>1)&0x5555555555555555);   // IBS0
               ibs0 = (ibs0&0x3333333333333333) + ((ibs0>>2)&0x3333333333333333);
               ibs0 = (ibs0+(ibs0>>4)) & 0x0F0F0F0F0F0F0F0F;
               word2 += (ibs0+(ibs0>>8)) & 0x00FF00FF00FF00FF;
            }
            word1 = (word1+(word1>>16)) & 0x0000FFFF0000FFFF;
            HomHomCount += (word1+(word1>>32)) & 0xFFFFFFFF;
            word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;
            IBS0Count += (word2+(word2>>32)) & 0xFFFFFFFF;
            if(IBS0Count >= HomHomCount * 0.475) continue; // C_Hom < 52.5%
            // 3rd quarter
            word1 = word2 = 0;
            m2 = (stop1>>2);
            for(int m = (stop1>>3); m < m2; m++){
               word = SLG[0][id1][m] & SLG[0][id2][m]; // HomHom
               ibs0 = word & (SLG[1][id1][m] ^ SLG[1][id2][m]);   // IBS0
               word = word - ((word>>1)&0x5555555555555555);
               word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
               word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
               word1 += (word+(word>>8)) & 0x00FF00FF00FF00FF;
               ibs0 = ibs0 - ((ibs0>>1)&0x5555555555555555);   // IBS0
               ibs0 = (ibs0&0x3333333333333333) + ((ibs0>>2)&0x3333333333333333);
               ibs0 = (ibs0+(ibs0>>4)) & 0x0F0F0F0F0F0F0F0F;
               word2 += (ibs0+(ibs0>>8)) & 0x00FF00FF00FF00FF;
            }
            word1 = (word1+(word1>>16)) & 0x0000FFFF0000FFFF;
            HomHomCount += (word1+(word1>>32)) & 0xFFFFFFFF;
            word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;
            IBS0Count += (word2+(word2>>32)) & 0xFFFFFFFF;
            if(IBS0Count >= HomHomCount * 0.45) continue; // C_Hom < 55%

            het1Count = hetInOnePersonCount[id1][0];
            het2Count = hetInOnePersonCount[id2][0];
            m3 = (het1Count < het2Count? het1Count: het2Count);
            HetHomCount = (stop1 << 5) - het1Count - het2Count -
               ((HomHomCount + missingInOnePersonCount[id1][0] + missingInOnePersonCount[id2][0])<<1);
            threshold = (0.5-kincutoff2)*m3 - HetHomCount*0.25;   // kinship < 0.0625
            if(IBS0Count >= threshold) continue;   // Kinship < 0.0625

            // 4th quarter
            word1 = word2 = 0;
            for(int m = (stop1>>2); m < stop1; m++){
               word = SLG[0][id1][m] & SLG[0][id2][m]; // HomHom
               ibs0 = word & (SLG[1][id1][m] ^ SLG[1][id2][m]);   // IBS0
               word = word - ((word>>1)&0x5555555555555555);
               word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
               word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
               word1 += (word+(word>>8)) & 0x00FF00FF00FF00FF;
               ibs0 = ibs0 - ((ibs0>>1)&0x5555555555555555);   // IBS0
               ibs0 = (ibs0&0x3333333333333333) + ((ibs0>>2)&0x3333333333333333);
               ibs0 = (ibs0+(ibs0>>4)) & 0x0F0F0F0F0F0F0F0F;
               word2 += (ibs0+(ibs0>>8)) & 0x00FF00FF00FF00FF;
            }
            word1 = (word1+(word1>>16)) & 0x0000FFFF0000FFFF;
            HomHomCount += (word1+(word1>>32)) & 0xFFFFFFFF;
            word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;
            IBS0Count += (word2+(word2>>32)) & 0xFFFFFFFF;
            if(IBS0Count >= HomHomCount * 0.375) continue; // C_Hom < 62.5%

            if(i==j && id2 <= id1) continue;
            HetHomCount = (stop1 << 7) - het1Count - het2Count -
               ((HomHomCount + missingInOnePersonCount[id1][0] + missingInOnePersonCount[id2][0])<<1);
            threshold = (0.5-kincutoff2)*m3 - HetHomCount*0.25;   // kinship < 0.0442
            if(IBS0Count >= threshold || IBS0Count >= HomHomCount*0.375)
               continue;   // Kinship < 0.0442 or CHom < 62.5%
            m1 = m2 = m3 = 0;
            for(int m = 0; m < stop1; m++){
               word1 = ~(SLG[0][id1][m] | SLG[1][id1][m]);
               word2 = ~(SLG[0][id2][m] | SLG[1][id2][m]);
               for(word = word1 & word2; word; word &= (word-1), m3++);
               for(word = (~SLG[0][id1][m]) & SLG[1][id1][m] & word2; word; word &= (word-1), m1++);
               for(word = word1 & (~SLG[0][id2][m]) & SLG[1][id2][m]; word; word &= (word-1), m2++);
            }
            het1Count -= m1;
            het2Count -= m2;
            HetHomCount += m1 + m2 + (m3<<1);
            m3 = (het1Count < het2Count? het1Count: het2Count);
            threshold = (0.5-kincutoff2)*m3 - HetHomCount*0.25;  // kinship < 0.0442
            if(IBS0Count >= threshold || IBS0Count >= HomHomCount * 0.375) continue;
            rpList[thread].Push(id1);
            rpList[thread].Push(id2);
         }  // end of id2
      }  // end of id1
   }  // end of loop blocking
#ifdef _OPENMP
}
#endif
}
   delete []loopIndex;
   for(int i = 0; i < idCount; i++){
      delete []missingInOnePersonCount[i];
      delete []hetInOnePersonCount[i];
   }
   delete []missingInOnePersonCount;
   delete []hetInOnePersonCount;
}

void Engine::IntegratedRelationshipInference()
{
   if(shortCount==0) error("No genotype data");
   printf("Autosome genotypes stored in %d", Bit64==64? longCount:shortCount);
   printf(" words for each of %d individuals.\n", idCount);

   printf("\nOptions in effect:\n");
   printf("\t--related\n");
   if(relativedegree)
      printf("\t--degree %d\n", relativedegree);
   if(Bit64Flag)
      printf("\t--sysbit 64\n");
   if(CoreCount)
      printf("\t--cpus %d\n", CoreCount);
   if(lessmemFlag)
      printf("\t--lessmem\n");
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");

   // default degree of relatedness is 1
   int originaldegree = relativedegree;
   if(relativedegree == 0) relativedegree = 1;

   int stop1, stop2;
   if(Bit64==64){
      stop1 = 64;
      stop2 = (stop1<<3);
      if(relativedegree == 2) stop1 <<= 3;
      if(stop1 > longCount) stop1 = longCount;
      if(stop2 > longCount) stop2 = longCount;
   }else{
      stop1 = 256;
      stop2 = (stop1<<3);
      if(relativedegree == 2) stop1 <<= 3;
      if(stop1 > shortCount) stop1 = shortCount;
      if(stop2 > shortCount) stop2 = shortCount;
   }

   printf("Sorting autosomes...\n");
   if(Bit64==64){
      if(SLG[0]==NULL){
         if(relativedegree==1)
            ConvertLGtoSLG(LG, markerCount, SLG, (stop2 < longCount)? (stop2<<6): markerCount);
         else if(relativedegree==2)
            ConvertLGtoSLG(LG, markerCount, SLG, (stop1 < longCount)? (stop1<<6): markerCount);
      }
   }else{
      if(SG[0]==NULL){
         if(relativedegree==1)
            ConvertGGtoSG(GG, markerCount, SG, (stop2 < shortCount)? (stop2<<4): markerCount);
         else if(relativedegree==2)
            ConvertGGtoSG(GG, markerCount, SG, (stop1 < shortCount)? (stop1<<4): markerCount);
      }
   }

   bool IBDvalidFlag = false;
   if(Bit64==64)
      IBDvalidFlag = PreSegment(/*chrSeg, totalLength, segmessage*/);
   if(!IBDvalidFlag){
      printf("%s\n", (const char*)segmessage);
      printf("  Inference will be based on kinship estimation only.\n");
   }
   int notMissingCount;
   int afterCount[6];
   int id1, id2;
   Kinship kin;
   double kinship, smaller, CHet, ibdprop, ibd2prop;
   IntArray pairs(0);
   Vector phis(0), pi0s(0);
   for(int f = 0; f < ped.familyCount; f++){
      if(id[f].Length()<2) continue;   // no pairs in family f
      kin.Setup(*ped.families[f]);
      for(int i = 0; i < id[f].Length(); i++)
         for(int j = i+1; j < id[f].Length(); j++){
            pairs.Push(geno[id[f][i]]);
            pairs.Push(geno[id[f][j]]);
            double phi = kin(ped[id[f][i]], ped[id[f][j]]);
            phis.Push(phi);
            double pi0 = 0.0;
            if(phi < 0.2)
               pi0 = 1-4*phi;
            else if(phi < 0.3 && ped[id[f][i]].isSib(ped[id[f][j]]))
               pi0 = 0.25;
            pi0s.Push(pi0);
         }
   }
   int pairCount = pairs.Length()/2;
   if(pairCount){ // within-family inference
      IntArray HetHetCounts, IBS0Counts, het1Counts, het2Counts, HomHomCounts, IBSCounts;
      //****************************Compute kinship coefficient*****************************
      if(Bit64==64)
         KinshipInSubset64Bit(pairs, HetHetCounts, IBS0Counts, het1Counts, het2Counts, HomHomCounts, IBSCounts);
      else
         KinshipInSubset(pairs, HetHetCounts, IBS0Counts, het1Counts, het2Counts, HomHomCounts, IBSCounts);
      //****************************Compute kinship coefficient*****************************

      Vector ibdprops, maxLengths, ibd2props, maxLengths2;
      if(IBDvalidFlag && Bit64 == 64)
         IBDSegInSubset64Bit(pairs, ibdprops, maxLengths, ibd2props, maxLengths2);
      String outfile;
      outfile.Copy(prefix);
      outfile.Add(".kin");
      FILE *fp = fopen(outfile, "wt");
      fprintf(fp, "FID\tID1\tID2\tN_SNP\tZ0\tPhi\tHetHet\tIBS0\tHetConc\tHomIBS0\tKinship");
      if(IBDvalidFlag)
         fprintf(fp, "\tIBD1Seg\tIBD2Seg\tPropIBD\tInfType");
      fprintf(fp, "\tError\n");
      int degree;
      double errorFlag, inflation, ibs0;
      int beforeCount[6];
      Vector IBS0PO(0), IBS0FS(0), IBS0L1(0);
      String type;
      for(int i = 0; i < 6; i++) beforeCount[i] = afterCount[i] = 0;
      for(int p = 0; p < pairCount; p++){
         if(!het1Counts[p] && !het2Counts[p] && !HetHetCounts[p]) continue;
         id1 = pairs[p*2]; id2 = pairs[p*2+1];
         kinship = (HetHetCounts[p] - IBS0Counts[p]*2.0) / (HetHetCounts[p]*2+het1Counts[p]+het2Counts[p]);
         notMissingCount = HetHetCounts[p]+het1Counts[p]+het2Counts[p]+HomHomCounts[p];
         ibs0 = IBS0Counts[p]*1.0/notMissingCount;
         CHet = HetHetCounts[p] * 1.0 / (HetHetCounts[p]+het1Counts[p]+het2Counts[p]);
         degree = 4;
         if(phis[p] > 0.0442)
            degree = int(-log(phis[p])/log(2.0) - 0.5);
         if(degree < 4){
            beforeCount[degree] ++;
            if(degree == 1)
            if(pi0s[p]==0) beforeCount[5] ++;
         }else
            beforeCount[4] ++;
         if(IBDvalidFlag){ // IBD seg
            ibdprop = ibdprops[p];
            ibd2prop = ibd2props[p];
            double pi = ibd2prop + ibdprop * 0.5;
            if(CHet<0.8){  // not MZ/Dup
               if(pi > 0.3535534){  // 1st-degree
                  if(ibdprop + ibd2prop > 0.96 || (ibdprop + ibd2prop > 0.9 && ibd2prop <= 0.08))
                     type = "PO";
                  else if(ibd2prop > 0.08)
                     type = "FS";
                  else
                     type = "2nd";
               }else if(pi > 0.1767767){  // 2nd-degree
                  if(pi > 0.32 && ibd2prop > 0.15)
                     type = "FS";
                  else
                     type = "2nd";
               }else if(pi > 0.08838835)
                  type = "3rd";
               else if(pi > 0.04419417)
                  type = "4th";
               else
                  type = "UN";
            }else // Duplicate
               type="Dup/MZ";
            if(type=="Dup/MZ"){
               afterCount[0]++;
               errorFlag = phis[p]==0.5? 0: 1;
            }else if(type=="PO"){
               afterCount[1]++;
               IBS0PO.Push(ibs0);
               IBS0L1.Push(ibs0);
               afterCount[5]++;
               errorFlag = phis[p]==0.25 && pi0s[p] == 0? 0: 1;
            }else if(type=="FS"){
               afterCount[1]++;
               IBS0FS.Push(ibs0);
               IBS0L1.Push(ibs0);
               errorFlag = phis[p]==0.25 && pi0s[p] == 0.25? 0: 1;
            }else if(type=="2nd"){
               afterCount[2]++;
               inflation = (phis[p] > 0)? pi*0.5 / phis[p]: -1;
               if(inflation > 2 || inflation < 0.5)
                  errorFlag = 1;
               else if(inflation > 1.4142 || inflation < 0.70711)
                  errorFlag = 0.5;
               else
                  errorFlag = 0;
            }else if(type=="3rd"){
               afterCount[3]++;
               inflation = (phis[p] > 0)? pi*0.5 / phis[p]: -1;
               if(inflation > 2 || inflation < 0.5)
                  errorFlag = 1;
               else if(inflation > 1.4142 || inflation < 0.70711)
                  errorFlag = 0.5;
               else
                  errorFlag = 0;
            }else if(type=="4th"){
               afterCount[4]++;
               inflation = (phis[p] > 0)? pi*0.5 / phis[p]: 0.51; // flag=0.5 for UN
               if(inflation > 2 || inflation < 0.5)
                  errorFlag = 1;
               else if(inflation > 1.4142 || inflation < 0.70711)
                  errorFlag = 0.5;
               else
                  errorFlag = 0;
            }else{   // type="UN"
               afterCount[4]++;
               errorFlag = phis[p]>0? 1: 0;
            }
         }else{   // no IBD seg
            errorFlag = 0;
            inflation = (phis[p] > 0)? kinship / phis[p]: -1;
            if(phis[p] > 0.03){  // up to 4th-degree relative
               if(inflation > 2 || inflation < 0.5)
                  errorFlag = 1;
               else if(inflation > 1.4142 || inflation < 0.70711)
                  errorFlag = 0.5;
            }else if(phis[p] < 0.005){  // unrelated pair
               if(kinship > 0.0442) errorFlag = 1;
               else if(kinship > 0.0221) errorFlag = 0.5;
            }else{   // distant relatives
               if(kinship < -0.0221) errorFlag = 1;
               else errorFlag = 0.5;
            }
            degree = 4;
            if(kinship > 0.0442)
               degree = int(-log(kinship)/log(2.0) - 0.5);
            if(degree < 4){
               afterCount[degree] ++;
               if(degree == 1 && errorrateCutoff != _NAN_ && ibs0 < errorrateCutoff)
                  afterCount[5]++;
            }else
               afterCount[4] ++;
            if(degree==1 && phis[p]==0.25){
               if(pi0s[p]==0)
                  IBS0PO.Push(ibs0);
               else
                  IBS0FS.Push(ibs0);
            }
            if(degree==1)
               IBS0L1.Push(ibs0);
         }
         fprintf(fp, "%s\t%s\t%s\t%d\t%.3lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf",
            (const char*)ped[phenoid[id1]].famid,
            (const char*)ped[phenoid[id1]].pid,
            (const char*)ped[phenoid[id2]].pid,
            notMissingCount, pi0s[p], phis[p],
            HetHetCounts[p]*1.0/notMissingCount, ibs0,
            CHet, // CHet
            ibs0/(ibs0+(IBSCounts[p]-HetHetCounts[p])*1.0/notMissingCount), kinship);
         if(IBDvalidFlag)
            fprintf(fp, "\t%.4lf\t%.4lf\t%.4lf\t%s\t%G\n",
               ibdprops[p], ibd2props[p], ibd2props[p] + ibdprops[p]*0.5, (const char*)type, errorFlag);
         else
            fprintf(fp, "\t%G\n", errorFlag);
      }  // end of pairs
      bool pedigreeFlag = false;
      for(int i = 0; i < 6; i++)
         if(beforeCount[i]) pedigreeFlag = true;
      if(pedigreeFlag){
         if(errorrateCutoff==_NAN_ && !IBDvalidFlag){
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
         if(ped.familyCount < 2) {
            if(ped.familyCount==1 && ped.families[0]->famid=="0")
               warning("All individuals with family ID 0 are considered as relatives.\n");
            printf("There is only one family.\n");
            return;
         }
      }
      fclose(fp);
   }else
      printf("Each family consists of one individual.\n");

   if(relativedegree < 3)  // fast algorithm for inference of close relatives
      printf("A subset of informative SNPs will be used to screen close relatives.\n");
   printf("Relationship inference across families starts at %s", currentTime());
#ifdef _OPENMP
   printf("%d CPU cores are used...\n", defaultMaxCoreCount);
#endif
   IntArray *allpairs0 = new IntArray [defaultMaxCoreCount];
   IntArray *allpairs = new IntArray [defaultMaxCoreCount];
   // most computationally intensive here: SCREENING RELATIVES
   if(Bit64==64){
      if(relativedegree < 3)
         ScreenCloseRelativesInSubset64Bit(allpairs0);
      else // more distant relatives
         ComputeLongRobustKinship64BitWithFilter(allpairs0, false);
   }else
      ScreenCloseRelativesInSubset(allpairs0);
   //****************SCREENING RELATIVES end************************

   long long int midrelativeCount = 0;
   for(int t = 0; t < defaultMaxCoreCount; t++)
      allpairs[t].Dimension(0);
   for(int t = 0; t < defaultMaxCoreCount; t++){
      int pairCount = allpairs0[t].Length()/2;
      for(int i = 0; i < pairCount; i++){
         id1 = allpairs0[t][i*2];
         id2 = allpairs0[t][i*2+1];
         if(ped[phenoid[id1]].famid != ped[phenoid[id2]].famid){
            allpairs[t].Push(id1);
            allpairs[t].Push(id2);
         }
      }
      midrelativeCount += allpairs[t].Length()/2;
   }
   delete []allpairs0;
   if(midrelativeCount==0){
      printf("                                           ends at %s", currentTime());
      printf("No close relatives are inferred.\n\n");
      return;
   }
   double lowerbound = pow(2.0, -(relativedegree+1.5));
   double kincutoff2 = 0.125;
   double kincutoff1 = kincutoff2 * 0.7071068;  // 0.0884
   if(relativedegree == 2){
      kincutoff2 *= 0.5;
      kincutoff1 *= 0.5;
   }
   if(midrelativeCount > 0){
      if(relativedegree < 3){
         printf("  Stages 1&2 (with %d SNPs): %lli pairs of relatives are detected (with kinship > %.4lf)\n",
            markerCount < 32768? markerCount: 32768, midrelativeCount, kincutoff2);
         printf("                               Screening ends at %s", currentTime());
         fflush(stdout);
      }
   }else { // not screened
      for(int t = 0; t < defaultMaxCoreCount; t++)
         allpairs[t].Dimension(0);
      midrelativeCount = (idCount & 1)? (idCount-1)/2*idCount: idCount/2*(idCount-1);
      long long int blockSize = (midrelativeCount-1)/defaultMaxCoreCount+1;
      long long int index = 0;
      for(int i = 0; i < idCount; i++)
         for(int j = i+1; j < idCount; j++){
            int thread = index / blockSize;
            allpairs[thread].Push(i);
            allpairs[thread].Push(j);
            index ++;
         }
   }
   IntArray HetHetCounts, IBS0Counts, het1Counts, het2Counts, HomHomCounts, IBSCounts;
   Vector ibdprops, maxLengths, ibd2props, maxLengths2;

   String outfile(prefix);
   outfile.Add(".kin0");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "FID1\tID1\tFID2\tID2\tN_SNP\tHetHet\tIBS0\tHetConc\tHomIBS0\tKinship");
   if(IBDvalidFlag) fprintf(fp, "\tIBD1Seg\tIBD2Seg\tPropIBD\tInfType");
   fprintf(fp,"\n");
   for(int j = 0; j < 6; j++)
      afterCount[j] = 0;
   int fpCount = 0;

   for(int t = 0; t < defaultMaxCoreCount; t++){
      int pairCount = allpairs[t].Length()/2;
      if(pairCount == 0) continue;

      //****************************Compute kinship coefficient*****************************
      if(Bit64==64)
         KinshipInSubset64Bit(allpairs[t], HetHetCounts, IBS0Counts, het1Counts, het2Counts, HomHomCounts, IBSCounts);
      else
         KinshipInSubset(allpairs[t], HetHetCounts, IBS0Counts, het1Counts, het2Counts, HomHomCounts, IBSCounts);
      //****************************Compute kinship coefficient*****************************

      if(IBDvalidFlag && Bit64 == 64)
         IBDSegInSubset64Bit(allpairs[t], ibdprops, maxLengths, ibd2props, maxLengths2);

      for(int p = 0; p < pairCount; p++){
         if(!het1Counts[p] && !het2Counts[p] && !HetHetCounts[p]) continue;
         id1 = allpairs[t][p*2];
         id2 = allpairs[t][p*2+1];
         smaller = HetHetCounts[p] + (het1Counts[p] < het2Counts[p]? het1Counts[p]: het2Counts[p]);
         kinship = 0.5 - ((het1Counts[p]+het2Counts[p])*0.25+IBS0Counts[p])/smaller;
         if(IBDvalidFlag){
            ibdprop = ibdprops[p];
            ibd2prop = ibd2props[p];
         }
         if(relativedegree == 1){
            if(IBDvalidFlag){
               double pi = ibd2prop + ibdprop * 0.5;
               if(((kinship < lowerbound) && // kinship < 0.177 AND following criteria is not 1st-degree
                  ((pi < 0.3535534) || (ibd2prop <= 0.08 && ibdprop+ibd2prop <= 0.9))) ) {fpCount++; continue;}
            }else{
               if(kinship < lowerbound) {fpCount++; continue;}
            }
         }else{   // degree==2
            if(IBDvalidFlag){
               if((kinship < 0) || ((kinship < lowerbound) && (ibd2prop*0.5 + ibdprop*0.25 < lowerbound))) {fpCount++; continue;}
            }else{
               if(kinship < lowerbound) {fpCount++; continue;}
            }
         }
         CHet = HetHetCounts[p] * 1.0/ (HetHetCounts[p] + het1Counts[p] + het2Counts[p]);
         notMissingCount = HetHetCounts[p] + het1Counts[p] + het2Counts[p] + HomHomCounts[p];
         fprintf(fp, "%s\t%s\t%s\t%s\t%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf",
            (const char*)ped[phenoid[id1]].famid, (const char*)ped[phenoid[id1]].pid,
            (const char*)ped[phenoid[id2]].famid, (const char*)ped[phenoid[id2]].pid,
            notMissingCount, HetHetCounts[p]*1.0/notMissingCount, IBS0Counts[p]*1.0/notMissingCount,
            CHet, // CHet
            IBS0Counts[p]*1.0/(IBS0Counts[p]+IBSCounts[p]-HetHetCounts[p]), kinship);
         if(IBDvalidFlag){
            double pi = ibd2prop + ibdprop * 0.5;
            fprintf(fp, "\t%.4lf\t%.4lf\t%.4lf\t", ibdprop, ibd2prop, pi);
            if(CHet<0.8){
               if(pi > 0.3535534){  // 1st-degree
                  if(ibdprop + ibd2prop > 0.96 || (ibdprop + ibd2prop > 0.9 && ibd2prop <= 0.08)){
                     fprintf(fp, "PO");
                     afterCount[1]++;
                     afterCount[5]++;
                  }else if(ibd2prop > 0.08){
                     fprintf(fp, "FS");
                     afterCount[1]++;
                  }else{
                     fprintf(fp, "2nd");
                     afterCount[2]++;
                  }
               }else if(pi > 0.1767767){  // 2nd-degree
                  if(pi > 0.32 && ibd2prop > 0.15){
                     fprintf(fp, "FS");
                     afterCount[1]++;
                  }else{
                     fprintf(fp, "2nd");
                     afterCount[2]++;
                  }
               }else if(pi > 0.08838835){
                  fprintf(fp, "3rd");
                  afterCount[3]++;
               }else if(pi > 0.04419417){
                  fprintf(fp, "4th");
               }else
                  fprintf(fp, "UN");
            }else{ // Duplicate
               fprintf(fp, "Dup/MZ");
               afterCount[0]++;
            }
         }
         fprintf(fp, "\n");
      }
   }  // end of t loop
   fclose(fp);
   delete []allpairs;
   long long int relativeCount = midrelativeCount - fpCount;
   if(IBDvalidFlag){
      relativeCount = afterCount[0] + afterCount[1];
      if(relativedegree >= 2) relativeCount += afterCount[2];
      if(relativedegree >= 3) relativeCount += afterCount[3];
   }
   if(relativedegree < 3)  // fast algorithm for inference of close relatives
      printf("  Final Stage (with %d SNPs): %lli pairs of relatives (up to %d%s-degree) are confirmed\n",
      markerCount, relativeCount, relativedegree, relativedegree==1?"st":"nd");
   printf("                               Inference ends at %s", currentTime());
   if(relativeCount==0){
      printf("No cryptic relatedness (up to the %d-degree) is found.\n", relativedegree);
      return;
   }else if(IBDvalidFlag)
      printRelationship(NULL, afterCount);
   printf("\nBetween-family relatives (kinship >= %.5lf) saved in file %s\n",
      lowerbound, (const char*)outfile);
   if(relativedegree == 1){
      printf("\nNote only duplicates and 1st-degree relatives are included in the inference.\n");
      printf("  Specifying '--degree 2' if a higher degree relationship inference is needed.\n\n");
   }
   relativedegree = originaldegree;
}

void Engine::ScreenCloseRelativesInSubset(IntArray rpList[])
{
   if(shortCount < 256 || !bigdataFlag) return;  // screening is not needed
   if(idCount > 0x1FFFFF){
      printf("This version of KING supports up to %d samples.\n", 0x1FFFFF);
      printf("Please contact the KING authors to allow an even larger sample size.\n");
      return;
   }
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
   int stop1, start2, stop2;
   stop1 = short_prescan;
   if(stop1 > shortCount) stop1 = shortCount;
   start2 = stop1;
   stop2 = (short_prescan<<3);
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
   int m1, m2, m3;
   unsigned short int word, word2;
   int **missingInOnePersonCount = new int * [idCount];
   int **hetInOnePersonCount = new int * [idCount];
   unsigned short int **missingWordIndex = new unsigned short int *[idCount];
   for(int i = 0; i < idCount; i++){
      missingInOnePersonCount[i] = new int [2];
      hetInOnePersonCount[i] = new int [2];
      for(int j = 0; j < 2; j++)
         missingInOnePersonCount[i][j] = hetInOnePersonCount[i][j] = 0;
      missingWordIndex[i] = new unsigned short int [(stop2-1)/16+1];
      for(int m = 0; m < (stop2-1)/16+1; m++)
         missingWordIndex[i][m] = 0;
   }
   unsigned short int *masks = new unsigned short int[stop2];
   for(int i = 0; i < stop2;  i++)
      masks[i] = 0xFFFF;
   if(stop2 == shortCount)
      masks[stop2-1] = markerCount%16? ((1 << (markerCount % 16))-1): 0xFFFF;
   for(int i = 0; i < idCount; i++){
      for(int m = 0; m < stop1; m++){
         word = (~SG[0][i][m]) & (~SG[1][i][m]) & masks[m];
         if(word){   // not all missing
            missingInOnePersonCount[i][0] += oneoneCount[word];
            missingWordIndex[i][m/16] |= 1<<(m%16);
         }
      }
      if(relativedegree==1)
         for(int m = start2; m < stop2; m++){
            word = (~SG[0][i][m]) & (~SG[1][i][m]) & masks[m];
            if(word){   // some missing
               missingInOnePersonCount[i][1] += oneoneCount[word];
               missingWordIndex[i][m/16] |= 1<<(m%16);
            }
         }
      for(int m = 0; m < stop1; m++)
         hetInOnePersonCount[i][0] += oneoneCount[(~SG[0][i][m]) & SG[1][i][m]];
      if(relativedegree==1)
         for(int m = start2; m < stop2; m++)
            hetInOnePersonCount[i][1] += oneoneCount[(~SG[0][i][m]) & SG[1][i][m]];
   }
   delete []masks;
   const int cutoffMissingCount = stop1*16-MINSNPCOUNT;
   for(int i = 0; i < defaultMaxCoreCount; i++)
      rpList[i].Dimension(0);
   double threshold, smaller, kinship;
   int id1, id2, IBS0Count, het1Count, het2Count, HetHomCount, HomHomCount, HetHetCount;
   const char LOOPBLOCKINGSIZE=32;
   const char SHIFTSIZE=5;
   unsigned int loopIndexLength = (unsigned int)((idCount-1)/LOOPBLOCKINGSIZE+1)*(((idCount-1)/LOOPBLOCKINGSIZE+1)+1)/2;
   unsigned int *loopIndex = new unsigned int [loopIndexLength];
   unsigned int index = 0;
   for(unsigned int i = 0; i < idCount; i += LOOPBLOCKINGSIZE){
      unsigned int iShifted = (i<<(16-SHIFTSIZE));
      for(int j = i; j < idCount; j += LOOPBLOCKINGSIZE)
         loopIndex[index++] = iShifted | (j>>SHIFTSIZE);
   }
   int thread = 0;
   char revbase[65536];
   for(int i = 0; i < 16; i++)
      revbase[shortbase[i]] = i;
   char rightmost[65536];
   for(unsigned int i = 0; i < 65536; i++)
      rightmost[i] = revbase[i&(-i)];
   if(relativedegree == 1)
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(HomHomCount, IBS0Count, HetHetCount, HetHomCount, het1Count, het2Count,\
      threshold, id1, id2, kinship, smaller, m1, m2, m3, thread, word, word2)
{
   thread = omp_get_thread_num();
   #pragma omp for
#endif
   for(unsigned long int k = 0; k < loopIndexLength; k++){
      int i = (loopIndex[k]>>16)<<SHIFTSIZE;
      int iMax = (i > idCount - LOOPBLOCKINGSIZE? idCount: i + LOOPBLOCKINGSIZE);
      int j = (loopIndex[k]&0xFFFF)<<SHIFTSIZE;
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
               HomHomCount += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1]];
               IBS0Count += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1] & (SG[1][id1][m1] ^ SG[1][id2][m1])];
            }
            if(IBS0Count >= HomHomCount * 0.5) continue; // C_Hom < 50%
            // 2nd quarter
            for(m1 = stop1/4, m2 = stop1/2; m1 < m2; m1++){
               HomHomCount += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1]];
               IBS0Count += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1] & (SG[1][id1][m1] ^ SG[1][id2][m1])];
            }
            if(IBS0Count >= HomHomCount * 0.4) continue; // C_Hom < 60%
            // 3rd quarter
            for(m1 = stop1/2, m2 = stop1*3/4; m1 < m2; m1++){
               HomHomCount += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1]];
               IBS0Count += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1] & (SG[1][id1][m1] ^ SG[1][id2][m1])];
            }
            if(IBS0Count >= HomHomCount * 0.375) continue; // C_Hom < 62.5%
            // 4th quarter
            for(m1 = stop1*3/4; m1 < stop1; m1++){
               HomHomCount += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1]];
               IBS0Count += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1] & (SG[1][id1][m1] ^ SG[1][id2][m1])];
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

            for(HetHomCount = het1Count = het2Count = m1 = 0; m1 < (stop1-1)/16+1; m1++){
               for(word = missingWordIndex[id1][m1] & missingWordIndex[id2][m1];
                  word; word &= (word-1)){   // MissMiss
                  m2 = m1*16+rightmost[word];   // m2 is the word with nonmissingness
                  for(word2 = ~(SG[0][id1][m2] | SG[0][id2][m2] | SG[1][id1][m2] | SG[1][id2][m2])&0xFFFF;
                       word2; word2 &= (word2-1))
                     HetHomCount ++;
               }
               for(word = missingWordIndex[id2][m1];  // HetMiss
                  word; word &= (word-1)){
                  m2 = m1*16+rightmost[word];
                  for(word2 = ~(SG[0][id1][m2] | SG[0][id2][m2] | SG[1][id2][m2]) & SG[1][id1][m2];
                      word2; word2 &= (word2-1))
                     het1Count ++;
               }
               for(word = missingWordIndex[id1][m1];  // MissHet
                  word; word &= (word-1)){
                  m2 = m1*16+rightmost[word];
                  for(word2 = ~(SG[0][id1][m2] | SG[0][id2][m2] | SG[1][id1][m2]) & SG[1][id2][m2];
                      word2; word2 &= (word2-1))
                     het2Count ++;
               }
            }
            het1Count = hetInOnePersonCount[id1][0] - het1Count;
            het2Count = hetInOnePersonCount[id2][0] - het2Count;
            HetHomCount = stop1 * 32 - HomHomCount * 2 - het1Count - het2Count +
               (-missingInOnePersonCount[id1][0] - missingInOnePersonCount[id2][0] + HetHomCount) * 2;
            m3 = (het1Count < het2Count? het1Count: het2Count);
            threshold = (0.5-kincutoff1*0.707107)*m3 - HetHomCount*0.25;  // kinship < 0.0625
            cutoffDist = (0.5-kincutoff1)*m3 - HetHomCount*0.25;  // C_Hom < 70% && kinship < 0.0884
            if(HomHomCount * 0.7 > cutoffDist) cutoffDist = HomHomCount * 0.7;
            if(cutoffDist < threshold) threshold = cutoffDist;
            if(IBS0Count >= threshold) continue;
//            rawrelativeCount[thread] ++;

            // Stage 2: few pairs
            for(HetHomCount -= (missingInOnePersonCount[id1][1]+missingInOnePersonCount[id2][1]),  // different non-missing genotypes
               m1 = start2; m1 < stop2; m1++) // different-genotpe Count
                  HetHomCount += oneoneCount[SG[0][id1][m1]^SG[0][id2][m1]];
            het1Count += hetInOnePersonCount[id1][1];
            het2Count += hetInOnePersonCount[id2][1];
            cutoffDist = (2-kincutoff2*4) * (het1Count < het2Count? het1Count: het2Count)-IBS0Count*4;
            if(HetHomCount >= cutoffDist) continue; // too big distance already
            cutoffDist = (cutoffDist+IBS0Count*4-HetHomCount)*0.25;
            for(m1 = start2; m1 < (start2+stop2)/2; m1++)
               IBS0Count += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1] & (SG[1][id1][m1] ^ SG[1][id2][m1])];
            if(IBS0Count >= cutoffDist * 1.1) continue;
            for(m1 = (start2+stop2)/2; (m1 < stop2) && (IBS0Count < cutoffDist); m1++)
               IBS0Count += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1] & (SG[1][id1][m1] ^ SG[1][id2][m1])];
            if(IBS0Count >= cutoffDist) continue;   // too big distance
            // now add the MissMissCount back

            for(m1 = start2/16, m3=0; m1 < (stop2-1)/16+1; m1++)
               for(word = missingWordIndex[id1][m1] & missingWordIndex[id2][m1];
                  word; word &= (word-1)){   // MissMiss
                  m2 = m1*16+rightmost[word];   // m2 is the word with nonmissingness
                  for(word2 = ~(SG[0][id1][m2] | SG[0][id2][m2] | SG[1][id1][m2] | SG[1][id2][m2])&0xFFFF;
                    word2; word2 &= (word2-1))
                     m3++;
               }
            HetHomCount += (m3<<1);

            for(m1 = start2/16, m3=0; m1 < (stop2-1)/16+1; m1++)
               for(word = missingWordIndex[id2][m1];  // HetMiss
                  word; word &= (word-1)){
                  m2 = m1*16+rightmost[word];
                  for(word2 = ~(SG[0][id1][m2] | SG[0][id2][m2] | SG[1][id2][m2]) & SG[1][id1][m2];
                    word2; word2 &= (word2-1))
                     m3 ++;
               }
            het1Count -= m3;
            HetHomCount += m3;

            for(m1 = start2/16, m3=0; m1 < (stop2-1)/16+1; m1++)
               for(word = missingWordIndex[id1][m1];  // MissHet
                  word; word &= (word-1)){
                  m2 = m1*16+rightmost[word];
                  for(word2 = ~(SG[0][id1][m2] | SG[0][id2][m2] | SG[1][id1][m2]) & SG[1][id2][m2];
                      word2; word2 &= (word2-1))
                     m3 ++;
               }
            het2Count -= m3;
            HetHomCount += m3;

            HetHetCount = (het1Count + het2Count - HetHomCount)/2;
            het1Count -= HetHetCount;
            het2Count -= HetHetCount;
            smaller = HetHetCount + (het1Count < het2Count? het1Count: het2Count);
            kinship = 0.5 - ((het1Count+het2Count)*0.25+IBS0Count)/smaller;
            if(kinship < kincutoff2) continue;
            rpList[thread].Push(id1);
            rpList[thread].Push(id2);
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
   for(unsigned long int k = 0; k < loopIndexLength; k++){
      int i = (loopIndex[k]>>16)<<SHIFTSIZE;
      int iMax = (i > idCount - LOOPBLOCKINGSIZE? idCount: i + LOOPBLOCKINGSIZE);
      int j = (loopIndex[k]&0xFFFF)<<SHIFTSIZE;
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
                  HomHomCount += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1]];
                  IBS0Count += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1] & (SG[1][id1][m1] ^ SG[1][id2][m1])];
               }
               if(IBS0Count >= HomHomCount * 0.5) continue; // C_Hom < 50%
               // 2nd quarter
               for(m1 = stop1/16, m2 = stop1/8; m1 < m2; m1++){
                  HomHomCount += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1]];
                  IBS0Count += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1] & (SG[1][id1][m1] ^ SG[1][id2][m1])];
               }
               if(IBS0Count >= HomHomCount * 0.45) continue; // C_Hom < 55%
               // 3rd quarter
               for(m1 = stop1/8, m2 = stop1/4; m1 < m2; m1++){
                  HomHomCount += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1]];
                  IBS0Count += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1] & (SG[1][id1][m1] ^ SG[1][id2][m1])];
               }
               if(IBS0Count >= HomHomCount * 0.42) continue; // C_Hom < 58%

               // 4th quarter
               for(m1 = stop1/4; m1 < stop1; m1++){
                  HomHomCount += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1]];
                  IBS0Count += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1] & (SG[1][id1][m1] ^ SG[1][id2][m1])];
               }
               if(IBS0Count >= HomHomCount * 0.375) continue; // C_Hom < 62.5%
            }else{   // more distant
               // 1st quarter
               for(HomHomCount = IBS0Count = m1 = 0, m2 = stop1/16; m1 < m2; m1++){
                  HomHomCount += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1]];
                  IBS0Count += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1] & (SG[1][id1][m1] ^ SG[1][id2][m1])];
               }
               if(IBS0Count >= HomHomCount * 0.55) continue; // C_Hom < 45%
               // 2nd quarter
               for(m1 = stop1/16, m2 = stop1/8; m1 < m2; m1++){
                  HomHomCount += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1]];
                  IBS0Count += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1] & (SG[1][id1][m1] ^ SG[1][id2][m1])];
               }
               if(IBS0Count >= HomHomCount * 0.5) continue; // C_Hom < 50%
               // 3rd quarter
               for(m1 = stop1/8, m2 = stop1/4; m1 < m2; m1++){
                  HomHomCount += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1]];
                  IBS0Count += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1] & (SG[1][id1][m1] ^ SG[1][id2][m1])];
               }
               if(IBS0Count >= HomHomCount * 0.45) continue; // C_Hom < 55%

               // 4th quarter
               for(m1 = stop1/4; m1 < stop1; m1++){
                  HomHomCount += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1]];
                  IBS0Count += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1] & (SG[1][id1][m1] ^ SG[1][id2][m1])];
               }
               if(IBS0Count >= HomHomCount * 0.45) continue; // C_Hom < 55%
            }

            m3 = (hetInOnePersonCount[id1][0] < hetInOnePersonCount[id2][0])?
                  hetInOnePersonCount[id1][0]: hetInOnePersonCount[id2][0];
            HetHomCount = stop1 * 32 - HomHomCount * 2 - hetInOnePersonCount[id1][0] - hetInOnePersonCount[id2][0]
               - missingInOnePersonCount[id1][0] - missingInOnePersonCount[id2][0];
            threshold = (0.5-kincutoff1)*m3 - HetHomCount*0.25;   // kinship < 0.0625
            if(IBS0Count >= threshold) continue;

            for(HetHomCount = het1Count = het2Count = m1 = 0; m1 < (stop1-1)/16+1; m1++){
               for(word = missingWordIndex[id1][m1] & missingWordIndex[id2][m1];
                  word; word &= (word-1)){   // MissMiss
                  m2 = m1*16+rightmost[word];   // m2 is the word with nonmissingness
                  for(word2 = ~(SG[0][id1][m2] | SG[0][id2][m2] | SG[1][id1][m2] | SG[1][id2][m2])&0xFFFF;
                       word2; word2 &= (word2-1))
                     HetHomCount ++;
               }
               for(word = missingWordIndex[id2][m1];  // HetMiss
                  word; word &= (word-1)){
                  m2 = m1*16+rightmost[word];
                  for(word2 = ~(SG[0][id1][m2] | SG[0][id2][m2] | SG[1][id2][m2]) & SG[1][id1][m2];
                      word2; word2 &= (word2-1))
                     het1Count ++;
               }
               for(word = missingWordIndex[id1][m1];  // MissHet
                  word; word &= (word-1)){
                  m2 = m1*16+rightmost[word];
                  for(word2 = ~(SG[0][id1][m2] | SG[0][id2][m2] | SG[1][id1][m2]) & SG[1][id2][m2];
                      word2; word2 &= (word2-1))
                     het2Count ++;
               }
            }
            het1Count = hetInOnePersonCount[id1][0] - het1Count;
            het2Count = hetInOnePersonCount[id2][0] - het2Count;
            HetHomCount = stop1 * 32 - HomHomCount * 2 - het1Count - het2Count +
               (-missingInOnePersonCount[id1][0] - missingInOnePersonCount[id2][0] + HetHomCount) * 2;
            m3 = (het1Count < het2Count? het1Count: het2Count);
            threshold = (0.5-kincutoff1)*m3 - HetHomCount*0.25;  // kinship < 0.0625
            if(IBS0Count >= threshold) continue;
//            rawrelativeCount[thread] ++;
           /*
            // Stage 2: few pairs
            for(HetHomCount -= (missingInOnePersonCount[id1][1]+missingInOnePersonCount[id2][1]),  // different non-missing genotypes
               m1 = start2; m1 < stop2; m1++) // different-genotpe Count
                  HetHomCount += oneoneCount[SG[0][id1][m1]^SG[0][id2][m1]];
            het1Count += hetInOnePersonCount[id1][1];
            het2Count += hetInOnePersonCount[id2][1];
            double cutoffDist = (2-kincutoff2*4) * (het1Count < het2Count? het1Count: het2Count)-IBS0Count*4;
            if(HetHomCount >= cutoffDist) continue; // too big distance already
            cutoffDist = (cutoffDist+IBS0Count*4-HetHomCount)*0.25;
            for(m1 = start2; m1 < (start2+stop2)/2; m1++)
               IBS0Count += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1] & (SG[1][id1][m1] ^ SG[1][id2][m1])];
            if(IBS0Count >= cutoffDist * 1.1) continue;
            for(m1 = (start2+stop2)/2; (m1 < stop2) && (IBS0Count < cutoffDist); m1++)
               IBS0Count += oneoneCount[SG[0][id1][m1] & SG[0][id2][m1] & (SG[1][id1][m1] ^ SG[1][id2][m1])];
            if(IBS0Count >= cutoffDist) continue;   // too big distance
            // now add the MissMissCount back

            for(m1 = start2/16, m3=0; m1 < (stop2-1)/16+1; m1++)
               for(word = missingWordIndex[id1][m1] & missingWordIndex[id2][m1];
                  word; word &= (word-1)){   // MissMiss
                  m2 = m1*16+rightmost[word];   // m2 is the word with nonmissingness
                  for(word2 = ~(SG[0][id1][m2] | SG[0][id2][m2] | SG[1][id1][m2] | SG[1][id2][m2])&0xFFFF;
                    word2; word2 &= (word2-1))
                     m3++;
               }
            HetHomCount += (m3<<1);

            for(m1 = start2/16, m3=0; m1 < (stop2-1)/16+1; m1++)
               for(word = missingWordIndex[id2][m1];  // HetMiss
                  word; word &= (word-1)){
                  m2 = m1*16+rightmost[word];
                  for(word2 = ~(SG[0][id1][m2] | SG[0][id2][m2] | SG[1][id2][m2]) & SG[1][id1][m2];
                    word2; word2 &= (word2-1))
                     m3 ++;
               }
            het1Count -= m3;
            HetHomCount += m3;

            for(m1 = start2/16, m3=0; m1 < (stop2-1)/16+1; m1++)
               for(word = missingWordIndex[id1][m1];  // MissHet
                  word; word &= (word-1)){
                  m2 = m1*16+rightmost[word];
                  for(word2 = ~(SG[0][id1][m2] | SG[0][id2][m2] | SG[1][id1][m2]) & SG[1][id2][m2];
                      word2; word2 &= (word2-1))
                     m3 ++;
               }
            het2Count -= m3;
            HetHomCount += m3;
            HetHetCount = (het1Count + het2Count - HetHomCount)/2;
            het1Count -= HetHetCount;
            het2Count -= HetHetCount;
            smaller = HetHetCount + (het1Count < het2Count? het1Count: het2Count);
            kinship = 0.5 - ((het1Count+het2Count)*0.25+IBS0Count)/smaller;
            if(kinship < kincutoff2) continue;
            */
            rpList[thread].Push(id1);
            rpList[thread].Push(id2);
         }  // end of id2
      }  // end of id1
   }  // end of loop blocking
#ifdef _OPENMP
}
#endif
   delete []loopIndex;
   for(int i = 0; i < idCount; i++){
      delete []missingInOnePersonCount[i];
      delete []hetInOnePersonCount[i];
      delete []missingWordIndex[i];
   }
   delete []missingInOnePersonCount;
   delete []hetInOnePersonCount;
   delete []missingWordIndex;
}

