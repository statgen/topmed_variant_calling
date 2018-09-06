//////////////////////////////////////////////////////////////////////
// ibdsegment.cpp
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
// August 23, 2018

#include <math.h>
#include "analysis.h"
#ifdef __ZLIB_AVAILABLE__
  #include <zlib.h>
#endif
#ifdef _OPENMP
  #include <omp.h>
#endif

void Engine::ComputeIBDSegment64Bit()
{
   if(idCount > 700000){
      printf("This version of KING --ibdseg supports up to %d samples.\n", 700000);
      printf("Please contact the KING author to allow an even larger sample size.\n");
      return;
   }
   printf("\nOptions in effect:\n");
   printf("\t--ibdseg\n");
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

   bool IBDvalidFlag = PreSegment(/*chrSeg, totalLength, segmessage*/);
   if(!IBDvalidFlag){
      printf("%s\n", (const char*)segmessage);
      printf("  Note chromosomal positions can be sorted conveniently using other tools such as PLINK.\n");
      return;
   }

   printf("IBD segment analysis starts at %s", currentTime());
   int segCount = chrSeg.Length()>>2;
   int id1, id2, segstart, segstop;
   double maxLength[2], prop[2];
   const char BLOCKSIZE=8;
   const char SHIFTSIZE=3;
   unsigned int blockCount = (idCount-1)/BLOCKSIZE+1;
   unsigned int loopIndexLength = (blockCount & 1) ? blockCount * ((blockCount+1)/2):
      (blockCount/2) * (blockCount+1);
   unsigned short int **loopIndex = new unsigned short int *[loopIndexLength];
   for(unsigned int k = 0; k < loopIndexLength; k++)
      loopIndex[k] = new unsigned short int [2];
   unsigned int index = 0;
   for(int i = 0; i < idCount; i += BLOCKSIZE)
      for(int j = i; j < idCount; j += BLOCKSIZE){
         loopIndex[index][0] = (i>>SHIFTSIZE);
         loopIndex[index][1] = (j>>SHIFTSIZE);
         index ++;
      }
#ifdef _OPENMP
   printf("%d CPU cores are used...\n", defaultMaxCoreCount);
#endif
#ifdef __ZLIB_AVAILABLE__
   const int coreCount = defaultMaxCoreCount;
   gzFile fps[coreCount];
   String outfile2(prefix);
   outfile2.Add(".segments.gz");
   fps[0] = gzopen((const char*)outfile2, "wb");
#ifdef _OPENMP
   StringArray outfiles(defaultMaxCoreCount);
   for(int c = 1; c < defaultMaxCoreCount; c++){
      outfiles[c].Copy(prefix);
      outfiles[c] += (c+1);
      outfiles[c].Add("$$$.segments.gz");
      fps[c] = gzopen((const char*)outfiles[c], "wb");
   }
#endif
#endif
   IntArray *allsegs = new IntArray [defaultMaxCoreCount];
   for(int c = 0; c < defaultMaxCoreCount; c++)
      allsegs[c].Dimension(0);
   int thread = 0;
   unsigned long long int word;//, word2, nomiss;
   int pbuffer = 0;
   char buffer[0x20000];
   char tempbuffer[1024];
   bool skipFlag;
   int icCount, cCount, localMin, localMax, minExtraBit, maxExtraBit;
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(id1, id2, icCount, cCount, minExtraBit, maxExtraBit, \
      segstart, segstop, word, localMin, localMax, skipFlag, \
      thread, maxLength, prop, pbuffer, buffer, tempbuffer)
{
#endif
   IntArray startPos[2][BLOCKSIZE][BLOCKSIZE];
   IntArray stopPos[2][BLOCKSIZE][BLOCKSIZE];
   IntArray startExtraBit[2][BLOCKSIZE][BLOCKSIZE];
   IntArray stopExtraBit[2][BLOCKSIZE][BLOCKSIZE];
   IntArray tempStart, tempStop, mergedStart, mergedStop, mergedBit, newchrSeg, newExtraBit;
#ifdef __ZLIB_AVAILABLE__
   IntArray interval[2];
#endif
#ifdef _OPENMP
   thread = omp_get_thread_num();
   pbuffer = 0;
   if(thread==0)
      pbuffer += sprintf(&buffer[pbuffer], "FID1\tID1\tFID2\tID2\tIBDType\tChr\tStartMB\tStopMB\tStartSNP\tStopSNP\t\tN_SNP\tLength\n");
   #pragma omp for
#endif
   for(unsigned int p = 0; p < loopIndexLength; p++){
      int i = (int)loopIndex[p][0]<<SHIFTSIZE;
      int iMax = i<idCount-BLOCKSIZE? i+BLOCKSIZE: idCount;
      int j = (int)loopIndex[p][1]<<SHIFTSIZE;
      int jMax = j<idCount-BLOCKSIZE? j+BLOCKSIZE: idCount;
      int jMin = j;
      for(int w = 0; w < 2; w++)
         for(int s = 0; s < BLOCKSIZE; s++)
            for(int t = 0; t < BLOCKSIZE; t++){
               startPos[w][s][t].Dimension(0);
               stopPos[w][s][t].Dimension(0);
               startExtraBit[w][s][t].Dimension(0);
               stopExtraBit[w][s][t].Dimension(0);
            }
      for(int seg = 0; seg < segCount; seg++){
         int chrsegMin = chrSeg[seg<<2];
         int chrsegMax = chrSeg[(seg<<2)|1];
         for(id1 = i; id1 < iMax; id1++){
            char ii = id1 - i;
            if(i == j) jMin = id1 + 1;
            for(id2 = jMin; id2 < jMax; id2++){
               minExtraBit = chrSeg[(seg<<2)|2];
               maxExtraBit = chrSeg[(seg<<2)|3];
               for(localMin = chrsegMin; localMin <= chrsegMax; localMin++){
                  for(; (localMin <= chrsegMax) &&
                  (((LG[0][id1][localMin] & LG[0][id2][localMin] &
                  (LG[1][id1][localMin] ^ LG[1][id2][localMin])) || //AA x aa
                  ((LG[1][id1][localMin] | LG[1][id2][localMin])==0))); // Neither is AA or Aa
                  localMin++); // keep passing if localMin includes AA x aa
                  if((localMin < chrsegMax) &&
                  (LG[0][id1][localMin+1] & LG[0][id2][localMin+1] &
                  (LG[1][id1][localMin+1] ^ LG[1][id2][localMin+1]))==0)
                     break;   // get out of the loop only when two 0 AA x aa words in a row
               }
               if(localMin >= chrsegMax ||
                  positions[(chrsegMax<<6)|0x3F] - positions[localMin<<6] < 2.5)
                  continue;// pass if the segment is shorter than 2.5MB
               if(localMin > chrsegMin){
                  word = LG[0][id1][localMin-1] & LG[0][id2][localMin-1] & (LG[1][id1][localMin-1] ^ LG[1][id2][localMin-1]);
                  if(word == 0){
                     localMin --;
                     minExtraBit = 0;
                  }else
                     for(minExtraBit = 0; (word & (1<<(63-minExtraBit))) == 0; minExtraBit++);
               }
               for(localMax = chrsegMax; localMax >= localMin; localMax--){
                  for(; (LG[0][id1][localMax] & LG[0][id2][localMax] &
                     (LG[1][id1][localMax] ^ LG[1][id2][localMax])) ||  //AA x aa
                     ((LG[1][id1][localMax] | LG[1][id2][localMax])==0);//None is AA or Aa
                     localMax--); // keep passing if localMax includes AA x aa
                  if((localMax > localMin) &&
                     (LG[0][id1][localMax-1] & LG[0][id2][localMax-1] &
                     (LG[1][id1][localMax-1] ^ LG[1][id2][localMax-1]))==0)
                     break;   // get out of the loop only when two 0 AA x aa words in a row
               }
               if(localMax < localMin + 2 ||   // 1 or 2 words are not enough for an IBD segment
                  positions[(localMax<<6)|0x3F] - positions[localMin<<6] < 2.5)
                  continue;// pass if the segment is shorter than 2.5MB
               if(localMax < chrsegMax){
                  word = LG[0][id1][localMax+1] & LG[0][id2][localMax+1] & (LG[1][id1][localMax+1] ^ LG[1][id2][localMax+1]);
                  if(word == 0){
                     localMax ++;
                     maxExtraBit = 0;
                  }else
                     for(maxExtraBit = 0; (word & (1<<maxExtraBit)) == 0; maxExtraBit++);
               }
               // start the IBD2 segment analysis
               char jj = id2 - j;
               tempStart.Dimension(0); tempStop.Dimension(0);
               for(int m = localMin; m <= localMax; m++){
                  for(; (m <= localMax) &&  // keep passing if m word includes IC or does not include C
                  ((((LG[0][id1][m] ^ LG[0][id2][m]) | (LG[1][id1][m] ^ LG[1][id2][m])) &
                  (LG[0][id1][m] | LG[1][id1][m]) & (LG[0][id2][m] | LG[1][id2][m])) ||
                  (LG[1][id1][m] & LG[1][id2][m] & (~(LG[0][id1][m]^LG[0][id2][m])))==0);
                  m++); //includes IC (Het x Hom or AA x aa) OR does not include any C (AA x AA or Het x Het)
                  for(cCount = 0, segstart = m; m <= localMax &&  // no IC
                  ((((LG[0][id1][m] ^ LG[0][id2][m]) | (LG[1][id1][m] ^ LG[1][id2][m])) &
                  (LG[0][id1][m] | LG[1][id1][m]) & (LG[0][id2][m] | LG[1][id2][m]))) == 0; m++)
                     cCount += popcount(LG[1][id1][m] & LG[1][id2][m] & (~(LG[0][id1][m]^LG[0][id2][m])));
                  if(cCount < 10 || positions[((m-1)<<6)|0x3F] - positions[segstart<<6] < 0.1)
                     continue;   // continue only when cCount>=10 AND >0.2Mb
                  for(word=0; (segstart >= localMin) && (word==0 || (word & (word-1)) == 0); segstart--)
                     word = (((LG[0][id1][segstart] ^ LG[0][id2][segstart]) |
                        (LG[1][id1][segstart] ^ LG[1][id2][segstart])) &
                        (LG[0][id1][segstart] | LG[1][id1][segstart]) &
                        (LG[0][id2][segstart] | LG[1][id2][segstart]));// IC
                  if(segstart < localMin && (word==0 || (word & (word-1))==0))
                     segstart ++;
                  else
                     segstart += 2; // now 0 or 1 inconsistency
                  tempStart.Push(segstart);
                  for(word = 0, m--; (m <= localMax) && (word == 0 || (word & (word-1)) == 0); m++)
                     word = (((LG[0][id1][m] ^ LG[0][id2][m]) |
                        (LG[1][id1][m] ^ LG[1][id2][m])) &
                        (LG[0][id1][m] | LG[1][id1][m]) &
                        (LG[0][id2][m] | LG[1][id2][m]));// IC
                  if(m > localMax && (word == 0 || (word & (word-1)) == 0))
                     m--;
                  else
                     m -= 2; // 0 or 1 inconsistency
                  tempStop.Push(m);
               }  // end of scan indexed by m
               skipFlag = false;
               newchrSeg.Dimension(0);
               newExtraBit.Dimension(0);
               int tempcount = tempStart.Length();
               if(tempcount == 0){  // no IBD2 segments
                  newchrSeg.Push(localMin);
                  newchrSeg.Push(localMax);
                  newExtraBit.Push(minExtraBit);
                  newExtraBit.Push(maxExtraBit);
                  skipFlag = true;
               }
               if(!skipFlag){
                  mergedStart.Dimension(0);
                  mergedStop.Dimension(0);
                  mergedBit.Dimension(0);
                  for(int t = 0; t < tempcount-1; t++){
                     double gap = positions[tempStart[t+1]<<6]-positions[(tempStop[t]<<6)|0x3F];
                     if(tempStart[t+1] - tempStop[t] < 3){
                        tempStart[t+1] = tempStart[t];// merge if 1 word in-between
                     }else if(gap < 5.0 && tempStart[t+1] - tempStop[t] < 100){
                        cCount = 0; // consistency C (AA x AA or AA x Aa) count
                        icCount = -2;   // inconsistency IC (AA x aa) count
                        for(int m = tempStop[t]+1; m < tempStart[t+1]; m++){
                           icCount += popcount((((LG[0][id1][m] ^ LG[0][id2][m]) | (LG[1][id1][m] ^ LG[1][id2][m])) &
                              (LG[0][id1][m] | LG[1][id1][m]) & (LG[0][id2][m] | LG[1][id2][m])));
                           cCount += popcount(LG[1][id1][m] & LG[1][id2][m] & (~(LG[0][id1][m]^LG[0][id2][m])));
                        }
                        if(cCount > icCount*3) // merge if C_IBD2 >
                           tempStart[t+1] = tempStart[t];
                        else if(positions[(tempStop[t]<<6)|0x3F] - positions[tempStart[t]<<6] >= 2.5){
                           mergedStart.Push(tempStart[t]);  // IBD2 segments need to be > 2.5MB
                           mergedStop.Push(tempStop[t]);
                        } // else discard the left interval
                     }else if(positions[(tempStop[t]<<6)|0x3F] - positions[tempStart[t]<<6] >= 2.5){
                        mergedStart.Push(tempStart[t]);  // IBD2 segments need to be > 2.5MB
                        mergedStop.Push(tempStop[t]);    // No gap to consider
                     } // else discard the left interval
                  }
                  if(positions[(tempStop[tempcount-1]<<6)|0x3F] - positions[tempStart[tempcount-1]<<6] >= 2.5){
                     mergedStart.Push(tempStart[tempcount-1]);  // IBD2 segments need to be > 2.5MB
                     mergedStop.Push(tempStop[tempcount-1]);
                  }
                  tempcount = mergedStart.Length();
                  if(tempcount == 0){  // no IBD2 segments
                     newchrSeg.Push(localMin);
                     newchrSeg.Push(localMax);
                     newExtraBit.Push(minExtraBit);
                     newExtraBit.Push(maxExtraBit);
                     skipFlag = true;
                  }else
                     for(int t = 0; t < tempcount; t++){
                        if(mergedStart[t] == localMin)
                           mergedBit.Push(minExtraBit);
                        else{
                           int m = mergedStart[t] - 1;
                           word = ((LG[0][id1][m] ^ LG[0][id2][m]) | (LG[1][id1][m] ^ LG[1][id2][m]))
                              & (LG[0][id1][m] | LG[1][id1][m]) & (LG[0][id2][m] | LG[1][id2][m]);
                           int bit = 0;
                           for(; bit < 63 && (word & (1<<(63-bit))) == 0; bit++);
                           mergedBit.Push(bit);
                        }
                        if(mergedStop[t] == localMax)
                           mergedBit.Push(maxExtraBit);
                        else{
                           int m = mergedStop[t] + 1;
                           word = ((LG[0][id1][m] ^ LG[0][id2][m]) | (LG[1][id1][m] ^ LG[1][id2][m]))
                              & (LG[0][id1][m] | LG[1][id1][m]) & (LG[0][id2][m] | LG[1][id2][m]);
                           int bit = 0;
                           for(; bit < 63 && (word & (1<<bit)) == 0; bit++);
                           mergedBit.Push(bit);
                        }
                     }
               }
               if(!skipFlag){
                  for(int t = 0; t < tempcount; t++){
                     startPos[1][ii][jj].Push(mergedStart[t]);
                     stopPos[1][ii][jj].Push(mergedStop[t]);
                     startExtraBit[1][ii][jj].Push(mergedStart[t]>chrsegMin? mergedBit[t*2]+1: mergedBit[t*2]);
                     stopExtraBit[1][ii][jj].Push(mergedStop[t]<chrsegMax? mergedBit[t*2+1]+1: mergedBit[t*2+1]);
                  }
                  if(mergedStart[0] > localMin){
                     newchrSeg.Push(localMin);
                     newExtraBit.Push(minExtraBit);
                     if(mergedBit[0]){
                        newchrSeg.Push(mergedStart[0]-2);
                        newExtraBit.Push(64-mergedBit[0]);
                     }else{
                        newchrSeg.Push(mergedStart[0]-1);
                        newExtraBit.Push(0);
                     }
                  }
                  for(int t = 0; t < tempcount-1; t++){
                     if(mergedBit[t*2+1]){
                        newchrSeg.Push(mergedStop[t]+2);
                        newExtraBit.Push(64 - mergedBit[t*2+1]);
                     }else{
                        newchrSeg.Push(mergedStop[t]+1);
                        newExtraBit.Push(0);
                     }
                     if(mergedBit[(t+1)*2]){
                        newchrSeg.Push(mergedStart[t+1]-2);
                        newExtraBit.Push(64 - mergedBit[(t+1)*2]);
                     }else{
                        newchrSeg.Push(mergedStart[t+1]-1);
                        newExtraBit.Push(0);
                     }
                  }
                  if(mergedStop[tempcount-1] < localMax){
                     if(mergedBit[tempcount*2-1]){
                        newchrSeg.Push(mergedStop[tempcount-1]+2);
                        newExtraBit.Push(64 - mergedBit[tempcount*2-1]);
                     }else{
                        newchrSeg.Push(mergedStop[tempcount-1]+1);
                        newExtraBit.Push(0);
                     }
                     newchrSeg.Push(localMax);
                     newExtraBit.Push(maxExtraBit);
                  }
               }
               // IBD1 starts here
               int newcount = newchrSeg.Length()/2;
               for(int s = 0; s < newcount; s++){
                  int chrsegMin = newchrSeg[s<<1];
                  int chrsegMax = newchrSeg[(s<<1)|1];
                  tempStart.Dimension(0);
                  tempStop.Dimension(0);
                  for(int m = chrsegMin; m <= chrsegMax; m++){
                     for(; (m <= chrsegMax) &&  // keep passing if m word
                     ((LG[0][id1][m] & LG[0][id2][m] & (LG[1][id1][m] ^ LG[1][id2][m])) ||
                     ((LG[1][id1][m] & LG[1][id2][m] & (LG[0][id1][m] | LG[0][id2][m]))==0));
                     m++); //includes IC (AA x aa) OR does not include any C (AA x AA or AA x Aa)
                     for(cCount = 0, segstart = m; m <= chrsegMax &&
                        (LG[0][id1][m] & LG[0][id2][m] & (LG[1][id1][m] ^ LG[1][id2][m]))==0; m++)
                        cCount += popcount(LG[1][id1][m] & LG[1][id2][m] & (LG[0][id1][m] | LG[0][id2][m]));
                     if(cCount < 10 || positions[((m-1)<<6)|0x3F] - positions[segstart<<6] < 0.2)
                        continue;   // continue only when cCount>=10 AND >0.2Mb
                     for(segstart--; segstart >= chrsegMin &&  // icCount==0
                        (LG[0][id1][segstart] & LG[0][id2][segstart] & (LG[1][id1][segstart] ^ LG[1][id2][segstart]))==0;
                        segstart--); // keep passing if free of IC (AA x aa)
                     tempStart.Push(segstart+1);
                     tempStop.Push(m-1);
                  }  // end of scan indexed by m
                  tempcount = tempStart.Length();
                  if(tempcount == 0) continue;
                  mergedStart.Dimension(0);
                  mergedStop.Dimension(0);
                  for(int t = 0; t < tempcount-1; t++){
                     double gap = positions[tempStart[t+1]<<6]-positions[(tempStop[t]<<6)|0x3F];
                     if(tempStart[t+1] - tempStop[t] < 3){
                        tempStart[t+1] = tempStart[t];// merge if 1 word in-between
                     }else if(gap < 5.0 && tempStart[t+1] - tempStop[t] < 100){
                        cCount = 0; // consistency C (AA x AA or AA x Aa) count
                        icCount = -2;   // inconsistency IC (AA x aa) count
                        for(int m = tempStop[t]+1; m < tempStart[t+1]; m++){
                           icCount += popcount(LG[0][id1][m] & LG[0][id2][m] & (LG[1][id1][m] ^ LG[1][id2][m]));
                           cCount += popcount(LG[1][id1][m] & LG[1][id2][m] & (LG[0][id1][m] | LG[0][id2][m]));
                        }
                        if(cCount > icCount*6) // merge if C_IBD1 > 85.7%
                           tempStart[t+1] = tempStart[t];
                        else if(positions[(tempStop[t]<<6)|0x3F] - positions[tempStart[t]<<6] > 2.5){
                           mergedStart.Push(tempStart[t]);  // IBD1 segments need to be > 2.5MB
                           mergedStop.Push(tempStop[t]);
                        } // else discard the left interval
                     }else if(positions[(tempStop[t]<<6)|0x3F] - positions[tempStart[t]<<6] > 2.5){
                        mergedStart.Push(tempStart[t]);  // IBD1 segments need to be > 2.5MB
                        mergedStop.Push(tempStop[t]);    // No gap to consider
                     } // else discard the left interval
                  }
                  if(positions[(tempStop[tempcount-1]<<6)|0x3F] - positions[tempStart[tempcount-1]<<6] > 2.5){
                     mergedStart.Push(tempStart[tempcount-1]);  // IBD1 segments need to be > 2.5MB
                     mergedStop.Push(tempStop[tempcount-1]);
                  }
                  tempcount = mergedStart.Length();
                  for(int t = 0; t < tempcount; t++){
                     startPos[0][ii][jj].Push(mergedStart[t]);
                     stopPos[0][ii][jj].Push(mergedStop[t]);
                  }
                  for(int t = 0; t < tempcount; t++){
                     int m = mergedStart[t];
                     if(m > chrsegMin){
                        m --; // AA x aa
                        word = LG[0][id1][m] & LG[0][id2][m] & (LG[1][id1][m] ^ LG[1][id2][m]);
                        int bit = 0;
                        for(; bit < 63 && (word & (1<<(63-bit))) == 0; bit++);
                        startExtraBit[0][ii][jj].Push(bit);
                     }else
                        startExtraBit[0][ii][jj].Push(newExtraBit[s*2]);
                     m = mergedStop[t];
                     if(m < chrsegMax){
                        m ++; // AA x aa
                        word = LG[0][id1][m] & LG[0][id2][m] & (LG[1][id1][m] ^ LG[1][id2][m]);
                        int bit = 0;
                        for(; bit < 63 && (word & (1<<bit)) == 0; bit++);
                        stopExtraBit[0][ii][jj].Push(bit);
                     }else
                        stopExtraBit[0][ii][jj].Push(newExtraBit[s*2+1]);
                  }// end of t
               }// end of s
            }  // end of id2
         }  // end of id1
      }//end of seg
      for(id1 = i; id1 < iMax; id1++){
         char ii = id1 - i;
         if(i == j) jMin = id1 + 1;
         for(id2 = jMin; id2 < jMax; id2++){
            char jj = id2 - j;
#ifdef __ZLIB_AVAILABLE__
            for(int k = 0; k < 2; k++)
               interval[k].Dimension(0);
#endif
            for(int k = 0; k < 2; k++){
               int newsegCount = startPos[k][ii][jj].Length();
               maxLength[k] = prop[k] = 0;
               if(newsegCount==0) continue;
               for(int seg = 0; seg < newsegCount; seg++){
                  segstart = (startPos[k][ii][jj][seg]<<6)-startExtraBit[k][ii][jj][seg];
                  if(stopPos[k][ii][jj][seg] == longCount-1)
                     segstop = markerCount - 1;
                  else
                     segstop = ((stopPos[k][ii][jj][seg]<<6)|0x3F)+stopExtraBit[k][ii][jj][seg];
                  double length = positions[segstop] - positions[segstart];
                  if(length <= 2.5) continue;
                  prop[k] += length;
                  if(length > maxLength[k]) maxLength[k] = length;
#ifdef __ZLIB_AVAILABLE__
                  interval[k].Push(segstart);
                  interval[k].Push(segstop);
#endif
               }
               prop[k] /= totalLength;
            }

            if(maxLength[0] < 10.0 && maxLength[1] < 10.0) continue;
            if(relativedegree == 1 &&  // pi < 0.355 && no IBD2
               int(-log(prop[1]*0.5+prop[0]*0.25)/log(2.0) - 0.5) > 1 && prop[1] < 0.08) continue;
            else if(relativedegree >1 && int(-log(prop[1]*0.5+prop[0]*0.25)/log(2.0) - 0.5) > relativedegree) continue;
            else if(relativedegree < 0 && int(-log(prop[1]*0.5+prop[0]*0.25)/log(2.0) - 0.5) <= -relativedegree) continue;
            allsegs[thread].Push(id1);
            allsegs[thread].Push(id2);
            allsegs[thread].Push((int(maxLength[0]*10+0.5)<<16) | int(prop[0]*10000+0.5));
            allsegs[thread].Push((int(maxLength[1]*10+0.5)<<16) | int(prop[1]*10000+0.5));
#ifdef __ZLIB_AVAILABLE__
            sprintf(tempbuffer, "%s\t%s\t%s\t%s",
               (const char*)ped[phenoid[id1]].famid, (const char*)ped[phenoid[id1]].pid,
               (const char*)ped[phenoid[id2]].famid, (const char*)ped[phenoid[id2]].pid);
            for(int k = 0; k < 2; k++)
               for(int s = 0; s < interval[k].Length()/2; s++){
                  pbuffer += sprintf(&buffer[pbuffer], "%s\t%s",
                     tempbuffer, k==0? "IBD1": "IBD2");
                  segstart = interval[k][s*2];
                  segstop = interval[k][s*2+1];
                  pbuffer += sprintf(&buffer[pbuffer], "\t%d\t%.3lf\t%.3lf\t%s\t%s\t%d\t%.1lf\n",
                     chromosomes[segstart], positions[segstart], positions[segstop],
                     (const char*)snpName[segstart], (const char*)snpName[segstop],
                     segstop-segstart+1, positions[segstop]-positions[segstart]);
               }
            if(pbuffer > 0xFFFF){  // buffer big enough for writing
               gzwrite(fps[thread], buffer, pbuffer);
               pbuffer = 0;
            }
#endif
         }  // end of id2
      }  // end of id1
   }  // end of loop blocking
#ifdef __ZLIB_AVAILABLE__
   if(pbuffer)
      gzwrite(fps[thread], buffer, pbuffer);
#endif
#ifdef _OPENMP
}  // extra bracket for omp
#endif                
   for(unsigned int k = 0; k < loopIndexLength; k++)
      delete []loopIndex[k];
   delete []loopIndex;
   pbuffer = 0;
   int items[2];
   String outfile(prefix);
   outfile.Add(".seg");
   FILE *fp = fopen(outfile, "wb");
   pbuffer += sprintf(&buffer[pbuffer], "FID1\tID1\tFID2\tID2\tMaxIBD1\tMaxIBD2\tIBD1Seg\tIBD2Seg\tPropIBD\tInfType\n");
   for(int c = 0; c < defaultMaxCoreCount; c++){
      int allsegCount = allsegs[c].Length()/4;
      for(int i = 0; i < allsegCount; i++){
         id1 = allsegs[c][i*4];
         id2 = allsegs[c][i*4+1];
         items[0] = allsegs[c][i*4+2];
         items[1] = allsegs[c][i*4+3];
         prop[0] = (items[0]&0xFFFF)*0.0001; // IBD1Seg
         prop[1] = (items[1]&0xFFFF)*0.0001; // IBD2Seg
         double pi = prop[1] + prop[0]*0.5;
         pbuffer += sprintf(&buffer[pbuffer], "%s\t%s\t%s\t%s\t%.1lf\t%.1lf\t%.4lf\t%.4lf\t%.4lf\t",
            (const char*)ped[phenoid[id1]].famid, (const char*)ped[phenoid[id1]].pid,
            (const char*)ped[phenoid[id2]].famid, (const char*)ped[phenoid[id2]].pid,
            (items[0]>>16)*0.1,  (items[1]>>16)*0.1,
            prop[0], prop[1], pi);
         if(prop[1] > 0.7)
            pbuffer += sprintf(&buffer[pbuffer], "Dup/MZ");
         else if(prop[0]+prop[1] > 0.96 || (prop[0]+prop[1] > 0.9 && prop[1] <= 0.08))
            pbuffer += sprintf(&buffer[pbuffer], "PO");
         else if(pi > 0.3535534 && prop[1] > 0.08)
            pbuffer += sprintf(&buffer[pbuffer], "FS");
         else if(pi > 0.1767767){
            if(pi > 0.32 && prop[1] > 0.15)
               pbuffer += sprintf(&buffer[pbuffer], "FS");
            else
               pbuffer += sprintf(&buffer[pbuffer], "2nd");
         }else if(pi > 0.08838835)
            pbuffer += sprintf(&buffer[pbuffer], "3rd");
         else if(pi > 0.04419417)
            pbuffer += sprintf(&buffer[pbuffer], "4th");
         else
            pbuffer += sprintf(&buffer[pbuffer], "UN");
         pbuffer += sprintf(&buffer[pbuffer], "\n");
         if(pbuffer > 0xFFFF){  // buffer big enough for writing
            fwrite(buffer, 1, pbuffer, fp);
            pbuffer = 0;
         }
      }
   }
   if(pbuffer>0)
      fwrite(buffer, 1, pbuffer, fp);
   fclose(fp);
   delete []allsegs;

#ifdef __ZLIB_AVAILABLE__
   for(int c = 0; c < defaultMaxCoreCount; c++)
      gzclose(fps[c]);
#ifdef _OPENMP
   FILE *fp0 = fopen(outfile2, "ab");
   for(int c = 1; c < defaultMaxCoreCount; c++){
      fp = fopen(outfiles[c], "rb");
      int count = fread(buffer, 1, 0x20000, fp);
      for(; count == 0x20000; count = fread(buffer, 1, 0x20000, fp))
         fwrite(buffer, 1, 0x20000, fp0);
      if(count)
         fwrite(buffer, 1, count, fp0);
      fclose(fp);
      remove(outfiles[c]);
   }
   fclose(fp0);
#endif
#endif
   printf("                       ends at %s", currentTime());
   printf("\nNote with relationship inference as the primary goal, the following filters are implemented:\n");
   printf("  Individual pairs with no long IBD segments (>10Mb) are excluded.\n");
   printf("  Super short IBD segments (<2.5Mb) are not counted here.\n");
   printf("Summary statistics of IBD segments for individual pairs saved in file %s\n", (const char*)outfile);
#ifdef __ZLIB_AVAILABLE__
   printf("IBD segments saved in a gzipped file %s\n", (const char*)outfile2);
#endif
}

void Engine::ROH()
{
   printf("\nOptions in effect:\n");
   printf("\t--roh\n");
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");

   bool IBDvalidFlag = PreSegment(/*chrSeg, totalLength, segmessage*/);
   if(!IBDvalidFlag){
      printf("%s\n", (const char*)segmessage);
      printf("  Note chromosomal positions can be sorted conveniently using tools such as PLINK.\n");
      return;
   }
   printf("Run of homozygosity analysis starts at %s", currentTime());
   int segCount = chrSeg.Length()/4;
   int cCount, icCount, segstart, segstop;
   unsigned long long int word;
   Vector rohprops(idCount);
   rohprops.Zero();
   Vector maxLengths(idCount);
   maxLengths.Zero();
   IntArray *allsegments = new IntArray[idCount];
   for(int id = 0; id < idCount; id++)
      allsegments[id].Dimension(0);
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
   private(word, segstart, segstop, cCount, icCount)
{
#endif
   IntArray tempStart, tempStop, mergedStart, mergedStop;
   IntArray startPos, stopPos, startExtraBit, stopExtraBit, interval;
#ifdef _OPENMP
   #pragma omp for
#endif
   for(int id = 0; id < idCount; id++){
      startPos.Dimension(0);
      stopPos.Dimension(0);
      startExtraBit.Dimension(0);
      stopExtraBit.Dimension(0);
      for(int seg = 0; seg < segCount; seg++){  // go through each segment
         int chrsegMin = chrSeg[seg<<2];
         int chrsegMax = chrSeg[(seg<<2)|1];
         int minExtraBit = chrSeg[(seg<<2)|2];
         int maxExtraBit = chrSeg[(seg<<2)|3];
         tempStart.Dimension(0);
         tempStop.Dimension(0);
         for(int m = chrsegMin; m <= chrsegMax; m++){
            for(; (m <= chrsegMax) &&  // keep passing if m word includes IC (Aa) or does not include any C (AA)
            (((~LG[0][id][m]) & LG[1][id][m]) || (LG[0][id][m] & LG[1][id][m])==0);
            m++); //includes IC (Aa) OR does not include any C (AA)
            for(cCount = 0, segstart = m; m <= chrsegMax &&
               ((~LG[0][id][m]) & LG[1][id][m])==0; m++)
               cCount += popcount(LG[0][id][m] & LG[1][id][m]);
            if(cCount < 10 || positions[((m-1)<<6)|0x3F] - positions[segstart<<6] < 0.2)
               continue;   // continue only when cCount>=10 AND >0.2Mb
            for(segstart--; segstart >= chrsegMin &&  // icCount==0
               ((~LG[0][id][segstart]) & LG[1][id][segstart])==0; segstart--);
            tempStart.Push(segstart+1);
            tempStop.Push(m-1);
         }  // end of scan indexed by m
         int tempcount = tempStart.Length();
         if(tempcount == 0) continue;
         mergedStart.Dimension(0);
         mergedStop.Dimension(0);
         for(int t = 0; t < tempcount-1; t++){
            double gap = positions[tempStart[t+1]<<6]-positions[(tempStop[t]<<6)|0x3F];
            if(tempStart[t+1] - tempStop[t] < 3){
               tempStart[t+1] = tempStart[t];// merge if 1 word in-between
            }else if(gap < 5.0 && tempStart[t+1] - tempStop[t] < 100){
               cCount = 0; // consistency C (AA) count
               icCount = -2;   // inconsistency IC (Aa) count
               for(int m = tempStop[t]+1; m < tempStart[t+1]; m++){
                  icCount += popcount((~LG[0][id][m]) & LG[1][id][m]);
                  cCount += popcount(LG[0][id][m] & LG[1][id][m]);
               }
               if(cCount > icCount*3) // merge if C_roh > 75%
                  tempStart[t+1] = tempStart[t];
               else if(positions[(tempStop[t]<<6)|0x3F] - positions[tempStart[t]<<6] > 2.5){
                  mergedStart.Push(tempStart[t]);  // ROH segments need to be > 2.5MB
                  mergedStop.Push(tempStop[t]);
               } // else discard the left interval
            }else if(positions[(tempStop[t]<<6)|0x3F] - positions[tempStart[t]<<6] > 2.5){
               mergedStart.Push(tempStart[t]);  // ROH segments need to be > 2.5MB
               mergedStop.Push(tempStop[t]);    // No gap to consider
            } // else discard the left interval
         }
         if(positions[(tempStop[tempcount-1]<<6)|0x3F] - positions[tempStart[tempcount-1]<<6] > 2.5){
            mergedStart.Push(tempStart[tempcount-1]);  // ROH segments need to be > 2.5MB
            mergedStop.Push(tempStop[tempcount-1]);
         }
         tempcount = mergedStart.Length();
         for(int t = 0; t < tempcount; t++){
            startPos.Push(mergedStart[t]);
            stopPos.Push(mergedStop[t]);
         }
         for(int t = 0; t < tempcount; t++){
            int m = mergedStart[t];
            if(m > chrsegMin){
               m --; // Aa
               word = LG[1][id][m] & (~LG[0][id][m]);
               int bit = 0;
               for(; bit < 63 && (word & (1<<(63-bit))) == 0; bit++);
               startExtraBit.Push(bit);
            }else
               startExtraBit.Push(minExtraBit);
            m = mergedStop[t];
            if(m < chrsegMax){
               m ++; // Aa
               word = LG[1][id][m] & (~LG[0][id][m]);
               int bit = 0;
               for(; bit < 63 && (word & (1<<bit)) == 0; bit++);
                  stopExtraBit.Push(bit);
            }else
               stopExtraBit.Push(maxExtraBit);
         }// end of a ROH segment indexed by t
      }  // end of a segment indexed by seg
      int newsegCount = startPos.Length();
      double prop = 0, maxLength = 0, length;
      interval.Dimension(0);
      for(int seg = 0; seg < newsegCount; seg++){
         segstart = (startPos[seg]<<6)-startExtraBit[seg];
         if(stopPos[seg] == longCount-1)
            segstop = markerCount - 1;
         else
            segstop = ((stopPos[seg]<<6)|0x3F)+stopExtraBit[seg];
         length = positions[segstop] - positions[segstart];
         if(length > maxLength)
            maxLength = length;
         if(length > 2.5) {
            prop += length;
            interval.Push(segstart);
            interval.Push(segstop);
         }
      }
      if(maxLength < 10.0) prop=0.0;
      else {
         prop /= totalLength;
         rohprops[id] = prop;
         maxLengths[id] = maxLength;
      }
      if(prop > 0.005)
         for(int s = 0; s < interval.Length()/2; s++){
            segstart = interval[s*2];
            segstop = interval[s*2+1];
            allsegments[id].Push(segstart);
            allsegments[id].Push(segstop);
         }
   }  // end of id
#ifdef _OPENMP
}
#endif
   IntArray excessiveROH(0);
   String outfile;
   outfile.Copy(prefix);
   outfile.Add(".roh");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "FID\tID\tFA\tMO\tSEX");
   if(ped.affectionCount)
      fprintf(fp, "\t%s",
      ped.affectionNames[0]=="DISEASEKING"?"DISEASE":(const char*)ped.affectionNames[0]);
   fprintf(fp, "\tMaxROH\tF_ROH\n");
   for(int id = 0; id < idCount; id++){
      fprintf(fp, "%s\t%s\t%s\t%s\t%d",
         (const char*)ped[phenoid[id]].famid,
         (const char*)ped[phenoid[id]].pid,
         ped[phenoid[id]].father?(const char*)ped[phenoid[id]].fatid:"0",
         ped[phenoid[id]].father?(const char*)ped[phenoid[id]].motid:"0",
         ped[phenoid[id]].sex);
      if(ped.affectionCount)
         fprintf(fp, "\t%d", ped[phenoid[id]].affections[0]);
      fprintf(fp,"\t%.1lf\t%.4lf\n", maxLengths[id], rohprops[id]);
      if(rohprops[id] > 0.0884) // parents are 2nd-degree or closer
         excessiveROH.Push(id);
   }
   fclose(fp);
   printf("Run of homozygosity analysis   ends at %s", currentTime());
   printf("\nNote with relationship inference as the primary goal, the following filters are implemented:\n");
   printf("  Individuals with no long ROH (>10Mb) are considered as outbred (with F_ROH=0).\n");
   printf("  Super short ROH (<2.5Mb) are not counted here.\n");
   if(excessiveROH.Length()){
      printf("\nThe following %d persons have excessive run of homozygosity (F>0.0884):\n", excessiveROH.Length());
      for(int i = 0; i < excessiveROH.Length(); i++){
         int id = excessiveROH[i];
         printf("  Person (%s, %s) has inbreeding coefficient = %.4lf\n",
            (const char*)ped[phenoid[id]].famid,
            (const char*)ped[phenoid[id]].pid, rohprops[id]);
      }
      printf("\n");
   }
   // print distribution of ROH
   printf("Run of homozygosity summary saved in file %s\n", (const char*)outfile);

#ifdef __ZLIB_AVAILABLE__
   gzFile fps;
   outfile = prefix;
   outfile.Add(".rohseg.gz");
   fps = gzopen((const char*)outfile, "wb");
   char buffer[0x20000], tempbuffer[1024];
   int pbuffer = 0;
   pbuffer += sprintf(&buffer[pbuffer], "FID\tID\tChr\tStartMB\tStopMB\tStartSNP\tStopSNP\t\tN_SNP\tLength\n");
   for(int id = 0; id < idCount; id++)
      if(allsegments[id].Length()){
         sprintf(tempbuffer, "%s\t%s",
         (const char*)ped[phenoid[id]].famid, (const char*)ped[phenoid[id]].pid);
         for(int s = 0; s < allsegments[id].Length()/2; s++){
            segstart = allsegments[id][s*2];
            segstop = allsegments[id][s*2+1];
            pbuffer += sprintf(&buffer[pbuffer], "%s\t%d\t%.3lf\t%.3lf\t%s\t%s\t%d\t%.1lf\n",
               tempbuffer,
               chromosomes[segstart], positions[segstart], positions[segstop],
               (const char*)snpName[segstart], (const char*)snpName[segstop],
               segstop-segstart+1, positions[segstop]-positions[segstart]);
         }
         if(pbuffer > 0xFFFF){  // buffer big enough for writing
            gzwrite(fps, buffer, pbuffer);
            pbuffer = 0;
         }
      }
   if(pbuffer)
      gzwrite(fps, buffer, pbuffer);
   gzclose(fps);
   printf("Run of homozygosity segments saved in file %s\n\n", (const char*)outfile);
#endif
   delete []allsegments;
}


void Engine::IBD2SegInOnePair64Bit(int id1, int id2, IntArray &chrSeg, double totalLength, double &ibd2prop, double &maxLength)
{
   int segCount = chrSeg.Length()/4;
   unsigned long long int word;
   int HetHetCount, HetHomCount, segstart, segstop, blockstart;
   IntArray startPos(0), stopPos(0);
   for(int seg = 0; seg < segCount; seg++){
      segstart = -1;
      int chrsegMin = chrSeg[seg<<2];
      int chrsegMax = chrSeg[(seg<<2)|1];
      for(int m = chrsegMin; m <= chrsegMax; m++){
         blockstart = m;
         for(HetHomCount = 0; m <= chrsegMax && HetHomCount < 5; m++){   // scan until 5 inconsistencies
            word = (LG[0][id1][m] ^ LG[0][id2][m]) & (LG[0][id1][m] | LG[1][id1][m]) & (LG[0][id2][m] | LG[1][id2][m]);
            HetHomCount += popcount(word);
         }
         m--;  // m is now blockstop
         if(segstart == -1 && m < blockstart+2)  // already too many inconsistencies with 1 or 2 words
            continue;
         HetHetCount = 0;
         for(int bm = blockstart; (bm <= m) && (HetHetCount < 95); bm++){
            word = ~(LG[0][id1][bm] | LG[0][id2][bm]) & LG[1][id1][bm] & LG[1][id2][bm];
            HetHetCount += popcount(word);
         }
         if(HetHetCount >= 95){   // IBD2 seg identified
            if(segstart == -1){   // start a new IBD2 seg
               segstart = blockstart;
               if(segstart > chrsegMin){
                  segstart --;
                  word = (LG[0][id1][segstart] ^ LG[0][id2][segstart]) &
                        (LG[0][id1][segstart] | LG[1][id1][segstart]) &
                        (LG[0][id2][segstart] | LG[1][id2][segstart]);
                  if(word==0 || ((word & (word-1))==0))   // one inconsistency only
                     for(segstart--; segstart>=chrsegMin &&
                           ((LG[0][id1][segstart] ^ LG[0][id2][segstart]) &
                           (LG[0][id1][segstart] | LG[1][id1][segstart]) &
                           (LG[0][id2][segstart] | LG[1][id2][segstart]))==0;
                           segstart--);
                  segstart++;
               }
            }
            segstop = m;
            if(m == chrsegMax){  // the end of chrseg
                  startPos.Push(segstart);
                  stopPos.Push(chrsegMax);
            }
         }else{   // HetHetCount < 95: not an IBD2 seg
            if(m == chrsegMax){
               if(segstart > -1)
                  startPos.Push(segstart);
               else
                  startPos.Push(blockstart);
               stopPos.Push(chrsegMax);
            }else if(segstart > -1){   // complete an existing IBD2 seg
                  startPos.Push(segstart);
                  for(segstop++; ((LG[0][id1][segstop] ^ LG[0][id2][segstop]) &
                     (LG[0][id1][segstop] | LG[1][id1][segstop]) &
                     (LG[0][id2][segstop] | LG[1][id2][segstop]))==0; segstop++);
                  word = (LG[0][id1][segstop] ^ LG[0][id2][segstop]) &
                     (LG[0][id1][segstop] | LG[1][id1][segstop]) &
                     (LG[0][id2][segstop] | LG[1][id2][segstop]);
                  if((word&(word-1))==0)  // one inconsistency only
                     for(segstop++; ((LG[0][id1][segstop] ^ LG[0][id2][segstop]) &
                     (LG[0][id1][segstop] | LG[1][id1][segstop]) &
                     (LG[0][id2][segstop] | LG[1][id2][segstop]))==0; segstop++);
                  segstop--;
                  stopPos.Push(segstop);
                  segstart = -1;
            }
         }  // end of deciding IBD2seg or not
      }  // end of one segment
   }  // end of all segments
   ibd2prop = maxLength = 0.0;
   int newsegCount = startPos.Length();
   if(newsegCount){
      double length;
      for(int seg = 0; seg < newsegCount; seg++){
         if(seg+1 < newsegCount && // two IBD2 Segments can be merged
               positions[startPos[seg+1]<<6]-positions[(stopPos[seg]<<6)|0x3F]<1.0 &&
               (positions[(stopPos[seg]<<6)|0x3F] - positions[stopPos[seg]<<6] > 1.0) &&
               (positions[(stopPos[seg+1]<<6)|0x3F] - positions[stopPos[seg+1]<<6] > 1.0) &&
               chromosomes[startPos[seg+1]<<6] == chromosomes[(stopPos[seg]<<6)|0x3F]){
               startPos[seg+1] = startPos[seg];
               continue;
         }
         if(stopPos[seg] == longCount-1)
               length = positions[markerCount-1] - positions[startPos[seg]<<6];
         else
               length = positions[(stopPos[seg]<<6)|0x3F] - positions[startPos[seg]<<6];
         ibd2prop += length;
         if(length > maxLength)
            maxLength = length;
      }
      if(maxLength < 10.0) ibd2prop=0.0;
      else
         ibd2prop /= totalLength;
   }
}

void Engine::ComputeIBD2Segment64Bit()
{
   printf("\nOptions in effect:\n");
   printf("\t--ibd2seg\n");
   if(Bit64Flag)
      printf("\t--sysbit 64\n");
   if(CoreCount)
      printf("\t--cpus %d\n", CoreCount);
   if(lessmemFlag)
      printf("\t--lessmem\n");
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");

//   IntArray chrSeg;
//   double totalLength;
//   String segmessage;
   bool IBDvalidFlag = PreSegment(/*chrSeg, totalLength, segmessage*/);
   if(!IBDvalidFlag){
      printf("%s\n", (const char*)segmessage);
      printf("  Note chromosomal positions can be sorted conveniently using other tools such as PLINK.\n");
      return;
   }
   printf("IBD2 segment analysis starts at %s", currentTime());
   int segCount = chrSeg.Length()>>2;
   int HetHetCount, HetHomCount, segstart, segstop, blockstart, id1, id2;
   FILE **fps;
   String outfile;
   outfile.Copy(prefix);
   outfile.Add(".seg2");
   const int BLOCKSIZE=4;
   IntArray loopIndex[2];
   loopIndex[0].Dimension(0);
   loopIndex[1].Dimension(0);
   for(int i = 0; i < idCount; i += BLOCKSIZE)
      for(int j = i; j < idCount; j += BLOCKSIZE){
         loopIndex[0].Push(i);
         loopIndex[1].Push(j);
      }
   int loopIndexLength = loopIndex[0].Length();
#ifdef _OPENMP
   printf("%d CPU cores are used.\n", defaultMaxCoreCount);
   fps = new FILE *[defaultMaxCoreCount];
   StringArray outfiles(defaultMaxCoreCount);
   for(int c = 1; c < defaultMaxCoreCount; c++){
      outfiles[c].Copy(prefix);
      outfiles[c] += (c+1);
      outfiles[c].Add("$$$.seg2");
      fps[c] = fopen(outfiles[c], "wt");
   }
#else
   fps = new FILE *[1];
#endif
   fps[0] = fopen(outfile, "wt");
   fprintf(fps[0], "FID1\tID1\tFID2\tID2\tN_IBD2\tMaxIBD2\tPr_IBD2\n");
   int thread = 0;
   unsigned long long int word;
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(HetHetCount, HetHomCount, id1, id2, segstart, segstop, word, \
      blockstart, thread)
{
#endif
   IntArray startPos[BLOCKSIZE][BLOCKSIZE];
   IntArray stopPos[BLOCKSIZE][BLOCKSIZE];
#ifdef _OPENMP
   thread = omp_get_thread_num();
   #pragma omp for
#endif
   for(int k = 0; k < loopIndexLength; k++){
      int i = loopIndex[0][k];
      int iMax = i<idCount-BLOCKSIZE? i+BLOCKSIZE: idCount;
      int j = loopIndex[1][k];
      int jMax = j<idCount-BLOCKSIZE? j+BLOCKSIZE: idCount;
      int jMin = j;
      for(int s = 0; s < BLOCKSIZE; s++)
         for(int t = 0; t < BLOCKSIZE; t++){
            startPos[s][t].Dimension(0);
            stopPos[s][t].Dimension(0);
         }
      for(int seg = 0; seg < segCount; seg++){
         for(id1 = i; id1 < iMax; id1++){
            char ii = id1 - i;
            if(i == j) jMin = id1 + 1;
            for(id2 = jMin; id2 < jMax; id2++){
               char jj = id2 - j;
               int chrsegMin = chrSeg[seg<<2];
               int chrsegMax = chrSeg[(seg<<2)|1];

               segstart = -1;
               for(int m = chrsegMin; m <= chrsegMax; m++){
                  blockstart = m;
                  for(HetHomCount = 0; m <= chrsegMax && HetHomCount < 5; m++){   // scan until 5 inconsistencies
                     word = (LG[0][id1][m] ^ LG[0][id2][m]) & (LG[0][id1][m] | LG[1][id1][m]) & (LG[0][id2][m] | LG[1][id2][m]);
                     HetHomCount += popcount(word);
                  }
                  m--;  // m is now blockstop
                  if(segstart == -1 && m < blockstart+2)  // already too many inconsistencies with 1 or 2 words
                     continue;
                  HetHetCount = 0;
                  for(int bm = blockstart; (bm <= m) && (HetHetCount < 95); bm++){
                     word = ~(LG[0][id1][bm] | LG[0][id2][bm]) & LG[1][id1][bm] & LG[1][id2][bm];
                     HetHetCount += popcount(word);
                  }
                  if(HetHetCount >= 95){   // IBD2 seg identified
                     if(segstart == -1){   // start a new IBD2 seg
                        segstart = blockstart;
                        if(segstart > chrsegMin){
                           segstart --;
                           word = (LG[0][id1][segstart] ^ LG[0][id2][segstart]) &
                              (LG[0][id1][segstart] | LG[1][id1][segstart]) &
                              (LG[0][id2][segstart] | LG[1][id2][segstart]);
                           if(word==0 || ((word & (word-1))==0))   // one inconsistency only
                              for(segstart--; segstart>=chrsegMin &&
                                 ((LG[0][id1][segstart] ^ LG[0][id2][segstart]) &
                                 (LG[0][id1][segstart] | LG[1][id1][segstart]) &
                                 (LG[0][id2][segstart] | LG[1][id2][segstart]))==0;
                                 segstart--);
                           segstart++;
                        }
                     }
                     segstop = m;
                     if(m == chrsegMax){  // the end of chrseg
                        startPos[ii][jj].Push(segstart);
                        stopPos[ii][jj].Push(chrsegMax);
                     }
                  }else{   // HetHetCount < 95: not an IBD2 seg
                     if(segstart > -1 && m == chrsegMax){
                        startPos[ii][jj].Push(segstart);
                        stopPos[ii][jj].Push(chrsegMax);
                     }else if(segstart > -1){   // complete an existing IBD2 seg
                        startPos[ii][jj].Push(segstart);
                        for(segstop++; ((LG[0][id1][segstop] ^ LG[0][id2][segstop]) &
                           (LG[0][id1][segstop] | LG[1][id1][segstop]) &
                           (LG[0][id2][segstop] | LG[1][id2][segstop]))==0; segstop++);
                        word = (LG[0][id1][segstop] ^ LG[0][id2][segstop]) &
                           (LG[0][id1][segstop] | LG[1][id1][segstop]) &
                           (LG[0][id2][segstop] | LG[1][id2][segstop]);
                        if((word&(word-1))==0)  // one inconsistency only
                           for(segstop++; ((LG[0][id1][segstop] ^ LG[0][id2][segstop]) &
                           (LG[0][id1][segstop] | LG[1][id1][segstop]) &
                           (LG[0][id2][segstop] | LG[1][id2][segstop]))==0; segstop++);
                        segstop--;
                        stopPos[ii][jj].Push(segstop);
                        segstart = -1;
                     }
                  }  // end of deciding IBD2seg or not
               }  // end of one segment
            }  // end of id2
         }  // end of id1
      }// end of seg
      for(id1 = i; id1 < iMax; id1++){
         char ii = id1 - i;
         if(i == j) jMin = id1 + 1;
         for(id2 = jMin; id2 < jMax; id2++){
            char jj = id2 - j;
            int newsegCount = startPos[ii][jj].Length();
            if(newsegCount==0) continue;
            double length, maxLength=0, prop=0;
            for(int seg = 0; seg < newsegCount; seg++){
               if(seg+1 < newsegCount && // two IBD2 Segments can be merged
                  positions[startPos[ii][jj][seg+1]<<6]-positions[(stopPos[ii][jj][seg]<<6)|0x3F]<1.0 &&
                  (positions[(stopPos[ii][jj][seg]<<6)|0x3F] - positions[stopPos[ii][jj][seg]<<6] > 1.0) &&
                  (positions[(stopPos[ii][jj][seg+1]<<6)|0x3F] - positions[stopPos[ii][jj][seg+1]<<6] > 1.0) &&
                  chromosomes[startPos[ii][jj][seg+1]<<6] == chromosomes[(stopPos[ii][jj][seg]<<6)|0x3F]){
                  startPos[ii][jj][seg+1] = startPos[ii][jj][seg];
                  continue;
               }
               if(stopPos[ii][jj][seg] == longCount-1)
                  length = positions[markerCount-1] - positions[startPos[ii][jj][seg]<<6];
               else
                  length = positions[(stopPos[ii][jj][seg]<<6)|0x3F] - positions[(startPos[ii][jj][seg]<<6)];
               if(length > 0.1) prop += length;
               if(length > maxLength) maxLength = length;
            }
            if(maxLength < 10.0) continue;
            prop /= totalLength;
            fprintf(fps[thread], "%s\t%s\t%s\t%s\t%d\t%.1lf\t%.4lf\n",
               (const char*)ped[phenoid[id1]].famid, (const char*)ped[phenoid[id1]].pid,
               (const char*)ped[phenoid[id2]].famid, (const char*)ped[phenoid[id2]].pid,
               newsegCount, maxLength, prop);
         }  // end of id2
      }  // end of id1
   }  // end of loop blocking
#ifdef _OPENMP
}  // extra bracket for omp
#endif
#ifdef _OPENMP
   for(int c = 0; c < defaultMaxCoreCount; c++)
      fclose(fps[c]);
   fps[0] = fopen(outfile, "at");
   char buffer[1024];
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
   printf("                        ends at %s", currentTime());
   printf("Summary of IBD2 segments (longer than 1MB) saved in file %s\n\n", (const char*)outfile);
}


void Engine::ComputeIBD2Segment()
{
//   IntArray chrSeg;
//   double totalLength;
//   String segmessage;
   bool IBDvalidFlag = PreSegment(/*chrSeg, totalLength, segmessage*/);
   if(!IBDvalidFlag){
      printf("%s\n", (const char*)segmessage);
      printf("  Note chromosomal positions can be sorted conveniently using tools such as PLINK.\n");
      return;
   }

   printf("\nOptions in effect:\n");
   printf("\t--ibd2seg\n");
   if(CoreCount)
      printf("\t--cpus %d\n", CoreCount);
   if(lessmemFlag)
      printf("\t--lessmem\n");
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");

   printf("\nIBD2 segment analysis starts at %s", currentTime());
   int segCount = chrSeg.Length()/4;
   // initilize segCounts here
   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];

   double prop, maxLength;
   int HetHetCount, HetCount;
   int start, stop, blockstart;
   int id1, id2;
   FILE **fps;
   String outfile;
   outfile.Copy(prefix);
   outfile.Add(".seg2");
   const int BLOCKSIZE=4;
   IntArray loopIndex[2];
   loopIndex[0].Dimension(0); loopIndex[1].Dimension(0);
   for(int i = 0; i < idCount; i += BLOCKSIZE)
      for(int j = i; j < idCount; j += BLOCKSIZE){
         loopIndex[0].Push(i);
         loopIndex[1].Push(j);
      }
#ifdef _OPENMP
   printf("%d CPU cores are used.\n", defaultMaxCoreCount);
   fps = new FILE *[defaultMaxCoreCount];
   StringArray outfiles(defaultMaxCoreCount);
   for(int c = 1; c < defaultMaxCoreCount; c++){
      outfiles[c].Copy(prefix);
      outfiles[c] += (c+1);
      outfiles[c].Add("$$$.seg2");
      fps[c] = fopen(outfiles[c], "wt");
   }
#else
   fps = new FILE *[1];
#endif
   fps[0] = fopen(outfile, "wt");
   fprintf(fps[0], "FID1\tID1\tFID2\tID2\tNHetHet\tNHet\tHetConc\tN_IBD2\tMaxIBD2\tPr_IBD2\n");
   int loopIndexLength = loopIndex[0].Length();
   int thread = 0;
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(HetHetCount, HetCount, id1, id2, start, stop, \
      blockstart, prop, maxLength, thread)
{
#endif
   IntArray startPos[BLOCKSIZE][BLOCKSIZE];
   IntArray stopPos[BLOCKSIZE][BLOCKSIZE];
   int totalHetHetCounts[BLOCKSIZE][BLOCKSIZE];
   int totalHetCounts[BLOCKSIZE][BLOCKSIZE];
#ifdef _OPENMP
   thread = omp_get_thread_num();
   #pragma omp for
#endif
   for(int k = 0; k < loopIndexLength; k++){
      int i = loopIndex[0][k];
      int iMax = i<idCount-BLOCKSIZE? i+BLOCKSIZE: idCount;
      int j = loopIndex[1][k];
      int jMax = j<idCount-BLOCKSIZE? j+BLOCKSIZE: idCount;
      int jMin = j;
      for(int s = 0; s < BLOCKSIZE; s++)
         for(int t = 0; t < BLOCKSIZE; t++){
            totalHetHetCounts[s][t] = totalHetCounts[s][t] = 0;
            startPos[s][t].Dimension(0);
            stopPos[s][t].Dimension(0);
         }
      for(int seg = 0; seg < segCount; seg++){
         for(id1 = i; id1 < iMax; id1++){
            char ii = id1 - i;
            if(i == j) jMin = id1 + 1;
            for(id2 = jMin; id2 < jMax; id2++){
               char jj = id2 - j;
               HetHetCount = HetCount = 0;
               start = -1;
               blockstart = chrSeg[seg*4];
               for(int m = chrSeg[seg*4]; m <= chrSeg[seg*4+1]; m++){
                  HetHetCount += oneoneCount[~GG[0][id1][m] & GG[1][id1][m] & ~GG[0][id2][m] & GG[1][id2][m]];
                  HetCount += oneoneCount[(~GG[0][id1][m] & GG[1][id1][m] & (GG[0][id2][m] | GG[1][id2][m])) | // HetNomiss
                     (~GG[0][id2][m] & GG[1][id2][m] & (GG[0][id1][m] | GG[1][id1][m])) ];   // NomissHet
                  if(HetCount >= 100){ // start a new block
                     if(HetHetCount >= 95){  // IBD2 segment
                        stop = m;
                        if(start == -1)  // a new IBD2 segment
                           start = blockstart;
                     }else if(start != -1){ // complete an old IBD2 segment
                        startPos[ii][jj].Push(start);
                        stopPos[ii][jj].Push(stop);
                        start = -1;
                     }
                     totalHetHetCounts[ii][jj] += HetHetCount;
                     totalHetCounts[ii][jj] += HetCount;
                     HetHetCount = HetCount = 0;
                     blockstart = m+1;
                  }else if(m==chrSeg[seg*4+1] && start != -1){   // complete an old IBD2 segment
                     startPos[ii][jj].Push(start);
                     if(HetHetCount > HetCount-2) // add this block to the IBD2 segment
                        stopPos[ii][jj].Push(m);
                     else  // without this block
                        stopPos[ii][jj].Push(stop);
                     totalHetHetCounts[ii][jj] += HetHetCount;
                     totalHetCounts[ii][jj] += HetCount;
                  }
               }  // end of word m
            }  // end of id2
         }  // end of id1
      }// end of seg
      for(id1 = i; id1 < iMax; id1++){
         char ii = id1 - i;
         if(i == j) jMin = id1 + 1;
         for(id2 = jMin; id2 < jMax; id2++){
            char jj = id2 - j;
            int newsegCount = startPos[ii][jj].Length();
            prop = 0;
            maxLength = 0;
            if(newsegCount==0) continue;
            double length;
            for(int seg = 0; seg < newsegCount; seg++){
               if(stopPos[ii][jj][seg] == shortCount-1)
                  length = positions[markerCount-1] - positions[startPos[ii][jj][seg]*16];
               else
                  length = positions[stopPos[ii][jj][seg]*16+15] - positions[startPos[ii][jj][seg]*16];
               prop += length;
               if(length > maxLength)
                  maxLength = length;
            }
            if(maxLength < 10.0) continue;
            prop /= totalLength;
            fprintf(fps[thread], "%s\t%s\t%s\t%s\t%d\t%d\t%.4lf\t%d\t%.1lf\t%.4lf\n",
               (const char*)ped[phenoid[id1]].famid, (const char*)ped[phenoid[id1]].pid,
               (const char*)ped[phenoid[id2]].famid, (const char*)ped[phenoid[id2]].pid,
               totalHetHetCounts[ii][jj], totalHetCounts[ii][jj], totalHetHetCounts[ii][jj]*1.0/totalHetCounts[ii][jj],
               newsegCount, maxLength, prop);
         }  // end of id2
      }  // end of id1
   }  // end of loop blocking
#ifdef _OPENMP
}  // extra bracket for omp
#endif
#ifdef _OPENMP
   for(int c = 0; c < defaultMaxCoreCount; c++)
      fclose(fps[c]);

   fps[0] = fopen(outfile, "at");
   char buffer[1024];
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
   printf("                        ends at %s", currentTime());
   printf("Summary of IBD2 segments (longer than 1MB) saved in file %s\n\n", (const char*)outfile);
}

bool Engine::PreSegment(/*IntArray & chrSeg, double & totalLength, String &segmessage, */bool printFlag)
{
   char buffer[256];
   chrSeg.Dimension(0);
   int startChr, startPos, stopPos;
   for(startPos = 0; startPos < markerCount && chromosomes[startPos]<1; startPos++);
   if(startPos == markerCount) return false;
   for(int m = 0; m < markerCount - 1; m++){
      if(positions[m+1] < positions[m] && chromosomes[m+1] == chromosomes[m]){
         sprintf(buffer, "Positions unsorted: %s at %d, %s at %d.",
            (const char*)snpName[m], int(positions[m]*1000000+0.5),
            (const char*)snpName[m+1], int(positions[m+1]*1000000+0.5) );
         segmessage = buffer;
         return false;
      }else if(chromosomes[m+1] < chromosomes[m]){ // chromosome not ordered
         sprintf(buffer, "Chromosomes unsorted: %s on chr %d, %s on chr %d.",
            (const char*)snpName[m], chromosomes[m],
            (const char*)snpName[m+1], chromosomes[m+1] );
         segmessage = buffer;
         return false;
      }
   }
   chrSeg.Push(startPos);
   startChr = chromosomes[startPos];
   for(int m = startPos+1; m < markerCount; m++)
      if(chromosomes[m] != startChr){ // chromosome change
         chrSeg.Push(m-1);
         chrSeg.Push(m);
         startChr = chromosomes[m];
      }else if(positions[m] - positions[m-1] > 1){  // gap too big
         chrSeg.Push(m-1);
         chrSeg.Push(m);
      }
   chrSeg.Push(markerCount-1);
   const int bit = Bit64==64? 64: 16;
   const int bitshift = Bit64==64? 6: 4;
   IntArray tmpchrSeg(0);
   int chrsegCount = chrSeg.Length()/2;
   for(int seg = 0; seg < chrsegCount; seg++){
      int extraBit = chrSeg[seg*2]%bit;
      if(extraBit){
         startPos = (chrSeg[seg*2]>>bitshift)+1;
         extraBit = bit - extraBit;
      }else
         startPos = (chrSeg[seg*2]>>bitshift);
      tmpchrSeg.Push(startPos);
      tmpchrSeg.Push(((chrSeg[seg*2+1]+1)>>bitshift)-1);
      tmpchrSeg.Push(extraBit);
      extraBit = (chrSeg[seg*2+1]+1)%bit;
      tmpchrSeg.Push(extraBit);
   }
   chrSeg.Dimension(0);
   chrsegCount = tmpchrSeg.Length()/4;
   int wordcount = 256>>bitshift;   // 4
   for(int seg = 0; seg < chrsegCount; seg++){
      if(tmpchrSeg[seg*4+1] >= tmpchrSeg[seg*4] + wordcount){
         for(startPos = tmpchrSeg[seg*4]; positions[(startPos+1)<<bitshift] - positions[startPos<<bitshift] > 10.0; startPos++);
         for(stopPos = tmpchrSeg[seg*4+1]; positions[(stopPos<<bitshift)-1] - positions[((stopPos-1)<<bitshift)-1] > 10.0; stopPos--);
         if(stopPos - startPos + 1 < wordcount) continue;
         if(positions[((stopPos+1)<<bitshift)-1] - positions[startPos<<bitshift] < 10.0) continue; // only deal with large segments
         chrSeg.Push(startPos);
         chrSeg.Push(stopPos);
         if(startPos!=tmpchrSeg[seg*4]) chrSeg.Push(0);
         else chrSeg.Push(tmpchrSeg[seg*4+2]);
         if(stopPos!=tmpchrSeg[seg*4+1]) chrSeg.Push(0);
         else chrSeg.Push(tmpchrSeg[seg*4+3]);
      }
   }
   if(chrSeg.Length()==0) {
      sprintf(buffer, "No informative IBD segments.");
      segmessage = buffer;
      return false;
   }

   String outfile;
   outfile.Copy(prefix);
   outfile.Add("allsegs.txt");
   FILE *fp = fopen(outfile, "wt");
   sprintf(buffer, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
      "Segment", "Chr", "StartMB", "StopMB", "Length", "N_SNP", "StartSNP", "StopSNP");
   fprintf(fp, "%s", buffer);

   totalLength = 0.0;
   chrsegCount = chrSeg.Length()/4;
   for(int seg = 0; seg < chrsegCount; seg++){
      startPos = (chrSeg[seg*4]<<bitshift)-chrSeg[seg*4+2];
      stopPos = ((chrSeg[seg*4+1]+1)<<bitshift)-1+chrSeg[seg*4+3];
      sprintf(buffer, "%d\t%d\t%.3lf\t%.3lf\t%.3lf\t%d\t%s\t%s\n",
         seg+1, chromosomes[startPos],
         positions[startPos], positions[stopPos], positions[stopPos] - positions[startPos],
         stopPos - startPos + 1,
         (const char*)snpName[startPos], (const char*)snpName[stopPos]);
      fprintf(fp, "%s", buffer);
      totalLength += positions[stopPos] - positions[startPos];
   }
   fclose(fp);

   IntArray AACounts, AaCounts, missingCounts;
   ComputeAlleleFrequency64Bit(AACounts, AaCounts, missingCounts, 64);
   int tempcount = AACounts.Length();
   int majorCount = 0;
   for(int m = 0; m < tempcount; m++)  // NAA/Naa > 3, i.e., MAF > 0.635
      if(AACounts[m] > (idCount - AACounts[m] - AaCounts[m] - missingCounts[m])*3)
         majorCount++;
   if(majorCount > 409.6)  // >10% of SNPs need to switch allele labels
      error("\nToo many first alleles as the major allele (~%.1lf%%). Please use plink1.9 --make-bed to regenerate the genotype data again.\n", majorCount/40.96);

   if(printFlag){
      printf("Total length of chromosomal segments usable for IBD segment analysis is %.1lf MB\n", totalLength);
      printf("  Information of these chromosomal segments can be found in file %s\n\n", (const char*)outfile);
   }
   if(totalLength < 100){  // total segment length shorter than 100MB not allowed
      sprintf(buffer, "Segments too short.");
      segmessage = buffer;
      return false;
   }
   segmessage.Clear();
   return true;
}


