//////////////////////////////////////////////////////////////////////
// SubsetKING.cpp
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

#include <math.h>
#include "analysis.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

void Engine::DuplicateInSubset64Bit(IntArray & pairList, IntArray & HetHetCounts,
   IntArray & DiffHomCounts, IntArray & HomHomCounts, IntArray & notMissingHetCounts)
{
   int pairCount = pairList.Length()/2;
   if(pairCount==0) return;
   HetHetCounts.Dimension(pairCount);
   DiffHomCounts.Dimension(pairCount);
   HomHomCounts.Dimension(pairCount);
   notMissingHetCounts.Dimension(pairCount);
   int HetHetCount, DiffHomCount, m1;
   unsigned long long int word, word1, word2;
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
      private(HetHetCount, DiffHomCount, m1, word, word1, word2)
#endif
   for(int p = 0; p < pairCount; p++){
      int id1 = pairList[p*2];
      int id2 = pairList[p*2+1];
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
      HetHetCounts[p] = HetHetCount;
      HomHomCounts[p] = HomHomCount;
      DiffHomCounts[p] = DiffHomCount;
      notMissingHetCounts[p] = notMissingHetCount;
   }  // end of p loop for pairs
}

void Engine::IBDSegInSubset64Bit(IntArray & pairList, Vector & ibdprops, Vector & maxLengths, Vector & ibd2props, Vector & maxLengths2, IntArray *ibd1segs, IntArray *ibd2segs)
{
   double ibdprop, maxLength;
   int segCount = (chrSeg.Length()>>2);
   int pairCount = pairList.Length()/2;
   ibdprops.Dimension(pairCount);
   maxLengths.Dimension(pairCount);
   ibd2props.Dimension(pairCount);
   maxLengths2.Dimension(pairCount);
   unsigned long long int word;
   int id1, id2, cCount, icCount, segstart, localMin, localMax, minExtraBit, maxExtraBit;
   bool skipFlag;
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(id1, id2, ibdprop, maxLength, word, segstart, \
      skipFlag, localMin, localMax, minExtraBit, maxExtraBit, cCount, icCount)
{
#endif
   IntArray startPos[2], stopPos[2], startExtraBit[2], stopExtraBit[2];
   IntArray tempStart, tempStop, mergedStart, mergedStop, mergedBit, newchrSeg, newExtraBit;
#ifdef _OPENMP
   #pragma omp for
#endif
   for(int p = 0; p < pairCount; p++){
      id1 = pairList[p<<1];
      id2 = pairList[(p<<1)|1];
      for(int k = 0; k < 2; k++){
         startPos[k].Dimension(0);
         stopPos[k].Dimension(0);
         startExtraBit[k].Dimension(0);
         stopExtraBit[k].Dimension(0);
      }
      for(int seg = 0; seg < segCount; seg++){
         int chrsegMin = localMin = chrSeg[seg<<2];
         int chrsegMax = localMax = chrSeg[(seg<<2)|1];
         minExtraBit = chrSeg[(seg<<2)|2];
         maxExtraBit = chrSeg[(seg<<2)|3];
         for(; localMin <= chrsegMax; localMin++){
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
         if(localMin >= chrsegMax || bp[(chrsegMax<<6)|0x3F] - bp[localMin<<6] < 2500000)
            continue;// pass if the segment is shorter than 2.5MB
         if(localMin > chrsegMin){
            word = LG[0][id1][localMin-1] & LG[0][id2][localMin-1] & (LG[1][id1][localMin-1] ^ LG[1][id2][localMin-1]);
            if(word == 0){
               localMin --;
               minExtraBit = 0;
            }else
               for(minExtraBit = 0; (word & (1<<(63-minExtraBit))) == 0; minExtraBit++);
         }
         for(; localMax >= localMin; localMax--){
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
            bp[(localMax<<6)|0x3F] - bp[localMin<<6] < 2500000)
            continue;// pass if the segment is shorter than 2.5MB
         if(localMax < chrsegMax){
            word = LG[0][id1][localMax+1] & LG[0][id2][localMax+1] & (LG[1][id1][localMax+1] ^ LG[1][id2][localMax+1]);
            if(word == 0){
               localMax ++;
               maxExtraBit = 0;
            }else
               for(maxExtraBit = 0; (word & (1<<maxExtraBit)) == 0; maxExtraBit++);
         }
         // IBD2 segment analysis starts here
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
            if(cCount < 10 || bp[((m-1)<<6)|0x3F] - bp[segstart<<6] < 100000)
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
               int gap = bp[tempStart[t+1]<<6]-bp[(tempStop[t]<<6)|0x3F];
               if(tempStart[t+1] - tempStop[t] < 3){
                  tempStart[t+1] = tempStart[t];// merge if 1 word in-between
               }else if(gap < 5000000 && tempStart[t+1] - tempStop[t] < 100){
                  cCount = 0; // consistency C (AA x AA or AA x Aa) count
                  icCount = -2;   // inconsistency IC (AA x aa) count
                  for(int m = tempStop[t]+1; m < tempStart[t+1]; m++){
                     icCount += popcount((((LG[0][id1][m] ^ LG[0][id2][m]) | (LG[1][id1][m] ^ LG[1][id2][m])) &
                        (LG[0][id1][m] | LG[1][id1][m]) & (LG[0][id2][m] | LG[1][id2][m])));
                     cCount += popcount(LG[1][id1][m] & LG[1][id2][m] & (~(LG[0][id1][m]^LG[0][id2][m])));
                  }
                  if(cCount > icCount*3) // merge if C_IBD2 >
                     tempStart[t+1] = tempStart[t];
                  else if(bp[(tempStop[t]<<6)|0x3F] - bp[tempStart[t]<<6] >= 2500000){
                     mergedStart.Push(tempStart[t]);  // IBD2 segments need to be > 2.5MB
                     mergedStop.Push(tempStop[t]);
                  } // else discard the left interval
               }else if(bp[(tempStop[t]<<6)|0x3F] - bp[tempStart[t]<<6] >= 2500000){
                  mergedStart.Push(tempStart[t]);  // IBD2 segments need to be > 2.5MB
                  mergedStop.Push(tempStop[t]);    // No gap to consider
               } // else discard the left interval
            }
            if(bp[(tempStop[tempcount-1]<<6)|0x3F] - bp[tempStart[tempcount-1]<<6] >= 2500000){
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
               startPos[1].Push(mergedStart[t]);
               stopPos[1].Push(mergedStop[t]);
               startExtraBit[1].Push(mergedStart[t]>chrsegMin? mergedBit[t*2]+1: mergedBit[t*2]);
               stopExtraBit[1].Push(mergedStop[t]<chrsegMax? mergedBit[t*2+1]+1: mergedBit[t*2+1]);
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

         // IBD1 segment analysis starts here
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
               if(cCount < 10 || bp[((m-1)<<6)|0x3F] - bp[segstart<<6] < 200000)
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
               int gap = bp[tempStart[t+1]<<6]-bp[(tempStop[t]<<6)|0x3F];
               if(tempStart[t+1] - tempStop[t] < 3){
                  tempStart[t+1] = tempStart[t];// merge if 1 word in-between
               }else if(gap < 5000000 && tempStart[t+1] - tempStop[t] < 100){
                  cCount = 0; // consistency C (AA x AA or AA x Aa) count
                  icCount = -2;   // inconsistency IC (AA x aa) count
                  for(int m = tempStop[t]+1; m < tempStart[t+1]; m++){
                     icCount += popcount(LG[0][id1][m] & LG[0][id2][m] & (LG[1][id1][m] ^ LG[1][id2][m]));
                     cCount += popcount(LG[1][id1][m] & LG[1][id2][m] & (LG[0][id1][m] | LG[0][id2][m]));
                  }
                  if(cCount > icCount*6) // merge if C_IBD1 > 85.7%
                     tempStart[t+1] = tempStart[t];
                  else if(bp[(tempStop[t]<<6)|0x3F] - bp[tempStart[t]<<6] > 2500000){
                     mergedStart.Push(tempStart[t]);  // IBD1 segments need to be > 2.5MB
                     mergedStop.Push(tempStop[t]);
                  } // else discard the left interval
               }else if(bp[(tempStop[t]<<6)|0x3F] - bp[tempStart[t]<<6] > 2500000){
                  mergedStart.Push(tempStart[t]);  // IBD1 segments need to be > 2.5MB
                  mergedStop.Push(tempStop[t]);    // No gap to consider
               } // else discard the left interval
            }
            if(bp[(tempStop[tempcount-1]<<6)|0x3F] - bp[tempStart[tempcount-1]<<6] > 2500000){
               mergedStart.Push(tempStart[tempcount-1]);  // IBD1 segments need to be > 2.5MB
               mergedStop.Push(tempStop[tempcount-1]);
            }
            tempcount = mergedStart.Length();
            for(int t = 0; t < tempcount; t++){
               startPos[0].Push(mergedStart[t]);
               stopPos[0].Push(mergedStop[t]);
            }
            for(int t = 0; t < tempcount; t++){
               int m = mergedStart[t];
               if(m > chrsegMin){
                  m --; // AA x aa
                  word = LG[0][id1][m] & LG[0][id2][m] & (LG[1][id1][m] ^ LG[1][id2][m]);
                  int bit = 0;
                  for(; bit < 63 && (word & (1<<(63-bit))) == 0; bit++);
                  startExtraBit[0].Push(bit);
               }else
                  startExtraBit[0].Push(newExtraBit[s*2]);
               m = mergedStop[t];
               if(m < chrsegMax){
                  m ++; // AA x aa
                  word = LG[0][id1][m] & LG[0][id2][m] & (LG[1][id1][m] ^ LG[1][id2][m]);
                  int bit = 0;
                  for(; bit < 63 && (word & (1<<bit)) == 0; bit++);
                  stopExtraBit[0].Push(bit);
               }else
                  stopExtraBit[0].Push(newExtraBit[s*2+1]);
            }// end of t
         }// end of s
      }//end of seg

      maxLength = 0.0;
      for(int k = 0; k < 2; k++){
         int newsegCount = startPos[k].Length();
         double localmaxLength = 0.0;
         ibdprop = 0.0;
         if(newsegCount){
            double length;
            for(int seg = 0; seg < newsegCount; seg++){
               localMin = (startPos[k][seg]<<6)-startExtraBit[k][seg];
               if(stopPos[k][seg] == longCount-1)
                  localMax = markerCount - 1;
               else
                  localMax = ((stopPos[k][seg]<<6)|0x3F)+stopExtraBit[k][seg];
               length = bp[localMax] - bp[localMin];
               if(length > 2500000) ibdprop += length;
               if(length > localmaxLength)
                  localmaxLength = length;
               if(ibd1segs && k==0){// IBD1
                  ibd1segs[p].Push(localMin);
                  ibd1segs[p].Push(localMax);
               }else if(ibd2segs && k==1){// IBD2
                  ibd2segs[p].Push(localMin);
                  ibd2segs[p].Push(localMax);
               }
            }
            ibdprop /= totalLength;
         }
         if(localmaxLength > maxLength) maxLength = localmaxLength;
         if(k==0){
            ibdprops[p] = ibdprop;
            maxLengths[p] = localmaxLength;
         }else{
            ibd2props[p] = ibdprop;
            maxLengths2[p] = localmaxLength;
         }
      }  // end of k
      if(maxLength < 10.0)
         ibdprops[p] = ibd2props[p] = maxLengths[p] = maxLengths2[p] = 0.0;
   }  // end of all pairs
#ifdef _OPENMP
}
#endif
}

void Engine::KinshipInSubset64Bit(IntArray & pairList, IntArray & HetHetCounts,
   IntArray & IBS0Counts, IntArray & het1Counts, IntArray & het2Counts, IntArray & HomHomCounts, IntArray & IBSCounts){
   int pairCount = pairList.Length()/2;
   if(pairCount==0) return;
   int id1, id2;
   HetHetCounts.Dimension(pairCount);
   HetHetCounts.Zero();
   IBS0Counts.Dimension(pairCount);
   IBS0Counts.Zero();
   het1Counts.Dimension(pairCount);
   het1Counts.Zero();
   het2Counts.Dimension(pairCount);
   het2Counts.Zero();
   HomHomCounts.Dimension(pairCount);
   HomHomCounts.Zero();
   IBSCounts.Dimension(pairCount);
   IBSCounts.Zero();
   unsigned long long int word, words[6], ibs0;
   const int BLOCKSIZE=8;   // cache blocking size for individual pairs
   const int CACHESIZE=512;  // cache blocking size for SNP words
   int HetHetCount[BLOCKSIZE], IBS0Count[BLOCKSIZE], het1Count[BLOCKSIZE];
   int het2Count[BLOCKSIZE], HomHomCount[BLOCKSIZE], IBSCount[BLOCKSIZE];
   if(pairCount < (idCount<<10))   // two few pairs, most commonly seen
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
   private(HetHetCount, IBS0Count, het1Count, het2Count, HomHomCount, IBSCount, id1, id2, word, words, ibs0)
#endif
   for(int blockp = 0; blockp < pairCount; blockp += BLOCKSIZE){
      int pMax = blockp > pairCount-BLOCKSIZE? pairCount: blockp+BLOCKSIZE;
      for(int p = blockp; p < pMax; p++){
         int pp = p - blockp;
         HetHetCount[pp] = IBS0Count[pp] = het1Count[pp] = het2Count[pp] = HomHomCount[pp] = IBSCount[pp]=0;
      }
      for(int blockm = 0; blockm < longCount; blockm += CACHESIZE){
         int mMax = (blockm > longCount-CACHESIZE) ? longCount: blockm+CACHESIZE;
         for(int p = blockp; p < pMax; p++){
            id1 = pairList[p*2];
            id2 = pairList[p*2+1];
            int pp = p - blockp;
            for(int m = 0; m < 6; m++) words[m] = 0;
            for(int m = blockm; m < mMax; m++){
               word = LG[1][id1][m] & LG[1][id2][m];  // AAAA, AaAa, or AAAa
               ibs0 = word & ~(LG[0][id1][m] | LG[0][id2][m]); // HetHet
               word = word - ((word>>1)&0x5555555555555555);
               word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
               word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
               words[5] += (word+(word>>8)) & 0x00FF00FF00FF00FF;
               word = ibs0 - ((ibs0>>1)&0x5555555555555555);   // HetHet
               word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
               word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
               words[0] += (word+(word>>8)) & 0x00FF00FF00FF00FF;
               word = (~LG[0][id1][m]) & LG[1][id1][m] & LG[0][id2][m];
               word = word - ((word>>1)&0x5555555555555555);   // HetHom
               word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
               word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
               words[1] += (word+(word>>8)) & 0x00FF00FF00FF00FF;
               word = LG[0][id1][m] & (~LG[0][id2][m]) & LG[1][id2][m];
               word = word - ((word>>1)&0x5555555555555555);   // HomHet
               word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
               word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
               words[2] += (word+(word>>8)) & 0x00FF00FF00FF00FF;
               word = LG[0][id1][m] & LG[0][id2][m];  // HomHom
               ibs0 = word & (LG[1][id1][m] ^ LG[1][id2][m]);  // IBS0
               word = word - ((word>>1)&0x5555555555555555);   // HomHom
               word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
               word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
               words[3] += (word+(word>>8)) & 0x00FF00FF00FF00FF;
               ibs0 = ibs0 - ((ibs0>>1)&0x5555555555555555);   // IBS0
               ibs0 = (ibs0&0x3333333333333333) + ((ibs0>>2)&0x3333333333333333);
               ibs0 = (ibs0+(ibs0>>4)) & 0x0F0F0F0F0F0F0F0F;
               words[4] += (ibs0+(ibs0>>8)) & 0x00FF00FF00FF00FF;
            }
            for(int m = 0; m < 6; m++){
               words[m] = (words[m]+(words[m]>>16)) & 0x0000FFFF0000FFFF;
               words[m] = (words[m]+(words[m]>>32)) & 0xFFFFFFFF;
            }
            HetHetCount[pp] += words[0];
            het1Count[pp] += words[1];
            het2Count[pp] += words[2];
            HomHomCount[pp] += words[3];
            IBS0Count[pp] += words[4];
            IBSCount[pp] += words[5];
         }  // end of p
      }  // end of blockm
      for(int p = blockp; p < pMax; p++){
         int pp = p - blockp;
         HetHetCounts[p] += HetHetCount[pp];
         IBS0Counts[p] += IBS0Count[pp];
         het1Counts[p] += het1Count[pp];
         het2Counts[p] += het2Count[pp];
         HomHomCounts[p] += HomHomCount[pp];
         IBSCounts[p] += IBSCount[pp];
      }  // end of p
   }  // end of blockp
   else{ // too many pairs, worth pre-computing miss & het, RARE
      int m1, m2, m3;
      unsigned long long int word1, word2;
      int *missingInOnePersonCount = new int [idCount];
      int *hetInOnePersonCount = new int [idCount];
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
      private(m1, m2, m3, word)
#endif
      for(int i = 0; i < idCount; i++){
         for(m1 = m2 = m3 = 0; m1 < longCount; m1++){   // not all non-missing
            m2 += popcount((~LG[0][i][m1]) & LG[1][i][m1]);  // Het
            for(word = ~(LG[0][i][m1] | LG[1][i][m1]); word; word &= (word-1), m3++);
         }
         hetInOnePersonCount[i] = m2;
         missingInOnePersonCount[i] = m3;
      }
      for(int p = 0; p < pairCount; p++){
         id1 = pairList[p*2];
         id2 = pairList[p*2+1];
         het1Counts[p] = hetInOnePersonCount[id1];
         het2Counts[p] = hetInOnePersonCount[id2];
         HetHetCounts[p] = missingInOnePersonCount[id1] + missingInOnePersonCount[id2] - (longCount<<6);
      }
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
      private(HetHetCount, IBS0Count, het1Count, het2Count, HomHomCount, \
      m1, m2, m3, id1, id2, word, word1, word2, ibs0)
#endif
      for(int blockp = 0; blockp < pairCount; blockp += BLOCKSIZE){
         int pMax = blockp > pairCount-BLOCKSIZE? pairCount: blockp+BLOCKSIZE;
         for(int p = blockp; p < pMax; p++){
            int pp = p - blockp;
            HetHetCount[pp] = IBS0Count[pp] = het1Count[pp] = het2Count[pp] = HomHomCount[pp] = 0;
         }
         for(int blockm = 0; blockm < longCount; blockm += CACHESIZE){
            int mMax = (blockm > longCount-CACHESIZE) ? longCount: blockm+CACHESIZE;
            for(int p = blockp; p < pMax; p++){
               id1 = pairList[p*2];
               id2 = pairList[p*2+1];
               int pp = p - blockp;
               word1 = word2 = 0;
               for(int m = blockm; m < mMax; m++){
                  word = LG[0][id1][m] & LG[0][id2][m]; // HomHom
                  ibs0 = word & (LG[1][id1][m] ^ LG[1][id2][m]);   // IBS0
                  word = word - ((word>>1)&0x5555555555555555);   // HomHom
                  word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
                  word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
                  word1 += (word+(word>>8)) & 0x00FF00FF00FF00FF;
                  ibs0 = ibs0 - ((ibs0>>1)&0x5555555555555555);   // IBS0
                  ibs0 = (ibs0&0x3333333333333333) + ((ibs0>>2)&0x3333333333333333);
                  ibs0 = (ibs0+(ibs0>>4)) & 0x0F0F0F0F0F0F0F0F;
                  word2 += (ibs0+(ibs0>>8)) & 0x00FF00FF00FF00FF;
               }
               word1 = (word1+(word1>>16)) & 0x0000FFFF0000FFFF;
               HomHomCount[pp] += (word1+(word1>>32)) & 0xFFFFFFFF;  // HomHom
               word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;
               IBS0Count[pp] += (word2+(word2>>32)) & 0xFFFFFFFF;    // IBS0
               m1 = m2 = m3 = 0;
               for(int m = blockm; m < mMax; m++){
                  word1 = ~(LG[0][id1][m] | LG[1][id1][m]);
                  word2 = ~(LG[0][id2][m] | LG[1][id2][m]);
                  for(word = word1 & word2; word; word &= (word-1), m3++);
                  for(word = (~LG[0][id1][m]) & LG[1][id1][m] & word2; word; word &= (word-1), m1++);
                  for(word = word1 & (~LG[0][id2][m]) & LG[1][id2][m]; word; word &= (word-1), m2++);
               }
               het1Count[pp] -= m1;
               het2Count[pp] -= m2;
               HetHetCount[pp] -= m3;
            }  // end of p
         }  // end of blockm
         for(int p = blockp; p < pMax; p++){
            int pp = p - blockp;
            HomHomCounts[p] += HomHomCount[pp];
            IBS0Counts[p] += IBS0Count[pp];
            het1Counts[p] += het1Count[pp];
            het2Counts[p] += het2Count[pp];
            HetHetCounts[p] += HetHetCount[pp];
         }  // end of p
      }  // end of blockp
// N = HomHomCount+het1Count+het2Count-HetHetCount+miss1+miss2-MissMiss
      for(int p = 0; p < pairCount; p++){
         HetHetCounts[p] += het1Counts[p] + het2Counts[p] + HomHomCounts[p];
         het1Counts[p] -= HetHetCounts[p];
         het2Counts[p] -= HetHetCounts[p];
      }
      delete []missingInOnePersonCount;
      delete []hetInOnePersonCount;
   }
}

void Engine::IBD2SegInSubset64Bit(IntArray & pairList, Vector & ibd2props, Vector & maxLengths)
{
   double ibd2prop, maxLength;
   int segCount = chrSeg.Length()/4;
   int pairCount = pairList.Length()/2;
   ibd2props.Dimension(pairCount);
   maxLengths.Dimension(pairCount);
   unsigned long long int word;
   int id1, id2, HetHetCount, HetHomCount, segstart, segstop, blockstart;
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(id1, id2, ibd2prop, maxLength, word, HetHetCount, HetHomCount, segstart, segstop, blockstart)
{
#endif
   IntArray startPos, stopPos;
#ifdef _OPENMP
   #pragma omp for
#endif
   for(int p = 0; p < pairCount; p++){
      id1 = pairList[p<<1];
      id2 = pairList[(p<<1)|1];
      startPos.Dimension(0);
      stopPos.Dimension(0);
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
      int newsegCount = startPos.Length();
      ibd2prop = maxLength = 0.0;
      if(newsegCount){
         double length;
         for(int seg = 0; seg < newsegCount; seg++){
            if(seg+1 < newsegCount && // two IBD2 Segments can be merged
               bp[startPos[seg+1]<<6]-bp[(stopPos[seg]<<6)|0x3F]<1000000 &&
               (bp[(stopPos[seg]<<6)|0x3F] - bp[stopPos[seg]<<6] > 1000000) &&
               (bp[(stopPos[seg+1]<<6)|0x3F] - bp[stopPos[seg+1]<<6] > 1000000) &&
               chromosomes[startPos[seg+1]<<6] == chromosomes[(stopPos[seg]<<6)|0x3F]){
               startPos[seg+1] = startPos[seg];
               continue;
            }
            if(stopPos[seg] == longCount-1)
               length = bp[markerCount-1] - bp[startPos[seg]<<6];
            else
               length = bp[(stopPos[seg]<<6)|0x3F] - bp[startPos[seg]<<6];
            ibd2prop += length;
            if(length > maxLength)
               maxLength = length;
         }
         if(maxLength < 10000000) ibd2prop=0.0;
         else
            ibd2prop /= totalLength;
      }
      ibd2props[p] = ibd2prop;
      maxLengths[p] = maxLength;
   }  // end of all pairs
#ifdef _OPENMP
}
#endif
}

void Engine::KinshipInSubset(IntArray & pairList, IntArray & HetHetCounts,
   IntArray & IBS0Counts, IntArray & het1Counts, IntArray & het2Counts, IntArray & HomHomCounts, IntArray & IBSCounts)
{
   int pairCount = pairList.Length()/2;
   if(pairCount==0) return;
   int id1, id2;
   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   HetHetCounts.Dimension(pairCount);
   HetHetCounts.Zero();
   IBS0Counts.Dimension(pairCount);
   IBS0Counts.Zero();
   het1Counts.Dimension(pairCount);
   het1Counts.Zero();
   het2Counts.Dimension(pairCount);
   het2Counts.Zero();
   HomHomCounts.Dimension(pairCount);
   HomHomCounts.Zero();
   IBSCounts.Dimension(pairCount);
   IBSCounts.Zero();
   const int BLOCKSIZE=64;   // cache blocking size for individual pairs
   const int CACHESIZE=1024;  // cache blocking size for SNP words
   int HetHetCount[BLOCKSIZE], IBS0Count[BLOCKSIZE], het1Count[BLOCKSIZE],
      het2Count[BLOCKSIZE], HomHomCount[BLOCKSIZE], IBSCount[BLOCKSIZE];
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
   private(HetHetCount, IBS0Count, het1Count, het2Count, HomHomCount, IBSCount, id1, id2)
#endif
   for(int blockp = 0; blockp < pairCount; blockp += BLOCKSIZE){
      int pMax = blockp > pairCount-BLOCKSIZE? pairCount: blockp+BLOCKSIZE;
      for(int p = blockp; p < pMax; p++){
         int pp = p - blockp;
         HetHetCount[pp] = IBS0Count[pp] = het1Count[pp] = het2Count[pp] = HomHomCount[pp] = IBSCount[pp] = 0;
      }
      for(int blockm = 0; blockm < shortCount; blockm += CACHESIZE){
         int mMax = (blockm > shortCount-CACHESIZE) ? shortCount: blockm+CACHESIZE;
         for(int p = blockp; p < pMax; p++){
            id1 = pairList[p*2];
            id2 = pairList[p*2+1];
            int pp = p - blockp;
            for(int m = blockm; m < mMax; m++){
               HetHetCount[pp] += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
               IBS0Count[pp] += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
               het1Count[pp] += oneoneCount[(~GG[0][id1][m]) & GG[1][id1][m] & GG[0][id2][m]]; // Het1Hom2
               het2Count[pp] += oneoneCount[GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]]; // Hom1Het2
               HomHomCount[pp] += oneoneCount[GG[0][id1][m] & GG[0][id2][m]]; // HomHom
               IBSCount[pp] += oneoneCount[GG[1][id1][m] & GG[1][id2][m]];
            }  // end of m
         }  // end of p
      }  // end of blockm
      for(int p = blockp; p < pMax; p++){
         int pp = p - blockp;
         HetHetCounts[p] += HetHetCount[pp];
         IBS0Counts[p] += IBS0Count[pp];
         het1Counts[p] += het1Count[pp];
         het2Counts[p] += het2Count[pp];
         HomHomCounts[p] += HomHomCount[pp];
         IBSCounts[p] += IBSCount[pp];
      }  // end of p
   }  // end of blockp
}

void Engine::IBD2SegInSubset(IntArray & pairList, Vector & ibd2props, Vector & maxLengths)
{
   int id1, id2;
   double ibd2prop, maxLength;
   int segCount = chrSeg.Length()/4;
   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   int pairCount = pairList.Length()/2;
   ibd2props.Dimension(pairCount);
   maxLengths.Dimension(pairCount);
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(id1, id2, ibd2prop, maxLength)
{
#endif
   IntArray startPos, stopPos;
   int localHetHetCount, localHetCount, totalHetHetCount, totalHetCount;
   int start, stop, blockstart;
#ifdef _OPENMP
   #pragma omp for
#endif
   for(int p = 0; p < pairCount; p++){
      id1 = pairList[p*2];
      id2 = pairList[p*2+1];
      totalHetHetCount = totalHetCount = 0;
      startPos.Dimension(0);
      stopPos.Dimension(0);
      for(int seg = 0; seg < segCount; seg++){
         localHetHetCount = localHetCount = 0;
         start = -1;
         blockstart = chrSeg[seg*4];
         for(int m = chrSeg[seg*4]; m <= chrSeg[seg*4+1]; m++){
               localHetHetCount += oneoneCount[~GG[0][id1][m] & GG[1][id1][m] & ~GG[0][id2][m] & GG[1][id2][m]];
               localHetCount += oneoneCount[(~GG[0][id1][m] & GG[1][id1][m] & (GG[0][id2][m] | GG[1][id2][m])) | // HetNomiss
                  (~GG[0][id2][m] & GG[1][id2][m] & (GG[0][id1][m] | GG[1][id1][m])) ];   // NomissHet
               if(localHetCount >= 100){ // start a new block
                  if(localHetHetCount >= 95){  // IBD2 segment
                     stop = m;
                     if(start == -1){  // a new IBD2 segment
                        start = blockstart;
                        for(start--; start>=chrSeg[seg*4] &&
                        ((~GG[0][id1][start] & GG[1][id1][start] & GG[0][id2][start])
                        | (~GG[0][id2][start] & GG[1][id2][start] & GG[0][id1][start]))==0; start--);
                        if(start>=chrSeg[seg*4] && // one error only
                        oneoneCount[(~GG[0][id1][start] & GG[1][id1][start] & GG[0][id2][start])
                        | (~GG[0][id2][start] & GG[1][id2][start] & GG[0][id1][start])]==1)
                           for(start--; start>=chrSeg[seg*4] &&
                           ((~GG[0][id1][start] & GG[1][id1][start] & GG[0][id2][start])
                           | (~GG[0][id2][start] & GG[1][id2][start] & GG[0][id1][start]))==0; start--);
                        start++;
                     }
                  }else if(start != -1){ // complete an old IBD2 segment
                     startPos.Push(start);
                     for(stop++; ((~GG[0][id1][stop] & GG[1][id1][stop] & GG[0][id2][stop])
                     | (~GG[0][id2][stop] & GG[1][id2][stop] & GG[0][id1][stop]))==0; stop++);
                     if(oneoneCount[(~GG[0][id1][stop] & GG[1][id1][stop] & GG[0][id2][stop])
                     | (~GG[0][id2][stop] & GG[1][id2][stop] & GG[0][id1][stop])]==1)  // one error only
                        for(stop++; ((~GG[0][id1][stop] & GG[1][id1][stop] & GG[0][id2][stop])
                        | (~GG[0][id2][stop] & GG[1][id2][stop] & GG[0][id1][stop]))==0; stop++);
                     stop--;
                     stopPos.Push(stop);
                     start = -1;
                  }
                  totalHetHetCount += localHetHetCount;
                  totalHetCount += localHetCount;
                  localHetHetCount = localHetCount = 0;
                  blockstart = m+1;
               }else if(m==chrSeg[seg*4+1] && start != -1){   // complete an old IBD2 segment
                  startPos.Push(start);
                  if(localHetHetCount > localHetCount-2) // add this block to the IBD2 segment
                     stopPos.Push(m);
                  else{  // without this block
                     for(stop++; ((~GG[0][id1][stop] & GG[1][id1][stop] & GG[0][id2][stop])
                     | (~GG[0][id2][stop] & GG[1][id2][stop] & GG[0][id1][stop]))==0; stop++);
                     stop--;
                     stopPos.Push(stop);
                  }
                  totalHetHetCount += localHetHetCount;
                  totalHetCount += localHetCount;
               }
         }  // end of one segment
      }// end of all segments

      int newsegCount = startPos.Length();
      ibd2prop = maxLength = 0.0;
      if(newsegCount){
            double length;
            for(int seg = 0; seg < newsegCount; seg++){
               if(seg+1 < newsegCount && // two IBD2 Segments can be merged
                  chromosomes[startPos[seg+1]*16] == chromosomes[stopPos[seg]*16+15] &&
                  bp[startPos[seg+1]*16]-bp[stopPos[seg]*16+15]<1000000){
                  startPos[seg+1] = startPos[seg];
                  continue;
               }
               if(stopPos[seg] == shortCount-1)
                  length = bp[markerCount-1] - bp[startPos[seg]*16];
               else
                  length = bp[stopPos[seg]*16+15] - bp[startPos[seg]*16];
               if(length > 100000) ibd2prop += length;
               if(length > maxLength) maxLength = length;
            }
            if(maxLength < 10000000) ibd2prop=0.0;
            else
               ibd2prop /= totalLength;
      }
      ibd2props[p] = ibd2prop;
      maxLengths[p] = maxLength;
   }  // end of all pairs
#ifdef _OPENMP
}
#endif
}

