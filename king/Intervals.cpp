//////////////////////////////////////////////////////////////////////
// Intervals.cpp
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
// March 13, 2019

#include "Intervals.h"
#include "QuickIndex.h"

void SegmentUnion(IntArray &A, IntArray &B, IntArray &C)
{
   if((A.Length()&1) || (B.Length()&1)) return;
   IntArray AB=A; AB.Append(B);
   QuickIndex Idx(AB);
   C.Dimension(0);
   int intervalCount = AB.Length()/2;
   for(int i = 0; i < intervalCount; i++){
      for(; (i < intervalCount) && (Idx[i*2+1] & 1); i++){C.Push(Idx[i*2]); C.Push(Idx[i*2+1]);}  // 1010...10
      if(i < intervalCount){
         C.Push(Idx[i*2]);  // 00
         for(i++; (Idx[i*2+1] & 1)==0; i++);// 0101...01
         C.Push(Idx[i*2]+1);    // 11
      }
   }
   int index = 0;
   intervalCount = C.Length()/2;
   int stop;
   for(int i = 0; i < intervalCount; i++){
      C[index*2] = AB[C[i*2]];
      stop = i;
      for(int j = i+1; j < intervalCount && AB[C[j*2]] == AB[C[j*2-1]]; j++) stop = j;
      C[index*2+1] = AB[C[stop*2+1]];
      index++;
      i += (stop-i);
   }
   if(index < intervalCount) C.Dimension(index*2);
   printf("First time use of this function only...\n");
   A.Print();
   B.Print();
   C.Print();
}

void SegmentIntersect(IntArray &A, IntArray &B, IntArray &C)
{
   if((A.Length()&1) || (B.Length()&1)) return;
   IntArray AB=A; AB.Append(B);
   QuickIndex Idx(AB);
   C.Dimension(0);
   int intervalCount = AB.Length()/2;
   for(int i = 0; i < intervalCount; i++){
      for(; (i < intervalCount) && (Idx[i*2+1] & 1); i++);  // 1010...10
      if(i < intervalCount){
         C.Push(Idx[i*2+1]);  // 00
         for(i++; (Idx[i*2+1] & 1)==0; i++){C.Push(Idx[i*2]); C.Push(Idx[i*2+1]);}// 0101...01
         C.Push(Idx[i*2]);    // 11
      }
   }
   int index = 0;
   intervalCount = C.Length()/2;
   for(int i = 0; i < intervalCount; i++)
      if(AB[C[i*2]] < AB[C[i*2+1]]){
         C[index*2] = AB[C[i*2]];
         C[index*2+1] = AB[C[i*2+1]];
         index++;
      }
   if(index < intervalCount) C.Dimension(index*2);
}

double RoRP(IntArray &RP1, IntArray &RP2, IntArray &R1R2, IntArray &positionBP)
{
   IntArray join2; SegmentIntersect(RP1, RP2, join2);
   double join2Length = SegmentLength(join2, positionBP); // Length of PR1 & PR2
   if(join2Length < 50) return -1; // pi1 < 0.02, cannot be related
   IntArray join3; SegmentIntersect(join2, R1R2, join3);
   double join3Length = SegmentLength(join3, positionBP); // Length of PR1 & PR2 & R1R2
   return join3Length / join2Length;
}

double SegmentLength(IntArray &segs, IntArray &positionBP)
{
   double length = 0.0;
   int segCount = segs.Length()/2;
   for(int i = 0; i < segCount; i++)
      length += positionBP[segs[i*2+1]] - positionBP[segs[i*2]];
   return length;
}

double JoinLength(IntArray &A, IntArray &B, IntArray &positionBP)
{
   IntArray C;
   SegmentIntersect(A, B, C);
   return SegmentLength(C, positionBP);
}

/*
void SegmentUnion(IntArray &A, IntArray &B, IntArray &C)
{
   C.Dimension(0);
   IntArray AB=A; AB.Append(B);
   int ABCount = AB.Length();
   if(ABCount==0) return;
   QuickIndex Idx(AB);
   int c = 0;
   for(int i = 0; i < ABCount; i++){
      bool IsEnd = Idx[i] & 1;   // (Even, Odd) where Odd for end of an interval
      if(!IsEnd){
         if((c++)==0) C.Push(AB[Idx[i]]); // Start ...
      }else{
         if((c--)==1) C.Push(AB[Idx[i]]); // Start Stop ...
      }
   }
}
*/

/*
   int ABCount = AB.Length();
   if(ABCount==0) return;
   int c = 0;
   for(int i = 0; i < ABCount; i++){
      bool IsEnd = Idx[i] & 1;   // (Even, Odd) where Odd for end of an interval
      if(!IsEnd){
         if((c++)==1) C.Push(AB[Idx[i]]); // Start Start ...
      }else{
         if((c--)==2) C.Push(AB[Idx[i]]); // Start Start Stop ...
      }
   }
*/
/*
void SegmentIntersect(IntArray &A, IntArray &B, IntArray &C)
{
   C.Dimension(0);
   if((A.Length()&1) || (B.Length()&1)) return;
   IntArray AB=A; AB.Append(B);
   QuickIndex Idx(AB);
   int intervalCount = AB.Length()/2;
   int start, stop;
   for(int i = 0; i < intervalCount; i++){
      for(; (i < intervalCount) && (Idx[i*2+1] & 1); i++);  // 001010...10
      if(i < intervalCount){
         start = AB[Idx[i*2+1]];
         for(i++; (Idx[i*2+1] & 1)==0; i++){  // 0101...01
            stop = AB[Idx[i*2]];
            if(stop > start){C.Push(start); C.Push(stop);}
            start = AB[Idx[i*2+1]];
         }
         stop = AB[Idx[i*2]];                // 11
         if(stop > start){C.Push(start); C.Push(stop);}
      }
   }
}
*/

