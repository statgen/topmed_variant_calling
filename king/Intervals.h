#ifndef __Intervals_h__
#define __Intervals_h__
#include "IntArray.h"

double RoRP(IntArray &RP1, IntArray &RP2, IntArray &R1R2, IntArray &positionBP);
double SegmentLength(IntArray &segs, IntArray &positionBP);
double JoinLength(IntArray &A, IntArray &B, IntArray &positionBP);
void SegmentIntersect(IntArray &A, IntArray &B, IntArray &C);
void SegmentUnion(IntArray &A, IntArray &B, IntArray &C);

#endif
