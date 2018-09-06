//////////////////////////////////////////////////////////////////////
// diseaseGEE.h
// Author: Wei-Min Chen
// March 16, 2005

#ifndef __diseaseGEE_H__
#define __diseaseGEE_H__

#include "Pedigree.h"
#include "IntArray.h"
#include "MathMatrix.h"
#include "MathVector.h"
#include "MathCholesky.h"
#include "VCGEE.h"

class GEE_DIS: public GEE{
//   void constraint(void){}
   void RefreshD(int f);
//   void RefreshOD(int f){}
public:
//   double OR[6];
//   double rho[6]; // correlation between relative pair
   IntArray * diseases;
   int disease;
   IntArray mCovariate;
   Vector *resid;
   void solve();
   GEE_DIS(Pedigree & pedigree);
   ~GEE_DIS();
   void InitCoef();
   void summary(){}
   void print();
};

#endif


