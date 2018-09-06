//////////////////////////////////////////////////////////////////////
// GEE.h
// Author: Wei-Min Chen
// Oct 10, 2004

#ifndef __VCGEE_H__
#define __VCGEE_H__

#include "Pedigree.h"
#include "IntArray.h"
#include "MathMatrix.h"
#include "MathVector.h"
#include "MathCholesky.h"

class GEE{
protected:
   Vector delta;
   Vector delta2;
   Matrix W2, B[5];

   virtual void GetPhi(int f);

   Vector SEvariances, SEvariances_R;

   Matrix DVD;
   Vector SEcoef, SEcoef_R;
   Matrix CovCoef;
   Matrix CovCoef_R;
   Cholesky chol;
   Matrix D;            // D for GEE
   Matrix Omega;        // variace-covariance matrix of trait
   Matrix OmegaInv;
   Matrix *OD;
   Matrix Phi;
   Matrix Delta;
   int parCount;        // number of variance components
   int coefCount;       // number of regression coefficients
   int size;            // size of score

   inline int Index(int u, int v){
      if(u==v) return u;
      else if(u<v) return size + v*(v-1)/2 + u;
      else return size + u*(u-1)/2 + v;
   }
   virtual void RefreshO(int f){}
   virtual void RefreshOD(int f){}
   void Refresh(int f);
   virtual void InitCoef(){}
   virtual void summary(){}
   virtual int constraint();
   virtual int StopRule();
public:
   Matrix covariance;   // variance-covariance matrix of vc
   Matrix covariance_R;

   FILE *polyfp;
   String prefix;
   bool moreFlag;
   bool polyFlag;
   bool saturatedMean;
   Vector meanPerFamily;
   int ValidFamilies;
   int ValidPersons;
   IntArray isNuclear;
   IntArray nuclearP1;
   IntArray nuclearP2;

   Vector *traits;
   Matrix *covariates;

   double deltaScale;
   int AtBorder;
   IntArray borderIndex;
   Matrix newDVD;

   Vector coef;         // regression coefficients
   Vector variances;    // variance components
   IntArray * pheno;

   int LoopCount;
   double Epsilon;
   double loglik;

//   inline Matrix OmegaInverse(Matrix & M);
//   Matrix BlockInverse(Matrix & M, int blockCount, int extra);
//   void Block2Inverse(Matrix & M);

//   Matrix sProductPhi(Matrix & M);
   Pedigree & ped;
   GEE(Pedigree & pedigree);
   virtual void print(){}
   virtual void solve();
   virtual ~GEE();
};


#endif
