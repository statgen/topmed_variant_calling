//////////////////////////////////////////////////////////////////////
// VCLINEAR.h
// Author: Wei-Min Chen
// May 13, 2005

#ifndef _VC_LINEAR_H_
#define _VC_LINEAR_H_

#include "Pedigree.h"
#include "IntArray.h"
#include "MathMatrix.h"
#include "MathVector.h"
#include "MathCholesky.h"
#include "VCGEE.h"

class GEEVC_LINEAR:public GEE{
protected:
   virtual void RefreshOD(int f);
public:
   double H2;
   double seH2;
   double totalVariance;
   double stat;
   double LOD;
   double pvalue;
   IntArray personValid;

   Matrix * varianceComponents;
   Matrix PhiX;
   Matrix PhiM;
   GEEVC_LINEAR(Pedigree & pedigree);
   ~GEEVC_LINEAR();

   void init();
   void InitCoef();
   virtual void summary();
   void print();
   double residual(int p);

   int trait;
   IntArray mCovariate;
};

class POLY:public GEEVC_LINEAR{
protected:
   void RefreshO(int f);
   void RefreshOD(int f);
   int StopRule();
public:
   void InitCoef();
   POLY(Pedigree & pedigree):GEEVC_LINEAR(pedigree){}
   ~POLY(){}
};

class GEEVC_LINKAGE:public GEEVC_LINEAR{
protected:
   void RefreshO(int f);
public:
   double h2;
   Vector *ibd;
   void InitCoef();
   void summary();
   GEEVC_LINKAGE(Pedigree & pedigree):GEEVC_LINEAR(pedigree){ibd=NULL;}
   ~GEEVC_LINKAGE(){if(ibd) delete []ibd;}
};

class GEEVC_ASSOC:public GEEVC_LINKAGE{
public:
   Vector IBS;
   void InitCoef();
   void summary();
   GEEVC_ASSOC(Pedigree & pedigree):GEEVC_LINKAGE(pedigree){/*IBS=new Vector[ped.familyCount];*/}
   ~GEEVC_ASSOC(){/*if(IBS) delete []IBS;*/}
};

#endif
