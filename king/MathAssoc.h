#ifndef __MATH_ASSOC_H__
#define __MATH_ASSOC_H__

#include "MathVector.h"

// Measures of association based on chi-sq
//

class AssocChi
   {
   public:
      double sum;          // values scored N
      double chisq;        // chi-square value
      double df;           // degrees freedom df
      double prob;         // significance level p
      double lop;          // -log10 of significance
      double cramrv;       // between 0 and 1 - Cramer's V
      double ccc;          // measure of association - depends on I and J

      int isValid;

   AssocChi();

   void Calc(int ** nn, int ni, int nj);
   };

// Measures of Association based on entropy
//

class AssocEntropy
   {
   public:
      double   sum;     // values scored N
      double   h;       // entropy of whole table
      double   hx;      // entropy of the x distribution
      double   hy;      // entropy of the y distribution
      double   hygx;    // entropy of y given x
      double   hxgy;    // entropy of x given y
      double   uygx;    // dependency of x on y
      double   uxgy;    // dependency of y on x
      double   uxy;     // symmetrical dependency of x and y

      int   isValid;

   AssocEntropy();

   void Calc(int **nn, int ni, int nj);
   };

// Spearman's Rank Correlation
void Spearman(Vector & v1, Vector & v2,
              double & rankD, double & zD, double & probD,
              double & spearmanR, double & probR);

#endif


