#include "MathAssoc.h"
#include "MathVector.h"
#include "MathConstant.h"
#include "MathStats.h"

// Measures of association based on chi-sq
//

AssocChi::AssocChi()
   {
   sum = chisq = df = prob = lop = cramrv = ccc = 0.0;
   isValid = 0;
   }

void AssocChi::Calc(int ** nn, int ni, int nj)
   {
   int nnj, nni, j, i, minij;
   double expected, temp;
   Vector sumi(ni), sumj(nj);

   sumi.Zero();
   sumj.Zero();
   sum = 0.0;

   nni = ni;
   nnj = nj;

   // Get the row totals
   for (i = 0; i < ni; i++)
      {
      for ( j = 0; j < nj; j++)
         {
         sumi[i] += nn[i][j];
         sum += nn[i][j];
         }
      if ( sumi[i] < FPMIN) --nni;            // eliminate zero rows by reducing number
      }

   // Get the column totals
   for (j = 0; j < nj; j++)
      {
      for ( i = 0; i < ni; i++ )
         sumj[j] += nn[i][j];
      if ( sumj[j] < FPMIN) --nnj;        // eliminated any zero columns
      }

   df = nni * nnj - nni - nnj + 1;           // corrected degrees of freedom

   chisq = 0.0;

   for (i = 0; i < ni; i++)
      {
      for (j = 0; j < nj; j++)
         {
         expected = sumj[j] * sumi[i] / sum;
         temp = nn[i][j] - expected;
         chisq += temp * temp / (expected + TINY);
         }
      }

   prob = df ? gammq ( 0.5 * df, 0.5 * chisq ) : 1.0;
   lop = prob > 1e-100 ? -log10(prob) : 99.999;
   minij = nni < nnj ? nni - 1 : nnj - 1;
   cramrv = minij ? sqrt ( chisq / ( sum * minij ) ) : 0.0;
   ccc = sum ? sqrt ( chisq / ( chisq + sum) ) : 0.0;

   isValid = 1;
   }

// Measures of association based on entropy
//

AssocEntropy::AssocEntropy()
   {
   sum = h= hx = hy = hygx = hxgy = uygx= uxgy = uxy = 0.0;
   isValid = 0;
   }

void AssocEntropy::Calc(int ** nn, int ni, int nj)
   {
   int   i, j;
   double p;
   Vector sumi(ni), sumj(nj);

   sumi.Zero();
   sumj.Zero();
   sum = 0.0;

   // Get the row totals
   for (i = 0; i < ni; i++)
      for ( j = 0; j < nj; j++)
         {
         sumi[i] += nn[i][j];
         sum += nn[i][j];
         }

   // Get the column totals
   for (j = 0; j < nj; j++)
      for ( i = 0; i < ni; i++ )
         sumj[j] += nn[i][j];

   // Entropy of the x distribution
   hx = 0.0;
   for (i = 0; i < ni; i++)
      if (sumi[i] > FPMIN)
         {
         p = sumi[i] / sum;
         hx -= p * log(p);
         }

   // Entropy of the y distribution
   hy = 0.0;
   for (j = 0; j < nj; j++)
      if (sumj[j] > FPMIN)
         {
         p = sumj[j] / sum;
         hx -= p * log(p);
         }

   // Total entropy for the table
   h = 0.0;
   for (i = 0; i < ni; i++)
      for (j = 0; j < nj; j++)
         if (nn[i][j] > 0)
            {
            p = nn[i][j] / sum;
            h -= p * log ( p );
            }

   hygx = h - hx;
   hxgy = h - hy;

   uygx = (hy - hygx) / ( hy + TINY );
   uxgy = (hx - hxgy) / ( hx + TINY );
   uxy = 2.0 * ( hx + hy - h ) / ( hx + hy + TINY );

   isValid = 1;
   }

// Spearman's Rank Correlation
//

double crank(Vector & v)
// Replaces elements of a sorted array by their rank,
// and returns the sum t^3 - t, where t is the number
// of elements in each tie
{
   int j = 0, ji, jt;
   double t, rank, result = 0;

   while (j < v.dim - 1)
      {
      if (v[j] != v[j+1])
         // no tie
         {
         v[j] = j + 1;
         j++;
         }
      else
         {
         // how far does the tie go?
         for (jt = j + 1; jt < v.dim && v[jt] == v[j]; jt++)
            ;
         rank = 0.5 * (j + jt + 1);
         for (ji = j; ji < jt; ji++)
            v[ji] = rank;
         t = jt - j;
         result += t*t*t - t;
         j = jt;
         }
      }
   if (j == v.dim - 1) v[j] = j;  // rank for last element

   return result;
   }

void Spearman(Vector & v1, Vector & v2,
              double & rankD, double & zD, double & probD,
              double & spearmanR, double & probSR)
   {
   double varD, sg, sf, fac, en3n, en, df, aveD, t;
   Vector wksp1, wksp2;

   wksp1.Copy(v1);
   wksp2.Copy(v2);

   wksp1.Sort(wksp2);
   sf = crank(wksp1);
   wksp2.Sort(wksp1);
   sg = crank(wksp2);

   rankD = 0;
   for (int j = 0; j < v1.dim; j++)
      // sum the square difference of ranks
      {
      double temp = wksp1[j] - wksp2[j];
      rankD += temp  * temp;
      }

   en = v1.dim;
   en3n = en*en*en - en;
   aveD = en3n / 6.0 - (sf + sg) / 12.0;
   fac = (1.0 - sf/en3n) * (1.0 - sg/en3n);
   varD = ((en - 1.0)*en*en*(en + 1.0)*(en + 1.0)/36.0)*fac;
   zD = (rankD - aveD) / sqrt(varD);
   probD = erfcc(fabs(zD)/1.4142136);

   spearmanR = (1.0 - (6.0/en3n)*(rankD+(sf+sg)/12.0))/sqrt(fac);
   fac = (spearmanR + 1.0) * (1.0 - spearmanR);
   if (fac)
      {
      t = (spearmanR) * sqrt((en - 2.0)/fac);
      df = en - 2.0;
      probSR = betai(0.5 * df, 0.5, df/(df + t*t));
      }
   else
      probSR = 0.0;
   }


