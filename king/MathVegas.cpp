#include "MathVegas.h"
#include "MathConstant.h"

#include <math.h>

Random Vegas::rand;

Vegas::Vegas() : d(NDMX, MXDIM), di(NDMX, MXDIM), dt(MXDIM),
   dx(MXDIM), r(NDMX), x(MXDIM), xi(MXDIM, NDMX), xin(NDMX)
   {
   itmx = 5;
   ncall = 1000;
   tgral = sd = chi2a;

   d.Zero(); di.Zero(); dt.Zero(); dx.Zero();
   r.Zero(); x.Zero(); xi.Zero(); xin.Zero();

   ia = new int [MXDIM];
   kg = new int [MXDIM];

   nd = -1;

   vfunc = NULL;
   }

Vegas::~Vegas()
   {
   delete [] ia;
   delete [] kg;
   }

void Vegas::Init(Matrix & Region, int level)
{
   x.Dimension(Region.cols);

   switch (level)
      {
      case 0 :
         {
         mds = ndo = 1;       // disable stratified sampling with mds = 1
         for (int j = 1; j <= Region.cols; j++) xi[j][1] = 1.0;
         }
      case 1 :
         si = swgt = schi = 0.0;
      case 2 :
         {
         nd = NDMX;
         ng = 1;
         if (mds)             // setup stratification
            {
            ng = (int) pow(ncall / 2.0 + 0.25, 1.0 / Region.cols);
            mds = 1;
            if (2 * ng - NDMX >= 0)
               {
               mds = -1;
               npg = ng / NDMX + 1;
               nd = ng / npg;
               ng = npg * nd;
               }
            }
         int k = 1;
         for  (int i = 1; i <= Region.cols; i++) k *= ng;
         npg = max( (int) ncall/k, 2 );
         calls = double (npg) * double (k);
         dxg = 1.0 / ng;
         dv2g = 1.0;
         for (int i = 1; i <= Region.cols; i++) dv2g *= dxg;
         dv2g *= calls;
         dv2g *= dv2g / npg / npg / (npg - 1.0);
         xnd = nd;
         dxg *= xnd;
         xjac = 1.0 / calls;
         for (int j = 1; j <= Region.cols; j++)
            {
            dx[j] = Region[2][j] - Region[1][j];
            xjac *= dx[j];
            }
         if (nd != ndo)
            {
            for (int i = 1; i <= max(nd, ndo); i++)
               r[i] = 1.0;
            for (int j = 1; j <= Region.cols; j++)
               Rebin(ndo / xnd, xi[j]);
            ndo = nd;
            }
         }
      default :
         ;
      }
   }

double Vegas::Integrate(Matrix & Region)
   {
   if (nd == -1) Init(Region);

   for (int it = 1; it <= itmx; it++)
      {
      int    k;
      double ti = 0.0, tsi = 0.0;
      for (int j = 1; j <= Region.cols; j++)
         {
         kg[j] = 1;
         for (int i = 1; i <= nd; i++)
            d[i][j] = di[i][j] = 0.0;
         }
      while (true)
         {
         double fb = 0.0, f2b = 0.0;
         for (k = 1; k <= npg; k++)
            {
            wgt = xjac;
            for (int j = 1; j <= Region.cols; j++)
               {
               xn = (kg[j] - rand.Next())*dxg + 1.0;
               ia[j] = min( (int) xn, NDMX );
               if (ia[j] <= 1)
                  {
                  ia[j] = 1;
                  xo = xi[j][1];
                  rc = (xn - 1) * xo;
                  }
               else
                  {
                  xo = xi[j][ia[j]] - xi[j][ia[j] - 1];
                  rc = xi[j][ia[j] - 1] + (xn - ia[j]) * xo;
                  }
               x[j] = Region[1][j] + rc * dx[j];
               wgt *= xo * xnd;
               }
            double f = wgt * func(x);
            double f2 = f * f;
            fb += f;
            f2b += f2;
            for (int j = 1; j <= Region.cols; j++)
               {
               di[ia[j]][j] += f;
               if (mds >= 0) d[ia[j]][j] += f2;
               }
            }
         f2b = sqrt(f2b * npg);
         f2b = (f2b - fb) * (f2b + fb);
         if (f2b <= 0.0) f2b = TINY;
         ti += fb;
         tsi += f2b;
         if (mds < 0)      // Use stratified sampling
            for (int j = 1; j <= Region.cols; j++)
               d[ia[j]][j] += f2b;
         for (k = Region.cols; k >= 1; k--)
            {
            kg[k] %= ng;
            if (++kg[k] != 1) break;
            }
         if (k < 1) break;
         }
      // Final results for this iteration
      tsi *= dv2g;
      wgt = 1.0 / tsi;
      si += wgt * ti;
      schi += wgt * ti * ti;
      swgt += wgt;
      tgral = si / swgt;
      chi2a = (schi - si * (tgral)) / ( it - 0.9999 );
      if (chi2a < 0.0) chi2a = 0.0;
      sd = sqrt(1/swgt);

      // If the function looks very close to constant
      // skip the refinement step...
      if (tsi < 1e-8) continue;

      // Grid refinement is damped, to avoid rapid
      // destabilizing changes, and compressed in
      // range by ALPH
      for (int j = 1; j <= Region.cols; j++)
         {
         xo = d[1][j];
         xn = d[2][j];
         d[1][j] = (xo + xn) / 2.0;
         dt[j] = d[1][j];
         for (int i = 2; i < nd; i++)
            {
            rc = xo + xn;
            xo = xn;
            xn = d [i+1][j];
            d[i][j] = (rc + xn) / 3.0;
            dt[j] += d[i][j];
            }
         d[nd][j] = (xo + xn) / 2.0;
         dt[j] += d[nd][j];
         }
      for (int j = 1; j <= Region.cols; j++)
         {
         rc = 0.0;
         for (int i = 1; i <= nd; i++)
            {
            if (d[i][j] < TINY) d[i][j] = TINY;
            r[i] = pow((1.0 - d[i][j]/dt[j]) /
                   (log(dt[j]) - log(d[i][j])), ALPH);
            rc += r[i];
            }
         Rebin(rc / xnd, xi[j]);
         }
      }
   return tgral;
   }

void Vegas::Rebin(double rc, Vector & xi)
   {
   double dr = 0, xo = 0;

   for (int k = 0, i = 1; i < nd; i++)
      {
      while (rc > dr)
         dr += r[++k];
      if (k > 1) xo = xi[k - 1];
      xn = xi[k];
      dr -= rc;
      xin[i] = xn - (xn - xo) * dr / r[k];
      }
   for (int i = 1; i < nd; i++)
      xi[i] = xin[i];
   xi[nd] = 1.0;
   }

