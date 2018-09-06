#include "MathMiser.h"
#include "MathConstant.h"
#include "Random.h"

#define  MINNODE    60       // minimum points for bisection
#define  MINLEAF    15       // minimum points in a region
#define  R_D_MIN    15       // minimum points used for research
#define  R_D_MAX    9999     // maximum number of points used for research
#define  R_D_FRAC   0.1      // proportion of total points used for research

double MathMiser::Integrate(Matrix & volume)
   {
   double average = 0.0, variance = 0.0;
   unsigned long seed = 0;

   sobol.Init(volume.cols);

   midpoint.Dimension(volume.cols);
   point.Dimension(volume.cols);
   minl.Dimension(volume.cols);
   minr.Dimension(volume.cols);
   maxl.Dimension(volume.cols);
   maxr.Dimension(volume.cols);

   stack[0].points = ncall;
   stack[0].weight = 1.0;
   stack[0].region = volume;

   for (int ptr = 1; ptr--; )
      {
      MiserStack & top = stack[ptr];

      // If very few points are left do straight MC over the region
      if (top.points < MINNODE)
         {
         double sum = 0.0, sum2 = 0.0;
         for (int n = 0; n < top.points; n++)
            {
            RandomPoint(top.region, point);
            double fval = func(point);
            sum += fval;
            sum2 += fval * fval;
            }
         average += top.weight * sum / top.points;
         variance += (top.weight * top.weight) *
            max(0.0, (sum2 - sum*sum/top.points) / (top.points * top.points));
         continue;
         }

      // Caculate midpoint for the current region
      for (int j = 0; j < midpoint.dim; j++)
         midpoint[j] = 0.5 * (top.region[0][j] + top.region[1][j]);

      // Explore the variance of the function in this region...
      int quota = min(max((int) (top.points * R_D_FRAC), R_D_MIN), R_D_MAX);
      top.points -= quota;

      minr.Set(FPMAX);
      minl.Set(FPMAX);
      maxr.Set(-FPMAX);
      maxl.Set(-FPMAX);

      for (int n = 0; n < quota; n++)
         {
         RandomPoint(top.region, point);
         double fval = func(point);

         for (int j = 0; j < point.dim; j++)
            if (point[j] <= midpoint[j])
               {
               minl[j] = min(minl[j], fval);
               maxl[j] = max(maxl[j], fval);
               }
            else
               {
               minr[j] = min(minr[j], fval);
               maxr[j] = max(maxr[j], fval);
               }
         }

      // Choose the region giving biggest reduction in the variance
      double bestvar = FPMAX, varl;
      int    bestdir = -1;

      for (int j = 0; j < point.dim; j++)
         {
         // have we got two points on each half of the region?
         if (maxl[j] > minl[j] && maxr[j] > minr[j])
            {
            double sigmal = max (TINY, pow(maxl[j] - minl[j], 2.0 / 3.0));
            double sigmar = max (TINY, pow(maxr[j] - minr[j], 2.0 / 3.0));
            double sigma = sigmal + sigmar;
            if (sigma < bestvar)
               {
               bestvar = sigma;
               bestdir = j;
               varl = sigmal;
               }
            }
         }

      if (bestdir == -1) { varl = 1.0; bestdir = RAND(seed) % point.dim; }

      // Divide the remaining points among the two regions
      int pointsl = MINLEAF + int ((top.points - 2*MINLEAF) * varl / bestvar);
      int pointsr = top.points - pointsl;

      // Integrate the two sub regions, starting with the smallest one
      stack[ptr + 1].region = stack[ptr].region;
      stack[ptr + 1].weight = stack[ptr].weight *= 0.5;

      if (pointsl > pointsr)
         {
         stack[ptr].points = pointsl;
         stack[ptr++].region[1][bestdir] = midpoint[bestdir];
         stack[ptr].points = pointsr;
         stack[ptr++].region[0][bestdir] = midpoint[bestdir];
         }
      else
         {
         stack[ptr].points = pointsr;
         stack[ptr++].region[0][bestdir] = midpoint[bestdir];
         stack[ptr].points = pointsl;
         stack[ptr++].region[1][bestdir] = midpoint[bestdir];
         }
      }

   tgral = 1.0;
   for (int j = 0; j < point.dim; j++)
      tgral *= (volume[1][j] - volume[0][j]);

   stdev  = tgral * sqrt(variance);
   tgral *= average;

   return tgral;
   }

Random tr;

void MathMiser::RandomPoint(Matrix & region, Vector & point)
   {
   sobol.Next(point);

   for (int i = 0; i < region.cols; i++)
      point[i] = region[0][i] + ( region[1][i] - region[0][i] ) * point[i];
   }


