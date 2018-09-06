#include "MapFunction.h"
#include "MathConstant.h"

#include <math.h>

double DistanceToRecombination(double distance)
   {
   return (1.0 - exp(-2.0 * distance)) * 0.5;
   }

double RecombinationToDistance(double recombination)
   {
   return (log(max(1.0 - 2 * recombination, 1e-7)) * -0.5);
   }

double KosambiDistanceToRecombination(double distance)
   {
   double e_to_4x = exp(4.0 * distance);

   return (0.5 * (e_to_4x - 1.0) / (e_to_4x + 1.0));
   }

double RecombinationToKosambiDistance(double theta)
   {
   return 0.25 * log((1.0 + 2*theta) / max(1.0 - 2.0*theta, 1e-7));
   }
