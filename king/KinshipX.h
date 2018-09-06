#ifndef __KINSHIPX_H__
#define __KINSHIPX_H__

#include "Pedigree.h"
#include "MathMatrix.h"

class KinshipX
   {
   public:
      Matrix    allPairs;
      Family *  fam;

      KinshipX() : allPairs()
         { fam = NULL; }

      void Setup(Family & f);

      double operator () (Person & p1, Person & p2);

   };

#endif

