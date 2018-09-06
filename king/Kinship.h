#ifndef __KINSHIP_H__
#define __KINSHIP_H__

#include "Pedigree.h"
#include "MathMatrix.h"

class Kinship
   {
   public:
      Matrix    allPairs;
      Family *  fam;

      Kinship() : allPairs()
         { fam = NULL; }

      void Setup(Family & f);

      bool isInbred();

      double operator () (Person & p1, Person & p2);

   };

#endif

