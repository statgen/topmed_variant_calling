#ifndef __MATINGS_H__
#define __MATINGS_H__

#include "Pedigree.h"

class Matings
   {
   public:
      // Number of distinct matings in the pedigree
      int      matingCount;
      int      founders;

      // Map linking each non-founder to a mating
      IntArray matingMap;

      // Index all the matings in a family
      void ListMatings(Family * family);

      // Lookup the mating index for a specific offspring
      int  LookupMating(int serial);
      int  LookupMating(Person & p);

   private:
      void InitializeHash(int size);
      int  LookupMating(int father, int mother);

      IntArray hash;
      IntArray hashId;
   };



#endif

