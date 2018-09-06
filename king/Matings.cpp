#include "Matings.h"

#define MATING_HASH_PADDING      2     /* Number of empty slots per mating,
                                          which are used to speed up searching
                                          at the cost of increased memory use */

void Matings::ListMatings(Family * family)
   {
   founders = family->founders;
   matingCount = 0;

   InitializeHash(family->count);

   for (int i = family->founders; i < family->nonFounders; i++)
      {
      Person & p = family->ped[family->path[i]];

      matingMap[i - family->founders] = LookupMating(p.father->serial, p.mother->serial);
      }
   }

void Matings::InitializeHash(int size)
   {
   size *= MATING_HASH_PADDING;

   hash.Dimension(size);
   hash.Set(-1);

   hashId.Dimension(size);
   hashId.Set(-1);
   }

int Matings::LookupMating(int father, int mother)
   {
   int id = father * hash.Length() + mother;
   int h  = father * MATING_HASH_PADDING;

   while (true)
      {
      if (hash[h] == -1)
         {
         hashId[h] = id;
         return hash[h] = matingCount++;
         }

      if (hashId[h] == id)
         return hash[h];

      h++;

      if (h == hash.Length()) h = 0;
      }
   }

int Matings::LookupMating(Person & p)
   {
   return LookupMating(p.serial - founders);
   }

int Matings::LookupMating(int serial)
   {
   return matingMap[serial - founders];
   }
