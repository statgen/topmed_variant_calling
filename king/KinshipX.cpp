#include "KinshipX.h"

void KinshipX::Setup(Family & f)
   {
   allPairs.Dimension(f.count, f.count);

   for (int i = 0; i < f.founders; i++)
      {
      bool isMale = f.ped[f.path[i]].sex == SEX_MALE;
      for (int j = 0; j < f.founders; j++)
         allPairs[i][j] = 0.0;
      allPairs[i][i] = isMale ? 1.0 : 0.5;
      }

   for (int i = f.founders; i < f.count; i++)
      {
      Person * p = &(f.ped[f.path[i]]);
      int k = p->father->traverse;
      int l = p->mother->traverse;

      bool isMale = f.ped[f.path[i]].sex == SEX_MALE;
      allPairs[i][i] = isMale ? 1.0 : 0.5 + allPairs[k][l] * 0.5;

      for (int j = 0; j < i; j++)
         if (!p->isMzTwin(f.ped[f.path[j]]))
            allPairs[i][j] = allPairs[j][i] = isMale ?
               allPairs[l][j] : (allPairs[k][j] + allPairs[l][j]) * 0.5;
         else
            allPairs[j][i] = allPairs[i][j] = allPairs[i][i];
      }

   fam = &f;
   }

double KinshipX::operator() (Person & p1, Person & p2)
   {
   int i = p1.traverse;
   int j = p2.traverse;

   return allPairs[i][j];
   }




