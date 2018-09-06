////////////////////////////////////////////////////////////////////// 
// merlin/MerlinSort.cpp 
// (c) 2000-2007 Goncalo Abecasis
// 
// This file is distributed as part of the MERLIN source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile MERLIN.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Tuesday December 18, 2007
// 
 
#include "MerlinSort.h"
#include "MathConstant.h"

void SortFamilies(Pedigree & ped)
   {
   IntArray scores;

   for (int f = 0; f < ped.familyCount; f++)
      {
      int & founders = ped.families[f]->founders;
      int & count = ped.families[f]->count;
      int * & path = ped.families[f]->path;

      // Define a score for each individual based on...
      //     * No. of genotyped markers (primarily)
      //     * Affection status (tie-breaker)
      //
      scores.Dimension(count);

      for (int i = founders; i < count; i++)
         scores[i] = ped[path[i]].ngeno * 2 + (ped.affectionCount == 0 ? 0 :
                     ped[path[i]].affections[0] == 2);

      // Optimize path so descendants with higher scores appear first
      for (int i = founders + 1; i < count; i++)
         {
         // Non-founders must always follow founders
         int new_pos = founders;

         // In addition they must follow their father and any of his MZ twins
         Person & father = *ped[path[i]].father;

         if (father.traverse >= new_pos)
            new_pos = father.traverse + 1;

         // If father is an MZ twin, can't move this individual above his co-twins
         if (father.zygosity & 1)
            for (int j = 0; j < father.sibCount; j++)
               if (father.sibs[j]->zygosity == father.zygosity &&
                   father.sibs[j]->traverse >= new_pos)
                  new_pos = father.sibs[j]->traverse + 1;

         // In addition they must follow their mother and any of his MZ twins
         Person & mother = *ped[path[i]].mother;

         if (mother.traverse >= new_pos)
            new_pos = mother.traverse + 1;

         // If mother is an MZ twin, can't move this individual above her co-twins
         if (mother.zygosity & 1)
            for (int j = 0; j < mother.sibCount; j++)
               if (mother.zygosity == mother.sibs[j]->zygosity &&
                   mother.sibs[j]->traverse >= new_pos)
                  new_pos = mother.sibs[j]->traverse + 1;

         // Subject to these constraints, place individual above any others
         // with lower informativeness scores
         while (scores[new_pos] > scores[i] && new_pos < i)
             new_pos++;

         if (new_pos != i)
            {
            int person_to_move = path[i];
            int saved_score = scores[i];

            for (int move = i; move > new_pos; move--)
               {
               scores[move] = scores[move-1];
               path[move] = path[move-1];
               ped[path[move]].traverse++;
               }

            ped[person_to_move].traverse = new_pos;
            path[new_pos] = person_to_move;
            scores[new_pos] = saved_score;
            }
         }
      }
   }

 
