#ifndef __PEELERNODES_H__
#define __PEELERNODES_H__

#include "Pedigree.h"

class PeelerNode
   {
   public:
      ~PeelerNode();

   protected:
      static Vector scratch;
   };

class MatingNode;
class PersonNode;

class PersonNode : public PeelerNode
   {
   public:
      Person * person;

      IntArray states;
      Vector   probabilities;

      void Clear();

      void PeelDescendants(MatingNode * mating, double (*trans) (int, int, int));
      void PeelAncestors(MatingNode * mating, double (*trans) (int, int, int));

      double Probability()   { return probabilities.Sum(); }
   };

class MatingNode : public PeelerNode
   {
   public:
      IntArray mstates, pstates;
      Vector   probabilities;

      PersonNode * father;
      PersonNode * mother;

      void Initialize(PersonNode * father, PersonNode * mother);

      void PeelFather();
      void PeelMother();
      void PeelOffspring(PersonNode * child, double (*trans) (int, int, int));

      double Probability()   { return probabilities.Sum(); }
   };

#endif

