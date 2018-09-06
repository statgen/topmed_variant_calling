#include "PeelerNodes.h"

Vector PeelerNode::scratch;

PeelerNode::~PeelerNode()
   {
   }

void PersonNode::PeelAncestors(MatingNode * mating, double (* trans)(int, int, int))
   {
   for (int j = 0; j < states.Length(); j++)
      {
      double p = 0.0;

      for (int i = 0; i < mating->mstates.Length(); i++)
         p += trans(mating->father->states[mating->pstates[i]],
                    mating->mother->states[mating->mstates[i]],
                    states[j]) * mating->probabilities[i];

      probabilities[j] *= p;
      }
   }

void PersonNode::PeelDescendants(MatingNode * mating, double (* trans)(int, int, int))
   {
   scratch.Dimension(states.Length());
   scratch.Zero();

   IntArray & index = person->sex == SEX_MALE ? mating->pstates : mating->mstates;

   for (int i = 0; i < index.Length(); i++)
      scratch[index[i]] += mating->probabilities[i];

   for (int j = 0; j < states.Length(); j++)
      probabilities[j] *= scratch[j];
   }

void PersonNode::Clear()
   {
   states.Clear();
   probabilities.Clear();
   }

void MatingNode::PeelFather()
   {
   for (int j = 0; j < pstates.Length(); j++)
      probabilities[j] *= father->probabilities[pstates[j]];
   }

void MatingNode::PeelMother()
   {
   for (int j = 0; j < mstates.Length(); j++)
      probabilities[j] *= mother->probabilities[mstates[j]];
   }

void MatingNode::PeelOffspring(PersonNode * child, double (*trans) (int, int, int))
   {
   for (int i = pstates.Length() - 1; i >= 0; i--)
      {
      double p = 0.0;

      for (int j = 0; j < child->states.Length(); j++)
         p += trans(pstates[i], mstates[i], child->states[j]) *
              child->probabilities[j];

      if (p > 0.0)
         probabilities[i] *= p;
      else
         mstates.Delete(i),
         pstates.Delete(i),
         probabilities.Delete(i);
      }
   }

void MatingNode::Initialize(PersonNode * father, PersonNode * mother)
   {
   mstates.Dimension(father->states.Length() * mother->states.Length());
   pstates.Dimension(mstates.Length());

   probabilities.Dimension(mstates.Length());
   probabilities.Set(1.0);

   for (int i = 0; i < father->states.Length(); i++)
      for (int j = 0; j < mother->states.Length(); j++)
         pstates[i] = father->states[i],
         mstates[j] = mother->states[j];
   }


