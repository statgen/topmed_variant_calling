#ifndef __IBD_H__
#define __IBD_H__

#include "Pedigree.h"

#include <stdio.h>

class IBD
   {
   public:
      double p0, p1, p2;

   IBD()
      { p0 = p1 = p2 = 0.0; }
   IBD(double zero, double one, double two)
      { p0 = zero; p1 = one; p2 = two; }

   void defaultSib()
      { p0 = p2 = 0.25; p1 = 0.5; }

   void defaultSelf()
      { p0 = p1 = 0; p2 = 1.0; }

   void defaultUnrelated()
      { p0 = 1.0; p1 = p2 = 0.0; }

   void defaultFounderOffspring()
      { p0 = p2 = 0.0; p1 = 1.0; }

   double expected()
      { return 0.5 * p1 + p2; }

   bool isValid()
      { return (p0 + p1 + p2) == 1.0; }

   IBD & operator = (IBD & rhs)
      { p0 = rhs.p0;
        p1 = rhs.p1;
        p2 = rhs.p2;
        return (*this); }

   bool operator == (IBD & rhs);
   bool operator != (IBD & rhs);

   IBD * SimpleIBD(int marker, Person & p1, Person & p2);
   };

struct IBDKey
   {
   int serialLo;
   int serialHi;

   void SelectPair(Person & p1, Person & p2);
   };

struct IBDPair
   {
   int   serialLo;
   int   serialHi;
   IBD   ibd;

   void Assign(IBDKey & key, IBD & i)
      {
      serialLo = key.serialLo;
      serialHi = key.serialHi;
      ibd = i;
      }
   };

class IBDList
   {
   public:
      IBDPair * list;
      int       size, count;

      IBDList();
      ~IBDList();

      IBD * Lookup(Person & p1, Person & p2);
      void  Append(Person & p1, Person & p2, IBD & ibd);
      void  Sort(Pedigree & ped);
      bool  IsRangeEmpty(int low, int high);

   private:
      void Grow();
   };

class IBDTable
   {
   public:
      IBDList * markers;

      IBDTable();
      ~IBDTable();

      void Load(Pedigree & ped, FILE * f);
      void Load(Pedigree & ped, const char * filename);
      void Load(Pedigree & ped, const char * filename, Vector & LocusMap);

      IBD * Lookup(int marker, Person & p1, Person & p2);

      bool  HaveFamily(int marker, Family * f);

      bool isEmpty()
         { return markers == NULL; }
   };


#endif
