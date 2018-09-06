////////////////////////////////////////////////////////////////////// 
// libsrc/PedigreeAlleles.h 
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
 
#ifndef __PEDALLELES_H__
#define __PEDALLELES_H__

#include "LongInt.h"

class Alleles{
   public:
      char geno;
      Alleles(){ geno = 0; }

      char operator [] (int i)
      { return (i == 1) ? (geno&15) : (geno>>4);}

      void AssignGenotype(int G1, int G2)
      {geno = char(G1 + (G2<<4));}

      void AssignGenotype(int G)
      {geno = char(G);}

      // is the genotype fully defined?
      bool isKnown(){ return geno != 0; }
      bool isHeterozygous()
         { return isKnown() && ((geno&15) != (geno>>4)); }
      bool isHomozygous()
         { return isKnown() && ((geno&15) == (geno>>4)); }
      bool hasAllele(int a)
         { return ((geno&15) == a) || ((geno>>4) == a); }

      // in a bi-allelic system (a, NOT a)
      bool isHeterozygousFor(int a){ return isHeterozygous() && hasAllele(a); }
      bool isHomozygousFor(int a){ return !(isHeterozygousFor(a)); }

      // how may alleles a in this genotype?
      int countAlleles(int a)
         { return (((geno&15) == a) ? 1 : 0) + (((geno>>4) == a) ? 1 : 0); }

      // what is the other allele, assuming genotype is (a, X)
      int otherAllele(int a)
         { return (((geno&15) == a) ? (geno>>4) : (geno&15)); }

      // are two unordered genotypes identical?
      int identicalTo(Alleles & al)
         { return (al.geno == geno) ||
                  ((al[2]>>4)+(al[1]<<4) == geno);}

      // how many alleles are identical by state
      int countIBS(Alleles & al)
         { return  ((geno&15) == al[1]) ?
                    (((geno>>4) == al[2]) ? 2 : 1) :
                  (  ((geno&15) == al[2]) ?
                    (((geno>>4) == al[1]) ? 2 : 1) :
                   ((((geno>>4) == al[1]) || ((geno>>4) == al[2])) ? 1 : 0));
         }

      int operator == (Alleles & rhs) { return identicalTo(rhs); }
      int operator != (Alleles & rhs) { return !identicalTo(rhs); }

      char Hi()
         { return (geno&15) > (geno>>4) ? (geno&15) : (geno>>4); }
      char Lo()
         { return (geno&15) > (geno>>4) ? (geno>>4) : (geno&15); }

      int SequenceCoded()
         { return isKnown() ? Hi() * (Hi() - 1) / 2 + Lo() : 0; }

      longint BinaryCoded()
         {
         if (isKnown())
            {
            longint allele1(1);
            longint allele2(1);

            allele1 <<= (geno&15) - 1;
            allele2 <<= (geno>>4) - 1;

            return allele1 | allele2;
            }
         else
            return NOTZERO;
         }

      void Intersect(Alleles & gen)
         {
         char a1 = Lo(), a2 = Hi();
         char b1 = gen.Lo(), b2 = gen.Hi();

         if (a1 == b1 && a2 == b2)
            return;
         if (a1 == b1 || a1 == b2)
            geno = (a1<<4) + a1;
         else if (a2 == b1 || a2 == b2)
            geno = (a2<<4) + a2;
         else
            geno = 0;
         }
                
      void Intersect(char allele)
         {
         if ((geno&15) != allele && (geno>>4) != allele)
            geno = 0;
         else
            geno = (allele << 4) + allele;
         }

      bool AddAllele(char allele)
         {
         if ((geno&15) == allele || (geno>>4) == allele)
            return true;

         if ((geno&15) != 0 && (geno>>4) != 0)
            return false;

         if ((geno&15) == 0) geno |= allele; else geno |= (allele>>4);
         return true;
         }

      void Wipe() {geno=0;}
   };

#endif


