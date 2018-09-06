#include "IBD.h"
#include "StringMap.h"
#include "Sort.h"
#include "Error.h"
#include "math.h"

#include <string.h>
#include <ctype.h>
#include <stdlib.h>

int CompareRangeToIBDPair(IBDKey * key, IBDPair * probe);
int CompareKeyToIBDPair(IBDKey * key, IBDPair * probe);
int CompareIBDPair(IBDPair * test, IBDPair * probe);

bool IBD::operator == (IBD & rhs)
   {
   return (p0 == rhs.p0 && p1 == rhs.p1 && p2 == rhs.p2);
   }

bool IBD::operator != (IBD & rhs)
   {
   return (p0 != rhs.p0 || p1 != rhs.p1 || p2 != rhs.p2);
   }

IBD * IBD::SimpleIBD(int m, Person & p, Person & q)
   {
   p0 = p1 = p2 = 0.0;

   // if comparing to self
   if (q.serial == p.serial || p.isMzTwin(q)) {
      p2 = 1.0;
      return this; }

   Person * person1 = &p;
   Person * person2 = &q;
   Person::Order(person1, person2);

   if (person1->isFounder())
      {
      if (person2->isFounder()) {
         p0 = 1.0;
         return this; }

      if (person2->father == person1 || person2->mother == person1) {
         p1 = 1.0;
         return this; }

      error("Your data set includes extended pedigrees (eg, family %s)\n"
            "Please run prelude, simwalk2 and finale to estimate IBD",
            (const char *) p.famid);
      }

   if (!p.isSib(q))
      error("Your data set includes extended pedigrees (eg, family %s)\n"
            "Please run prelude, simwalk2 and finale to estimate IBD",
            (const char *) p.famid);

   // we will only handle cases where both parents are typed
   if ( !p.father->isGenotyped(m) ||
        !p.mother->isGenotyped(m))
      error("Your data set includes untyped parents (eg, family %s)\n"
            "Please run prelude, simwalk2 and finale to estimate IBD",
            (const char *) p.famid);

   if (!p.isGenotyped(m) || !q.isGenotyped(m))
      {
      defaultSib();
      return this;
      }

   // if one parent is homozygous we can easily evaluate
   // ibd for the other meiosis
   int uninformative = 0;

   if ( p.father->markers[m].isHomozygous() )
        uninformative = p.father->markers[m][0];

   if ( p.mother->markers[m].isHomozygous() )
      {
      // no information if both are homozygous
      if (uninformative) {
         defaultSib();
         return this; }

      uninformative = p.mother->markers[m][0];
      }

   if ( uninformative )
      {
      if ( p.markers[m].otherAllele(uninformative) ==
           q.markers[m].otherAllele(uninformative) )
         p1 = p2 = 0.5;
      else
         p0 = p1 = 0.5;
      return this;
      }

   // If both parents are heterozygous
   // identity by state is a good starting point
   // (with two exceptions!)

   int ibs = p.markers[m].countIBS( q.markers[m] );

   switch (ibs) {
      case 0:
         p0 = 1;
         break;
      case 1:
         if (p.markers[m].isHomozygous() ||
             q.markers[m].isHomozygous())
            p1 = 1;
         else
            {
            int shared = q.markers[m].hasAllele(p.markers[m][0]) ?
                         p.markers[m][0] : p.markers[m][1];
            if (p.father->markers[m].hasAllele(shared) &&
                p.mother->markers[m].hasAllele(shared))
               p0 = 1;
            else
               p1 = 1;
            }
         break;
      case 2:
         if (p.markers[m].isHomozygous() ||
             !p.father->markers[m].identicalTo(p.mother->markers[m]))
            p2 = 1;
         else
            p0 = p2 = .5;
         break;
      }

   return this;
   }

int CompareRangeToIBDPair(IBDKey * key, IBDPair * probe)
   {
   return  (key->serialLo <= probe->serialLo) ?
          ((key->serialHi >= probe->serialHi) ? 0 : -1) : 1;
   }

int CompareKeyToIBDPair(IBDKey * key, IBDPair * probe)
   {
   int result = key->serialLo - probe->serialLo;

   if (result) return result;

   return key->serialHi - probe->serialHi;
   }

int CompareIBDPair(IBDPair * test, IBDPair * probe)
   {
   int result = test->serialLo - probe->serialLo;

   if (result) return result;

   return test->serialHi - probe->serialHi;
   }

void IBDKey::SelectPair(Person & p1, Person & p2)
   {
   serialLo = p1.serial;
   serialHi = p2.serial;

   if (serialLo > serialHi)
      {
      int temp = serialHi;
      serialHi = serialLo;
      serialLo = temp;
      }
   }

IBDList::IBDList()
   {
   size = 1024;
   count = 0;
   list = new IBDPair [size];
   }

IBDList::~IBDList()
   {
   delete [] list;
   }

void IBDList::Grow()
   {
   int newSize = size * 2;
   IBDPair * newList = new IBDPair [newSize];

   memcpy(newList, list, size * sizeof(IBDPair));

   delete [] list;

   list = newList;
   size = newSize;
   }

IBD * IBDList::Lookup(Person & p1, Person & p2)
   {
   static IBD selfIBD(0.0, 0.0, 1.0);
   static IBD founderIBD(1.0, 0.0, 0.0);
   static IBD founderOffspringIBD(0.0, 1.0, 0.0);

   if (&p1 == &p2 || p1.isMzTwin(p2))
      return &selfIBD;

   IBDKey key;

   key.SelectPair(p1, p2);

   IBDPair * result = (IBDPair *) BinarySearch
                      (&key, list, count, sizeof(IBDPair),
                       COMPAREFUNC CompareKeyToIBDPair);

   if (result == NULL)
      {
      if (p1.isFounder() && p2.isFounder())
         return & founderIBD;
      else if (p1.isFounder() && (p2.mother == &p1 || p2.father == &p1) ||
               p2.isFounder() && (p1.mother == &p2 || p1.father == &p2))
         return & founderOffspringIBD;
      else return NULL;
/*         error("IBDList.Lookup: Couldn't find IBD status for pair:\n"
               "     Family, Person : %s, %s\n"
               "     Family, Person : %s, %s\n\n"
               "This error can result from mismatched pedigree and ibd files\n",
               (const char *) p1.famid, (const char *) p1.pid,
               (const char *) p2.famid, (const char *) p2.pid);
*/      }

   return &(result->ibd);
   }

bool IBDList::IsRangeEmpty(int lo, int hi)
   {
   IBDKey key;

   key.serialLo = lo;
   key.serialHi = hi;

   return BinarySearch(&key, list, count, sizeof(IBDPair),
          COMPAREFUNC CompareRangeToIBDPair) == NULL;
   }

void IBDList::Append(Person & p1, Person & p2, IBD & ibd)
   {
   if (count == size)
      Grow();

   IBDKey key;
   key.SelectPair(p1, p2);

   list[count++].Assign(key, ibd);
   }

void IBDList::Sort(Pedigree & ped)
   {
   QuickSort(list, count, sizeof(IBDPair), COMPAREFUNC CompareIBDPair);

   for (int i = 1; i < count; i++)
      if (list[i].serialLo == list[i - 1].serialLo &&
          list[i].serialHi == list[i - 1].serialHi &&
          list[i].ibd != list[i - 1].ibd)
         {
         Person & p1 = ped[list[i].serialLo];
         Person & p2 = ped[list[i].serialHi];

         error("Conflicting IBD values listed for pair:\n"
               "   Family, Person : %s, %s\n"
               "   Family, Person : %s, %s\n",
               (const char *) p1.famid, (const char *) p1.pid,
               (const char *) p2.famid, (const char *) p2.pid);
         }
   }

IBDTable::IBDTable()
   {
   markers = NULL;
   }

IBDTable::~IBDTable()
   {
   if (markers != NULL)
      delete [] markers;
   }

IBD * IBDTable::Lookup(int marker, Person & p1, Person & p2)
   {
   if (!isEmpty())
      return markers[marker].Lookup(p1, p2);
   else
      {
      static IBD ibd;
      return ibd.SimpleIBD(marker, p1, p2);
      }
   }

bool IBDTable::HaveFamily(int marker, Family * f)
   {
   if ( isEmpty() ||                                      // use internal engine
        f->count <= f->mzTwins + 3 && f->founders == 2 || // a trivial family
        f->count == f->founders)                          // cases and controls
      return true;

   else
      return !(markers[marker].IsRangeEmpty(f->first, f->last));
   }


void IBDTable::Load(Pedigree & ped, const char * filename)
   {
   FILE * f = fopen(filename, "rb");
   if (f == NULL) return;
   Load(ped, f);
   fclose(f);
   }

void IBDTable::Load(Pedigree & ped, FILE * f)
{
   if (markers != NULL)
      delete [] markers;

   markers = new IBDList [ped.markerCount];

   char buffer[BUFSIZE];
   char famid[BUFSIZE], pid1[BUFSIZE], pid2[BUFSIZE], markerid[BUFSIZE];
   Person * p1;
   Person * p2;
   int      marker;
   IBD      ibd;
   bool isEmpty=true;

   while (fgets(buffer, BUFSIZE - 1, f) != NULL) {
      int check = sscanf(buffer, " %s %s %s %s %lf %lf %lf",
                        famid, pid1, pid2, markerid, &ibd.p0,
                        &ibd.p1, &ibd.p2);
      if (check < 7) continue;

      marker = ped.LookupMarker(markerid);

      if (marker < 0) continue;

      p1 = ped.FindPerson(famid, pid1);
      p2 = ped.FindPerson(famid, pid2);
      if (p1 == NULL || p2 == NULL) continue;
      isEmpty = false;
      markers[marker].Append(*p1, *p2, ibd);
   }
   if(isEmpty)
      warning("No IBD loaded. This error could be due to missing --marker in Merlin --ibd analysis\n");
   for (int i = 0; i < ped.markerCount; i++)
      markers[i].Sort(ped);
}

void IBDTable::Load(Pedigree & ped, const char *filename, Vector & locusMap)
{
   FILE * f = fopen(filename, "rb");
   if (f == NULL) return;
   if (markers != NULL) delete [] markers;
   markers = new IBDList [locusMap.Length()];

   char buffer[BUFSIZE], famid[BUFSIZE], pid1[BUFSIZE], pid2[BUFSIZE];
   Person * p1;
   Person * p2;
   int      marker;
   IBD      ibd;
   int      check;
   double   position;

   fgets(buffer, BUFSIZE - 1, f);
   while (fgets(buffer, BUFSIZE - 1, f) != NULL) {
      check = sscanf(buffer, " %s %s %s %lf %lf %lf %lf",
         famid, pid1, pid2, &position, &ibd.p0, &ibd.p1, &ibd.p2);
      if(check != 7) continue;
      p1 = ped.FindPerson(famid, pid1);
      p2 = ped.FindPerson(famid, pid2);
      if (p1 == NULL || p2 == NULL) continue;
      marker = -1;
      for(int i = 0; i < locusMap.Length(); i++)
         if(fabs(locusMap[i] - position) < 0.001){
            marker = i;
            break;
         }
      if (marker < 0) continue;
      markers[marker].Append(*p1, *p2, ibd);
   }
   for (int i = 0; i < locusMap.Length(); i++) markers[i].Sort(ped);
   fclose(f);
}

