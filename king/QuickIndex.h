#ifndef __QUICKINDEX_H__
#define __QUICKINDEX_H__

#include "MathVector.h"
#include "StringArray.h"
#include "IntArray.h"
#include "StringMap.h"

class QuickIndex : public IntArray
   {
   public:
      QuickIndex();
      QuickIndex(const IntArray & source_data)
         { Index(source_data); }
      QuickIndex(const StringArray & source_data)
         { Index(source_data); }
      QuickIndex(const Vector & source_data)
         { Index(source_data); }

      void Index(const IntArray & source_data);
      void Index(const StringArray & source_data);
      void Index(const Vector & source_data);
      void IndexCounts(const StringIntMap & source_data);

   private:
      const void * source;
      int    datatype;

      bool IsBefore(int i, int j);
      void Sort();
   };

#endif

