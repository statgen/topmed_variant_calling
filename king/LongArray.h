#ifndef __LONGINTARRAY_H__
#define __LONGINTARRAY_H__

#include "LongInt.h"

class LongArray
   {
   private:
      longint * items;
      int size, count;

      void Grow(int new_size);
      static int Compare(int * a, int * b);

   public:
      static int alloc;

      LongArray(int start_size = 0);
      LongArray(LongArray & source);
      ~LongArray();

      LongArray & operator = (const LongArray & rhs);

      longint & operator [] (int index) { return items[index]; }

      int  Append(longint value);
      void Push(longint value)      { Append(value); }
      longint Pop()                 { return items[--count]; }
      longint Peek() const          { return items[count - 1]; }
      longint &Last() const         { return items[count - 1]; }

      int  Delete(int index);
      void InsertAt(int index, longint value);

      int  Find(longint value) const;
      void Sort();

      void Zero();
      void Set(longint value);

      int  Length()                 { return count; }
      void Dimension(int new_count) { Grow(new_count); count = new_count; }
      void Clear()                  { count = 0; }

      void Swap(int i, int j)
           { longint tmp = items[i]; items[i] = items[j]; items[j] = tmp; }

      void Reverse();

      operator longint * () { return items; }

      bool operator == (const LongArray & rhs) const;
      bool operator != (const LongArray & rhs) const;

      int Hash(int initval);
   };

#endif /* __LONGINTARRAY_H */


