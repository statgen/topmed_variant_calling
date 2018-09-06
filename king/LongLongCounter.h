#ifndef __LONGLONGCOUNTER_H_
#define __LONGLONGCOUNTER_H_

#include "LongHash.h"

class LongCounter : public LongHash<int>
   {
   public:
      LongCounter();

      void IncrementCount(long long key);
      void DecrementCount(long long key);
      int  GetCount(long long key);
   };

#endif


