////////////////////////////////////////////////////////////////////// 
// libsrc/StringArray.cpp 
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
 
#include "StringArray.h"
#include "Sort.h"
#include "Error.h"

int StringArray::alloc = 32;

StringArray::StringArray(int startsize)
   {
   count = startsize;
   size = (startsize + alloc) / alloc * alloc;
   strings = new String * [size];
   for (int i = 0; i < count; i++)
      strings[i] = new String;
   };

StringArray::StringArray(StringArray & rhs)
   {
   count = rhs.count;
   size = (rhs.count + alloc) / alloc * alloc;
   strings = new String * [size];

   for (int i = 0; i < count; i++)
      strings[i] = new String(rhs[i]);;
   }

StringArray::~StringArray()
   {
   for (int i = 0; i < count; i++)
      delete strings[i];
   delete [] strings;
   }

int StringArray::CharLength()
   {
   int charlen = 0;
   for (int i = 0; i < count; i++)
      charlen += strings[i]->Length();
   return charlen;
   }

void StringArray::Read(const char * filename)
   {
   IFILE f = ifopen(filename, "rb");
   if (f == NULL) return;
   Read(f);
   ifclose(f);
   }

void StringArray::Write(const char * filename)
   {
   FILE * f = fopen(filename, "wt");
   if (f == NULL) return;
   Write(f);
   fclose(f);
   }

void StringArray::WriteLine(const char * filename)
   {
   FILE * f = fopen(filename, "wt");
   if (f == NULL) return;
   WriteLine(f);
   fclose(f);
   }

void StringArray::Read(FILE * f)
   {
   while (!feof(f))
      {
      Grow(count + 1);
      strings[count] = new String;
      strings[count]->ReadLine(f);
      count++;
      }
   }

void StringArray::Write(FILE * f)
   {
   for(int i = 0; i < count; i++)
      strings[i]->WriteLine(f);
   }

void StringArray::WriteLine(FILE * f)
   {
   for (int i = 0; i < count; i++)
      fprintf(f, "%s%c", (const char *) (*strings[i]), i == count-1 ? '\n' : '\t');
   }

#ifdef __ZLIB_AVAILABLE__
void StringArray::Read(IFILE & f)
   {
   while (!ifeof(f))
      {
      Grow(count + 1);
      strings[count] = new String;
      strings[count]->ReadLine(f);
      count++;
      }
   }
#endif

void StringArray::Grow(int newsize)
   {
   if (newsize >= size)
      {
      if ((newsize >> 1) >= size)
         size = (newsize + alloc) / alloc * alloc;
      else
         {
         size = alloc;
         while (size <= newsize)
            size *= 2;
         }
      String ** tmp = new String * [size];
      for (int i = 0; i < count; i++) tmp[i] = strings[i];
      delete [] strings;
      strings = tmp;
      }
   }

void StringArray::Clear()
   {
   for (int i = 0; i < count; i++)
      delete strings[i];
   count = 0;
   }

int StringArray::AddColumns(const String & s, char ch)
   {
   for (int pos = 0; pos <= s.Length(); pos++)
      {
      int oldpos = pos;
      pos = s.FindChar(ch, pos);
      if (pos == -1) pos = s.Length();
      Grow(count + 1);
      strings[count++] = new String(s.Mid(oldpos, pos - 1));
      };

   return count;
   }

int StringArray::AddTokens(const String & s, char ch)
   {
   for (int pos = 0; pos < s.Length(); pos++)
      {
      while (pos < s.Length() && s[pos] == ch) pos++;
      int oldpos = pos;

      while (pos < s.Length() && s[pos] != ch) pos++;

      if (oldpos < s.Length())
         {
         Grow(count + 1);
         strings[count++] = new String(s.Mid(oldpos, pos - 1));
         }
      }

   return count;
   }

int StringArray::AddTokens(const String & s, const String & separators)
   {
   for (int pos = 0; pos < s.Length(); pos++)
      {
      while (pos < s.Length() && separators.FindChar(s[pos]) != -1) pos++;
      int oldpos = pos;

      while (pos < s.Length() && separators.FindChar(s[pos]) == -1) pos++;

      if (oldpos < s.Length())
         {
         Grow(count + 1);
         strings[count++] = new String(s.Mid(oldpos, pos - 1));
         }
      }

   return count;
   }

int StringArray::Dimension(int newcount)
   {
   if (newcount > count)
      {
      Grow(newcount);
      for (int i = count; i < newcount; i++)
         strings[i] = new String;
      count = newcount;
      }
   else if (newcount < count)
      {
      for (int i = newcount; i < count; i++)
         delete strings[i];
      count = newcount;
      }

   return count;
   }

int StringArray::Find(const String & s) const
   {
   for (int i = 0; i < count; i++)
      if (*(strings[i]) == s)
         return i;
   return -1;
   }

int StringArray::FastFind(const String & s) const
   {
   for (int i = 0; i < count; i++)
      if (strings[i]->FastCompare(s) == 0)
         return i;
   return -1;
   }

int StringArray::SlowFind(const String & s) const
   {
   for (int i = 0; i < count; i++)
      if (strings[i]->SlowCompare(s) == 0)
         return i;
   return -1;
   }

int StringArray::Add(const String & s)
   {
   Grow(count + 1);
   strings[count] = new String(s);
   return ++count;
   }

void StringArray::InsertAt(int position, const String & s)
   {
   Grow(count + 1);
   for (int i = count; i > position; i--)
      strings[i] = strings[i - 1];
   strings[position] = new String(s);
   count++;
   }

String & StringArray::Last() const
   {
   if (!count) error("StringArray: Null String Access");
   return *(strings[count - 1]);
   }

void StringArray::Delete(int index)
   {
   delete strings[index];
   count--;
   for ( ; index < count; index++)
      strings[index] = strings[index + 1];
   }

StringArray & StringArray::operator = (const StringArray & rhs)
   {
   Clear();
   for (int i = 0; i < rhs.count; i++)
      Push(*rhs.strings[i]);
   return *this;
   }

bool StringArray::operator == (const StringArray & rhs)
   {
   if (count != rhs.count) return false;
   for (int i = 0; i < rhs.count; i++)
      if (*strings[i] != *rhs.strings[i])
         return false;
   return true;
   }

void StringArray::Sort()
   {
   QuickSort(strings, count, sizeof(String *), ComparisonForSort);
   }

int StringArray::ComparisonForSort(const void * a, const void * b)
   {
   String * string1 = *(String **) a;
   String * string2 = *(String **) b;

   return Compare(*string1, *string2);
   }

String StringArray::Pop()
   {
   String result = *(strings[count - 1]);

   Dimension(count - 1);

   return result;
   }

void StringArray::Trim()
   {
   for (int i = 0; i < count; i++)
      strings[i]->Trim();
   }

void StringArray::Print()
   {
   for (int i = 0; i < count; i++)
      printf("%s\n", (const char *) (*strings[i]));
   }

void StringArray::PrintLine()
   {
   for (int i = 0; i < count; i++)
      printf("%s%c", (const char *) (*strings[i]), i == count - 1 ? '\n' : '\t');
   }
 
