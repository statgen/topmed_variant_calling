////////////////////////////////////////////////////////////////////// 
// libsrc/InputFile.h 
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
 
#ifndef __INPUTFILE_H__
#define __INPUTFILE_H__

#ifdef  __gnu_linux__
#ifndef __ZLIB_AVAILABLE__
#define __ZLIB_AVAILABLE__
#endif
#endif

#ifdef  __ZLIB_AVAILABLE__

#include <zlib.h>
#include <stdio.h>

class IFILE
   {
   public:
      bool gzMode;
      union
         {
         gzFile gzHandle;
         FILE * handle;
         };

   IFILE()
      {
      gzMode = false;
      handle = NULL;
      }

   IFILE(const char * filename, const char * mode);

   operator void * ()
      { return gzMode ? (void *) gzHandle : (void *) handle; }

   IFILE operator = (const IFILE & rhs)
      {
      if ((gzMode = rhs.gzMode) == true)
         gzHandle = rhs.gzHandle;
      else
         handle = rhs.handle;

      return *this;
      }

   IFILE operator = (FILE * rhs)
      {
      gzMode = false;
      handle = rhs;
      return *this;
      }

   IFILE operator = (gzFile & rhs)
      {
      gzMode = true;
      gzHandle = rhs;
      return *this;
      }

   bool operator == (void * rhs)
      {
      if (rhs != NULL)
         return false;
      return gzMode ? gzHandle == rhs : handle == rhs;
      }
   };

inline IFILE ifopen(const char * filename, const char * mode)
   { IFILE file(filename, mode); return file; }

inline int ifclose(IFILE & file)
   { return file.gzMode ? gzclose(file.gzHandle) : fclose(file.handle); }

inline int ifgetc(IFILE & file)
   { return file.gzMode ? gzgetc(file.gzHandle) : fgetc(file.handle); }

inline void ifrewind(IFILE & file)
   { if (file.gzMode) gzrewind(file.gzHandle); else rewind(file.handle); }

inline int ifeof(IFILE & file)
   { return file.gzMode ? gzeof(file.gzHandle) : feof(file.handle); }

#else

#include <stdio.h>

class IFILE
   {
   public:
      FILE * handle;

      IFILE()
         { handle = NULL; }
      IFILE(const char * filename, const char * mode)
         { handle = fopen(filename, mode); }

      operator FILE *()
         { return handle; }

      IFILE & operator = (FILE * rhs)
         { handle = rhs; return *this; }

      IFILE & operator = (const IFILE & rhs)
         { handle = rhs.handle; return * this; }

      bool operator == (void * rhs)
         {
         if (rhs != NULL)
            return false;
         return handle == rhs;
         }
   };

inline IFILE ifopen(const char * filename, const char * mode)
   { IFILE file(filename, mode); return file; }

inline int ifclose(IFILE & file)
   { return fclose(file.handle); }

inline int ifgetc(IFILE & file)
   { return fgetc(file.handle); }

inline void ifrewind(IFILE & file)
   { rewind(file.handle); }

inline int ifeof(IFILE & file)
   { return feof(file.handle); }

#endif

#endif

 
