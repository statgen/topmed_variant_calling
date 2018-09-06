#include "WindowsHelper.h"
#ifdef    __WIN32__
#ifndef   __GNUC__
#include <dir.h>

void WildCardArguments(int & argc, char ** & argv)
   {
   if (argc < 2) return;

   int  count = 0;
   for (int i = 1; i < argc; i++)
      {
      struct ffblk blk;

      int done = findfirst(argv[i], &blk, 0);
      while(!done)
         {
         done = findnext(&blk);
         count++;
         }
      }

   char ** new_argv = new char * [count + 1];
   int     new_argc = 1;

   new_argv[0] = argv[0];
   for (int i = 1; i < argc; i++)
      {
      struct ffblk blk;

      int done = findfirst(argv[i], &blk, 0);
      while (!done && new_argc <= count)
         {
         new_argv[new_argc++] = strdup(blk.ff_name);
         done = findnext(&blk);
         }
      }

   argc = new_argc;
   argv = new_argv;
   }

#endif
#endif

