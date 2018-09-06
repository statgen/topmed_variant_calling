#include "MemoryInfo.h"

String & MemoryInfo(double bytes)
   {
   static String info;

   if (bytes < 1024)
      return info = "<1.0 kb";

   if (bytes < 1024. * 1024.)
      info.printf("%.1f kb", (bytes + 1023) / 1024.);
   else if (bytes < 1024. * 1024. * 1024.)
      info.printf("%.1f mb", (bytes + 1024. * 1024. - 1) / (1024. * 1024.));
   else if (bytes < 1024. * 1024. * 1024. * 1024.)
      info.printf("%.1f gb", bytes / (1024. * 1024. * 1024.));
   else
      info.printf("%.1f tb", bytes / (1024. * 1024. * 1024. * 1024.));

   return info;
   }
