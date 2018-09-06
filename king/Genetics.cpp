#include "Genetics.h"
#include "Error.h"

#include <stdio.h>

void ImprintingParameter::Status()
   {
   char * msg;

   switch (* (int *) var)
      {
      case I_NONE: msg = "NOT MODELLED"; break;
      case I_MATERNAL: msg = "MATERNAL"; break;
      case I_PATERNAL: msg = "PATERNAL"; break;
      case I_FULL: msg = "FULLY MODELLED"; break;
      case I_IMPRINTING: msg = "TEST IMPRINTING"; break;
      }

   printf("%30s : %15s (-%c[+|-|f|i|m|p])\n", description, msg, ch);
   }

void ImprintingParameter::Translate(char * value)
   {
   switch (tolower(*value))
      {
      case '-' : * (int *) var = I_NONE; break;
      case 'm' : * (int *) var = I_MATERNAL; break;
      case 'p' : * (int *) var = I_PATERNAL; break;
      case 'f' :
      case '+' : 
      case  0  : * (int *) var = I_FULL; break;
      case 'i' : * (int *) var = I_IMPRINTING; break;
      default  : warning("unknown parameter %c%s\n", ch, value);
      };
   }

void GeneticModelParameter::Status()
   {
   char * msg;

   switch (* (int *) var)
      {
      case GM_FREE: msg = "FREE"; break;
      case GM_RECESSIVE: msg = "RECESSIVE"; break;
      case GM_ADDITIVE: msg = "ADDITIVE"; break;
      case GM_DOMINANT: msg = "DOMINANT"; break;
      }

   printf("%30s : %15s (-%c[a|d|f|r])\n", description, msg, ch);
   }

void GeneticModelParameter::Translate(char * value)
   {
   switch (tolower(*value))
      {
      case 'a' : * (int *) var = GM_ADDITIVE; break;
      case 'd' : * (int *) var = GM_DOMINANT; break;
      case 'f' : * (int *) var = GM_FREE; break;
      case 'r' : * (int *) var = GM_RECESSIVE; break;
      default  : warning("unknown parameter %c%s\n", ch, value);
      };
   }



