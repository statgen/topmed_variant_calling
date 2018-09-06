#ifndef  __GENETICS_H__
#define  __GENETICS_H__

#include "Parameters.h"

// Genetic models
#define  GM_FREE        0
#define  GM_RECESSIVE   1
#define  GM_ADDITIVE    2
#define  GM_DOMINANT    3

// Constants for imprinting analysis
#define  I_NONE         0
#define  I_PATERNAL     1
#define  I_MATERNAL     2
#define  I_FULL         3
#define  I_IMPRINTING   4

// Constants for special effects
#define SFX_NONE        0
#define SFX_PATERNAL    1
#define SFX_MATERNAL    2

class ImprintingParameter : public Parameter
   {
   public:
   ImprintingParameter(char c, char * desc, int & v)
      : Parameter(c, desc, &v)
      {}

   virtual void Status();

   protected:
      virtual void Translate(char * value);
   };

class GeneticModelParameter : public Parameter
   {
   public:
   GeneticModelParameter(char c, char * desc, int & v)
      : Parameter(c, desc, &v)
      {}

   virtual void Status();

   protected:
      virtual void Translate(char * value);
   };

#endif

