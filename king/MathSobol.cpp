#include "MathSobol.h"
#include "Random.h"
#include "Error.h"

#include "stdlib.h"

int SobolSequence::poly_degrees[POLY_COUNT] =
   { 1,  2,  3,  3,  4,  4,  5,  5,  5,  5,  5,  5,
     6,  6,  6,  6,  6,  6,  7,  7,  7,  7,  7,  7,
     7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7
   };

int SobolSequence::poly_integers[POLY_COUNT] =
   {  0, 1,   1,  2,  1,  4,  2,  4,  7, 11, 13, 14,
      1, 13, 16, 19, 22, 25,  1,  4,  7,  8, 14, 19,
     21, 28, 31, 32, 37, 41, 42, 50, 55, 56, 59, 62
   };

SobolSequence::SobolSequence()
   {
   bits = NULL;
   }

SobolSequence::~SobolSequence()
   {
   if (bits != NULL) delete [] bits;
   }

void SobolSequence::Init(int dimensions)
   {
   if (dimensions > POLY_COUNT)
      numerror("Sobol sequences of > %d dimensions not supported", POLY_COUNT);

   x.Dimension(dim = dimensions);
   x.Set(0);
   bits = new IntArray[SOBOL_BITS];

   for (int i = 0; i < SOBOL_BITS; i++)
      bits[i].Dimension(dim);

   unsigned long seed = 0;

   for (int k = 0; k < dim; k++)
      {
      int degrees = poly_degrees[k];

      for (int j = 0; j < degrees; j++)
         // initialize the 0 to kth bit as random odd number <= 2^j - 1
         // and apply a left shift by SOBOL_BITS - j - 1
         {
         bits[j][k] = (RAND(seed) % (1 << j) * 2 | 1);
         bits[j][k] <<= (SOBOL_BITS - j - 1);
         }

      for (int j = degrees; j < SOBOL_BITS; j++)
         // Fill in the remaining values using recurrence
         {
         long poly = poly_integers[k];

         long i = bits[j - degrees][k];
         i ^= (i >> poly_degrees[k]);

         for (int l = j - degrees + 1; l < j; l++)
            {
            if (poly & 1) i ^= bits[l][k];
            poly >>= 1;
            }

         bits[j][k] = i;
         }
      }
   counter = 0;
   }

Vector & SobolSequence::Next(Vector & point)
   {
   long  i = counter, bit;

   for (bit = 0; bit < SOBOL_BITS; bit++)
      {
      if (!(i & 1)) break;
      i >>= 1;
      }

   if (bit == SOBOL_BITS) numerror("SobolSequence is too short");

   for (int k = 0; k < dim; k++)
      {
      x[k] ^= bits[bit][k];
      point[k] = x[k] * SOBOL_FACTOR;
      }

   counter++;

   return point;
   }


