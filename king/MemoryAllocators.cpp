#include "MemoryAllocators.h"

#include <stdlib.h>

char *** AllocateCharCube(int n, int rows, int cols)
   {
   char *** cube = new char ** [n];

   // Stop early if we are out of memory
   if (cube == NULL)
      return NULL;

   for (int i = 0; i < n; i++)
      {
      cube[i] = AllocateCharMatrix(rows, cols);

      // Safely unravel allocation if we run out of memory
      if (cube[i] == NULL)
         {
         while (i--)
            FreeCharMatrix(cube[i], rows);

         delete [] cube;

         return NULL;
         }
      }

   return cube;
   }

int ** AllocateIntMatrix(int rows, int cols)
   {
   int ** matrix = new int * [rows];

   // Stop early if we are out of memory
   if (matrix == NULL)
      return NULL;

   for (int i = 0; i < rows; i++)
      {
      matrix[i] = new int [cols];

      // Safely unravel allocation if we run out of memory
      if (matrix[i] == NULL)
         {
         while (i--)
            delete [] matrix[i];

         delete [] matrix;

         return NULL;
         }
      }

   return matrix;
   }

char ** AllocateCharMatrix(int rows, int cols)
   {
   char ** matrix = new char * [rows];

   // Stop early if we are out of memory
   if (matrix == NULL)
      return NULL;

   for (int i = 0; i < rows; i++)
      {
      matrix[i] = new char [cols];

      // Safely unravel allocation if we run out of memory
      if (matrix[i] == NULL)
         {
         while (i--)
            delete [] matrix[i];

         delete [] matrix;

         return NULL;
         }
      }

   return matrix;
   }

float ** AllocateFloatMatrix(int rows, int cols)
   {
   float ** matrix = new float * [rows];

   // Stop early if we are out of memory
   if (matrix == NULL)
      return NULL;

   for (int i = 0; i < rows; i++)
      {
      matrix[i] = new float [cols];

      // Safely unravel allocation if we run out of memory
      if (matrix[i] == NULL)
         {
         while (i--)
            delete [] matrix[i];

         delete [] matrix;

         return NULL;
         }
      }

   return matrix;
   }

void FreeCharCube(char *** & cube, int n, int rows)
   {
   for (int i = 0; i < n; i++)
      FreeCharMatrix(cube[i], rows);

   delete [] cube;

   cube = NULL;
   }

void FreeCharMatrix(char ** & matrix, int rows)
   {
   for (int i = 0; i < rows; i++)
      delete [] matrix[i];

   delete [] matrix;

   matrix = NULL;
   }

void FreeFloatMatrix(float ** & matrix, int rows)
   {
   for (int i = 0; i < rows; i++)
      delete [] matrix[i];

   delete [] matrix;

   matrix = NULL;
   }

void FreeIntMatrix(int ** & matrix, int rows)
   {
   for (int i = 0; i < rows; i++)
      delete [] matrix[i];

   delete [] matrix;

   matrix = NULL;
   }



