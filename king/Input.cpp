#include "Input.h"
#include "Error.h"
#include "Constant.h"

#include <stdio.h>
#include <string.h>

int InputPromptWidth = 25;

void safe_gets(char * buffer, int n)
   {
   buffer[0] = 0;

   fgets(buffer, n, stdin);

   for (char * ptr = buffer; *ptr != 0; ptr++)
      if (*ptr == '\n') 
         *ptr = 0;
   }
      
void Input(const char * prompt, int & n, int _default)
   {
   char buffer[BUFSIZE];

   int success;
   do {
      printf("%*s [%8d]: ", InputPromptWidth, prompt, _default);
      safe_gets(buffer, BUFSIZE);
      success = sscanf(buffer, "%d", &n);
      if (success == EOF)
         n = _default;
   } while (success == 0);
   }

void Input(const char * prompt, char & ch, char _default)
   {
   char buffer[BUFSIZE];

   int success;
   do {
      printf("%*s [%8c]: ", InputPromptWidth, prompt, _default);
      safe_gets(buffer, BUFSIZE);
      success = sscanf(buffer, "%c", &ch);
      if (success == EOF)
         ch = _default;
   } while (success == 0);
   }

void Input(const char * prompt, double & d, double _default)
   {
   char buffer[BUFSIZE];

   int success;
   do {
      printf("%*s [%8.2f]: ", InputPromptWidth, prompt, _default);
      safe_gets(buffer, BUFSIZE);
      success = sscanf(buffer, "%lf", &d);
      if (success == EOF)
         d = _default;
   } while (success == 0);
   }

void Input(const char * prompt, bool & b, bool _default)
   {
   char buffer[BUFSIZE];
   int success;
   char c;

   do {
      printf("%*s [%8s]: ", InputPromptWidth, prompt, _default ? "Y/n" : "y/N");
      safe_gets(buffer, BUFSIZE);
      success = sscanf(buffer, "%c", &c);
      if (success == EOF)
         b = _default;
      else
         switch (c)
            {
            case 'y' :
            case 'Y' :
               b = true;
               break;
            case 'n' :
            case 'N' :
               b = false;
               break;
            default :
               success = 0;
            }
   } while (success == 0);
   }


void Input(const char * prompt, char * s, char * _default)
   {
   char buffer[BUFSIZE];

   int success;
   do {
      printf("%*s [%8s]: ", InputPromptWidth, prompt, _default);
      safe_gets(buffer, BUFSIZE);
      success = sscanf(buffer, " %[^\n]", s);
      if (success == EOF)
         strcpy(s, _default);
   } while (success == 0);
   }

void InputBounds(const char * prompt, int & n, int  min, int max,
                 int _default)
   {
   Input(prompt, n, _default);
   while ((n < min) || (n > max))
      {
      printf("\n*** Input value must be between %d and %d ***\n", min, max);
      Input(prompt, n, _default);
      }
   }

void InputBounds(const char * prompt, double & d, double min, double max,
                 double _default)
   {
   Input(prompt, d, _default);
   while ((d < min) || (d > max))
      {
      printf("\n*** Input value must be between %.2f and %.2f ***\n", min, max);
      Input(prompt, d, _default);
      }
   }


