#include "GenotypeCompressor.h"

#include <stdio.h>

#define       BLOCK_SIZE (4 * 1024  * 1024)

#define       ENCODE_CONSTANT       (0 * 32)
#define       ENCODE_RUNLENGTH      (1 * 32)
#define       ENCODE_ONEBIT         (2 * 32)
#define       ENCODE_TWOBITS        (3 * 32)
#define       ENCODE_HUFFMAN3       (4 * 32)
#define       ENCODE_HUFFMAN4       (5 * 32)
#define       ENCODE_MASK            31

uchar * GenotypeCompressor::memoryBlocks[1024];
int    GenotypeCompressor::blockIndex = -1;
int    GenotypeCompressor::blockByte = BLOCK_SIZE;

void GenotypeCompressor::AllocateBlock()
   {
   memoryBlocks[++blockIndex] = new uchar [BLOCK_SIZE];
   blockByte = 0;
   }

void GenotypeCompressor::AllocateMemory(int n)
   {
   // printf("Allocating %d bytes\n", n);

   if (blockByte + n >= BLOCK_SIZE)
      AllocateBlock();
   }

uchar GenotypeCompressor::OddOneOut(uchar a, uchar b, uchar c)
   {
   return a ^ b ^ c;
   }

uchar GenotypeCompressor::EncodeTriplet(uchar a, uchar b, uchar c)
   {
   if (b > a) b--;
   if (c > a) c--;
   if (c > b) c--;

   return a * 8 + b * 2 + c;
   }

void GenotypeCompressor::DecodeTriplet(uchar triplet, uchar & a, uchar & b, uchar & c)
   {
   a = triplet / 8;
   b = (triplet & 7) / 2;
   c = triplet & 1;

   if (c >= b) c++;
   if (c >= a) c++;
   if (b >= a) b++;
   }

uchar * GenotypeCompressor::CompressGenotypes(uchar * genotypes, int n)
   {
   if (n == 0)
      return 0;

   // First we count the number of runs of identical genotypes
   // and the number of occurences of each of the four possible
   // genotypes ...

   int counts[4] = {0, 0, 0, 0};
   int runs = 0, length = 64;

   counts[int(genotypes[0])]++;
   for (int i = 1; i < n; i++)
      {
      counts[int(genotypes[i])]++;

      if (genotypes[i] == genotypes[i - 1])
         length++;
      else
         {
         runs += length / 64;
         length = 64;
         }
      }
   runs += length / 64;

   // Find out the genotype with highest and lowest counts ...
   int ihi, ilo, inlo;

   if (counts[1] > counts[0])
      ihi = inlo = 1, ilo = 0;
   else
      ihi = inlo = 0, ilo = 1;

   int zeros = (counts[0] == 0) + (counts[1] == 0) + (counts[2] == 0) + (counts[3] == 0);

   for (int i = 2; i < 4; i++)
      if (counts[i] > counts[ihi])
         ihi = i;
      else if (counts[i] <= counts[ilo])
         { inlo = ilo; ilo = i; }
      else if (counts[i] <= counts[inlo])
         inlo = i;

   int inhi = OddOneOut(ihi, ilo, inlo);

   // Intermediate variables for encoding
   uchar mask = 1;
   uchar byte = 0;

   // Finally, decide on encoding strategy and encode the data...
   if (zeros == 3)
      {
      AllocateMemory(1);
      uchar * start = memoryBlocks[blockIndex] + blockByte;
      uchar * block = memoryBlocks[blockIndex];

      block[blockByte++] = ENCODE_CONSTANT + ihi;

      return start;
      }
   else if (zeros == 2)
      {
      if (n / 8 < runs)
         {
         AllocateMemory((n + 15) / 8);
         uchar * start = memoryBlocks[blockIndex] + blockByte;
         uchar * block = memoryBlocks[blockIndex];

         block[blockByte++] = ENCODE_ONEBIT + ihi * 4 + inhi;

         for (int i = 0; i < n; i++)
            if (genotypes[i] == ihi)
               WRITEBIT(block, byte, mask, 0);
            else
               WRITEBIT(block, byte, mask, 1);

         if (mask > 1)
            block[blockByte++] = byte;

         return start;
         }
      }
   else if (zeros == 1)
      {
      int packed = (n * 2 - counts[ihi] + 7) / 8;

      if (packed < runs)
         {
         AllocateMemory(packed + 1);
         uchar * start = memoryBlocks[blockIndex] + blockByte;
         uchar * block = memoryBlocks[blockIndex];

         block[blockByte++] = ENCODE_HUFFMAN3 + EncodeTriplet(ihi, inhi, inlo);

         for (int i = 0; i < n; i++)
            if (genotypes[i] == ihi)
               WRITEBIT(block, byte, mask, 0);
             else if (genotypes[i] == inhi)
               {
               WRITEBIT(block, byte, mask, 1);
               WRITEBIT(block, byte, mask, 0);
               }
             else
               {
               WRITEBIT(block, byte, mask, 1);
               WRITEBIT(block, byte, mask, 1);
               }

         if (mask > 1)
            block[blockByte++] = byte;

         return start;
         }
      }
   else if (counts[ihi] > counts[ilo] + counts[inlo])
      {
      int packed = (n * 3 - counts[ihi] * 2 - counts[inhi] + 7) / 8;

      if (packed < runs)
         {
         AllocateMemory(packed + 1);
         uchar * start = memoryBlocks[blockIndex] + blockByte;
         uchar * block = memoryBlocks[blockIndex];

         block[blockByte++] = ENCODE_HUFFMAN4 + EncodeTriplet(ihi, inhi, inlo);

         for (int i = 0; i < n; i++)
            if (genotypes[i] == ihi)
               WRITEBIT(block, byte, mask, 0);
            else if (genotypes[i] == inhi)
               {
               WRITEBIT(block, byte, mask, 1);
               WRITEBIT(block, byte, mask, 0);
               }
            else if (genotypes[i] == inlo)
               {
               WRITEBIT(block, byte, mask, 1);
               WRITEBIT(block, byte, mask, 1);
               WRITEBIT(block, byte, mask, 0);
               }
            else
               {
               WRITEBIT(block, byte, mask, 1);
               WRITEBIT(block, byte, mask, 1);
               WRITEBIT(block, byte, mask, 1);
               }

         if (mask > 1)
            block[blockByte++] = byte;

         return start;
         }
      }
   else
      {
      int packed = (n + 3) / 4;

      if (packed < runs)
         {
         AllocateMemory(packed + 1);
         uchar * start = memoryBlocks[blockIndex] + blockByte;
         uchar * block = memoryBlocks[blockIndex];

         block[blockByte++] = ENCODE_TWOBITS;

         for (int i = 0; i < n; i++)
            {
            WRITEBIT(block, byte, mask, genotypes[i] & 1);
            WRITEBIT(block, byte, mask, genotypes[i] & 2);
            }

         if (mask > 1)
            block[blockByte++] = byte;

         return start;
         }
      }

   AllocateMemory(runs + 1);
   uchar * start = memoryBlocks[blockIndex] + blockByte;
   uchar * block = memoryBlocks[blockIndex];

   // If we get here, then use run length encoding ...
   block[blockByte++] = ENCODE_RUNLENGTH;
   length = 0;

   for (int i = 0; i < n; i++)
      {
      if (genotypes[i] == genotypes[i - 1])
         length++;
      else
         {
         if (length)
            block[blockByte++] = (length << 2) + genotypes[i - 1];
         length = 1;
         }

      if (length == 64)
         {
         block[blockByte++] = (length << 2) + genotypes[i - 1];
         length = 0;
         }
      }

   if (length)
      block[blockByte++] = (length << 2) + genotypes[n - 1];

   return start;
   }

char * GenotypeCompressor::Describe(uchar * compressed)
   {
   switch (*compressed & ~ENCODE_MASK)
      {
      case ENCODE_CONSTANT:
         return "ENCODE_CONSTANT";
      case ENCODE_ONEBIT:
         return "ENCODE_ONEBIT";
      case ENCODE_TWOBITS:
         return "ENCODE_TWOBITS";
      case ENCODE_HUFFMAN3:
         return "ENCODE_HUFFMAN3";
      case ENCODE_HUFFMAN4:
         return "ENCODE_HUFFMAN4";
      default:
         return "ENCODE_RUNLENGTH";
      }

   }

void GenotypeCompressor::RetrieveGenotypes(uchar * compressed, uchar * genotypes, int n)
   {
   uchar mask = 0;

   switch (*compressed & ~ENCODE_MASK)
      {
      case ENCODE_CONSTANT:
         {
         uchar out = *compressed & ENCODE_MASK;

         for (int i = 0; i < n; i++)
            genotypes[i] = out;
         } break;
      case ENCODE_ONEBIT:
         {
         uchar a = (*compressed & ENCODE_MASK) / 4;
         uchar b = (*compressed & ENCODE_MASK) & 3;

         for (int i = 0; i < n; i++)
            *genotypes++ = READBIT(compressed, mask) ? b : a;
         } break;
      case ENCODE_TWOBITS:
         {
         for (int i = 0; i < n; i++)
            *genotypes++ = READBIT(compressed, mask) + READBIT(compressed, mask) * 2;
         } break;
      case ENCODE_HUFFMAN3:
         {
         uchar a, b, c;
         DecodeTriplet(*compressed & ENCODE_MASK, a, b, c);

         for (int i = 0; i < n; i++)
            if (!READBIT(compressed, mask))
               *genotypes++ = a;
            else if (!READBIT(compressed, mask))
               *genotypes++ = b;
            else
               *genotypes++ = c;
         } break;
      case ENCODE_HUFFMAN4:
         {
         uchar a, b, c;
         DecodeTriplet(*compressed & ENCODE_MASK, a, b, c);
         uchar d = OddOneOut(a, b, c);

         for (int i = 0; i < n; i++)
            if (!READBIT(compressed, mask))
               *genotypes++ = a;
            else if (!READBIT(compressed, mask))
               *genotypes++ = b;
            else if (!READBIT(compressed, mask))
               *genotypes++ = c;
            else
               *genotypes++ = d;
         } break;
      default:
         {
         while (n > 0)
            {
            compressed++;
            int a = *compressed & 3;
            int l = *compressed >> 2;

            if (l == 0) l = 64;
            if (l > n)  l = n;

            n -= l;

            while (l--) *genotypes++ = a;
            }
         }
      }
   }

int GenotypeCompressor::MemoryAllocated()
   {
   return (blockIndex + 1) * BLOCK_SIZE;
   }

int GenotypeCompressor::MemoryInUse()
   {
   if (blockIndex == -1)
      return 0;

   return blockIndex * BLOCK_SIZE + blockByte;
   }
