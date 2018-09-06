#ifndef __GENOTYPE_COMPRESSOR_H__
#define __GENOTYPE_COMPRESSOR_H__

#ifndef  uchar
#define uchar        unsigned char
#endif

class GenotypeCompressor
   {
   public:
      static uchar * CompressGenotypes(uchar * genotypes, int n);
      static void    RetrieveGenotypes(uchar * compressed, uchar * genotypes, int n);
      static char *  Describe(uchar * compressed);

      static int     MemoryAllocated();
      static int     MemoryInUse();

   private:
      static uchar * memoryBlocks[1024];
      static int     blockIndex;
      static int     blockByte;

      static void    AllocateBlock();
      static void    AllocateMemory(int size);

      static uchar   OddOneOut(uchar a, uchar b, uchar c);
      static uchar   EncodeTriplet(uchar a, uchar b, uchar c);
      static void    DecodeTriplet(uchar triplet, uchar & a, uchar & b, uchar & c);

      static void    WRITEBIT(uchar * block, uchar & byte, uchar & mask, int bit)
         {
         if (bit) byte |= mask;
         mask *= 2;
         if (mask == 0)
            {
            block[blockByte++] = byte;
            mask = 1;
            byte = 0;
            }
         }

      static bool READBIT(uchar * & input, uchar & mask)
         {
         mask *= 2;

         if (mask == 0) mask = 1, input++;

         return *input & mask;
         }
   };

#endif


