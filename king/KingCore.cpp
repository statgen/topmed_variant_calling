//////////////////////////////////////////////////////////////////////
// KingCore.cpp
// (c) 2010-2019 Wei-Min Chen
//
// This file is distributed as part of the KING source code package
// and may not be redistributed in any form, without prior written
// permission from the author. Permission is granted for you to
// modify this file for your own personal use, but modified versions
// must retain this copyright notice and must not be distributed.
//
// Permission is granted for you to use this file to compile KING.
//
// All computer programs have bugs. Use this file at your own risk.
//
// Feb 25, 2019

#include <math.h>
#include "KingCore.h"
#include "Kinship.h"
#include "MathStats.h"
#include "MathSVD.h"
#include "QuickIndex.h"
#include "MathCholesky.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

KingEngine::KingEngine(Pedigree & pedigree):ped(pedigree)
{
   minMAF = 0.01;
   freqLoaded = false;
   individualInfo = false;
   GG[0] = NULL;
   LG[0] = NULL;
   SG[0] = NULL;
   SLG[0] = NULL;
   id = NULL;
   errorrateCutoff = _NAN_;
   for(int i = 0; i < 8; i++)
      base[i] = 1 << i;
   for(int i = 0; i < 16; i++)
      shortbase[i] = 1 << i;
   for(int i = 0; i < 256; i++){
      oneCount[i] = 0;
      for(int j = 0; j < 8; j++)
         if(i & base[j]) oneCount[i]++;
   }
   allflags = 0;
   QCwipe = false;
   genoFilterFlag = false;
   shortFlag = true;
   shortFlip = NULL;
   diagFlag = false;
   callrate = _NAN_;
   exclusionList.Dimension(0);
   for(int i = 0; i < 5; i++) inclusionList[i].Dimension(0);
   xsnpName.Dimension(0);
   markerCount = xmarkerCount = ymarkerCount = mtmarkerCount = shortCount = xshortCount = yshortCount = mtshortCount = 0;
   ysnpName.Dimension(0);
   mtsnpName.Dimension(0);
   start = stop = _NAN_;
   autoflipFlag = false;
   flipFlag = false;
   chrList.Dimension(0);
   SaveFormat = "KING";
   pedstatFlag = false;
   Bit64=0;
   Bit64Flag = false;
   defaultMaxCoreCount=defaultEfficientCoreCount=1;
   SEXCHR=23;
   lessmemFlag=false;
   G_SNP = NULL;
}

KingEngine::~KingEngine()
{
   if(GG[0]){
      delete []GG[0];
      delete []GG[1];
   }
   if(LG[0]){
      delete []LG[0];
      delete []LG[1];
   }
   if(SG[0]){
      delete []SG[0];
      delete []SG[1];
   }
   if(SLG[0]){
      delete []SLG[0];
      delete []SLG[1];
   }
   if(id) delete []id;
   if(shortFlip) delete []shortFlip;
}

void KingEngine::ComputeMZBySNP64Bit(IntArray &nonmissingMZCounts, IntArray &ibs1MZCounts, IntArray &HetHetMZCounts, IntArray &ibs0MZCounts)
{
   const int BLOCKSIZE=255;
   const int CACHESIZE=16; // total cache size: 2^(6 + 4 - 2 + 8 + 1) = 2^17 = 128KB
   const int CACHESIZE_BITPAR=CACHESIZE<<3;  // 128
   unsigned long long int word, genobit[4][CACHESIZE_BITPAR];
   unsigned char *pchar;
   int id1, id2;
   int thread = 0;
   int L0Count = (L0.Length()>>1);
   HetHetMZCounts.Dimension(markerCount);
   ibs0MZCounts.Dimension(markerCount);
   ibs1MZCounts.Dimension(markerCount);
   nonmissingMZCounts.Dimension(markerCount);
   HetHetMZCounts.Zero();
   ibs0MZCounts.Zero();
   ibs1MZCounts.Zero();
   nonmissingMZCounts.Zero();
   int **genobitCount[4];
   for(int k = 0; k < 4; k++){
      genobitCount[k] = new int *[defaultMaxCoreCount];
      for(int i = 0; i < defaultMaxCoreCount; i++)
         genobitCount[k][i] = NULL;
   }
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(word, genobit, pchar, thread, id1, id2)
{
   thread = omp_get_thread_num();
   for(int k = 0; k < 4; k++){
      genobitCount[k][thread] = new int [(longCount<<6)];
      for(int m = 0; m < markerCount; m++)
         genobitCount[k][thread][m] = 0;
   }
   #pragma omp for
#else
   for(int k = 0; k < 4; k++){
      genobitCount[k][thread] = new int [longCount<<6];
      for(int m = 0; m < markerCount; m++)
         genobitCount[k][thread][m] = 0;
   }
#endif
   for(int bi = 0; bi < L0Count; bi+=BLOCKSIZE){   // blocked individuals
      int iMax = (bi > L0Count-BLOCKSIZE) ? L0Count: bi+BLOCKSIZE;
      for(int blockb = 0; blockb < longCount; blockb += CACHESIZE){  // blocked SNP words
         int bMax = (blockb >= longCount-CACHESIZE) ? longCount: blockb+CACHESIZE;
         for(int  k = 0;  k < 4; k++)
            for(int j = 0; j < CACHESIZE_BITPAR; j++)
               genobit[k][j] = 0;
         for(int i = bi; i < iMax; i++){
            id1 = geno[L0[i*2]];
            id2 = geno[L0[i*2+1]];
            for(int b = blockb; b < bMax; b++){
               int bbb = ((b-blockb)<<3);
               word = (LG[0][id1][b] | LG[1][id1][b]) & (LG[0][id2][b] | LG[1][id2][b]); // nonmissing
               genobit[0][bbb] += (word&0x0101010101010101);
               for(int k = 1; k < 8; k++)
                  genobit[0][bbb|k] += ((word>>k)&0x0101010101010101);
               word &= (LG[0][id1][b] ^ LG[0][id2][b]);  // IBS1
               genobit[1][bbb] += (word&0x0101010101010101);
               for(int k = 1; k < 8; k++)
                  genobit[1][bbb|k] += ((word>>k)&0x0101010101010101);
               word = (~LG[0][id1][b]) & LG[1][id1][b] & (~LG[0][id2][b]) & LG[1][id2][b];   // HetHet
               genobit[2][bbb] += (word&0x0101010101010101);
               for(int k = 1; k < 8; k++)
                  genobit[2][bbb|k] += ((word>>k)&0x0101010101010101);
               word = LG[0][id1][b] & LG[0][id2][b] & (LG[1][id1][b]^LG[1][id2][b]);   // IBS0
               genobit[3][bbb] += (word&0x0101010101010101);
               for(int k = 1; k < 8; k++)
                  genobit[3][bbb|k] += ((word>>k)&0x0101010101010101);
            }  // end of blocked marker
         }  // end of blocked individual
         for(int b = blockb; b < bMax; b++){
            int bbb = ((b-blockb)<<3);
            for(int k = 0; k < 8; k++){   // shift
               int pos = (b<<6) | k;
               for(int bit = 0; bit < 3; bit++){
                  pchar = (unsigned char *)&genobit[bit][bbb|k];
                  for(int j = 0; j < 8; j++) // bit parallel
                     genobitCount[bit][thread][pos | (j<<3)] += pchar[j];
               }  // end of bit
            }  // end of 8 shifts
         }  // end of b
      }  // end of blockb
   }  // end of bi
#ifdef _OPENMP
}
#endif
   for(int i = 0; i < defaultMaxCoreCount; i++)
      if(genobitCount[0][i]){  // thread effectively used
         for(int m = 0; m < markerCount; m++){
            nonmissingMZCounts[m] += genobitCount[0][i][m];
            ibs1MZCounts[m] += genobitCount[1][i][m];
            HetHetMZCounts[m] += genobitCount[2][i][m];
            ibs0MZCounts[m] += genobitCount[3][i][m];
         }
         for(int k = 0; k < 4; k++)
            delete []genobitCount[k][i];
      }
   for(int k = 0; k < 4; k++)
      delete []genobitCount[k];
}

void KingEngine::ComputeTrioBySNP64Bit(IntArray &MItrioCounts, IntArray &HetInOffspringCounts, IntArray &nonmissingtrioCounts)
{
   const int BLOCKSIZE=255;
   const int CACHESIZE=8; // total cache size: 2^(6 + 4 - 2 + 8 + 1) = 2^17 = 128KB
   const int CACHESIZE_BITPAR=CACHESIZE<<3;
   unsigned long long int word, genobit[3][CACHESIZE_BITPAR];
   unsigned char *pchar;
   int id1, id2, id3;
   int thread = 0;
   int LtrioCount = Ltrio.Length()/3;
   MItrioCounts.Dimension(markerCount);
   HetInOffspringCounts.Dimension(markerCount);
   nonmissingtrioCounts.Dimension(markerCount);
   MItrioCounts.Zero();
   HetInOffspringCounts.Zero();
   nonmissingtrioCounts.Zero();
   int **genobitCount[3];
   for(int k = 0; k < 3; k++){
      genobitCount[k] = new int *[defaultMaxCoreCount];
      for(int i = 0; i < defaultMaxCoreCount; i++)
         genobitCount[k][i] = NULL;
   }
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(word, genobit, pchar, thread, id1, id2, id3)
{
   thread = omp_get_thread_num();
   for(int k = 0; k < 3; k++){
      genobitCount[k][thread] = new int [(longCount<<6)];
      for(int m = 0; m < markerCount; m++)
         genobitCount[k][thread][m] = 0;
   }
   #pragma omp for
#else
   for(int k = 0; k < 3; k++){
      genobitCount[k][thread] = new int [longCount<<6];
      for(int m = 0; m < markerCount; m++)
         genobitCount[k][thread][m] = 0;
   }
#endif
   for(int bi = 0; bi < LtrioCount; bi+=BLOCKSIZE){   // blocked individuals
      int iMax = (bi > LtrioCount-BLOCKSIZE) ? LtrioCount: bi+BLOCKSIZE;
      for(int blockb = 0; blockb < longCount; blockb += CACHESIZE){  // blocked SNP words
         int bMax = (blockb >= longCount-CACHESIZE) ? longCount: blockb+CACHESIZE;
         for(int  k = 0;  k < 3; k++)
            for(int j = 0; j < CACHESIZE_BITPAR; j++)
               genobit[k][j] = 0;
         for(int i = bi; i < iMax; i++){
            id1 = geno[Ltrio[i*3+1]];
            id2 = geno[Ltrio[i*3+2]];
            id3 = geno[Ltrio[i*3]];   // id3 is the child
            for(int b = blockb; b < bMax; b++){
               int bbb = ((b-blockb)<<3);
               word = LG[0][id1][b] & LG[0][id2][b] & (~(LG[1][id1][b] ^ LG[1][id2][b]))
               & (~LG[0][id3][b]) & LG[1][id3][b];   // AA x AA -> Aa, or aa x aa -> Aa
               genobit[0][bbb] += (word&0x0101010101010101);
               for(int k = 1; k < 8; k++)
                  genobit[0][bbb|k] += ((word>>k)&0x0101010101010101);
               word = (LG[0][id1][b] | LG[1][id1][b]) & (LG[0][id2][b] | LG[1][id2][b])
                  & (LG[0][id3][b] | LG[1][id3][b]);  // nonmissing
               genobit[1][bbb] += (word&0x0101010101010101);
               for(int k = 1; k < 8; k++)
                  genobit[1][bbb|k] += ((word>>k)&0x0101010101010101);
               word &= ((~LG[0][id3][b]) & LG[1][id3][b]); // HetInOffspring
               genobit[2][bbb] += (word&0x0101010101010101);
               for(int k = 1; k < 8; k++)
                  genobit[2][bbb|k] += ((word>>k)&0x0101010101010101);
            }  // end of blocked marker
         }  // end of blocked individual
         for(int b = blockb; b < bMax; b++){
            int bbb = ((b-blockb)<<3);
            for(int k = 0; k < 8; k++){   // shift
               int pos = (b<<6) | k;
               for(int bit = 0; bit < 3; bit++){
                  pchar = (unsigned char *)&genobit[bit][bbb|k];
                  for(int j = 0; j < 8; j++) // bit parallel
                     genobitCount[bit][thread][pos | (j<<3)] += pchar[j];
               }  // end of bit
            }  // end of 8 shifts
         }  // end of b
      }  // end of blockb
   }  // end of bi
#ifdef _OPENMP
}
#endif
   for(int i = 0; i < defaultMaxCoreCount; i++)
      if(genobitCount[0][i]){  // thread effectively used
         for(int m = 0; m < markerCount; m++){
            MItrioCounts[m] += genobitCount[0][i][m];
            nonmissingtrioCounts[m] += genobitCount[1][i][m];
            HetInOffspringCounts[m] += genobitCount[2][i][m];
         }
         for(int k = 0; k < 3; k++)
            delete []genobitCount[k][i];
      }
   for(int k = 0; k < 3; k++)
      delete []genobitCount[k];
}

void KingEngine::ComputePOBySNP64Bit(IntArray &HomHomCounts, IntArray &ibs0Counts, IntArray &nonmissingCounts)
{
   const int BLOCKSIZE=255;
   const int CACHESIZE=16; // total cache size: 2^(6 + 4 - 2 + 8 + 1) = 2^17 = 128KB
   const int CACHESIZE_BITPAR=CACHESIZE<<3;  // 128
   unsigned long long int word, genobit[3][CACHESIZE_BITPAR];
   unsigned char *pchar;
   int id1, id2;
   int thread = 0;
   int LpoCount = (Lpo.Length()>>1);
   HomHomCounts.Dimension(markerCount);
   ibs0Counts.Dimension(markerCount);
   nonmissingCounts.Dimension(markerCount);
   HomHomCounts.Zero();
   ibs0Counts.Zero();
   nonmissingCounts.Zero();    
   int **genobitCount[3];
   for(int k = 0; k < 3; k++){
      genobitCount[k] = new int *[defaultMaxCoreCount];
      for(int i = 0; i < defaultMaxCoreCount; i++)
         genobitCount[k][i] = NULL;
   }
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(word, genobit, pchar, thread, id1, id2)
{
   thread = omp_get_thread_num();
   for(int k = 0; k < 3; k++){
      genobitCount[k][thread] = new int [(longCount<<6)];
      for(int m = 0; m < markerCount; m++)
         genobitCount[k][thread][m] = 0;
   }
   #pragma omp for
#else
   for(int k = 0; k < 3; k++){
      genobitCount[k][thread] = new int [longCount<<6];
      for(int m = 0; m < markerCount; m++)
         genobitCount[k][thread][m] = 0;
   }
#endif
   for(int bi = 0; bi < LpoCount; bi+=BLOCKSIZE){   // blocked individuals
      int iMax = (bi > LpoCount-BLOCKSIZE) ? LpoCount: bi+BLOCKSIZE;
      for(int blockb = 0; blockb < longCount; blockb += CACHESIZE){  // blocked SNP words
         int bMax = (blockb >= longCount-CACHESIZE) ? longCount: blockb+CACHESIZE;
         for(int  k = 0;  k < 3; k++)
            for(int j = 0; j < CACHESIZE_BITPAR; j++)
               genobit[k][j] = 0;
         for(int i = bi; i < iMax; i++){
            id1 = geno[Lpo[i*2]];
            id2 = geno[Lpo[i*2+1]];
            for(int b = blockb; b < bMax; b++){
               int bbb = ((b-blockb)<<3);
               word = LG[0][id1][b] & LG[0][id2][b];   // HomHom
               genobit[0][bbb] += (word&0x0101010101010101);
               for(int k = 1; k < 8; k++)
                  genobit[0][bbb|k] += ((word>>k)&0x0101010101010101);
               word &= (LG[1][id1][b]^LG[1][id2][b]);   // IBS0
               genobit[1][bbb] += (word&0x0101010101010101);
               for(int k = 1; k < 8; k++)
                  genobit[1][bbb|k] += ((word>>k)&0x0101010101010101);
               word = (LG[0][id1][b] | LG[1][id1][b]) & (LG[0][id2][b] | LG[1][id2][b]); // nonmissing
               genobit[2][bbb] += (word&0x0101010101010101);
               for(int k = 1; k < 8; k++)
                  genobit[2][bbb|k] += ((word>>k)&0x0101010101010101);
            }  // end of blocked marker
         }  // end of blocked individual
         for(int b = blockb; b < bMax; b++){
            int bbb = ((b-blockb)<<3);
            for(int k = 0; k < 8; k++){   // shift
               int pos = (b<<6) | k;
               for(int bit = 0; bit < 3; bit++){
                  pchar = (unsigned char *)&genobit[bit][bbb|k];
                  for(int j = 0; j < 8; j++) // bit parallel
                     genobitCount[bit][thread][pos | (j<<3)] += pchar[j];
               }  // end of bit
            }  // end of 8 shifts
         }  // end of b
      }  // end of blockb
   }  // end of bi
#ifdef _OPENMP
}
#endif
   for(int i = 0; i < defaultMaxCoreCount; i++)
      if(genobitCount[0][i]){  // thread effectively used
         for(int m = 0; m < markerCount; m++){
            HomHomCounts[m] += genobitCount[0][i][m];
            ibs0Counts[m] += genobitCount[1][i][m];
            nonmissingCounts[m] += genobitCount[2][i][m];
         }
         for(int k = 0; k < 3; k++)
            delete []genobitCount[k][i];
      }
   for(int k = 0; k < 3; k++)
      delete []genobitCount[k];
}

void KingEngine::ComputeAlleleFrequency64Bit(IntArray &AACounts, IntArray &AaCounts, IntArray &missingCounts, int wordcount)
{
   const int BLOCKSIZE=255;
   const int CACHESIZE=16; // total cache size: 2^(6 + 4 - 2 + 8) = 2^16 = 64KB
   const int CACHESIZE_BITPAR=CACHESIZE<<3;  // 512
   unsigned long long int word, genobit[2][CACHESIZE_BITPAR];
   unsigned char *pchar;
   unsigned char byte;

   char revbase[256];
   for(int i = 0; i < 8; i++)
      revbase[base[i]] = i;
   char rightmost[256];
   for(int i = 0; i < 256; i++)
      rightmost[i] = revbase[i&(-i)];
   int thread = 0;

   int **genobitCount[3];
   for(int k = 0; k < 3; k++){
      genobitCount[k] = new int *[defaultMaxCoreCount];
      for(int i = 0; i < defaultMaxCoreCount; i++)
         genobitCount[k][i] = NULL;
   }

   int mylongCount=(wordcount && (wordcount<=longCount))?wordcount:longCount;
   int mymarkerCount=(mylongCount==longCount)? markerCount:(mylongCount<<6);
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(word, genobit, pchar, byte, thread)
{
   thread = omp_get_thread_num();
   for(int k = 0; k < 3; k++){
      genobitCount[k][thread] = new int [(mylongCount<<6)];
      for(int m = 0; m < mymarkerCount; m++)
         genobitCount[k][thread][m] = 0;
   }
   #pragma omp for
#else
   for(int k = 0; k < 3; k++){
      genobitCount[k][thread] = new int [mylongCount<<6];
      for(int m = 0; m < mymarkerCount; m++)
         genobitCount[k][thread][m] = 0;
   }
#endif
   for(int bi = 0; bi < idCount; bi+=BLOCKSIZE){   // blocked individuals
      int iMax = (bi > idCount-BLOCKSIZE) ? idCount: bi+BLOCKSIZE;
      for(int blockb = 0; blockb < mylongCount; blockb += CACHESIZE){  // blocked SNP words
         int bMax = (blockb >= mylongCount-CACHESIZE) ? mylongCount: blockb+CACHESIZE;
         for(int j = 0; j < CACHESIZE_BITPAR; j++)
            genobit[0][j] = genobit[1][j] = 0;
         for(int i = bi; i < iMax; i++){
            for(int b = blockb; b < bMax; b++){
               int bbb = ((b-blockb)<<3);
               for(int j = 0; j < 2; j++){
                  word = LG[j][i][b];  // jth bit word
                  genobit[j][bbb] += (word&0x0101010101010101);
                  for(int k = 1; k < 8; k++)
                     genobit[j][bbb|k] += ((word>>k)&0x0101010101010101);
               }  // end of bit
               word = ~(LG[0][i][b] | LG[1][i][b]);  // Miss
               pchar = (unsigned char*)&word;
               for(int k = 0; k < 8; k++)
                  for(byte = pchar[k]; byte; byte &= (byte-1))
                     genobitCount[2][thread][(b<<6)|(k<<3)|rightmost[byte]] ++;
            }  // end of blocked marker
         }  // end of blocked individual
         for(int b = blockb; b < bMax; b++){
            int bbb = ((b-blockb)<<3);
            for(int k = 0; k < 8; k++){   // shift
               int pos = (b<<6) | k;
               for(int bit = 0; bit < 2; bit++){
                  pchar = (unsigned char *)&genobit[bit][bbb|k];
                  for(int j = 0; j < 8; j++) // bit parallel
                     genobitCount[bit][thread][pos | (j<<3)] += pchar[j];
               }
            }  // end of 8 shifts
         }  // end of b
      }  // end of blockb
   }  // end of bi
#ifdef _OPENMP
}
#endif
   AACounts.Dimension(mymarkerCount);
   AACounts.Zero();
   AaCounts.Dimension(mymarkerCount);
   AaCounts.Zero();
   missingCounts.Dimension(mymarkerCount);
   missingCounts.Zero();
   for(int i = 0; i < defaultMaxCoreCount; i++)
      if(genobitCount[0][i]){  // thread effectively used
         for(int m = 0; m < mymarkerCount; m++){
            AACounts[m] += genobitCount[0][i][m] + genobitCount[1][i][m] + genobitCount[2][i][m];
            AaCounts[m] += genobitCount[1][i][m];
            missingCounts[m] += genobitCount[2][i][m];
         }
         for(int k = 0; k < 3; k++)
            delete []genobitCount[k][i];
      }
   for(int m = 0; m < mymarkerCount; m++){
      AACounts[m] -= idCount;
      AaCounts[m] -= AACounts[m];
   }
   for(int k = 0; k < 3; k++)
      delete []genobitCount[k];
}

void KingEngine::GenotypeCount64Bit(char genobit, int wordstart, int wordstop, int *genocount, unsigned long long int *words)
{
   unsigned long long int word;
   int wordslength = (wordstop-wordstart)<<3;
   int genocountlength = wordslength<<3;
   for(int m = 0; m < genocountlength; m++)
      genocount[m] = 0;
   unsigned char *pchar;
   const unsigned char BLOCKSIZE=255;
   for(int bi = 0; bi < idCount; bi += BLOCKSIZE){
      int iMax = (bi > idCount-BLOCKSIZE) ? idCount: bi+BLOCKSIZE;
      for(int m = 0; m < wordslength; m++)
         words[m] = 0;
      for(int m = wordstart; m < wordstop; m++){
         int mm = (m-wordstart)<<3;
         for(int id = bi; id < iMax; id++){
            word = LG[genobit][id][m];
            words[mm] += word&0x0101010101010101;
            for(int k = 1; k < 8; k++) // N*M*8: most computationally intensive
               words[mm|k] += (word>>k)&0x0101010101010101;
         }  // end of a block of individuals
      }  // end of m
      for(int m = wordstart; m < wordstop; m++){
         int mm = (m-wordstart)<<3;
         for(int k = 0; k < 8; k++){   // shift
            pchar = (unsigned char*)&words[mm|k];
            for(int j = 0; j < 8; j++) // bit parallel
               genocount[(mm<<3) | k | (j<<3)] += pchar[j];
         }  // end of 8 shifts
      }  // end of m
   }
}

void KingEngine::ComputeAlleleFrequency(Vector & frequencies)
{
   const int CACHESIZE = 1024;
   const int CACHESIZE_SNP = CACHESIZE*16;
   double freq[CACHESIZE_SNP];
   int missingCount[CACHESIZE_SNP];
   frequencies.Dimension(markerCount);
   unsigned short int word;
   char revbase[65536];
   for(int i = 0; i < 16; i++)
      revbase[shortbase[i]] = i;
   char rightmost[65536];
   for(int i = 0; i < 65536; i++)
      rightmost[i] = revbase[i&(-i)];
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
      private(freq, missingCount, word)
#endif
   for(int blockb = 0; blockb < shortCount; blockb += CACHESIZE){
      int bMax = (blockb > shortCount-CACHESIZE) ? shortCount: blockb+CACHESIZE;
      int bbMax = bMax - blockb;
      int mMax = (blockb > shortCount-CACHESIZE)? markerCount-blockb*16: (bMax-blockb)*16;
      for(int j = 0; j < CACHESIZE_SNP; j++){
         freq[j] = 0.0;
         missingCount[j] = 0;
      }
      for(int id = 0; id < idCount; id++)
         for(int b = blockb; b < bMax; b++){
            int bb = (b-blockb)*16;
            for(word = GG[0][id][b] & GG[1][id][b]; word; word &= (word-1))
               freq[bb+rightmost[word]] += 1;   // AA
            for(word = (~GG[0][id][b]) & GG[1][id][b]; word; word &= (word-1))
               freq[bb+rightmost[word]] += 0.5; // Aa
            for(word = (~GG[0][id][b]) & (~GG[1][id][b]) & 65535; word; word &= (word-1))
               missingCount[bb+rightmost[word]] ++;   // missing
         }  // end of b
      for(int m = 0; m < mMax; m++)
         freq[m] /= (idCount - missingCount[m]);
      for(int m = 0; m < mMax; m++)
         frequencies[blockb*16+m] = freq[m];
   }                                       
}

void KingEngine::ComputeAlleleFrequencyInX(Vector & frequencies)
{
   const int CACHESIZE = 1024;
   const int CACHESIZE_SNP = CACHESIZE*16;
   double freq[CACHESIZE_SNP];
   int missingCount[CACHESIZE_SNP];
   frequencies.Dimension(xmarkerCount);
   unsigned short int word;
   char revbase[65536];
   for(int i = 0; i < 16; i++)
      revbase[shortbase[i]] = i;
   char rightmost[65536];
   for(int i = 0; i < 65536; i++)
      rightmost[i] = revbase[i&(-i)];
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
      private(freq, missingCount, word)
#endif
   for(int blockb = 0; blockb < xshortCount; blockb += CACHESIZE){
      int bMax = (blockb > xshortCount-CACHESIZE) ? xshortCount: blockb+CACHESIZE;
      int bbMax = bMax - blockb;
      int mMax = (blockb > xshortCount-CACHESIZE)? xmarkerCount-blockb*16: (bMax-blockb)*16;
      for(int j = 0; j < CACHESIZE_SNP; j++){
         freq[j] = 0.0;
         missingCount[j] = 0;
      }
      for(int id = 0; id < idCount; id++)
         for(int b = blockb; b < bMax; b++){
            int bb = (b-blockb)*16;
            for(word = XG[0][id][b] & XG[1][id][b]; word; word &= (word-1))
               freq[bb+rightmost[word]] += 1;   // AA
            for(word = (~XG[0][id][b]) & XG[1][id][b]; word; word &= (word-1))
               freq[bb+rightmost[word]] += 0.5; // Aa
            for(word = (~XG[0][id][b]) & (~XG[1][id][b]) & 65535; word; word &= (word-1))
               missingCount[bb+rightmost[word]] ++;   // missing
         }  // end of b
      for(int m = 0; m < mMax; m++)
         freq[m] /= (idCount - missingCount[m]);
      for(int m = 0; m < mMax; m++)
         frequencies[blockb*16+m] = freq[m];
   }
}


void KingEngine::runKING()
{
   if(geno.Length()==0) {
      individualInfo = true;
      if(shortFlag)
         BuildShortBinary();
      else
         BuildBinary();
   }
   printf("Genotypes stored in %d integers for each of %d individuals.\n",
      shortCount, idCount);

   double smaller;
   int id1, id2;
   double kinship;
   int HetHetCount, IBS0Count, het1Count, het2Count, notMissingCount;
   int deg;

   char oneoneCount[65536];
   if(shortFlag)
      for(int i = 0; i < 65536; i++)
         oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   L0.Dimension(0);
   Lfs.Dimension(0);
   Lpo.Dimension(0);
   L2.Dimension(0);
   Vector IBS0L1(0);
   double ibs0;
   double tempErrorRate = (errorrateCutoff==_NAN_? 0.008: errorrateCutoff);

   for(int i = 0; i < ped.count; i++)
      for(int j = i+1; j < ped.count; j++){
         if(ped[i].famid == ped[j].famid) continue;
         id1 = geno[i]; id2 = geno[j];
         if(id1 < 0 || id2 < 0) continue;
         HetHetCount = IBS0Count = het1Count = het2Count = 0;
         for(int m = 0; m < shortCount; m++){
            HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
            IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
         }
         if(HetHetCount <= 2*IBS0Count) continue;
         for(int m = 0; m < shortCount; m++){
            het1Count += oneoneCount[(GG[0][id2][m] | GG[1][id2][m]) & (~GG[0][id1][m]) & GG[1][id1][m]];
            het2Count += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
         }
         smaller = het1Count < het2Count? het1Count: het2Count;
         kinship = 0.5 - (het1Count+het2Count)*0.25/smaller + (HetHetCount*0.5-IBS0Count)/smaller;
         if(kinship <= 0.044194) continue;
         if(markerCount < 6000 && kinship <= 0.088388) continue;
         deg = int(-log(kinship)/0.6931472-0.5);
         if(deg == 0) {L0.Push(i); L0.Push(j);}
         else if(deg == 1){
            notMissingCount = 0;
            for(int m = 0; m < shortCount; m++)
               notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
            ibs0 = IBS0Count*1.0/notMissingCount;
            IBS0L1.Push(ibs0);
            if(ibs0 > tempErrorRate) {Lfs.Push(i); Lfs.Push(j); }
            else if(ibs0 > 0.005){
               if(kinship > 0.24){  // FS
                  Lfs.Push(i); Lfs.Push(j);
               }else{   // PO
                  Lpo.Push(i); Lpo.Push(j);
               }
            }else{   // PO
               Lpo.Push(i); Lpo.Push(j);
            }
         }else if(deg == 2){L2.Push(i); L2.Push(j);}
      }

      if(errorrateCutoff == _NAN_){
         IntArray L1Count(10);
         L1Count.Zero();
         for(int i = 0; i < IBS0L1.Length(); i++)
            if(IBS0L1[i] < 0.01)
               L1Count[int(IBS0L1[i]*1000)] ++;
         int T = 9;
         for(; L1Count[T] && T >= 0; T--);
         if(T < 0)  // find the one with minimum counts
            for(T = 9; L1Count[T] > L1Count.Min(); T--);
         if(T > 0 && int(errorrateCutoff*1000)!=T) {
            errorrateCutoff = (T+0.5)*0.001;
            printf("Cutoff value between full siblings and parent-offspring is set at %.4f\n",
               errorrateCutoff);
         }
      }
}

void KingEngine::WritePlinkBinary(const char *prefix)
{
   IntArray removeflags(ped.count);
   removeflags.Set(0);
   int extraCount = 0;
   int removeCount = 0;
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
         if(geno[i] > -1){ // genotype data stored for the ith person
            if(ped[i].ngeno==0 && ped[i].isFounder()){// founder's genotype stored but all missing
               removeflags[i] = 1;
               removeCount ++;
            }else
               continue;
         }
         // untyped
         if(ped[i].isFounder())
            removeflags[i] = 1;
         else
            extraCount ++;
      }
   if(extraCount)
      printf("%d individuals with missing genotypes will be kept to connect pedigrees.\n",
         extraCount);

   String temp;
   String pedfile = prefix;
   pedfile.Add(".fam");
   FILE *fp = fopen(pedfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = 0; i < id[f].Length(); i++){
            if(removeflags[id[f][i]]) continue;
            fprintf(fp, "%s %s %s %s %d %d\n",
               (const char*)ped[id[f][i]].famid,
               (const char*)ped[id[f][i]].pid,
               (ped[id[f][i]].fatid.SubStr(0, 4) == "KING" && ped[id[f][i]].father && ped[id[f][i]].father->ngeno==0)? "0": (const char*)ped[id[f][i]].fatid,
               (ped[id[f][i]].motid.SubStr(0, 4) == "KING" && ped[id[f][i]].mother && ped[id[f][i]].mother->ngeno==0)? "0": (const char*)ped[id[f][i]].motid,
               ped[id[f][i]].sex,
               (ped.affectionCount && ped[id[f][i]].affections[0]> 0)?
                  ped[id[f][i]].affections[0]: -9);
   }
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(!removeflags[i] && geno[i] == -1)
            fprintf(fp, "%s %s %s %s %d %d\n",
               (const char*)ped[i].famid,
               (const char*)ped[i].pid,
               (ped[i].fatid.SubStr(0, 4) == "KING" && ped[i].father && ped[i].father->ngeno==0)? "0": (const char*)ped[i].fatid,
               (ped[i].motid.SubStr(0, 4) == "KING" && ped[i].mother && ped[i].mother->ngeno==0)? "0": (const char*)ped[i].motid,
               ped[i].sex,
               (ped.affectionCount && ped[i].affections[0]> 0)?
                  ped[i].affections[0]: -9);
   int inclusionCount = inclusionList[1].Length();
   for(int i = 0; i < inclusionCount; i++)
      if(inclusionList[3][i]!="" || inclusionList[4][i]!=""){
         fprintf(fp, "%s %s %s %s %d %d\n",
            (const char*)inclusionList[0][i],
            (const char*)inclusionList[1][i],
            inclusionList[3][i]!=""? (const char*)inclusionList[3][i]: 0,
            inclusionList[4][i]!=""? (const char*)inclusionList[4][i]: 0,
            (int)inclusionList[2][i], -9);
         extraCount ++;
      }
   fclose(fp);
   printf("  Family file saved in %s\n", (const char*)pedfile);

   String covfile = prefix;
   covfile.Add(".cov");
   if(ped.covariateCount || ped.haveTwins){
      bool AddTwinFlag = false;
      if(ped.haveTwins && ped.covariateNames.SlowFind("MZTWIN")==-1)
         AddTwinFlag = true;
      fp = fopen(covfile, "wt");
      if(fp == NULL) error("Cannot open %s to write.", (const char*)covfile);
      fprintf(fp, "FID IID");
      for(int i = 0; i < ped.covariateCount; i++)
         fprintf(fp, " %s", (const char*)ped.covariateNames[i]);
      if(AddTwinFlag)
         fprintf(fp, " MZTWIN");
      fprintf(fp, "\n");
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
            if(removeflags[i]) continue;
            fprintf(fp, "%s %s",
               (const char*)ped[i].famid,
               (const char*)ped[i].pid);
            for(int t = 0; t < ped.covariateCount; t++)
               if(ped[i].covariates[t] == _NAN_)
                  fprintf(fp, " -9");
               else
                  fprintf(fp, " %G", ped[i].covariates[t]);
            if(AddTwinFlag)
               fprintf(fp, " %d", ped[i].zygosity);
            fprintf(fp, "\n");
         }
      fclose(fp);
   }
   String phefile = prefix;
   phefile.Add(".phe");
   if(ped.traitCount || (ped.affectionCount > 1) ){
      int headerFlag = 0;
      fp = fopen(phefile, "wt");
      if(fp == NULL) error("Cannot open %s to write.", (const char*)phefile);
      fprintf(fp, "FID IID");
      for(int i = 0; i < ped.traitCount; i++)
         fprintf(fp, " %s", (const char*)ped.traitNames[i]);
      if(ped.affectionCount){
         if(ped.affectionNames[0] == "DISEASEKING")
            headerFlag = 1;
         else
            fprintf(fp, " DISEASEKING");
      }
      for(int i = 1; i < ped.affectionCount; i++)
         fprintf(fp, " %s", (const char*)ped.affectionNames[i]);
      fprintf(fp, "\n");
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
            if(removeflags[i]) continue;
            fprintf(fp, "%s %s",
               (const char*)ped[i].famid,
               (const char*)ped[i].pid);
            for(int t = 0; t < ped.traitCount; t++)
               if(ped[i].traits[t] == _NAN_)
                  fprintf(fp, " -9");
               else
                  fprintf(fp, " %G", ped[i].traits[t]);
            for(int t = headerFlag; t < ped.affectionCount; t++)
               if(ped[i].affections[t] == _NAN_)
                  fprintf(fp, " 0");
               else
                  fprintf(fp, " %d", ped[i].affections[t]);
            fprintf(fp, "\n");
         }
      fclose(fp);
   }

   String mapfile = prefix;
   mapfile.Add(".bim");
   fp = fopen(mapfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)mapfile);
   if(chromosomes.Length()){
      for(int m = 0; m < markerCount; m++)
         fprintf(fp, "%d\t%s\t0\t%d\t%s\t%s\n",
            chromosomes[m],
            (const char*)snpName[m],
            bp[m], //int(positions[m] * 1000000+0.5),
            (const char*)alleleLabel[0][m],
            (const char*)alleleLabel[1][m]);
   }else
      for(int m = 0; m < markerCount; m++)
         fprintf(fp, "0\tSNP%d\t0\t0\t1\t2\n", m+1);
   if(xbp.Length())
      for(int m = 0; m < xmarkerCount; m++)
         fprintf(fp, "%d\t%s\t0\t%d\t%s\t%s\n",
            SEXCHR,
            (const char*)xsnpName[m],
            xbp[m], //int(xpositions[m] * 1000000+0.5),
            (const char*)xalleleLabel[0][m],
            (const char*)xalleleLabel[1][m]);
   else
      for(int m = 0; m < xmarkerCount; m++)
         fprintf(fp, "%d\tXSNP%d\t0\t0\t1\t2\n", SEXCHR, m+1);

   if(ybp.Length())
      for(int m = 0; m < ymarkerCount; m++)
         fprintf(fp, "%d\t%s\t0\t%d\t%s\t%s\n",
            SEXCHR+1, (const char*)ysnpName[m],
            ybp[m], //int(ypositions[m] * 1000000 + 0.5),
            (const char*)yalleleLabel[0][m],
            (const char*)yalleleLabel[1][m]);
   else
      for(int m = 0; m < ymarkerCount; m++)
         fprintf(fp, "%d\tYSNP%d\t0\t0\t1\t2\n", SEXCHR+1, m+1);

   if(mtbp.Length())
      for(int m = 0; m < mtmarkerCount; m++)
         fprintf(fp, "%d\t%s\t0\t%d\t%s\t%s\n",
            SEXCHR+3, (const char*)mtsnpName[m],
            mtbp[m], //int(mtpositions[m] * 1000000 + 0.5),
            (const char*)mtalleleLabel[0][m],
            (const char*)mtalleleLabel[1][m]);
   else
      for(int m = 0; m < mtmarkerCount; m++)
         fprintf(fp, "%d\tMTSNP%d\t0\t0\t1\t2\n", SEXCHR+3, m+1);
   fclose(fp);
   printf("  Map file saved in %s\n", (const char*)mapfile);

   String bedfile = prefix;
   bedfile.Add(".bed");
   fp = fopen(bedfile, "wb");
   if(fp == NULL)
      error("Cannot open %s to write.", (const char*)bedfile);
   unsigned char alocus[60000];
   int byte, bit;
   unsigned char g;
   char buffer[0x20000];
   int thread = 0;
   int pbuffer = 0;
   pbuffer += sprintf(&buffer[pbuffer], "%c%c%c", 108, 27, 1);
   if(Bit64==64){
      for(int m = 0; m < markerCount; m++){
         int b = m / 64;
         int k = m % 64;
         for(int i = 0; i < (idCount-1-removeCount+extraCount)/4+1; i++)
            alocus[i] = 0;
         byte = bit = 0;
         for(int f = 0; f < ped.familyCount; f++)
            for(int i = 0; i < id[f].Length(); i++){
               if(removeflags[id[f][i]]) continue;
               unsigned long long int base = (unsigned long long int)1<<k;
               if(LG[0][geno[id[f][i]]][b] & base) // homozygote
                  g = (LG[1][geno[id[f][i]]][b] & base) ? 0: 3;
               else  // heterozygote
                  g = (LG[1][geno[id[f][i]]][b] & base)? 2: 1;
               if(g) alocus[byte] |= (g << bit);
               if(bit == 6){
                  bit = 0;
                  byte ++;
               }else
                  bit += 2;
            }
         if(extraCount){
            for(int f = 0; f < ped.familyCount; f++)
               for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
                  if(!removeflags[i] && geno[i] == -1){ // untyped but cannot be removed
                     alocus[byte] |= (1 << bit);
                     if(bit == 6){
                        bit = 0;
                        byte ++;
                     }else
                        bit += 2;
                  }
            for(int i = 0; i < inclusionCount; i++){
               if(inclusionList[3][i]!="" || inclusionList[4][i]!=""){
                  alocus[byte] |= (1 << bit);
                  if(bit == 6){
                     bit = 0;
                     byte ++;
                  }else
                     bit += 2;
               }
            }
         }
         for(int i = 0; i < (idCount-removeCount+extraCount-1)/4+1; i++){
            pbuffer += sprintf(&buffer[pbuffer], "%c", alocus[i]);
            if(pbuffer > 0xFFFF){  // buffer big enough for writing
               fwrite(buffer, 1, pbuffer, fp);
               pbuffer = 0;
            }
         }
      }
   }else{  // short genotypes for 32 bit system
      for(int m = 0; m < markerCount; m++){
         int b = m / 16;
         int k = m % 16;
         for(int i = 0; i < (idCount-1-removeCount+extraCount)/4+1; i++)
            alocus[i] = 0;
         byte = bit = 0;
         for(int f = 0; f < ped.familyCount; f++)
            for(int i = 0; i < id[f].Length(); i++){
               if(removeflags[id[f][i]]) continue;
               if(GG[0][geno[id[f][i]]][b] & shortbase[k]) // homozygote
                  g = (GG[1][geno[id[f][i]]][b] & shortbase[k]) ? 0: 3;
               else  // heterozygote
                  g = (GG[1][geno[id[f][i]]][b] & shortbase[k])? 2: 1;
               if(g) alocus[byte] |= (g << bit);
               if(bit == 6){
                  bit = 0;
                  byte ++;
               }else
                  bit += 2;
            }
         if(extraCount)
            for(int f = 0; f < ped.familyCount; f++)
               for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
                  if(!removeflags[i] && geno[i] == -1){ // untyped but cannot be removed
                     alocus[byte] |= (1 << bit);
                     if(bit == 6){
                        bit = 0;
                        byte ++;
                     }else
                        bit += 2;
                  }
         for(int i = 0; i < (idCount-removeCount+extraCount-1)/4+1; i++){
            pbuffer += sprintf(&buffer[pbuffer], "%c", alocus[i]);
            if(pbuffer > 0xFFFF){  // buffer big enough for writing
               fwrite(buffer, 1, pbuffer, fp);
               pbuffer = 0;
            }
         }
      }
   }
   for(int m = 0; m < xmarkerCount; m++){
      int b = m / 16;
      int k = m % 16;
      for(int i = 0; i < (idCount-removeCount+extraCount-1)/4+1; i++)
         alocus[i] = 0;
      byte = bit = 0;
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = 0; i < id[f].Length(); i++){
            if(removeflags[id[f][i]]) continue;
            if(XG[0][geno[id[f][i]]][b] & shortbase[k])
               g = (XG[1][geno[id[f][i]]][b] & shortbase[k]) ? 0: 3;
            else
               g = (XG[1][geno[id[f][i]]][b] & shortbase[k])? 2: 1;
            if(g) alocus[byte] |= (g << bit);
            if(bit == 6){
               bit = 0;
               byte ++;
            }else
               bit += 2;
         }
      if(extraCount)
         for(int f = 0; f < ped.familyCount; f++)
            for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
               if(!removeflags[i] && geno[i] == -1){
                  alocus[byte] |= (1 << bit);
                  if(bit == 6){
                     bit = 0;
                     byte ++;
                  }else
                     bit += 2;
               }
      for(int i = 0; i < (idCount-removeCount-1+extraCount)/4+1; i++){
         pbuffer += sprintf(&buffer[pbuffer], "%c", alocus[i]);
         if(pbuffer > 0xFFFF){  // buffer big enough for writing
            fwrite(buffer, 1, pbuffer, fp);
            pbuffer = 0;
         }
      }
   }
   for(int m = 0; m < ymarkerCount; m++){
      int b = m / 16;
      int k = m % 16;
      for(int i = 0; i < (idCount-removeCount-1+extraCount)/4+1; i++)
         alocus[i] = 0;
      byte = bit = 0;
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = 0; i < id[f].Length(); i++){
            if(removeflags[id[f][i]]) continue;
            if(YG[0][geno[id[f][i]]][b] & shortbase[k])
               g = (YG[1][geno[id[f][i]]][b] & shortbase[k]) ? 0: 3;
            else
               g = (YG[1][geno[id[f][i]]][b] & shortbase[k])? 2: 1;
            if(g) alocus[byte] |= (g << bit);
            if(bit == 6){
               bit = 0;
               byte ++;
            }else
               bit += 2;
         }
      if(extraCount)
         for(int f = 0; f < ped.familyCount; f++)
            for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
               if(!removeflags[i] && geno[i] == -1){
                  alocus[byte] |= (1 << bit);
                  if(bit == 6){
                     bit = 0;
                     byte ++;
                  }else
                     bit += 2;
               }
      for(int i = 0; i < (idCount-removeCount-1+extraCount)/4+1; i++){
         pbuffer += sprintf(&buffer[pbuffer], "%c", alocus[i]);
            if(pbuffer > 0xFFFF){  // buffer big enough for writing
               fwrite(buffer, 1, pbuffer, fp);
               pbuffer = 0;
            }
      }
   }
   for(int m = 0; m < mtmarkerCount; m++){
      int b = m / 16;
      int k = m % 16;
      for(int i = 0; i < (idCount-removeCount-1+extraCount)/4+1; i++)
         alocus[i] = 0;
      byte = bit = 0;
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = 0; i < id[f].Length(); i++){
            if(removeflags[id[f][i]]) continue;
            if(MG[0][geno[id[f][i]]][b] & shortbase[k])
               g = (MG[1][geno[id[f][i]]][b] & shortbase[k]) ? 0: 3;
            else
               g = (MG[1][geno[id[f][i]]][b] & shortbase[k])? 2: 1;
            if(g) alocus[byte] |= (g << bit);
            if(bit == 6){
               bit = 0;
               byte ++;
            }else
               bit += 2;
         }
      if(extraCount)
         for(int f = 0; f < ped.familyCount; f++)
            for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
               if(!removeflags[i] && geno[i] == -1){
                  alocus[byte] |= (1 << bit);
                  if(bit == 6){
                     bit = 0;
                     byte ++;
                  }else
                     bit += 2;
               }
      for(int i = 0; i < (idCount-removeCount-1+extraCount)/4+1; i++){
         pbuffer += sprintf(&buffer[pbuffer], "%c", alocus[i]);
         if(pbuffer > 0xFFFF){  // buffer big enough for writing
            fwrite(buffer, 1, pbuffer, fp);
            pbuffer = 0;
         }
      }
   }
   if(pbuffer>0)
      fwrite(buffer, 1, pbuffer, fp);
   fclose(fp);
   printf("  Binary genotypes saved in %s\n", (const char*)bedfile);
//   printf("Data saved as PLINK binary format in files %s, %s, %s",
//      (const char*)pedfile, (const char*)mapfile, (const char*)bedfile);
//   if(ped.covariateCount || ped.haveTwins) printf(", %s", (const char*)covfile);
//   if(ped.traitCount) printf(", %s", (const char*)phefile);
//   printf("\n");
}




