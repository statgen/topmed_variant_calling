#ifndef __KingCore_h__
#define __KingCore_h__

#include "Pedigree.h"
#include "MathMatrix.h"
#include "StringArray.h"
#define MINSNPCOUNT 512
#include "QuickIndex.h"

#define RelatedFLAG 0
#define KinshipFLAG 1
#define AutoQCFLAG 2
#define MtscoreFLAG 3
#define RiskFLAG 4
#define IBSFLAG 5
#define HomoFLAG 6
#define IBDSegFLAG 7
#define MDSFLAG 8
#define PCAFLAG 9
#define ClusterFLAG 10
#define BuildFLAG 11
#define BysampleFLAG 12
#define BysnpFLAG 13
#define TDTFLAG 14
#define UnrelatedFLAG 15
#define DuplicateFLAG 16
#define ROHFLAG 17
#define TOTALFLAGCOUNT 18

#define PLINKFLAG 31

class KingEngine{
   protected:
      unsigned char **G_SNP;
      unsigned short int **GG[2];
      unsigned short int **SG[2];
      unsigned long long int **LG[2];
      unsigned long long int **SLG[2];
      unsigned char base[8];
      unsigned short int shortbase[16];
      unsigned short int *shortFlip;
      char oneCount[256];
      IntArray *id;
      int idCount;
      IntArray geno;
      IntArray phenoid;
      int shortCount;
      int longCount;
      int markerCount;
      StringArray sampleName;
      StringArray traitNames;
      StringArray affectionNames;
      StringArray covariateNames;
      StringArray FID, PID, FA, MO, SEX;
      StringArray ZG;
      Matrix traits;
      Matrix affections;
      Matrix covariates;

      StringArray snpName;
      IntArray chromosomes;
      Vector positions;
      StringArray alleleLabel[2];
      IntArray sortedSNP;

      unsigned short int **XG[2];
      StringArray xalleleLabel[2];
      Vector xpositions;
      StringArray xsnpName;
      int xshortCount;
      int xmarkerCount;

      int ymarkerCount;
      int yshortCount;
      StringArray ysnpName;
      unsigned short int **YG[2];
      StringArray yalleleLabel[2];
      Vector ypositions;

      int mtmarkerCount;
      int mtshortCount;
      StringArray mtsnpName;
      unsigned short int **MG[2];
      StringArray mtalleleLabel[2];
      Vector mtpositions;

      void MakePed(void);
      void autoflip(void);
//      void FlipInDuplicate(void);
      void AfterBinaryLoaded(void);
   public:
      unsigned int allflags;
      StringArray exclusionList;
      StringArray inclusionList[3];
      bool diagFlag;
      bool detailFlag;
      bool shortFlag;
      bool runlengthFlag;
      bool lessmemFlag;
      String prefix;
      String plinkinput;
      Pedigree & ped;
      bool mafFilterFlag;
      bool genoFilterFlag;
      double minMAF;
      double callrate;
      bool freqLoaded;
      bool individualInfo;
      bool QCwipe;
      IntArray L0, Lpo, Lfs, L2/*, L3*/, Ltrio;
      double errorrateCutoff;
      bool autoflipFlag;
      bool flipFlag;
      String SaveFormat;
      bool familyFlag;
      bool pedstatFlag;
      bool bigdataFlag;
//      Vector HetBySample;
      QuickIndex bigdataIdx;
      int veryrareCount;
      bool autosomeOnly;
      bool genotypeOnly;
      int SEXCHR;

      IntArray chrList;
      double start, stop;
      int CoreCount;
      int defaultMaxCoreCount;
      int defaultEfficientCoreCount;

      char Bit64;
      bool Bit64Flag;

      KingEngine(Pedigree &);
      ~KingEngine();

      void BuildBinary();
//      void ReadBinaryData(const char *);
      void runKING();
      void WriteKingBinary(const char*);
      void WritePlinkBinary(const char*);
      void DumpBinary(const char*);

      void BuildShortBinary();
//      void ReadBinaryShortData(const char *);
//      void ReadPlinkBinaryData(const char *);
      void ReadKingBinaryData(const char *);

      void ConvertGGtoSG(unsigned short int **oldG[2], int oldCount, unsigned short int **newG[2], int newCount);
      void ConvertLGtoSLG(unsigned long long int **oldG[2], int oldCount, unsigned long long int **newG[2], int newCount);
      void ConvertMajorFromSNPtoIndividual();
      void ConvertMajorFromSNPtoIndividual64Bit();

      void ReadPlinkBinaryBigData(const char *filename, const char *famfile=NULL, const char *bimfile=NULL, const char *covfile=NULL, const char *phefile=NULL);
      void ReadPlinkBinaryBigDataSNPMajor(const char *, const char*, const char *famfile=NULL, const char *bimfile=NULL, const char *covfile=NULL, const char *phefile=NULL);
      void ReadPlinkBinaryBigDataWithLessMemory(const char *, const char *famfile=NULL, const char *bimfile=NULL, const char *covfile=NULL, const char *phefile=NULL);
      void ReadPlinkPedFile(const char *, const char *famfile=NULL);
      void ReadPlinkCovFile(const char *, const char *covfile=NULL);
      void ReadPlinkTraitFile(const char *, const char *phefile=NULL);
      void ReadPlinkMapFile(const char *, IntArray &, const char *bimfile=NULL);

      void ComputeAlleleFrequency64Bit(IntArray &AACounts, IntArray &AaCounts, IntArray &missingCounts, int wordcount=0);
      void ComputeMZBySNP64Bit(IntArray &nonmissingMZCounts, IntArray &ibs1MZCounts, IntArray &HetHetMZCounts, IntArray &ibs0MZCounts);
      void ComputePOBySNP64Bit(IntArray &HomHomCounts, IntArray &ibs0Counts, IntArray &nonmissingCounts);
      void ComputeTrioBySNP64Bit(IntArray & MITrioCounts, IntArray & HetInOffspringCounts, IntArray & nonmissingtrioCounts);
      void ComputeAlleleFrequency(Vector &frequencies);
      void ComputeAlleleFrequencyInX(Vector &frequencies);
      void GenotypeCount64Bit(char genobit, int wordstart, int wordstop, int *genocount, unsigned long long int *words);

      void ReadMultiplePlinkBinaryBigData(const char *);
      char *currentTime();
};

#endif


