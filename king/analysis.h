#ifndef __analysis_h__
#define __analysis_h__
#include "KingCore.h"

inline unsigned char popcount(unsigned long long int word)
{
   word = word - ((word>>1)&0x5555555555555555);
   word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
   word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
   word = (word+(word>>8)) & 0x00FF00FF00FF00FF;
   word = (word+(word>>16)) & 0x0000FFFF0000FFFF;
   return((word+(word>>32)) & 0xFF);
}

inline unsigned char popcount(unsigned short int word)
{
   word = word - ((word>>1)&0x5555);
   word = (word&0x3333) + ((word>>2)&0x3333);
   word = (word+(word>>4)) & 0x0F0F;
   return((word+(word>>8)) & 0xFF);
}

class Engine:public KingEngine{
   protected:
      void printRelationship(int *beforeCount, int *afterCount);
      int *pAACount, *pAaCount, *paaCount;
      int *pxAACount, *pxAaCount, *pxaaCount;
      int *pyAACount, *pyAaCount, *pyaaCount;
      int *pmtAACount, *pmtAaCount, *pmtaaCount;
      int *quality;
      Vector BetaSum;
      Matrix BetaSquareSum;
      Matrix freqBeta;

      IntArray covariatePC;
      int diagnosis(void);
      int qualityT;
      int cAge;
      int missingBase;
      bool uniqueIID;
      Matrix *pedKin;
      void SemifamilyKinship(void);
      double EmpP(double, Vector &);
      void ComputeDistanceMatrix(Matrix & Dist);
      void ComputeDistanceMatrix64Bit(Matrix & Dist);
      double ComputeAUC(Vector & risks, IntArray & diseaseStatus, int printFlag, int TOTALCUT);

   public:
      Engine(Pedigree &ped);
      ~Engine();

      int xflag;
// Relationship Inference
      int faster;
      int slower;
      double kinFilter;
      int relativedegree;
      int mincons;
      bool homogeneity;
      bool adjustFamily;
      void ComputeRunLength();
      void ComputeRunLengthInOneTrio(String &);
      void internalKING(int);
      void IBD2SegInOnePair64Bit(int id1, int id2, IntArray &chrSeg, double totalLength, double &prop, double &maxLength);
      bool PreSegment(bool printFlag=true);
      void ROH();
      void ROHOnly(IntArray & idList, int segment, IntArray & rohStorage, IntArray & rohIndex, bool LengthOnly=false);
      void PopulationROH();
      void ComputeIBD2Segment();
      void ComputeIBD2Segment64Bit();
      void ComputeIBDSegment64Bit();

      void adjustPC(int pcCount);

      void ComputeBigDataDuplicate();
      void ComputeBigDataDuplicate64Bit();
      void ComputeBigDataSecondDegree();
      void ComputeBigDataDistant();
      void ComputeBigDataPO();
      void ComputeBigDataKinship();
      void IntegratedRelationshipInference();
      void ComputeBigDataKinshipAdjustPC();

      int ScreenDuplicates(IntArray & rpList);
      int ScreenDuplicates64Bit(IntArray & rpList);
      long int ScreenCloseRelatives(IntArray & rpList);
      void ScreenCloseRelativesInSubset(IntArray & rpList);
      void ScreenCloseRelativesInSubset64Bit(IntArray & rpList);

//      void ComputeLongHomoKinship64Bit();
      void ComputeShortFastHomoKinship();
      void ComputeShortFastXHomoKinship();
      void ComputeShortExtendedIBS();
      void ComputeExtendedIBS64Bit();
      void ComputeShortSimpleIBS();
      void ComputeShortRobustKinship();
      void ComputeLongRobustKinship64Bit();
      void ComputeLongRobustKinship64BitWithFilter(IntArray & ids, bool WriteFlag=true);
      void ComputeShortRobustKinshipAdjustPC();
      void ComputeShortRobustXKinship();
      void ComputeShortSibIBS();
      void ComputeShortDuplicate();
      void KinshipInSubset(IntArray & pairList, IntArray & HetHetCounts,
         IntArray & IBS0Counts, IntArray & het1Counts, IntArray & het2Counts, IntArray & HomHomCounts);
      void KinshipInSubset64Bit(IntArray & pairList, IntArray & HetHetCounts,
         IntArray & IBS0Counts, IntArray & het1Counts, IntArray & het2Counts, IntArray & HomHomCounts, IntArray & IBS1Counts);
      void IBD2SegInSubset(IntArray & pairList, Vector & ibd2props, Vector & maxLengths);
      void IBD2SegInSubset64Bit(IntArray & pairList, Vector & ibd2props, Vector & maxLengths);
      void IBDSegInSubset64Bit(IntArray & pairList, Vector & ibdprops, Vector & maxLengths, Vector & ibd2props, Vector & maxLengths2);
      void IBDSegOnly(IntArray & pairList, int segment, IntArray & ibdsegStorage1, IntArray & ibdsegIndex1, IntArray & ibdsegStorage2, IntArray & ibdsegIndex2, bool LengthOnly=false, double MINSEGLENGTH=2.5, int MINCCOUNT=200);
      void NPL();
      void IBDmapping(int nperm);
      void HomozygosityMapping();
      void HomozygosityMappingMH(const char *popName="Reference");
      void HomozygosityMappingForQT();
      void HomozygosityMappingForQTMH(const char *popName="Reference");
      void IBDGDT();
      void IBDMI();
      void IBDVC();
      void AUCmapping();
      void AncestryInference();
      void PopulationDistance();
      void PopulationIBD();
      void AUCpredicting(IntArray &allchr, IntArray &allpos);
      void LocalH2();
      void ComputeLongRobustKinship();
      void ComputeFilteredRobustKinship();

// IBD Analysis
      IntArray chrSeg;
      double totalLength;
      String segmessage;

// Population Structure Analysis
      void pca(int every);
      void pca_family(int every);
      void pca_family_check(int every);
      void pca_projection(int every);
      void mds();
      void mds_moving();
      void mds_family();
      void mds_family_projection();
      void mds_kin_family();
      void mds_kin();

      void OneWindowProjection(int mv_chr, double mv_start, double mv_stop, Vector & ancestry);
      void SlidingWindows();
      Matrix localAncestry;

//      void mds_family_check();
      bool semifamilyFlag;
//      bool clusterFlag;
      bool projectFlag;
      bool unrelatedExtraction;

// Pedigree Reconstruction
      char specialChar;
      void rebuild();
      int BuildOneFamily(int f, IntArray, double, String &);
      int ClusterFamily(int pedrebuildFlag, int degree);
      void rebuild_semifamily();
      void IBDLength_3Pairs(IntArray &trio, int order, Vector &lengths);
      void IBDLength_trioMI(IntArray &trio, int order, Vector &lengths);

// QC
//      IntArray Ltrio;
      void MakeFamilyForMI(void);
      void countGenotype(void);
      void OutputIndividualInfo(void);
      void QC_By_SNP(void);
      void QC_By_SNP64Bit(void);
      void QC_By_Sample(void);
      void QC_By_Sample64Bit(void);
      void QC_WipeMI(void);

// autoQC & autoplots
      char *idMask;
      unsigned short int *gMask, *xMask, *yMask, *mtMask;
      void Monomorphic_SNP(IntArray & removeList);
      void xHeterozygosity_SNP(IntArray & removeList, double xHeterozygosity);
      void CallRate_SNP(double rateFilter, IntArray & removeList);
      void CallRate_xSNP(double rateFilter, IntArray & removeList);
      void CallRate_ySNP(double rateFilter, IntArray & removeList, bool lessthanFlag=true);
      void CallRate_Sample(double rateFilter, IntArray & removeList);
      void autoQC(double samplecallrate, double snpcallrate/*, double xHeterozygosity*/);
      void plotGenderError(IntArray & plotx, Vector & ploty, IntArray & plotz, double xHeterozygosity, int gendererrorCount);
      void plotRelationship();
      void plotIBDSeg();
      void plotIBD1vsIBD2(FILE *);
      void plotIBD2();
      void plotHetConcvsIBD2();
      void plotPopStructure();
      void plotAUCmapping();
      void plotNPL();
      void plotAncestry();
      void plotIBDmapping();
      void plotPopDist();
      void plotROHmapping(const char *stratName);
      void plotROHforQT();
      void plotPopROH();

// Tools
      int WriteFile;
//      void WritePeking();
      void WriteMerlin();
      void WritePlink();

// Binary Trait Association
      void ResidualQuantitative(int trait, int * vid, Vector & adjResid);
      void ResidualDisease(int * vid, Vector & adjResid);
      double ComputeDaviesP(double Q, Matrix & kernel);
      char KernelShort;
      void InvNorm(int tr);
      void SKAT();
      void SKAT_Batch_WeightedLinear(const char* skatfile);
      void SKAT_WeightedLinear();
      void SKAT_family_quantitative(int trait, StringArray &, IntArray *, Vector *, FILE *fp);
      void SKAT_family_quantitatives(const char* skatfile);
      void VCR_family_quantitative(int trait, StringArray &, IntArray *, Vector *, FILE *fp);
      void VCR_family_quantitatives(const char* vcrfile);
      void AdmixtureMapping(double, double);
      void MakeTrioForTDT(void);
      void TDT();
      void CountTDT(IntArray &trioList, IntArray &ncounts, IntArray & ntcounts);
      void CountTDT64Bit(IntArray &trioList, IntArray &ncounts, IntArray & ntcounts);
      void CountTDTinX(IntArray &trioList, IntArray &ncounts, IntArray & ntcounts);
      void rareTDT();
      void permuteTDT();
      int permuteCount;
      double OnePermuteTDT(IntArray &);
      void permuteRareTDT();
      double OnePermuteRareTDT(IntArray &, int);
      int TDT_T[2];
      int TDT_PeakPos;
      IntArray rmarkers;
      double rareMAF;
      double chisqFilter;
      void HetSharing();
      void ExtractUnrelated(void);
      void Haplotyping(bool SaveFlag=false);
      IntArray *hap[2];

// Quantitative Trait Association
      static double MACHEPS;
      static double MACHEPS_SQRT;
      static double cbrent;
      double minimize(double a, double b, double eps, double &funcx, int &numiter, int maxIter, bool quiet);

      int normalization;
      bool HeritFlag;
      bool slowFlag;
      void ROADTRIPS();

      int FixedEff;
      bool noiterFlag;
      bool svdoutFlag;
      String svdinfile;
      IntArray ID;
      double **UT;
      double **VR;
      //Matrix D;
      double **D;
      double *varVR;
      void ReadSVD();
      void ComputeSimilarity();
      void ComputeScoreSVD(bool quiet);
      void ComputeSVD();
      void ComputeMTSVD();

      IntArray ID0;
      IntArray validpheno;
      int **missingpheno;
      String dosagefile;
      String dfamfile;
      String dmapfile;
      Vector lambda0, tau0;
      void PreScan_FixedEff();
      void PreScan_MTSCORE();
      void PreScan_RSCORE();
      void GenomeScan();
      void GenomeScan64Bit();
      void ScoreScan(FILE *fp);
      void GenomeScanWithPermutation(FILE *fp);
      void DosageScan(FILE *fp);
      Vector *chisqs;
      void PreVC();
      void PostVC();
      void VC();
      void LMMSCORE();
      void FASTASSOC();
      void GRAMMAR();
      void RSCORE();
      void LMM();
      double fLL(double x);
      Matrix UX;
      Vector UY;
      Vector EV;
      int currentT;
      IntArray unrelatedList;
      void polygenic();
      void PrintPolygenic();
      void AssocFamHistory();
      void IBDGRM();
      void AllIBDSegments();

      Matrix means;
      Vector variances, heritabilities;
      IntArray traits, covariates, SampleSize, diseases;
      StringArray traitList;
      StringArray covariateList;
      int unrelatedFlag;
      bool effectFlag;
      bool CheckCovariates(Person & p);

//      void PrintGeneticRiskScore(const char* weightfile);
      void PrintGeneticRiskScoreSNPMajor(const char* weightfile, bool noflipFlag=false);
      double prevalence;
};

#endif

