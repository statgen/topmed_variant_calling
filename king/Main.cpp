////////////////////////////////////////////////////////////////////////
// Main.cpp
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
// March 28, 2019

#define VERSION "2.2"
//#define VERSION ""

#include <time.h>
#include "Pedigree.h"
#include "Parameters.h"
#include "analysis.h"
#include "MathStats.h"
#include "rplot.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

void ShowBanner()
{
#ifndef VERSION
   printf("KING - (c) 3/27/2019 Wei-Min Chen");
#else
   printf("KING " VERSION " - (c) 2010-2019 Wei-Min Chen");
#endif
   printf("\n");
}

int main(int argc, char **argv)
{
   ShowBanner();

   String binfile;
   ParameterList pl;
   bool flags[32];
   for(int i = 0; i < 32; i++)
      flags[i] = false;
   bool flipFlag = false;
   bool sibFlag = false;
   double kinFilter = _NAN_;
   bool xchrFlag = false;
   bool individual = false;
   bool adjustFamily = false;
   bool bpekingFlag = false;
   bool merlinFlag = false;
   bool distantFlag = false;
   bool poFlag = false;
   bool secondFlag = false;
   bool HeritFlag = false;
   bool detailFlag = false;
   bool autoflipFlag = false;
   bool rplotFlag = false;
   bool lessmemFlag = false;
   double callrateN=_NAN_;
   double callrateM=_NAN_;
   bool QCwipe = false;
   bool AssocVC = false;
   bool AssocFastAssoc = false;
   bool AssocRScoreAssoc = false;
   bool AssocGRAMMAR = false;
   bool AssocROADTRIPS = false;
   bool AssocADMIX = false;
   bool slowFlag = false;
   bool noiterFlag = false;
   String svdinfile("");
   bool svdoutFlag = false;
   bool AssocHet = false;
   bool pedstat = false;
   bool PhaseFlag = false;
   bool projectionFlag = false;
   double rareMAF = _NAN_;
   double prevalence = _NAN_;
   bool noflipFlag = false;
   int permuteCount= 0;
   int adjustPC = 0;
   double windowSize = 20.0;
   double moveSize = 10.0;
   String windowPar;
   int degree = 0;
   int normalization = 0;
   bool lmmFlag = false;
   int faster = 0;
   int slower = 0;
   int CoreCount = 0;
   int permCount = 1000;
//   double minMAF = _NAN_;
   double start = _NAN_;
   double stop = _NAN_;
   double pvalueFilter = _NAN_;
   bool AssocSKAT = false;
   bool AssocVCR = false;
   bool LinkageNPL = false;
   bool LinkageHE = false;
   bool LinkageIBD = false;
   bool LinkageIBDMI = false;
   bool LinkageGDT = false;
   bool LinkageAUC = false;
   bool LinkageAncestry = false;
   bool LinkagePOPDIST = false;
   bool LinkagePOPIBD = false;
   bool LinkageH2 = false;
   bool LinkageVC = false;
   bool LinkageGRM = false;
   bool LinkageIBDSEG = false;
   bool LinkagePOPROH = false;
   bool LinkageROHMAP = false;
   bool LinkageROHMAPQT = false;
   bool IBDAnalysis = false;
   bool LinkageMDS = false;
   int mincons = 0;
//   bool stratFLAG = false;
   String prefix("king");
   String weightfile("");
   String chr;
   String positions("");
   String stratName("");
//   int chr;
   String traitList, covariateList, diseaseList;
   String skatfile("");
   String vcrfile("");
   String dosagefile("");
   String Kernel("");
   String famfile("");
   String bimfile("");
   String phefile("");
   String covfile("");
   int sexchr=23;
   int IDadded = 0;
   bool splitPed = false;
   int Bit64 = 0;
   bool Bit64Flag = sizeof(void*)==8? true: false;

   BEGIN_LONG_PARAMETERS(longParameters)
      LONG_PARAMETER_GROUP("Close Relative Inference")
         LONG_PARAMETER("related", &flags[RelatedFLAG])
         LONG_PARAMETER("duplicate", &flags[DuplicateFLAG])
#ifndef VERSION
      LONG_PARAMETER_GROUP("Rare Variant Inference")
         LONG_PARAMETER("exact", &secondFlag) // 2nd-degree relationship
#endif
      LONG_PARAMETER_GROUP("Pairwise Relatedness Inference")
         LONG_PARAMETER("kinship", &flags[KinshipFLAG])
         LONG_PARAMETER("ibdseg", &flags[IBDSegFLAG])
         LONG_PARAMETER("ibs", &flags[IBSFLAG])
         LONG_PARAMETER("homog", &flags[HomoFLAG])
      LONG_PARAMETER_GROUP("Inference Parameter")
         LONG_INTPARAMETER("degree", &degree)
#ifndef VERSION
      LONG_PARAMETER_GROUP("Rare Variant Inference")
         LONG_PARAMETER("distant", &distantFlag) // distant relative pairs through sharing of ancient segments with a third relative
         LONG_PARAMETER("porel", &poFlag) // relatives of parent-offspring
         LONG_INTPARAMETER("adjustPC", &adjustPC)
      LONG_PARAMETER_GROUP("Computational Speedup")
         LONG_DOUBLEPARAMETER("noloss", &kinFilter)
         LONG_INTPARAMETER("faster", &faster)
         LONG_INTPARAMETER("slower", &slower)
#endif
      LONG_PARAMETER_GROUP("Relationship Application")
         LONG_PARAMETER("unrelated", &flags[UnrelatedFLAG])
         LONG_PARAMETER("cluster", &flags[ClusterFLAG])
         LONG_PARAMETER("build", &flags[BuildFLAG])
#ifndef VERSION
         LONG_PARAMETER("splitped", &splitPed)
         LONG_INTPARAMETER("idadded", &IDadded)
#endif
      LONG_PARAMETER_GROUP("QC Report")
         LONG_PARAMETER("bysample", &flags[BysampleFLAG])
         LONG_PARAMETER("bySNP", &flags[BysnpFLAG])
         LONG_PARAMETER("roh", &flags[ROHFLAG])
         LONG_PARAMETER("autoQC", &flags[AutoQCFLAG])
      LONG_PARAMETER_GROUP("QC Parameter")
         LONG_DOUBLEPARAMETER("callrateN", &callrateN)
         LONG_DOUBLEPARAMETER("callrateM", &callrateM)
#ifndef VERSION
         LONG_PARAMETER("semiped", &adjustFamily)
#endif
      LONG_PARAMETER_GROUP("Population Structure")
         LONG_PARAMETER("pca", &flags[PCAFLAG])
         LONG_PARAMETER("mds", &flags[MDSFLAG])
#ifndef VERSION
         LONG_PARAMETER("individual", &individual)
#endif
      LONG_PARAMETER_GROUP("Structure Parameter")
         LONG_PARAMETER("projection", &projectionFlag)
      LONG_PARAMETER_GROUP("Disease Association")
         LONG_PARAMETER("tdt", &flags[TDTFLAG])
      LONG_PARAMETER_GROUP("Quantitative Trait Association")
         LONG_PARAMETER("mtscore", &flags[MtscoreFLAG])
#ifndef VERSION
         LONG_PARAMETER("vc", &AssocVC)
         LONG_PARAMETER("rscore", &AssocRScoreAssoc)
         LONG_PARAMETER("fasta", &AssocFastAssoc)
         LONG_PARAMETER("grammar", &AssocGRAMMAR)
         //LONG_PARAMETER("herit", &HeritFlag)
      LONG_PARAMETER_GROUP("LMM Parameter")
         LONG_PARAMETER("lmm", &lmmFlag)
         LONG_STRINGPARAMETER("svdin", &svdinfile)
         LONG_STRINGPARAMETER("dosage", &dosagefile)
//         LONG_INTPARAMETER("fixedeff", &fixedeffCount)
         LONG_PARAMETER("noiter", &noiterFlag)
      LONG_PARAMETER_GROUP("Population Inference")
         LONG_PARAMETER("popdist", &LinkagePOPDIST)
         LONG_PARAMETER("popibd", &LinkagePOPIBD)
         LONG_PARAMETER("poproh", &LinkagePOPROH)
      LONG_PARAMETER_GROUP("Linkage Analysis")
         LONG_PARAMETER("npl", &LinkageNPL)
         LONG_PARAMETER("HEreg", &LinkageHE)
      LONG_PARAMETER_GROUP("IBD Analysis")
         LONG_PARAMETER("ibdall", &LinkageIBDSEG)
         LONG_PARAMETER("ibdGRM", &LinkageGRM)
         LONG_PARAMETER("ibdmap", &LinkageIBD)
         LONG_PARAMETER("ibdmds", &LinkageMDS)
         LONG_PARAMETER("ibdgdt", &LinkageGDT)
         LONG_PARAMETER("ibdMI", &LinkageIBDMI)
         LONG_PARAMETER("ibdH2", &LinkageH2)
         LONG_PARAMETER("ibdvc", &LinkageVC)
         LONG_PARAMETER("aucmap", &LinkageAUC)
         LONG_PARAMETER("ancestry", &LinkageAncestry)
      LONG_PARAMETER_GROUP("Homozygosity Mapping")
         LONG_PARAMETER("homomap", &LinkageROHMAP)
         LONG_PARAMETER("mthomo", &LinkageROHMAPQT)
      LONG_PARAMETER_GROUP("IBD Parameter")
         LONG_INTPARAMETER("mincons", &mincons)
         LONG_INTPARAMETER("nperm", &permCount)
         LONG_STRINGPARAMETER("strat", &stratName)
      LONG_PARAMETER_GROUP("Association Analysis")
         LONG_PARAMETER("admix", &AssocADMIX)
         LONG_PARAMETER("roadtrips", &AssocROADTRIPS)
//      LONG_PARAMETER_GROUP("Local Ancestry")
//         LONG_INTPARAMETER("window", &windowSize)
//         LONG_INTPARAMETER("move", &moveSize)
//      LONG_PARAMETER_GROUP("Rare Variant Analysis")
//         LONG_PARAMETER("skat", &AssocSKAT)
//         LONG_STRINGPARAMETER("batch", &skatfile)
//      LONG_PARAMETER_GROUP("Kernel-Based Analysis")
         LONG_STRINGPARAMETER("skat", &skatfile)
//         LONG_STRINGPARAMETER("vcr", &vcrfile)
         LONG_STRINGPARAMETER("kernel", &Kernel)
         //LONG_PARAMETER("phase", &PhaseFlag)
         //LONG_PARAMETER("hetsharing", &AssocHet)
//      LONG_PARAMETER_GROUP("Rare Variant Parameter")
//         LONG_INTPARAMETER("permut", &permuteCount)
//         LONG_DOUBLEPARAMETER("rare", &rareMAF)
#endif
      LONG_PARAMETER_GROUP("Association Model")
         LONG_STRINGPARAMETER("trait", &traitList)
         LONG_STRINGPARAMETER("covariate", &covariateList)
      LONG_PARAMETER_GROUP("Association Parameter")
         LONG_PARAMETER("invnorm", &normalization)
         //LONG_PARAMETER("slow", &slowFlag)
         LONG_DOUBLEPARAMETER("maxP", &pvalueFilter)
      LONG_PARAMETER_GROUP("Genetic Risk Score")
         LONG_PARAMETER("risk", &flags[RiskFLAG])
         LONG_STRINGPARAMETER("model", &weightfile)
         LONG_DOUBLEPARAMETER("prevalence", &prevalence)
         LONG_PARAMETER("noflip", &noflipFlag)
#ifndef VERSION
      LONG_PARAMETER_GROUP("Chromosome Position")
         LONG_STRINGPARAMETER("position", &positions)
         LONG_STRINGPARAMETER("chr", &chr)
         LONG_DOUBLEPARAMETER("start", &start)
         LONG_DOUBLEPARAMETER("stop", &stop)
//         LONG_STRINGPARAMETER("window", &windowPar)
#endif
      LONG_PARAMETER_GROUP("Computing Parameter")
         LONG_INTPARAMETER("cpus", &CoreCount)
#ifndef VERSION
         LONG_PARAMETER("lessmem", &lessmemFlag)
         LONG_INTPARAMETER("sysbit", &Bit64)
         LONG_PARAMETER("wipeMI", &QCwipe)
         LONG_PARAMETER("pedstat", &pedstat)
         LONG_PARAMETER("autoflip", &autoflipFlag)
         LONG_PARAMETER("flip", &flipFlag)
#endif
      LONG_PARAMETER_GROUP("Optional Input")
         LONG_STRINGPARAMETER("fam", &famfile)
         LONG_STRINGPARAMETER("bim", &bimfile)
         LONG_INTPARAMETER("sexchr", &sexchr)
#ifndef VERSION
         LONG_STRINGPARAMETER("phefile", &phefile)
         LONG_STRINGPARAMETER("covfile", &covfile)
#endif
      LONG_PARAMETER_GROUP("Output")
         LONG_STRINGPARAMETER("prefix", &prefix)
         LONG_PARAMETER("rplot", &rplotFlag)
#ifndef VERSION
         LONG_PARAMETER("svdout", &svdoutFlag)
      LONG_PARAMETER_GROUP("Rewrite In Format")
         LONG_PARAMETER("plink", &flags[PLINKFLAG])
         LONG_PARAMETER("king", &bpekingFlag)
         LONG_PARAMETER("merlin", &merlinFlag)
#endif
   END_LONG_PARAMETERS()

   pl.Add(new StringParameter('b', "Binary File", binfile));
   pl.Add(new LongParameters("Additional Options", longParameters));
   pl.Read(argc, argv);
   pl.Status();

   if(binfile.IsEmpty())
      error("Genotype files are required. e.g.,\n  king -b ex.bed --related\n\nPlease check the reference paper Manichaikul et al. 2010 Bioinformatics,\n\t\t\t\t\tChen et al. 2019,\n          or the KING website at http://people.virginia.edu/~wc9c/KING");

   unsigned int allflags = 0;
   StringArray flagNames(TOTALFLAGCOUNT);
   for(int i = 0; i < TOTALFLAGCOUNT; i++)
      switch(i){
         case DuplicateFLAG: flagNames[i] = "duplicate"; break;
         case RelatedFLAG: flagNames[i] = "related"; break;
         case AutoQCFLAG: flagNames[i] = "autoQC"; break;
         case MtscoreFLAG: flagNames[i] = "mtscore"; break;
         case RiskFLAG: flagNames[i] = "risk"; break;
         case PLINKFLAG: flagNames[i] = "plink"; break;
         case IBSFLAG: flagNames[i] = "ibs"; break;
         case KinshipFLAG: flagNames[i] = "kinship"; break;
         case HomoFLAG: flagNames[i] = "homog"; break;
         case IBDSegFLAG: flagNames[i] = "ibdseg"; break;
         case MDSFLAG: flagNames[i] = "mds"; break;
         case PCAFLAG: flagNames[i] = "pca"; break;
         case ClusterFLAG: flagNames[i] = "cluster"; break;
         case BuildFLAG: flagNames[i] = "build"; break;
         case BysampleFLAG: flagNames[i] = "bysample"; break;
         case BysnpFLAG: flagNames[i] = "bysnp"; break;
         case TDTFLAG: flagNames[i] = "tdt"; break;
         case UnrelatedFLAG: flagNames[i] = "unrelated"; break;
         case ROHFLAG: flagNames[i] = "roh"; break;
      }

   unsigned int base[32];
   for(int i = 0; i < 32; i++)
      base[i] = 1 << i;
   for(int i = 0; i < TOTALFLAGCOUNT; i++)
      if(flags[i])   // ith command is specified
         allflags |= base[i];

#ifdef VERSION
   if(allflags==0 && !rplotFlag){
      printf("\nPlease specify one of the following %d options:", TOTALFLAGCOUNT);
      for(int i = 0; i < TOTALFLAGCOUNT; i++)
         printf(" --%s", (const char*)flagNames[i]);
      printf("\n\n");
      flags[RelatedFLAG]=true;
   }else if(allflags & (allflags-1))  // more than 1 command
      if(!flags[ClusterFLAG]){   // --cluster is the only exception
         printf("The following analyses will run separately: ");
         for(int i = 0; i < TOTALFLAGCOUNT; i++)
            if(flags[i])
               printf(" --%s", (const char*)flagNames[i]);
         printf("\n\n");
      }
#endif

   time_t timer;
   struct tm *tblock;
   timer = time(NULL);
   tblock = localtime(&timer);
   printf("KING starts at %s", asctime(tblock));

   Pedigree ped;
   Engine engine(ped);
   engine.allflags = allflags;
   if(ped.familyCount && autoflipFlag) {
      engine.autoflipFlag = true;
      printf("Allele flipping is enabled.\n");
   }
   if(ped.familyCount && flipFlag) {
      engine.flipFlag = true;
      printf("Flipping at SNPs with labels AT or CG is enabled.\n");
   }
   if(normalization) engine.normalization = normalization;
   if(degree) engine.relativedegree = degree;
   //   engine.relativedegree = degree > 1? degree: 1;
   engine.prefix = prefix;
   if(pvalueFilter!=_NAN_){
       double t = ninv(pvalueFilter*0.5);
       engine.chisqFilter = t*t;
   }
   if(projectionFlag) engine.projectFlag = true;
   if(CoreCount) engine.CoreCount = CoreCount;
   if(sexchr != 23 && sexchr > 1) {
      engine.SEXCHR = sexchr;
      printf("Non-human samples are analyzed, with %d pairs of chromosomes\n", sexchr);
   }else if(sexchr <= 1)
      error("Sex chromosome %d out of range.\n", sexchr);
#ifdef _OPENMP
   int availableCoreCount = omp_get_max_threads();
   if(engine.CoreCount > 0){
      engine.defaultEfficientCoreCount = engine.CoreCount > 8? 8: engine.CoreCount;
      engine.defaultMaxCoreCount = engine.CoreCount;
   }else{
      if(availableCoreCount == 1)
         engine.defaultEfficientCoreCount = engine.defaultMaxCoreCount = 1;
      else{
         engine.defaultEfficientCoreCount = availableCoreCount > 8? 8: availableCoreCount;
         engine.defaultMaxCoreCount = availableCoreCount / 2;
      }
   }
#endif
   if(LinkageNPL || LinkageHE || LinkageIBD || LinkageGDT || LinkageAUC || LinkageAncestry ||
       LinkagePOPDIST || LinkagePOPIBD || LinkageIBDMI || LinkageH2 || LinkageIBDSEG ||
       LinkagePOPROH || LinkageROHMAP || LinkageROHMAPQT || LinkageVC || LinkageGRM || LinkageMDS)
      IBDAnalysis = true;
   // NOT autosomeonly:
   if(flags[KinshipFLAG] || flags[IBSFLAG] || flags[HomoFLAG] || flags[RiskFLAG] ||
      flags[BysnpFLAG] || flags[BysampleFLAG] || flags[AutoQCFLAG] ||
      flags[PLINKFLAG] || flags[TDTFLAG] || flags[MtscoreFLAG] ||
      AssocVC || AssocRScoreAssoc || IBDAnalysis || flags[BuildFLAG] ||
      AssocFastAssoc || AssocGRAMMAR || AssocSKAT || merlinFlag )
      engine.autosomeOnly = false;

   // NOT genotypeOnly:
   if(flags[RelatedFLAG] || flags[IBSFLAG] || flags[IBDSegFLAG] || flags[PLINKFLAG] ||
      flags[AutoQCFLAG] || flags[BysnpFLAG] || flags[RiskFLAG] || flags[TDTFLAG] || flags[MtscoreFLAG] ||
      AssocVC || AssocRScoreAssoc || flags[ROHFLAG] || IBDAnalysis || secondFlag || 
      flags[UnrelatedFLAG] || flags[ClusterFLAG] || flags[BuildFLAG] ||
      AssocFastAssoc || AssocGRAMMAR || AssocSKAT || merlinFlag )
      engine.genotypeOnly = false;

   if(prevalence!=_NAN_)
      engine.prevalence = prevalence;
   if(rplotFlag) engine.rplotFlag = true;
#ifndef VERSION
   if(mincons) engine.mincons = mincons;
   if(kinFilter!=_NAN_) flags[KinshipFLAG] = true;
   if(!chr.IsEmpty()){
      StringArray chrList;
      chrList.AddTokens(chr, ',');
      for(int i = 0; i < chrList.Length(); i++){
         if(chrList[i].IsNumber()){
            int k = int(chrList[i]);
            if(k > 0 && k < 27) engine.chrList.Push(k);
         }else if(chrList[i] == "X" || chrList[i] == "x")
            engine.chrList.Push(engine.SEXCHR);
         else if(chrList[i] == "Y" || chrList[i] == "y")
            engine.chrList.Push(engine.SEXCHR+1);
         else if(chrList[i] == "XY" || chrList[i] == "xy")
            engine.chrList.Push(engine.SEXCHR+2);
         else if(chrList[i] == "MT" || chrList[i] == "mt")
            engine.chrList.Push(engine.SEXCHR+3);
         else{ // 1-22
            StringArray tempList;
            tempList.AddTokens(chrList[i], '-');
            if(tempList.Length()==2 && tempList[0].IsNumber() && tempList[1].IsNumber()){
               int j = int(tempList[0]);
               int k = int(tempList[1]);
               if(j < 1 || j > engine.SEXCHR+3 || k < 1 || k > engine.SEXCHR+3) continue;
               if(j <= k) for(; j<=k; j++) engine.chrList.Push(j);
               else for(; j>=k; j--) engine.chrList.Push(j);
            }
         }
      }
      if(engine.chrList.Length()==1 && engine.chrList[0]==engine.SEXCHR)
         engine.xflag = true;
   }
   if(!windowPar.IsEmpty()){
      StringArray windowList;
      windowList.AddTokens(windowPar, ',');
      if(windowList.Length()>0 && windowList[0].IsNumber())
         windowSize = double(windowList[0]);
      if(windowList.Length()>1 && windowList[1].IsNumber())
         moveSize = double(windowList[1]);
   }
   if(slowFlag) engine.slowFlag = true;
   if(detailFlag) engine.detailFlag = true;
//   if(flags[HomoFLAG]) engine.mafFilterFlag = true;
//   if(minMAF != _NAN_) engine.minMAF = minMAF;
   if(start != _NAN_) engine.start = start;
   if(stop != _NAN_) engine.stop = stop;
   if(adjustFamily) engine.adjustFamily = adjustFamily;
//   if(runlengthFlag) engine.runlengthFlag = true;
   if(individual) engine.individualInfo = individual;
   if(QCwipe) {engine.QCwipe = QCwipe;}
   if(adjustFamily) engine.semifamilyFlag = true;
   if(flags[PLINKFLAG]) engine.SaveFormat = "PLINK";
   else if(merlinFlag) engine.SaveFormat = "MERLIN";
   else if(bpekingFlag) engine.SaveFormat = "KING";
   if(pedstat) engine.pedstatFlag = true;
   if(rareMAF!=_NAN_) engine.rareMAF = rareMAF;
   if(permuteCount > 0) engine.permuteCount = permuteCount;
   if(svdinfile!="") engine.svdinfile = svdinfile;
   if(svdoutFlag) engine.svdoutFlag = true;
   if(noiterFlag) engine.noiterFlag = true;
   if(lmmFlag) engine.FixedEff = 0;
   if(dosagefile!=""){
      StringArray fileList;
      fileList.AddTokens(dosagefile, ',');
      engine.dosagefile = fileList[0];
      if(fileList.Length()>1)
         engine.dfamfile = fileList[1];
      if(fileList.Length()>2)
         engine.dmapfile=fileList[2];
   }
   if(skatfile!="") AssocSKAT = true;
   if(vcrfile!="") AssocVCR = true;
#endif

   if(flags[RelatedFLAG] || (kinFilter!=_NAN_ && kinFilter < 1.0) || flags[DuplicateFLAG] || distantFlag || poFlag || secondFlag \
   || flags[ClusterFLAG] || flags[BuildFLAG] || flags[UnrelatedFLAG] || IBDAnalysis || bpekingFlag || flags[AutoQCFLAG] || individual || adjustPC) {
      engine.bigdataFlag = true;
      if(faster>1) engine.faster = faster;
      if(slower>1) engine.slower = slower;
   }

   if(Bit64 == 64){  // --sysbit 64
      if(Bit64Flag && (allflags & (base[AutoQCFLAG] | base[HomoFLAG] | base[PCAFLAG]))){
         printf("64-bit algorithm to be implemented.\n");
         Bit64 = 32;
         engine.Bit64 = 32;
      }else if(Bit64Flag){
         printf("64-bit algorithms will be used.\n");
         Bit64 = 64;
         engine.Bit64 = 64;
         engine.Bit64Flag = true;
      }else if(!Bit64Flag){
         printf("This is not a 64-bit system and request is ignored.\n");
         Bit64 = 32;
         engine.Bit64 = 32;
      }
   }else if(Bit64 == 32){  // --sysbit 32
      printf("32-bit algorithms will be used.\n");
      Bit64 = 32;
      engine.Bit64 = 32;
      engine.Bit64Flag = false;
   }else{   // sysbit not (properly) specified
      if(Bit64)
         printf("--sysbit %d is ignored. Please select either 32 or 64 bit for your system.\n", Bit64);
      if((!Bit64Flag) || (allflags & (base[AutoQCFLAG] | base[HomoFLAG] | base[PCAFLAG]))){
            Bit64 = 32;
            engine.Bit64 = 32;
      }else{
            Bit64 = 64;
            engine.Bit64 = 64;
      }
   }
//***************************************************************************
   // Read binary data here
   if(flags[RiskFLAG]){   // SNP-major for internal format
      if(weightfile.IsEmpty()) error("Please use --model <file> to specify a risk model.");
      if(!famfile.IsEmpty() || !bimfile.IsEmpty() || !covfile.IsEmpty() || !phefile.IsEmpty())
         engine.ReadPlinkBinaryBigDataSNPMajor(binfile, weightfile, famfile, bimfile, covfile, phefile);
      else
         engine.ReadPlinkBinaryBigDataSNPMajor(binfile, weightfile);
      if(noflipFlag)
         engine.PrintGeneticRiskScoreSNPMajor(weightfile, true);
      else
         engine.PrintGeneticRiskScoreSNPMajor(weightfile);
      exit(0);
   }

   if(lessmemFlag) engine.lessmemFlag = true;
   String filenames(binfile);
   if(filenames.Find(",")>-1){
//      Bit64 = 32;
//      engine.Bit64 = 32;
      if(!famfile.IsEmpty() || !bimfile.IsEmpty() || !covfile.IsEmpty() || !phefile.IsEmpty())
         engine.ReadMultiplePlinkBinaryBigData(binfile, famfile, bimfile, covfile, phefile);
      else
         engine.ReadMultiplePlinkBinaryBigData(binfile);
   }else{
      String line;
      StringArray tokens;
      IFILE input = ifopen(binfile, "rb");
      if(input == NULL)
         error("Genotype file %s cannot be opened", (const char*) binfile);
      // line 1: idCount markerCount
      tokens.Clear();
      line.ReadLine(input);
      ifclose(input);
      if((line[0] == 108) && (line[1] == 27)){
         String pedfile(binfile);
         int m = pedfile.Find(".bed");
         if(m == -1) error("Please use either MERLIN or KING binary format as input.");
         String filename=pedfile.SubStr(0,m);
         if(!famfile.IsEmpty() || !bimfile.IsEmpty() || !covfile.IsEmpty() || !phefile.IsEmpty())
            engine.ReadPlinkBinaryBigData(filename, famfile, bimfile, covfile, phefile);
         else
            engine.ReadPlinkBinaryBigData(filename);
      }else{
         String pedfile(binfile);
         int m = pedfile.Find(".king");
         if(m == -1) error("Please use either PLINK or KING binary format as input.");
         String filename=pedfile.SubStr(0,m+5);
         engine.ReadKingBinaryData(filename);
      }
   }
   if(!traitList.IsEmpty()){
      engine.traitList.AddTokens(traitList, ',');
      engine.traits.Dimension(0);
      for(int i = 0; i < engine.traitList.Length(); i++){
         int k = ped.traitNames.Find(engine.traitList[i]);
         if(k==-1){
            if((ped.affectionNames.Find(engine.traitList[i]))!=-1){
               k = ped.affectionNames.Find(engine.traitList[i]);
               for(int p = 0; p < ped.count; p++)
                  if(ped[p].affections[k]>0)
                     ped[p].traits.Push(ped[p].affections[k]-1);
                  else
                     ped[p].traits.Push(_NAN_);
               ped.traitNames.Push(ped.affectionNames[k]);
               engine.traits.Push(ped.traitCount);
               ped.traitCount++;
               printf("Binary trait %s is now considered as a quantitative trait.",
                  (const char*)ped.affectionNames[k]);
            }else
               printf("Trait %s cannot be found.\n",
                  (const char*)engine.traitList[i]);
         }else
            engine.traits.Push(k);
      }
      if(engine.traits.Length()==0 && engine.diseases.Length()==0)
         error("No traits are specified.");
   }else{
      engine.traits.Dimension(ped.traitCount);
      for(int i = 0; i < ped.traitCount; i++)
         engine.traits[i] = i;
   }

   if(!covariateList.IsEmpty()){
      engine.covariateList.AddTokens(covariateList, ',');
      engine.covariates.Dimension(0);
      for(int i = 0; i < engine.covariateList.Length(); i++){
         int k = ped.covariateNames.Find(engine.covariateList[i]);
         if(k==-1){
            if(engine.covariateList[i].Compare("sex")==0){
               for(int p = 0; p < ped.count; p++)
                  if(ped[p].sex == 0)
                     ped[p].covariates.Push(_NAN_);
                  else
                     ped[p].covariates.Push(ped[p].sex-1);
               ped.covariateNames.Push("sex");
               ped.covariateCount++;
               engine.covariates.Push(ped.covariateCount-1);
            }else if((k = ped.markerNames.Find(engine.covariateList[i])) > -1){
               for(int p = 0; p < ped.count; p++)
                  if(!ped[p].markers[k].isKnown())
                     ped[p].covariates.Push(_NAN_);
                  else
                     ped[p].covariates.Push(ped[p].markers[k].countAlleles(2));
               engine.covariateList[i] += "_";
               engine.covariateList[i] += ped.markerInfo[k]->GetAlleleLabel(2);
               ped.covariateNames.Push(engine.covariateList[i]);
               ped.covariateCount ++;
               engine.covariates.Push(ped.covariateCount-1);
            }else
               printf("Covariate %s cannot be found.\n", (const char*)engine.covariateList[i]);
         }else{
            if(engine.covariates.Find(k)==-1)
               engine.covariates.Push(k);
            else
               printf("Covariate %s is duplicated.\n",
                  (const char*)engine.covariateList[i]);
         }
      }
      if(engine.covariates.Length()==0)
         printf("No covariates are included in the analysis.\n");
   }
   if(diseaseList.IsEmpty()){
      engine.diseases.Dimension(ped.affectionCount);
      for(int i = 0; i < ped.affectionCount; i++)
         engine.diseases[i] = i;
   }

   if(adjustPC)
      engine.adjustPC(adjustPC);
   bool familyData = false;
   for(int f = 0; f < ped.familyCount; f++){
      int count = 0;
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(ped[i].ngeno >= MINSNPCOUNT) count++;
      if(count > 1) {
         familyData = true;
         break;
      }
   }
   if(!allflags && rplotFlag){
      if(engine.SplitPedigree()) plotSplitped(prefix);
      exit(0);
   }
   for(int flag = 0; flag < TOTALFLAGCOUNT; flag++)
      if(flags[flag])
         switch(flag){
            case ClusterFLAG:
               engine.semifamilyFlag = true;
               if(flags[MDSFLAG] || flags[BysampleFLAG] || flags[BysnpFLAG] || flags[TDTFLAG]){
                  engine.ClusterFamily(2, degree);
                  if(flags[MDSFLAG]) {
                     engine.mds_family();
                     if(rplotFlag) plotPopStructure(prefix, projectionFlag);
                  }
                  if(flags[BysampleFLAG]) {
                     if(Bit64==64)
                        engine.QC_By_Sample64Bit();
                     else
                        engine.QC_By_Sample();
                  }
                  if(flags[BysnpFLAG]) {
                     if(Bit64==64)
                        engine.QC_By_SNP64Bit();
                     else
                        engine.QC_By_SNP();
                  }
                  if(flags[TDTFLAG]) engine.TDT();
               }else{   // --cluster by itself
                  engine.ClusterFamily(0, degree);
                  if(rplotFlag && Bit64==64 && engine.totalLength>10000000) plotCluster(prefix);
               }
               break;
            case BysampleFLAG:
               if(!flags[ClusterFLAG]){
                  if(Bit64==64)
                     engine.QC_By_Sample64Bit();
                  else
                     engine.QC_By_Sample();
               }
               if(rplotFlag) printf("R plot for --bysample is not available.\n");
               break;
            case BysnpFLAG:
               if(!flags[ClusterFLAG]){
                  if(Bit64==64)
                     engine.QC_By_SNP64Bit();
                  else
                     engine.QC_By_SNP();
               }
               if(rplotFlag) printf("R plot for --bySNP is not available.\n");
               break;
            case TDTFLAG:
               if(!flags[ClusterFLAG])
                  engine.TDT();
               if(rplotFlag) printf("R plot for --tdt is not available.\n");
               break;
            case MDSFLAG:
               if(!flags[ClusterFLAG]){
                  if(Bit64==64){
                     if(projectionFlag)
                        engine.mds_projection();
                     else
                        engine.mds();
                  }else{
                  /*
                  if(familyData){
                     if(flags[KinshipFLAG]){
                        printf("Hidden relatedness across families is incorporated in MDS.\n");
                        engine.mds_kin_family();
                     }else if(adjustFamily){
                        printf("Relatedness will be verified within each family.\n");
                        engine.semifamilyFlag = true;
                        engine.mds_family();
                     }else
                        engine.mds_family();
                  }else{
                     if(flags[KinshipFLAG]){
                        printf("Hidden relatedness is incorporated in MDS.\n");
                        engine.mds_kin();
                     }else
                        engine.mds();
                  } */
                        engine.mds();
                  }
                  if(rplotFlag) plotPopStructure(prefix, projectionFlag);
               }
               break;
            case RelatedFLAG:
               if(Bit64==32 && (allflags & (base[AutoQCFLAG] | base[HomoFLAG] | base[PCAFLAG]))){
                  printf("--related is skipped.\n");
                  if(allflags & base[AutoQCFLAG]) printf("Please do not run --related together with --autoQC\n");
                  if(allflags & base[HomoFLAG]) printf("Please do not run --related together with --homog\n");
                  if(allflags & base[PCAFLAG]) printf("Please do not run --related together with --pca\n");
               }else{
                  if(engine.idCount < 10){
                     if(allflags & base[KinshipFLAG])
                        printf("\n--related is skipped for a rather small sample size.\n");
                     else{
                        printf("\n--related is replaced with --kinship for a small sample size.\n");
                        engine.ComputeLongRobustKinship64Bit();
                     }
                  }else{   // integrated inference
                     engine.IntegratedRelationshipInference();
                     if(rplotFlag) {
                        plotRelationship(prefix);
                        if(engine.SplitPedigree()) plotMIerror(prefix);
                        plotUniqueFamily(prefix, degree==0?1:degree, "related");
                     }
                  }
               }
               break;
            case DuplicateFLAG:
               if(Bit64==64)
                  engine.ComputeBigDataDuplicate64Bit();
               else
                  engine.ComputeBigDataDuplicate();
               if(rplotFlag) plotDuplicate(prefix);
               break;
            case UnrelatedFLAG:
               engine.unrelatedExtraction = true;
               engine.allflags |= base[ClusterFLAG];
               engine.semifamilyFlag = true;
               engine.ClusterFamily(2, degree);
               for(int i = 0; i < engine.ped.count; i++)
                  engine.ped[i].affections[0] = 1; // make sure all individuals are analyzed
               engine.mds_family_projection();
               if(rplotFlag) printf("R plot for --unrelated is not available.\n");
               break;
            case AutoQCFLAG:
               engine.autoQC(callrateN==_NAN_? 0.95: callrateN, callrateM==_NAN_? 0.95: callrateM);
               break;
            case MtscoreFLAG:
               engine.LMMSCORE();
               if(rplotFlag) printf("R plot for --mtscore is not available.\n");
               break;
            case BuildFLAG:
               bool built;
               if(engine.ClusterFamily(1, degree)){
                  if(IDadded) built = engine.rebuild(IDadded);
                  else built = engine.rebuild();
               }
               if(rplotFlag && built) plotBuild(prefix);
               break;
            case PCAFLAG:
               if(familyData && !projectionFlag){
                  if(adjustFamily){
                     printf("Relatedness will be verified within each family.\n");
                     engine.pca_family_check(1);
                  }else{
                     printf("Families are incorporated in PCA.\n");
                     engine.pca_family(1);
                  }
               }else
                  engine.pca(1);
               if(rplotFlag) plotPopStructure(prefix, projectionFlag);
               break;
            case IBSFLAG:
               if(Bit64 == 64)
                  engine.ComputeExtendedIBS64Bit();
               else
                  engine.ComputeShortExtendedIBS();
               if(rplotFlag) printf("R plot for --ibs is not available.\n");
               break;
            case KinshipFLAG:
               if(Bit64==64){
                  if(degree){
                     IntArray *temp=new IntArray[engine.defaultMaxCoreCount];
                     engine.ComputeLongRobustKinship64BitWithFilter(temp);
                     delete []temp;
                  }else
                     engine.ComputeLongRobustKinship64Bit();
               }else
                  engine.ComputeShortRobustKinship();
               if(rplotFlag) printf("R plot for --kinship is not available.\n");
               break;
            case IBDSegFLAG:
               if(Bit64==64){
                  if(ped.count < 10){
                     printf("\n--kinship analysis carried out instead for such a small sample size.\n");
                     if(degree){
                        IntArray *temp=new IntArray[engine.defaultMaxCoreCount];
#ifndef VERSION
                        engine.ComputeLongRobustKinship64BitWithFilter(temp);
#else
                        engine.ComputeLongRobustKinship64Bit();
#endif
                        delete []temp;
                     }else
                        engine.ComputeLongRobustKinship64Bit();
                  }else{
                     IntArray temp1(0), temp2(0);
                     if(degree)
#ifndef VERSION
                        //engine.ComputeIBDSegment64BitWithFilter(temp1, temp2);
                        engine.ComputeIBDSegment64Bit(temp1, temp2);
#else
                        engine.ComputeIBDSegment64Bit(temp1, temp2);
#endif
                     else
                        engine.ComputeIBDSegment64Bit(temp1, temp2);
                  }
               }else{
                  printf("--ibdseg is skipped.\n");
                  if(allflags & base[AutoQCFLAG]) printf("Please do not run --ibdseg together with --autoQC\n");
                  if(allflags & base[HomoFLAG]) printf("Please do not run --ibdseg together with --homog\n");
                  if(allflags & base[PCAFLAG]) printf("Please do not run --ibdseg together with --pca\n");
               }
               if(rplotFlag) {
                  plotIBDSeg(prefix);
                  if(degree) plotUniqueFamily(prefix, degree, "ibdseg");
                  else
                     printf("  For additional relative pair plots please use --degree, --degree 2, or --degree 3.\n");
               }
               break;
            case HomoFLAG:
               if(xchrFlag)
                  engine.ComputeShortFastXHomoKinship();
               else
                  engine.ComputeShortFastHomoKinship();
               if(rplotFlag) printf("R plot for --homog is not available.\n");
               break;
            case ROHFLAG:
               if(Bit64==64)
                  engine.ROH();
               else{
                  printf("--roh is skipped.\n");
                  if(allflags & base[AutoQCFLAG]) printf("Please do not run --roh together with --autoQC\n");
                  if(allflags & base[HomoFLAG]) printf("Please do not run --roh together with --homog\n");
                  if(allflags & base[PCAFLAG]) printf("Please do not run --roh together with --pca\n");
               }
               if(rplotFlag) printf("R plot for --roh is not available.\n");
               break;
            default:
               printf("Not available\n");
         }

#ifndef VERSION
   if(IBDAnalysis){
      if(LinkageNPL){
         engine.NPL();
         if(rplotFlag) plotNPL(prefix, sexchr);
      }
      if(LinkageHE){
         engine.HEreg();
         if(rplotFlag) plotHEreg(prefix, sexchr);
      }
      if(LinkageIBD){
         engine.IBDmapping(permCount);
         if(rplotFlag) plotIBDmapping(prefix, sexchr);
      }
      if(LinkageGDT)
         engine.IBDGDT();
      if(LinkageAncestry){
         engine.AncestryInference();
         if(rplotFlag) plotAncestry(prefix);
      }
      if(LinkagePOPDIST){
         engine.PopulationDistance();
         if(rplotFlag) plotPopDist(prefix);
      }
      if(LinkageAUC){
         if(!positions.IsEmpty()){
            StringArray posList;
            StringArray onePosition;
            IntArray allchr(0);
            IntArray allpos(0);
            posList.AddTokens(positions, ',');
            for(int i = 0; i < posList.Length(); i++){
               onePosition.Clear();
               onePosition.AddTokens(posList[i], ':');
               if(onePosition.Length()<2) error("Format: chr:pos,chr:pos,...");
               if(onePosition[0].IsNumber()){
                  allchr.Push(int(onePosition[0]));
                  allpos.Push(int(onePosition[1]));
               }//else if(chrList[i] == "X" || chrList[i] == "x")
            }
            engine.AUCpredicting(allchr, allpos);
         }else{
            engine.AUCmapping();
            if(rplotFlag) plotAUCmapping(prefix, sexchr);
         }
      }
      if(LinkagePOPIBD){
         engine.PopulationIBD();
//         if(rplotFlag) plotPopDist(prefix);
      }
      if(LinkagePOPROH){
         engine.PopulationROH();
         if(rplotFlag) plotPopROH(prefix, sexchr);
      }
      if(LinkageIBDMI)
         engine.IBDMI();
      if(LinkageH2)
         engine.LocalH2();
      if(LinkageROHMAP){
         if(stratName!="")
            engine.HomozygosityMappingMH(stratName);
         else
            engine.HomozygosityMapping();
         if(rplotFlag) plotROHmapping(prefix, stratName, sexchr);
      }
      if(LinkageROHMAPQT){
         if(stratName!="")
            engine.HomozygosityMappingForQTMH(stratName);
         else
            engine.HomozygosityMappingForQT();
         if(rplotFlag) plotROHforQT(prefix, sexchr);
      }
      if(LinkageVC)
         engine.IBDVC();
      if(LinkageGRM)
         engine.IBDGRM();
      if(LinkageIBDSEG)
         engine.AllIBDSegments();
      if(LinkageMDS){
         if(projectionFlag)
            engine.IBDMDS_Projection();
         else
            engine.IBDMDS();
         if(rplotFlag) plotPopStructure(prefix, projectionFlag);
      }
      exit(0);
   }
   if(splitPed){
      engine.SplitPedigree();
      if(rplotFlag) plotSplitped(prefix);
      exit(0);
   }
   if(PhaseFlag){
      engine.Haplotyping(1);
      exit(0);
   }
   if(HeritFlag)
      engine.HeritFlag = true;
   if(AssocSKAT){
      if(Kernel=="A" || Kernel=="a") engine.KernelShort='A'; // admixture/local ancestry adjusted
      else if(Kernel=="G" || Kernel=="g") engine.KernelShort='G'; // glocal ancestry adjusted
      else if(Kernel=="L") engine.KernelShort='L'; // weighted linear
      else if(Kernel=="l") engine.KernelShort='l'; // weighted linear
      else if(Kernel=="D" || Kernel=="d") engine.KernelShort='D'; // Distance, unadjusted
      else {
         engine.KernelShort='L';
         //printf("Please select L, D, G, or A as kernel for SKAT.\n");
         printf("The default Weighted Linear Kernel is used.\n");
      }
      engine.SKAT_family_quantitatives(skatfile);
      exit(0);
   }
   if(AssocVCR){
      engine.VCR_family_quantitatives(vcrfile);
      exit(0);
   }

   if(!flags[KinshipFLAG] && !individual) flags[KinshipFLAG] = true;
/*
   if(merlinFlag){
      engine.WriteMerlin();
      exit(0);
   }else if(flags[PLINKFLAG]){
      engine.WritePlink();
      exit(0);
   }else if(bpekingFlag && (!engine.genotypeOnly || !engine.autosomeOnly)){
      String pedfile = prefix;
      pedfile.Add(".king");
      engine.WriteKingBinary(pedfile);
      exit(0);
   }
   */
   if(ped.count < 5 && flags[HomoFLAG]){
      printf("Sample size is too small. Robust KING is used.\n");
   }

//      if(runlengthFlag){
//        engine.ComputeRunLength();
//      }else if(!onetrio.IsEmpty()){
//         engine.ComputeRunLengthInOneTrio(onetrio);
//      }else
      if(distantFlag){
         engine.ComputeBigDataDistant();
      }else if(poFlag){
         engine.ComputeBigDataPO();
      }else if(secondFlag){
         engine.ComputeBigDataSecondDegree();
      }else if(sibFlag)
         engine.ComputeShortSibIBS();
/*else if(AssocLMM){
         engine.LMM();
      }*/else if(AssocVC){
         engine.VC();
      }else if(AssocFastAssoc){
         engine.FASTASSOC();
      }else if(AssocGRAMMAR){
         engine.GRAMMAR();
      }else if(AssocRScoreAssoc){
         engine.RSCORE();
      }else if(AssocROADTRIPS){
         engine.ROADTRIPS();
/*      }else if(HeritFlag){
         engine.kinFlag = false;
         engine.HeritFlag = true;
         engine.LMM();
*/
      }else if(AssocADMIX){
         engine.AdmixtureMapping(windowSize, moveSize);
      }else if(AssocHet)
         engine.HetSharing();


#endif

   timer = time(NULL);
   tblock = localtime(&timer);
   printf("KING ends at %s", asctime(tblock));
}



//#include <math.h>
//#include "Kinship.h"
//#include "string.h"
//#include "MerlinSort.h"
/*
         if(lessmemFlag){
            if(!famfile.IsEmpty() || !bimfile.IsEmpty() || !covfile.IsEmpty() || !phefile.IsEmpty())
               engine.ReadPlinkBinaryBigDataWithLessMemory(filename, famfile, bimfile, covfile, phefile);
            else
               engine.ReadPlinkBinaryBigDataWithLessMemory(filename);
         }else{
*/

