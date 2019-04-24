//////////////////////////////////////////////////////////////////////
// qc.cpp
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
// Feb 22, 2019

#include "analysis.h"
#include "Kinship.h"
#include "QuickIndex.h"
#include "MathCholesky.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

void Engine::QC_By_Sample64Bit()
{
   int id1, id2, id3;
   unsigned long long int ibs0, notMissing, MI_trio, notMissing_trio;

   printf("\nOptions in effect:\n");
   printf("\t--bysample\n");
   if(CoreCount)
      printf("\t--cpus %d\n", CoreCount);
   if(lessmemFlag)
      printf("\t--lessmem\n");
   if(Bit64Flag)
      printf("\t--sysbit 64\n");
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");

   printf("QC-by-sample starts at %s", currentTime());

   MakeFamilyForMI();

   IntArray ibs0Count(idCount);
   ibs0Count.Zero();
   IntArray nonmissingCount(idCount);
   nonmissingCount.Zero();
   IntArray MItrioCount(idCount);
   MItrioCount.Zero();
   IntArray nonmissingtrioCount(idCount);
   nonmissingtrioCount.Zero();

   int tempint;
   int LpoCount = (Lpo.Length()>>1);    
   if(LpoCount)
      for(int j = 0; j < idCount; j++)
         ibs0Count[j] = nonmissingCount[j] = 0;
   for(int i = 0; i < LpoCount; i++){
      id1 = geno[Lpo[i*2]];
      id2 = geno[Lpo[i*2+1]];
      for(int b = 0; b < longCount; b++){
         ibs0 = LG[0][id1][b] & LG[0][id2][b] & (LG[1][id1][b] ^ LG[1][id2][b]);
         notMissing = (LG[0][id1][b] | LG[1][id1][b]) & (LG[0][id2][b] | LG[1][id2][b]);
         tempint = popcount(ibs0);
         ibs0Count[id1] += tempint;
         ibs0Count[id2] += tempint;
         tempint = popcount(notMissing);
         nonmissingCount[id1] += tempint;
         nonmissingCount[id2] += tempint;
      }
   }

   Vector ER;
   IntArray removeflag;
   double cutoffER = 0.003; // at the moment
   /*
   if(LpoCount){
      if(errorrateCutoff != _NAN_ && errorrateCutoff > 0)
         cutoffER = errorrateCutoff;
      printf("Error rate is set at %.4lf to determine the MI removal list.\n", cutoffER);
   } */
   if(LpoCount){
      ER.Dimension(idCount);
      ER.Zero();
      for(int i = 0; i < idCount; i++)
         if(nonmissingCount[i])
            ER[i] = ibs0Count[i]*1.0/nonmissingCount[i];

      Vector allerrors(100);
      int allerrorCount = 0;
      removeflag.Dimension(idCount);
      removeflag.Zero();
      int k;
      for(int f = 0; f < ped.familyCount; f++)
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
         allerrorCount = ped[i].sibCount + 2;
         if(allerrorCount < 3 || ped[i].sibs[0]->serial != i) continue;
         allerrors[0] = -1;
         if(ped[i].father){
            k = geno[ped[i].father->serial];
            allerrors[0] = (k>-1)? ER[k]: -1;
         }
         allerrors[1] = -1;
         if(ped[i].mother){
            k = geno[ped[i].mother->serial];
            allerrors[1] = (k>-1)? ER[k]: -1;
         }
         for(int s = 0; s < allerrorCount-2; s++){
            k = geno[ped[i].sibs[s]->serial];
            allerrors[s+2] = (k>-1)? ER[k]: -1;
         }
         k = 0;
         for(int j = 0; j < allerrorCount; j++)
            if(allerrors[j] >= cutoffER) {k = 1; break;}
         if(k == 0) continue; // no errors found in the family
         k = 0;
         for(int j = 0; j < allerrorCount; j++)
            if(allerrors[j] > -1 && allerrors[j] < cutoffER) {k = 1; break;}

         if(k == 1){ // clean family members are found
            if(allerrors[0] > -1 && allerrors[0] < cutoffER){
               if(ped[i].mother && geno[ped[i].mother->serial]>-1)
                  removeflag[geno[ped[i].mother->serial]] = 1;
            }else if(allerrors[1] > -1 && allerrors[1] < cutoffER){
               if(ped[i].father && geno[ped[i].father->serial]>-1)
                  removeflag[geno[ped[i].father->serial]] = 1;
            }else
               for(int s = 0; s < allerrorCount-2; s++){
                  k = geno[ped[i].sibs[s]->serial];
                  if(allerrors[s+2] >= cutoffER)   // to be removed
                     removeflag[k] = 1;
               }
         }else{   // not sure about the errors
            k = 0;
            for(int j = 0; j < allerrorCount; j++)
               if(allerrors[j] >= cutoffER) {k = 1; break;}
            if(k==0) continue; // no errors
            if(ped[i].father){
               k = geno[ped[i].father->serial];
               if(k>-1) removeflag[k] = 2;
            }
            if(ped[i].mother){
               k = geno[ped[i].mother->serial];
               if(k>-1) removeflag[k] = 2;
            }
            for(int s = 0; s < allerrorCount-2; s++){
               k = geno[ped[i].sibs[s]->serial];
               if(k>-1) removeflag[k] = 2;
            }
         }
      }
   }
   int LtrioCount = Ltrio.Length()/3;
   for(int i = 0; i < LtrioCount; i++){
      id1 = geno[Ltrio[i*3+1]];
      id2 = geno[Ltrio[i*3+2]];
      id3 = geno[Ltrio[i*3]];   // id3 is the child
      for(int b = 0; b < longCount; b++){
         MI_trio = LG[0][id1][b] & LG[0][id2][b] & (~(LG[1][id1][b] ^ LG[1][id2][b]))
            & (~LG[0][id3][b]) & LG[1][id3][b];
         notMissing_trio = (LG[0][id1][b] | LG[1][id1][b]) & (LG[0][id2][b] | LG[1][id2][b])
            & (LG[0][id3][b] | LG[1][id3][b]);
         tempint = popcount(MI_trio);
         MItrioCount[id1] += tempint;
         MItrioCount[id2] += tempint;
         MItrioCount[id3] += tempint;
         tempint = popcount(notMissing_trio);
         nonmissingtrioCount[id1] += tempint;
         nonmissingtrioCount[id2] += tempint;
         nonmissingtrioCount[id3] += tempint;
      }
   }
   if(allflags&(1<<ClusterFLAG))
      for(int f = 0; f < ped.familyCount; f++){
         int k = ped[ped.families[f]->first].pid.Find("->");
         if(ped.families[f]->famid.SubStr(0, 4) == "KING" && k > -1)
            for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
               k = ped[i].pid.Find("->"); // a family could have diff family ID
                 ped[i].famid = ped[i].pid.SubStr(0,k);
               ped[i].pid = ped[i].pid.SubStr(k+2);
               if(ped[i].fatid != "0")
                  ped[i].fatid = ped[i].fatid.SubStr(k+2);
               if(ped[i].motid != "0")
                  ped[i].motid = ped[i].motid.SubStr(k+2);
            }
      }

   String filename = prefix;
   filename.Add("bySample.txt");
   FILE *fp = fopen(filename, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)filename);
   fprintf(fp, "FID IID FA MO SEX N_SNP Missing Heterozygosity");
   if(xmarkerCount)
      fprintf(fp, " N_xSNP xHeterozygosity");
   if(ymarkerCount)
      fprintf(fp, " N_ySNP N_yHetero");
   if(mtmarkerCount)
      fprintf(fp, " N_mtSNP N_mtHetero");

   if(LpoCount)
      fprintf(fp, " N_pair N_MIp Err_MIp");
   if(LtrioCount)
      fprintf(fp, " N_trio N_MIt Err_MIt");
   if(LpoCount)
      fprintf(fp, " MI_Removal");
   fprintf(fp, "\n");

   int AaCount, notMissingCount;
   for(int i = 0; i < idCount; i++){
      AaCount = notMissingCount = 0;
      for(int m = 0; m < longCount; m++){
         AaCount += popcount(~LG[0][i][m] & LG[1][i][m]);
         notMissingCount += popcount(LG[0][i][m] | LG[1][i][m]);
      }
      fprintf(fp, "%s %s %s %s %d %d %.4lf %.4lf",
         (const char*)ped[phenoid[i]].famid, (const char*)ped[phenoid[i]].pid,
         (const char*)ped[phenoid[i]].fatid, (const char*)ped[phenoid[i]].motid,
         ped[phenoid[i]].sex, notMissingCount,
         markerCount? (markerCount-notMissingCount) * 1.0 / markerCount: 0,
         notMissingCount? AaCount * 1.0 / notMissingCount: 0);
      if(xmarkerCount){
         AaCount = notMissingCount = 0;
         for(int m = 0; m < xshortCount; m++){
            AaCount += popcount((unsigned short int)(~XG[0][i][m] & XG[1][i][m]));
            notMissingCount += popcount((unsigned short int)(XG[0][i][m] | XG[1][i][m]));
         }
         fprintf(fp, " %d %.4lf",
            notMissingCount,
            notMissingCount ? AaCount * 1.0 / notMissingCount: 0);
      }
      if(ymarkerCount){
         AaCount = notMissingCount = 0;
         for(int m = 0; m < yshortCount; m++){
            AaCount += popcount((unsigned short int)(~YG[0][i][m] & YG[1][i][m]));
            notMissingCount += popcount((unsigned short int)(YG[0][i][m] | YG[1][i][m]));
         }
         fprintf(fp, " %d %d", notMissingCount, AaCount);
      }
      if(mtmarkerCount){
         AaCount = notMissingCount = 0;
         for(int m = 0; m < mtshortCount; m++){
            AaCount += popcount((unsigned short int)(~MG[0][i][m] & MG[1][i][m]));
            notMissingCount += popcount((unsigned short int)(MG[0][i][m] | MG[1][i][m]));
         }
         fprintf(fp, " %d %d", notMissingCount, AaCount);
      }
      // MI
      if(LpoCount)
         fprintf(fp, " %d %d %.4lf",
            nonmissingCount[i], ibs0Count[i], ER[i]);
      if(LtrioCount)
         fprintf(fp, " %d %d %.4lf",
            nonmissingtrioCount[i], MItrioCount[i],
            nonmissingtrioCount[i] ? MItrioCount[i]*1.0/nonmissingtrioCount[i]: 0);
      if(LpoCount) fprintf(fp, " %G", removeflag[i]==2?0.5:removeflag[i]);
      fprintf(fp, "\n");
   }

   fclose(fp);
   printf("  QC-by-sample ends at %s", currentTime());
   printf("QC statistics by samples saved in file %s\n\n", (const char*)filename);
}

void Engine::QC_By_SNP64Bit()
{
   String pedfile = prefix;
   printf("\nOptions in effect:\n");
   printf("\t--bysnp\n");
   if(CoreCount)
      printf("\t--cpus %d\n", CoreCount);
   if(lessmemFlag)
      printf("\t--lessmem\n");
   if(Bit64Flag)
      printf("\t--sysbit 64\n");
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");

   printf("QC-by-SNP starts at %s", currentTime());
   if(!(allflags & (1<<BysampleFLAG)))
      MakeFamilyForMI();
   printf("Scanning autosomes for QC-by-SNP with %d CPU cores...\n", defaultMaxCoreCount);

   int id1, id2, id3;
   int L0Count = (L0.Length()>>1);
   int LpoCount = (Lpo.Length()>>1);
   int LtrioCount = Ltrio.Length()/3;
   IntArray AACounts, AaCounts, missingCounts;
   ComputeAlleleFrequency64Bit(AACounts, AaCounts, missingCounts);
   IntArray nonmissingMZCounts, ibs1MZCounts, HetHetMZCounts, ibs0MZCounts;
   if(L0Count)
      ComputeMZBySNP64Bit(nonmissingMZCounts, ibs1MZCounts, HetHetMZCounts, ibs0MZCounts);
   IntArray ibs0Counts, nonmissingCounts, HomHomCounts;
   if(LpoCount)
      ComputePOBySNP64Bit(HomHomCounts, ibs0Counts, nonmissingCounts);
   IntArray MItrioCounts, nonmissingtrioCounts, HetInOffspringCounts;
   if(LtrioCount)
      ComputeTrioBySNP64Bit(MItrioCounts, HetInOffspringCounts, nonmissingtrioCounts);

   pedfile.Add("bySNP.txt");
   FILE *fp = fopen(pedfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
   fprintf(fp, "SNP");
   if(chromosomes.Length() || xbp.Length()) fprintf(fp, " Chr");
   if(bp.Length() || xbp.Length()) fprintf(fp, " Pos");
   fprintf(fp, " Label_A Label_a Freq_A N N_AA N_Aa N_aa CallRate");
   if(L0Count) fprintf(fp, " N_MZ N_HetMZ N_errMZ Err_InMZ Err_InHetMZ");
   if(LpoCount) fprintf(fp, " N_PO N_HomPO N_errPO Err_InPO Err_InHomPO");
   if(LtrioCount) fprintf(fp, " N_trio N_HetOff N_errTrio Err_InTrio Err_InHetTrio");
   fprintf(fp, "\n");
   for(int m = 0; m < markerCount; m++){
      fprintf(fp, "%s %d %d %s %s %.4lf %d %d %d %d %.4lf",
            (const char*)snpName[m], chromosomes[m], bp[m],
            (const char*)alleleLabel[0][m], (const char*)alleleLabel[1][m],
            missingCounts[m] == idCount? 0:
               (AACounts[m] + AaCounts[m]*0.5) / (idCount - missingCounts[m]),
            idCount - missingCounts[m], AACounts[m], AaCounts[m],
            idCount - missingCounts[m] - AACounts[m] - AaCounts[m],
            1 - missingCounts[m]*1.0/idCount);
      if(L0Count){
         int temp = nonmissingMZCounts[m] - ibs0MZCounts[m] - ibs1MZCounts[m];
         fprintf(fp, " %d %d %d %.4lf %.4lf",
            nonmissingMZCounts[m], ibs1MZCounts[m]+HetHetMZCounts[m], nonmissingMZCounts[m] - temp,
            nonmissingMZCounts[m]? 1 - temp * 1.0 / nonmissingMZCounts[m]: 0,
            ibs1MZCounts[m]+HetHetMZCounts[m]? ibs1MZCounts[m] * 1.0 / (ibs1MZCounts[m] + HetHetMZCounts[m]): 0);
      }
      if(LpoCount)
         fprintf(fp, " %d %d %d %.4lf %.4lf",
            nonmissingCounts[m], HomHomCounts[m], ibs0Counts[m],
            nonmissingCounts[m]? ibs0Counts[m] * 1.0 / nonmissingCounts[m]: 0,
            HomHomCounts[m]? ibs0Counts[m] * 1.0 / HomHomCounts[m]: 0);
      if(LtrioCount)
         fprintf(fp, " %d %d %d %.4lf %.4lf",
            nonmissingtrioCounts[m], HetInOffspringCounts[m], MItrioCounts[m],
            nonmissingtrioCounts[m]? MItrioCounts[m] * 1.0 / nonmissingtrioCounts[m]: 0,
            HetInOffspringCounts[m]? MItrioCounts[m] * 1.0 / HetInOffspringCounts[m]: 0);
      fprintf(fp, "\n");
   }
   if(xmarkerCount) printf("Scanning chromosome X for QC-by-SNP...\n");
   double frequencies[16];
   int missingCount[16], AACount[16], AaCount[16];
   int ibs0Count[16], nonmissingCount[16], MItrioCount[16], nonmissingtrioCount[16];
   int ibs0MZCount[16], ibs1MZCount[16], ibs2MZCount[16], nonmissingMZCount[16];
   int HomHomCount[16], HetHetMZCount[16], HetInOffspringCount[16];
   int AA, Aa, missing;
   unsigned short int HetHet, ibs0, ibs1, notMissing, HomHom, MI_trio, HetInOffspring, T, NT;
   for(int m = 0; m < xmarkerCount; m+=16){
      int b = m/16;
      for(int i = 0; i < 16; i++){
         AACount[i] = AaCount[i] = missingCount[i] = 0;
         frequencies[i] = 0.0;
      }
      for(int i = 0; i < idCount; i++){
         AA = XG[0][i][b] & XG[1][i][b];
         Aa = (~XG[0][i][b]) & XG[1][i][b];
         missing = (~XG[0][i][b]) & (~XG[1][i][b]);
         for(int j = 0; j < 16; j++){
            if(AA & shortbase[j])
               AACount[j] ++;
            else if(Aa & shortbase[j])
               AaCount[j] ++;
            else if(missing & shortbase[j])
               missingCount[j] ++;
         }
      }
      for(int j = 0; j < 16; j++)
         if(missingCount[j] < idCount)
            frequencies[j] = (AACount[j] + AaCount[j]*0.5) / (idCount - missingCount[j]);
      // MZ
      if(L0Count)
         for(int j = 0; j < 16; j++)
            ibs0MZCount[j] = ibs1MZCount[j] = ibs2MZCount[j]
            = nonmissingMZCount[j] = HetHetMZCount[j] = 0;
      for(int i = 0; i < L0Count; i++){
         id1 = geno[L0[i*2]];
         id2 = geno[L0[i*2+1]];
         HetHet = (~XG[0][id1][b]) & XG[1][id1][b] & (~XG[0][id2][b]) & XG[1][id2][b];
         ibs0 = XG[0][id1][b] & XG[0][id2][b] & (XG[1][id1][b]^XG[1][id2][b]);
         ibs1 = ((~XG[0][id1][b]) & XG[1][id1][b] & XG[0][id2][b])
            | ((~XG[0][id2][b]) & XG[1][id2][b] & XG[0][id1][b]);
         notMissing = (XG[0][id1][b] | XG[1][id1][b]) & (XG[0][id2][b] | XG[1][id2][b]);
         for(int j = 0; j < 16; j++)
            if(notMissing & shortbase[j]){
               nonmissingMZCount[j] ++;
               if(ibs0 & shortbase[j])
                  ibs0MZCount[j] ++;
               else if(ibs1 & shortbase[j])
                  ibs1MZCount[j] ++;
               else if(HetHet & shortbase[j])
                  HetHetMZCount[j] ++;
         }
      }
      for(int j = 0; j < 16; j++)
         ibs2MZCount[j] = nonmissingMZCount[j] - ibs0MZCount[j] - ibs1MZCount[j];
      // MI in PO
      if(LpoCount)
         for(int j = 0; j < 16; j++)
            ibs0Count[j] = nonmissingCount[j] = HomHomCount[j] = 0;
      for(int i = 0; i < LpoCount; i++){
         if(ped[Lpo[i*2]].sex==1 && ped[Lpo[i*2+1]].sex==1)
            continue;   // dad does not pass genes to son
         if(ped[Lpo[i*2]].sex + ped[Lpo[i*2+1]].sex == 1)
            continue;   // sex unknown, cannot check
         id1 = geno[Lpo[i*2]];
         id2 = geno[Lpo[i*2+1]];
         ibs0 = XG[0][id1][b] & XG[0][id2][b] & (XG[1][id1][b] ^ XG[1][id2][b]);
         notMissing = (XG[0][id1][b] | XG[1][id1][b]) & (XG[0][id2][b] | XG[1][id2][b]);
         HomHom = XG[0][id1][b] & XG[0][id2][b];
         for(int j = 0; j < 16; j++){
            if(ibs0 & shortbase[j])
               ibs0Count[j] ++;
            if(notMissing & shortbase[j])
               nonmissingCount[j] ++;
            if(HomHom & shortbase[j])
               HomHomCount[j] ++;
         }
      }
      if(LtrioCount)
         for(int j = 0; j < 16; j++)
            MItrioCount[j] = nonmissingtrioCount[j] = HetInOffspringCount[j] = 0;
      for(int i = 0; i < LtrioCount; i++){
         id1 = geno[Ltrio[i*3+1]];
         id2 = geno[Ltrio[i*3+2]];
         id3 = geno[Ltrio[i*3]];   // id3 is the child
         MI_trio = XG[0][id1][b] & XG[0][id2][b] & (~(XG[1][id1][b] ^ XG[1][id2][b]))
            & (~XG[0][id3][b]) & XG[1][id3][b];
         notMissing = (XG[0][id1][b] | XG[1][id1][b]) & (XG[0][id2][b] | XG[1][id2][b])
            & (XG[0][id3][b] | XG[1][id3][b]);
         HetInOffspring = (~XG[0][id3][b]) & XG[1][id3][b] & notMissing;
         for(int j = 0; j < 16; j++){
            if(MI_trio & shortbase[j])
               MItrioCount[j] ++;
            if(notMissing & shortbase[j])
               nonmissingtrioCount[j] ++;
            if(HetInOffspring & shortbase[j])
               HetInOffspringCount[j] ++;
         }
      }
      for(int j = 0; j < 16; j++){
         if(m+j >= xmarkerCount) continue;
         fprintf(fp, "%s X %d %s %s %.4lf %d %d %d %d %.4lf",
            (const char*)xsnpName[m+j], xbp[m+j],
            (const char*)xalleleLabel[0][m+j], (const char*)xalleleLabel[1][m+j],
            frequencies[j],
            idCount - missingCount[j], AACount[j], AaCount[j],
            idCount - missingCount[j] - AACount[j] - AaCount[j],
            1 - missingCount[j]*1.0/idCount);
         if(L0Count)
            fprintf(fp, " %d %d %d %.4lf %.4lf",
               nonmissingMZCount[j], ibs1MZCount[j]+HetHetMZCount[j], nonmissingMZCount[j] - ibs2MZCount[j],
               nonmissingMZCount[j]? 1-ibs2MZCount[j] * 1.0 / nonmissingMZCount[j]: 0,
               ibs1MZCount[j]+HetHetMZCount[j]? ibs1MZCount[j] * 1.0 / (ibs1MZCount[j] + HetHetMZCount[j]): 0);
         if(LpoCount)
            fprintf(fp, " %d %d %d %.4lf %.4lf",
               nonmissingCount[j], HomHomCount[j], ibs0Count[j],
               nonmissingCount[j]? ibs0Count[j] * 1.0 / nonmissingCount[j]: 0,
               HomHomCount[j]? ibs0Count[j] * 1.0 / HomHomCount[j]: 0);
         if(LtrioCount)
            fprintf(fp, " %d %d %d %.4lf %.4lf",
               nonmissingtrioCount[j], HetInOffspringCount[j], MItrioCount[j],
               nonmissingtrioCount[j]? MItrioCount[j] * 1.0 / nonmissingtrioCount[j]: 0,
               HetInOffspringCount[j]? MItrioCount[j] * 1.0 / HetInOffspringCount[j]: 0);
         fprintf(fp, "\n");
      }
   }
   if(ymarkerCount) printf("Scanning chromosome Y for QC-by-SNP...\n");
   for(int m = 0; m < ymarkerCount; m+=16){
      int b = m/16;
      for(int i = 0; i < 16; i++){
         AACount[i] = AaCount[i] = missingCount[i] = 0;
         frequencies[i] = 0.0;
      }
      for(int i = 0; i < idCount; i++){
         AA = YG[0][i][b] & YG[1][i][b];
         Aa = (~YG[0][i][b]) & YG[1][i][b];
         missing = (~YG[0][i][b]) & (~YG[1][i][b]);
         for(int j = 0; j < 16; j++){
            if(AA & shortbase[j])
               AACount[j] ++;
            else if(Aa & shortbase[j])
               AaCount[j] ++;
            else if(missing & shortbase[j])
               missingCount[j] ++;
         }
      }
      for(int j = 0; j < 16; j++)
         if(missingCount[j] < idCount)
            frequencies[j] = (AACount[j] + AaCount[j]*0.5) / (idCount - missingCount[j]);
      // MZ
      if(L0Count)
         for(int j = 0; j < 16; j++)
            ibs0MZCount[j] = ibs1MZCount[j] = ibs2MZCount[j]
            = nonmissingMZCount[j] = HetHetMZCount[j] = 0;
      for(int i = 0; i < L0Count; i++){
         id1 = geno[L0[i*2]];
         id2 = geno[L0[i*2+1]];
         HetHet = (~YG[0][id1][b]) & YG[1][id1][b] & (~YG[0][id2][b]) & YG[1][id2][b];
         ibs0 = YG[0][id1][b] & YG[0][id2][b] & (YG[1][id1][b]^YG[1][id2][b]);
         ibs1 = ((~YG[0][id1][b]) & YG[1][id1][b] & YG[0][id2][b])
            | ((~YG[0][id2][b]) & YG[1][id2][b] & YG[0][id1][b]);
         notMissing = (YG[0][id1][b] | YG[1][id1][b]) & (YG[0][id2][b] | YG[1][id2][b]);
         for(int j = 0; j < 16; j++)
            if(notMissing & shortbase[j]){
               nonmissingMZCount[j] ++;
               if(ibs0 & shortbase[j])
                  ibs0MZCount[j] ++;
               else if(ibs1 & shortbase[j])
                  ibs1MZCount[j] ++;
               else if(HetHet & shortbase[j])
                  HetHetMZCount[j] ++;
            }
      }
      for(int j = 0; j < 16; j++)
         ibs2MZCount[j] = nonmissingMZCount[j] - ibs0MZCount[j] - ibs1MZCount[j];
      // MI in PO
      if(LpoCount)
         for(int j = 0; j < 16; j++)
            ibs0Count[j] = nonmissingCount[j] = HomHomCount[j] = 0;
      for(int i = 0; i < LpoCount; i++){
         id1 = geno[Lpo[i*2]];
         id2 = geno[Lpo[i*2+1]];
         ibs0 = YG[0][id1][b] & YG[0][id2][b] & (YG[1][id1][b] ^ YG[1][id2][b]);
         notMissing = (YG[0][id1][b] | YG[1][id1][b]) & (YG[0][id2][b] | YG[1][id2][b]);
         HomHom = YG[0][id1][b] & YG[0][id2][b];
         for(int j = 0; j < 16; j++){
            if(ibs0 & shortbase[j])
               ibs0Count[j] ++;
            if(notMissing & shortbase[j])
               nonmissingCount[j] ++;
            if(HomHom & shortbase[j])
               HomHomCount[j] ++;
         }
      }
      for(int j = 0; j < 16; j++){
         if(m+j >= ymarkerCount) continue;
         if(ysnpName.Length())
            fprintf(fp, "%s", (const char*)ysnpName[m+j]);
         else
            fprintf(fp, "SNPY%d", m+j+1);
         fprintf(fp, " Y");
         if(ybp.Length())
            fprintf(fp, " %d", ybp[m+j]);
         fprintf(fp, " %s %s %.4lf %d %d %d %d %.4lf",
            (const char*)yalleleLabel[0][m+j], (const char*)yalleleLabel[1][m+j],
            frequencies[j],
            idCount - missingCount[j], AACount[j], AaCount[j],
            idCount - missingCount[j] - AACount[j] - AaCount[j],
            1 - missingCount[j]*1.0/idCount);
         if(L0Count)
            fprintf(fp, " %d %d %d %.4lf %.4lf",
               nonmissingMZCount[j], ibs1MZCount[j]+HetHetMZCount[j], nonmissingMZCount[j] - ibs2MZCount[j],
               nonmissingMZCount[j]? 1-ibs2MZCount[j] * 1.0 / nonmissingMZCount[j]: 0,
               ibs1MZCount[j]+HetHetMZCount[j]? ibs1MZCount[j] * 1.0 / (ibs1MZCount[j] + HetHetMZCount[j]): 0);
         if(LpoCount)
            fprintf(fp, " %d %d %d %.4lf %.4lf",
               nonmissingCount[j], HomHomCount[j], ibs0Count[j],
               nonmissingCount[j]? ibs0Count[j] * 1.0 / nonmissingCount[j]: 0,
               HomHomCount[j]? ibs0Count[j] * 1.0 / HomHomCount[j]: 0);
         if(LtrioCount)
            fprintf(fp, " 0 0 0 0 0");
         fprintf(fp, "\n");
      }
   }
   for(int m = 0; m < mtmarkerCount; m+=16){
      int b = m/16;
      for(int i = 0; i < 16; i++){
         AACount[i] = AaCount[i] = missingCount[i] = 0;
         frequencies[i] = 0.0;
      }
      for(int i = 0; i < idCount; i++){
         AA = MG[0][i][b] & MG[1][i][b];
         Aa = (~MG[0][i][b]) & MG[1][i][b];
         missing = (~MG[0][i][b]) & (~MG[1][i][b]);
         for(int j = 0; j < 16; j++){
            if(AA & shortbase[j])
               AACount[j] ++;
            else if(Aa & shortbase[j])
               AaCount[j] ++;
            else if(missing & shortbase[j])
               missingCount[j] ++;
         }
      }
      for(int j = 0; j < 16; j++)
         if(missingCount[j] < idCount)
            frequencies[j] = (AACount[j] + AaCount[j]*0.5) / (idCount - missingCount[j]);

      // MZ
      if(L0Count)
         for(int j = 0; j < 16; j++)
            ibs0MZCount[j] = ibs1MZCount[j] = ibs2MZCount[j]
            = nonmissingMZCount[j] = HetHetMZCount[j] = 0;
      for(int i = 0; i < L0Count; i++){
         id1 = geno[L0[i*2]];
         id2 = geno[L0[i*2+1]];
         HetHet = (~MG[0][id1][b]) & MG[1][id1][b] & (~MG[0][id2][b]) & MG[1][id2][b];
         ibs0 = MG[0][id1][b] & MG[0][id2][b] & (MG[1][id1][b]^MG[1][id2][b]);
         ibs1 = ((~MG[0][id1][b]) & MG[1][id1][b] & MG[0][id2][b])
            | ((~MG[0][id2][b]) & MG[1][id2][b] & MG[0][id1][b]);
         notMissing = (MG[0][id1][b] | MG[1][id1][b]) & (MG[0][id2][b] | MG[1][id2][b]);
         for(int j = 0; j < 16; j++)
            if(notMissing & shortbase[j]){
               nonmissingMZCount[j] ++;
               if(ibs0 & shortbase[j])
                  ibs0MZCount[j] ++;
               else if(ibs1 & shortbase[j])
                  ibs1MZCount[j] ++;
               else if(HetHet & shortbase[j])
                  HetHetMZCount[j] ++;
            }
      }
      for(int j = 0; j < 16; j++){
         if(m+j >= mtmarkerCount) continue;
         ibs2MZCount[j] = nonmissingMZCount[j] - ibs0MZCount[j] - ibs1MZCount[j];
         if(mtsnpName.Length())
            fprintf(fp, "%s", (const char*)mtsnpName[m+j]);
         else
            fprintf(fp, "SNPMT%d", m+j+1);
         fprintf(fp, " MT");
         if(mtbp.Length())
            fprintf(fp, " %d", mtbp[m+j]);
         fprintf(fp, " %s %s %.4lf %d %d %d %d %.4lf",
            (const char*)mtalleleLabel[0][m+j], (const char*)mtalleleLabel[1][m+j],
            frequencies[j],
            idCount - missingCount[j], AACount[j], AaCount[j],
            idCount - missingCount[j] - AACount[j] - AaCount[j],
            1 - missingCount[j]*1.0/idCount);
         if(L0Count)
            fprintf(fp, " %d %d %d %.4lf %.4lf",
               nonmissingMZCount[j], ibs1MZCount[j]+HetHetMZCount[j], nonmissingMZCount[j] - ibs2MZCount[j],
               nonmissingMZCount[j]? 1-ibs2MZCount[j] * 1.0 / nonmissingMZCount[j]: 0,
               ibs1MZCount[j]+HetHetMZCount[j]? ibs1MZCount[j] * 1.0 / (ibs1MZCount[j] + HetHetMZCount[j]): 0);
         if(LpoCount)
            fprintf(fp, " 0 0 0 0 0");
         if(LtrioCount)
            fprintf(fp, " 0 0 0 0 0");
         fprintf(fp, "\n");
      }
   }
   fclose(fp);
   printf("QC-by-SNP ends at %s", currentTime());
   printf("QC statistics by SNPs saved in file %s\n\n", (const char*)pedfile);
}

void Engine::adjustPC(int pcCount_specified)
{
   printf("Pre-computing allele frequencies...\n");
#ifdef _OPENMP
   printf("  %d CPU cores are used parallelly to pre-compute allele frequencies\n",
      defaultMaxCoreCount);
#endif
   covariatePC.Dimension(0);
   char word[10];
   for(int i = 0; i < pcCount_specified; i++){
      sprintf(word, "PC%d", i+1);
      int k = ped.covariateNames.Find(word);
      if(k>-1)
         covariatePC.Push(k);
      else
         printf("%s is not found\n", word);
   }
   int pcCount = covariatePC.Length();
   if(pcCount > 0)
      printf("  %d PCs are used to predict allele frequencies\n", pcCount);
   else
      error("No PCs are found. Please check the covariate data.");

   Matrix X(idCount, pcCount+1);
   for(int i = 0; i < idCount; i++){
      X[i][0] = 1.0;
      for(int j = 0; j < pcCount; j++)
         X[i][j+1] = ped[phenoid[i]].covariates[covariatePC[j]];
   }
   Matrix XX(pcCount+1, pcCount+1);
   XX.Zero();
   for(int i = 0; i < pcCount+1; i++)
      for(int j = 0; j < pcCount+1; j++)
         for(int k = 0; k < idCount; k++)
            XX[i][j] += X[k][i] * X[k][j];
   Cholesky chol;
   chol.Decompose(XX);
   chol.Invert();
   XX.Dimension(pcCount+1, idCount);
   XX.Zero();
   for(int i = 0; i < pcCount+1; i++)
      for(int j = 0; j < idCount; j++)
         for(int k = 0; k < pcCount+1; k++)
            XX[i][j] += chol.inv[i][k] * X[j][k];

   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];

   int AACount, AaCount, MissingCount;
   const double RAREFREQ = minMAF > 0.05? minMAF: 0.05;
   int wordCount = 1;
   for(; wordCount < shortCount; wordCount <<= 1){
      AACount = AaCount = MissingCount = 0;
      for(int i = 0; i < idCount; i++){
         AACount += oneoneCount[GG[0][i][wordCount] & GG[1][i][wordCount]];
         AaCount += oneoneCount[(~GG[0][i][wordCount]) & GG[1][i][wordCount]];
         MissingCount += oneoneCount[~(GG[0][i][wordCount] | GG[1][i][wordCount]) & 65535];
      }
      if((AACount + AaCount * 0.5) / (idCount*16 - MissingCount) < RAREFREQ) break;
   }
   wordCount >>= 1;
   for(; wordCount < shortCount; wordCount++){
      AACount = AaCount = MissingCount = 0;
      for(int i = 0; i < idCount; i++){
         AACount += oneoneCount[GG[0][i][wordCount] & GG[1][i][wordCount]];
         AaCount += oneoneCount[(~GG[0][i][wordCount]) & GG[1][i][wordCount]];
         MissingCount += oneoneCount[~(GG[0][i][wordCount] | GG[1][i][wordCount]) & 65535];
      }
      if((AACount + AaCount * 0.5) / (idCount*16 - MissingCount) < RAREFREQ) break;
   }
   wordCount--;
   if((wordCount == shortCount) && (markerCount%16!=0)) wordCount--;
   int commonCount = wordCount * 16;
   
   freqBeta.Dimension(commonCount, pcCount+1);
   freqBeta.Zero();
   double freq[16], missingCount[16], temp;
   int AA, Aa, missing;
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(freq, missingCount, AA, Aa, missing, temp)
{
#endif
   int *G2[16], *G1[16], *GM[16];
   for(int m = 0; m < 16; m++){
      G2[m] = new int[idCount+1];
      G1[m] = new int[idCount+1];
      GM[m] = new int[idCount+1];
   }
#ifdef _OPENMP
   #pragma omp for
#endif
   for(int b = 0; b < wordCount; b++){ // loop over short integer
      for(int m = 0; m < 16; m++){
         freq[m] = 0;
         missingCount[m] = 0;
         G2[m][0] = G1[m][0] = GM[m][0] = 0;
      }
      for(int i = 0; i < idCount; i++){
         AA = GG[0][i][b] & GG[1][i][b];
         Aa = (~GG[0][i][b]) & GG[1][i][b];
         missing = (~GG[0][i][b]) & (~GG[1][i][b]) & 65535;
         for(int m = 0; m < 16; m++){
            if(missing & shortbase[m]){
               GM[m][0] ++;
               GM[m][GM[m][0]] = i;
               missingCount[m] ++;
            }else if(AA & shortbase[m]){
               G2[m][0] ++;
               G2[m][G2[m][0]] = i;
               freq[m] += 1.0;
            }else if(Aa & shortbase[m]){
               G1[m][0] ++;
               G1[m][G1[m][0]] = i;
               freq[m] += 0.5;
            }
         }
      }
      for(int m = 0; m < 16; m++){
         int pos = b*16+m;
         if(missingCount[m] >= idCount) continue; // all missing
         freq[m] /= (idCount-missingCount[m]);
         for(int j = 0; j < pcCount+1; j++){
            for(int i = 0; i < G2[m][0]; i++)
               freqBeta[pos][j] += XX[j][G2[m][i+1]];
            temp = 0.0;
            for(int i = 0; i < G1[m][0]; i++)
               temp += XX[j][G1[m][i+1]];
            freqBeta[pos][j] += temp * 0.5;
            temp = 0.0;
            for(int i = 0; i < GM[m][0]; i++)
               temp += XX[j][GM[m][i+1]];
            freqBeta[pos][j] += temp * freq[m];
         }
      }
   }// end of marker word loop
   for(int m = 0; m < 16; m++){
      delete []G2[m];
      delete []G1[m];
      delete []GM[m];
   }
#ifdef _OPENMP
}  // extra bracket for omp
#endif
   BetaSum.Dimension(pcCount+1);
   BetaSum.Zero();
   for(int i = 0; i < pcCount+1; i++){
      double sum = 0.0;
      for(int m = 0; m < commonCount; m++)
         sum += freqBeta[m][i];
      BetaSum[i] = sum;
   }
   BetaSquareSum.Dimension(pcCount+1, pcCount+1);
   BetaSquareSum.Zero();
   for(int i = 0; i < pcCount+1; i++)
      for(int j = 0; j < pcCount+1; j++){
         double sum = 0;
         for(int m = 0; m < commonCount; m++)
            sum += freqBeta[m][i] * freqBeta[m][j];
         BetaSquareSum[i][j] = sum;
   }
}

/*
      for(int i = 0; i < idCount; i++)
         for(int b = blockb; b < bMax; b++){
            int bb = (b-blockb)<<6;
            word = LG[0][i][b] & LG[1][i][b];  // AA
            pchar = (unsigned char*)&word;
            for(int k = 0; k < 8; k++)
               for(byte = pchar[k]; byte; byte &= (byte-1))
                  AACount[bb+(k<<3)+rightmost[byte]] ++;
            word = (~LG[0][i][b]) & LG[1][i][b];  // Aa
            pchar = (unsigned char*)&word;
            for(int k = 0; k < 8; k++)
               for(byte = pchar[k]; byte; byte &= (byte-1))
                  AaCount[bb+(k<<3)+rightmost[byte]] ++;
            word = (~LG[0][i][b]) & (~LG[1][i][b]);  // Miss
            pchar = (unsigned char*)&word;
            for(int k = 0; k < 8; k++)
               for(byte = pchar[k]; byte; byte &= (byte-1))
                  missingCount[bb+(k<<3)+rightmost[byte]] ++;
         }
*/


void Engine::MakeFamilyForMI()
{
   L0.Dimension(0);
   Lpo.Dimension(0);
   Ltrio.Dimension(0);                
   if(ped.haveTwins){
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
            if(!ped[i].zygosity) continue;
            for(int j = i+1; j <= ped.families[f]->last; j++){
               if(ped[j].isTwin(ped[i])){
                  L0.Push(i); L0.Push(j);
                  break;
               }
            }
         }
   }
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
         if(geno[i]==-1) continue;
         if(ped[i].father && geno[ped[i].father->serial]!=-1){
            Lpo.Push(i); Lpo.Push(ped[i].father->serial);
         }
         if(ped[i].mother && geno[ped[i].mother->serial]!=-1){
            Lpo.Push(i); Lpo.Push(ped[i].mother->serial);
         }
         if(ped[i].father && geno[ped[i].father->serial]!=-1
            && ped[i].mother && geno[ped[i].mother->serial]!=-1){
            Ltrio.Push(i);
            Ltrio.Push(ped[i].father->serial);
            Ltrio.Push(ped[i].mother->serial);
         }
      }

   int L0Count = (L0.Length()>>1);
   int LpoCount = (Lpo.Length()>>1);
   int LtrioCount = Ltrio.Length()/3;
   if(L0Count)
      printf("There are %d pairs of duplicates / MZ twins with estimated kinship > 0.49.\n", L0Count);
   printf("There are %d parent-offspring pairs and %d trios according to the pedigree.\n",
      LpoCount, LtrioCount);
   IntArray Lpo1F;
   Kinship kin;
   int id1, id2;

   if(semifamilyFlag){
      if(Bit64==32){
         IntArray connections[5000];
         double kinship;
         char oneoneCount[65536];
         int HetHetCount, IBS0Count, het1Count, notMissingCount;
         for(int i = 0; i < 65536; i++)
            oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
         printf("1st-degree relatives are treated as parent-offspring if IBS0 < %.4lf\n",
            errorrateCutoff);
         for(int f = 0; f < ped.familyCount; f++){
            if(id[f].Length() < 2) continue;
            kin.Setup(*ped.families[f]);
            Lpo1F.Dimension(0);
            for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
               if(geno[i]==-1) continue;
               if(ped[i].father && geno[ped[i].father->serial]!=-1
                  && ped[i].mother && geno[ped[i].mother->serial]!=-1);
               else if(ped[i].father && geno[ped[i].father->serial]!=-1){
                  Lpo1F.Push(i); Lpo1F.Push(ped[i].father->serial);
               }else if(ped[i].mother && geno[ped[i].mother->serial]!=-1){
                  Lpo1F.Push(i); Lpo1F.Push(ped[i].mother->serial);
               }
            }
            for(int i = 0; i < id[f].Length(); i++){
               id1 = geno[id[f][i]];
               for(int j = i+1; j < id[f].Length(); j++){
                  id2 = geno[id[f][j]];
                  if(kin(ped[id[f][i]], ped[id[f][j]]) > 0.005) continue;
                  HetHetCount = IBS0Count = het1Count = 0;
                  for(int m = 0; m < shortCount; m++){
                     HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
                     IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
                  }
                  if(HetHetCount < IBS0Count*2) continue;   // unrelated
                  for(int m = 0; m < shortCount; m++)
                     het1Count += (oneoneCount[(GG[0][id2][m] | GG[1][id2][m]) & (~GG[0][id1][m]) & GG[1][id1][m]]
                        + oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]]);
                  kinship = het1Count? (HetHetCount - IBS0Count*2)*1.0/het1Count: 0;
                  if(kinship > 0.49){  // MZ
                     L0.Push(id[f][i]); L0.Push(id[f][j]);
                  }else if(kinship > 0.15){ // PO
                     notMissingCount = 0;
                     for(int m = 0; m < shortCount; m++)
                        notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
                     if(IBS0Count < notMissingCount * errorrateCutoff){
                        Lpo.Push(id[f][i]); Lpo.Push(id[f][j]);
                        Lpo1F.Push(id[f][i]); Lpo1F.Push(id[f][j]);
                     }
                  }
               }
            }
            for(int i = 0; i < ped.families[f]->count; i++) connections[i].Dimension(0);
            for(int i = 0; i < Lpo1F.Length()/2; i++){
               connections[Lpo1F[i*2]-ped.families[f]->first].Push(Lpo1F[i*2+1]);
               connections[Lpo1F[i*2+1]-ped.families[f]->first].Push(Lpo1F[i*2]);
            }
            for(int i = 0; i < ped.families[f]->count; i++)
               if(connections[i].Length() >= 2)
                  for(int j = 0; j < connections[i].Length(); j++){
                     id1 = geno[connections[i][j]];
                     for(int k = j+1; k < connections[i].Length(); k++){
                        id2 = geno[connections[i][k]];
                        HetHetCount = IBS0Count = het1Count = 0;
                        for(int m = 0; m < shortCount; m++){
                           HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
                           het1Count += (oneoneCount[(GG[0][id2][m] | GG[1][id2][m]) & (~GG[0][id1][m]) & GG[1][id1][m]]
                              + oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]]);
                           IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
                        }
                        if(het1Count){
                           kinship = (HetHetCount - IBS0Count*2)*1.0/het1Count;
                           if(kinship < 0.0625){   // not 2nd-degree, likely unrelated
                              Ltrio.Push(ped.families[f]->first+i);
                              Ltrio.Push(connections[i][j]);
                              Ltrio.Push(connections[i][k]);
                              j = k = connections[i].Length();
                              break;
                           }
                        }
                     }  // k loop ends
                  }  // i, j loop ends
         }  // f loop ends
      }else{// Bit64
//         IntArray chrSeg;
//         double totalLength;
         bool IBDvalidFlag = false;
//         String segmessage;
         IBDvalidFlag = PreSegment(/*chrSeg, totalLength, segmessage, */false);
         if(!IBDvalidFlag) return;

         int maxcount = 0;
         for(int f = 0; f < ped.familyCount; f++)
            if(id[f].Length()>maxcount) maxcount = id[f].Length();
         IntArray *connections = new IntArray[maxcount];
         IntArray allpairs;
         Matrix pi;
         IntArray inverse;
         Vector ibdprops, maxLengths, ibd2props, maxLengths2;
         for(int f = 0; f < ped.familyCount; f++){
            int idcount = id[f].Length(); // idcount is not idCount, temporary only
            if(idcount < 2) continue;
            allpairs.Dimension(0);
            inverse.Dimension(ped.families[f]->count);
            for(int i = 0; i < idcount; i++)
               inverse[id[f][i] - ped.families[f]->first] = i;
            for(int i = 0; i < idcount; i++){
               id1 = geno[id[f][i]];
               for(int j = i+1; j < idcount; j++){
                  id2 = geno[id[f][j]];
                  allpairs.Push(id1);
                  allpairs.Push(id2);
               }
            }
            IBDSegInSubset64Bit(allpairs, ibdprops, maxLengths, ibd2props, maxLengths2);
            pi.Dimension(idcount,idcount);
            int count = 0;
            for(int i = 0; i < idcount; i++)
               for(int j = i+1; j < idcount; j++){
                  pi[i][j] = ibdprops[count];
                  pi[j][i] = ibd2props[count];
                  count ++;
               }
            Lpo1F.Dimension(0);
            for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
               if(ped[i].father && geno[ped[i].father->serial]!=-1
                  && ped[i].mother && geno[ped[i].mother->serial]!=-1);
               else if(ped[i].father && geno[ped[i].father->serial]!=-1){
                  Lpo1F.Push(inverse[i-ped.families[f]->first]);
                  Lpo1F.Push(inverse[ped[i].father->serial-ped.families[f]->first]);
               }else if(ped[i].mother && geno[ped[i].mother->serial]!=-1){
                  Lpo1F.Push(inverse[i-ped.families[f]->first]);
                  Lpo1F.Push(inverse[ped[i].mother->serial-ped.families[f]->first]);
               }
            kin.Setup(*ped.families[f]);
            for(int i = 0; i < idcount; i++)
               for(int j = i+1; j < idcount; j++){
                  if(kin(ped[id[f][i]], ped[id[f][j]]) > 0.2) continue;
                  if(pi[j][i] > 0.7){ // MZ
                     L0.Push(id[f][i]); L0.Push(id[f][j]);
                  }else if((pi[i][j]+pi[j][i] > 0.96) || (pi[i][j]+pi[j][i]>0.9 && pi[j][i]<0.08)){ // PO
                     Lpo.Push(id[f][i]); Lpo.Push(id[f][j]);
                     Lpo1F.Push(i); Lpo1F.Push(j);
                  }
               }
            for(int i = 0; i < idcount; i++) connections[i].Dimension(0);
            for(int i = 0; i < Lpo1F.Length()/2; i++){
               connections[Lpo1F[i*2]].Push(Lpo1F[i*2+1]);
               connections[Lpo1F[i*2+1]].Push(Lpo1F[i*2]);
            }
            for(int i = 0; i < idcount; i++){
               if(connections[i].Length() < 2) continue;
               for(int j = 0; j < connections[i].Length(); j++)
                  for(int k = j+1; k < connections[i].Length(); k++){
                     if(connections[i][j] < connections[i][k]){
                        id1 = connections[i][j];
                        id2 = connections[i][k];
                     }else{
                        id1 = connections[i][k];
                        id2 = connections[i][j];
                     }
                     if(pi[id1][id2]*0.5 + pi[id2][id1] < 0.125){   // > 2nd-deg
                        Ltrio.Push(id[f][i]);
                        Ltrio.Push(id[f][id1]);
                        Ltrio.Push(id[f][id2]);
                        j = k = connections[i].Length();
                        break;
                     }
                  }
            }  // i loop ends
         }  // f loop ends
         delete []connections;
      }  // if Bit64 ends
      if(Lpo.Length() > LpoCount*2 || Ltrio.Length() > LtrioCount*3)
         printf("There are additional %d parent-offspring pairs and %d trios by inference.\n",
            Lpo.Length()/2 - LpoCount, Ltrio.Length()/3 - LtrioCount);
      if(L0.Length())
         printf("There are %d pairs of duplicates in total.\n", L0.Length()/2);
   }
}

void Engine::QC_WipeMI()
{
   int id1, id2, id3;
   int ibs0, MI_trio, wipe;
   char oneoneCount[65536];
   int trioWipeCount=0;
   int poWipeCount=0;

//   if(!QCbySample && !QCbySNP) MakeFamilyForMI();
   if(!(allflags&(1<<BysampleFLAG)) && !(allflags&(1<<BysnpFLAG)))
      MakeFamilyForMI();
   if(Lpo.Length() || Ltrio.Length())
      for(int i = 0; i < 65536; i++)
         oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   for(int i = 0; i < Ltrio.Length()/3; i++){
      id1 = geno[Ltrio[i*3+1]];
      id2 = geno[Ltrio[i*3+2]];
      id3 = geno[Ltrio[i*3]];   // id3 is the child
      for(int b = 0; b < shortCount; b++){
         MI_trio = GG[0][id1][b] & GG[0][id2][b] & (~(GG[1][id1][b] ^ GG[1][id2][b]))
            & (~GG[0][id3][b]) & GG[1][id3][b];
         if(MI_trio){// MI
            wipe = ~MI_trio;
            for(int j = 0; j < 2; j++){
               GG[j][id1][b] &= wipe;
               GG[j][id2][b] &= wipe;
               GG[j][id3][b] &= wipe;
            }
            trioWipeCount += oneoneCount[MI_trio];
         }
      }
   }
   for(int i = 0; i < Lpo.Length()/2; i++){
      id1 = geno[Lpo[i*2]];
      id2 = geno[Lpo[i*2+1]];
      for(int b = 0; b < shortCount; b++){
         ibs0 = GG[0][id1][b] & GG[0][id2][b] & (GG[1][id1][b] ^ GG[1][id2][b]);
         if(ibs0){
            wipe = ~ibs0;
            for(int j = 0; j < 2; j++){
               GG[j][id1][b] &= wipe;
               GG[j][id2][b] &= wipe;
            }
            poWipeCount += oneoneCount[ibs0];
         }
      }
   }
   if(trioWipeCount==0 && poWipeCount==0){
      if(shortCount) printf("No MI errors are identified on autosomes.\n");
//      if(xshortCount==0) return;
   }else{
      printf("%d SNPs are zeroed out in %d parent-offspring trios.\n",
         trioWipeCount, Ltrio.Length()/3);
      printf("%d SNPs are zeroed out in %d parent-offspring pairs.\n",
         poWipeCount, Lpo.Length()/2);
   }
   trioWipeCount=poWipeCount=0;

   for(int i = 0; i < Ltrio.Length()/3; i++){
      id1 = geno[Ltrio[i*3+1]];
      id2 = geno[Ltrio[i*3+2]];
      id3 = geno[Ltrio[i*3]];   // id3 is the child
      for(int b = 0; b < xshortCount; b++){
         MI_trio = XG[0][id1][b] & XG[0][id2][b] & (~(XG[1][id1][b] ^ XG[1][id2][b]))
            & (~XG[0][id3][b]) & XG[1][id3][b];
         if(MI_trio){// MI
            wipe = ~MI_trio;
            for(int j = 0; j < 2; j++){
               XG[j][id1][b] &= wipe;
               XG[j][id2][b] &= wipe;
               XG[j][id3][b] &= wipe;
            }
            trioWipeCount += oneoneCount[MI_trio];
         }
      }
   }
   for(int i = 0; i < Lpo.Length()/2; i++){
      if(ped[Lpo[i*2]].sex==1 && ped[Lpo[i*2+1]].sex==1)
         continue;   // dad does not pass genes to son
      if(ped[Lpo[i*2]].sex + ped[Lpo[i*2+1]].sex == 1)
         continue;   // sex unknown, cannot check
      id1 = geno[Lpo[i*2]];
      id2 = geno[Lpo[i*2+1]];
      for(int b = 0; b < xshortCount; b++){
         ibs0 = XG[0][id1][b] & XG[0][id2][b] & (XG[1][id1][b] ^ XG[1][id2][b]);
         if(ibs0){
            wipe = ~ibs0;
            for(int j = 0; j < 2; j++){
               XG[j][id1][b] &= wipe;
               XG[j][id2][b] &= wipe;
            }
            poWipeCount += oneoneCount[ibs0];
         }
      }
   }
   if(xshortCount){
      if(trioWipeCount==0 && poWipeCount==0){
         printf("No MI errors are identified at X-chromosome SNPs.\n");
//         return;
      }else{
         printf("%d X-chromosome SNPs are zeroed out in %d parent-offspring trios.\n",
            trioWipeCount, Ltrio.Length()/3);
         printf("%d X-chromosome SNPs are zeroed out in %d parent-offspring pairs.\n",
            poWipeCount, Lpo.Length()/2);
      }
   }

   poWipeCount=0;
   for(int i = 0; i < Lpo.Length()/2; i++){
      id1 = geno[Lpo[i*2]];
      id2 = geno[Lpo[i*2+1]];
      for(int b = 0; b < yshortCount; b++){
         ibs0 = YG[0][id1][b] & YG[0][id2][b] & (YG[1][id1][b] ^ YG[1][id2][b]);
         if(ibs0){
            wipe = ~ibs0;
            for(int j = 0; j < 2; j++){
               YG[j][id1][b] &= wipe;
               YG[j][id2][b] &= wipe;
            }
            poWipeCount += oneoneCount[ibs0];
         }
      }
   }
   if(yshortCount){
      if(poWipeCount==0)
         printf("No MI errors are identified at Y-chromosome SNPs.\n");
      else
         printf("%d Y-chromosome SNPs are zeroed out in %d parent-offspring pairs.\n",
            poWipeCount, Lpo.Length()/2);
   }

   String filename = prefix;
   if(SaveFormat == "PLINK"){
      filename.Add("wipe");
      WritePlinkBinary(filename);
   }else if(SaveFormat == "MERLIN")
      WriteMerlin();
   else{
      filename.Add("wipe.bgeno");
      DumpBinary(filename);
   }
}

void Engine::QC_By_Sample()
{
   int id1, id2, id3;
   int ibs0, notMissing, MI_trio, notMissing_trio;

   printf("\nOptions in effect:\n");
   printf("\t--bysample\n");
   if(CoreCount)
      printf("\t--cpus %d\n", CoreCount);
   if(lessmemFlag)
      printf("\t--lessmem\n");
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");

   printf("QC-by-sample starts at %s", currentTime());

   countGenotype();
   MakeFamilyForMI();

   char oneoneCount[65536];
   if(Lpo.Length() || Ltrio.Length())
      for(int i = 0; i < 65536; i++)
         oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   IntArray ibs0Count(idCount);
   ibs0Count.Zero();
   IntArray nonmissingCount(idCount);
   nonmissingCount.Zero();
   IntArray MItrioCount(idCount);
   MItrioCount.Zero();
   IntArray nonmissingtrioCount(idCount);
   nonmissingtrioCount.Zero();

   if(Lpo.Length())
      for(int j = 0; j < idCount; j++)
         ibs0Count[j] = nonmissingCount[j] = 0;
   for(int i = 0; i < Lpo.Length()/2; i++){
      id1 = geno[Lpo[i*2]];
      id2 = geno[Lpo[i*2+1]];
      for(int b = 0; b < shortCount; b++){
         ibs0 = GG[0][id1][b] & GG[0][id2][b] & (GG[1][id1][b] ^ GG[1][id2][b]);
         notMissing = (GG[0][id1][b] | GG[1][id1][b]) & (GG[0][id2][b] | GG[1][id2][b]);
         ibs0Count[id1] += oneoneCount[ibs0];
         ibs0Count[id2] += oneoneCount[ibs0];
         nonmissingCount[id1] += oneoneCount[notMissing];
         nonmissingCount[id2] += oneoneCount[notMissing];
      }
   }

   Vector ER;
   IntArray removeflag;
   double cutoffER = 0.01; // at the moment
   if(Lpo.Length()){
      if(errorrateCutoff != _NAN_ && errorrateCutoff > 0)
         cutoffER = errorrateCutoff;
      printf("Error rate is set at %.4lf to determine the MI removal list.\n", cutoffER);
   }
   if(Lpo.Length()){
      ER.Dimension(idCount);
      ER.Zero();
      for(int i = 0; i < idCount; i++)
         if(nonmissingCount[i])
            ER[i] = ibs0Count[i]*1.0/nonmissingCount[i];

      Vector allerrors(100);
      int allerrorCount = 0;
      removeflag.Dimension(idCount);
      removeflag.Zero();
      int k;
      for(int f = 0; f < ped.familyCount; f++)
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
         allerrorCount = ped[i].sibCount + 2;
         if(allerrorCount < 3 || ped[i].sibs[0]->serial != i) continue;
         allerrors[0] = -1;
         if(ped[i].father){
            k = geno[ped[i].father->serial];
            allerrors[0] = (k>-1)? ER[k]: -1;
         }
         allerrors[1] = -1;
         if(ped[i].mother){
            k = geno[ped[i].mother->serial];
            allerrors[1] = (k>-1)? ER[k]: -1;
         }
         for(int s = 0; s < allerrorCount-2; s++){
            k = geno[ped[i].sibs[s]->serial];
            allerrors[s+2] = (k>-1)? ER[k]: -1;
         }
         k = 0;
         for(int j = 0; j < allerrorCount; j++)
            if(allerrors[j] >= cutoffER) {k = 1; break;}
         if(k == 0) continue; // no errors found in the family
         k = 0;
         for(int j = 0; j < allerrorCount; j++)
            if(allerrors[j] > -1 && allerrors[j] < cutoffER) {k = 1; break;}

         if(k == 1){ // clean family members are found
            if(allerrors[0] > -1 && allerrors[0] < cutoffER){
               if(ped[i].mother && geno[ped[i].mother->serial]>-1)
                  removeflag[geno[ped[i].mother->serial]] = 1;
            }else if(allerrors[1] > -1 && allerrors[1] < cutoffER){
               if(ped[i].father && geno[ped[i].father->serial]>-1)
                  removeflag[geno[ped[i].father->serial]] = 1;
            }else
               for(int s = 0; s < allerrorCount-2; s++){
                  k = geno[ped[i].sibs[s]->serial];
                  if(allerrors[s+2] >= cutoffER)   // to be removed
                     removeflag[k] = 1;
               }
         }else{   // not sure about the errors
            k = 0;
            for(int j = 0; j < allerrorCount; j++)
               if(allerrors[j] >= cutoffER) {k = 1; break;}
            if(k==0) continue; // no errors
            if(ped[i].father){
               k = geno[ped[i].father->serial];
               if(k>-1) removeflag[k] = 2;
            }
            if(ped[i].mother){
               k = geno[ped[i].mother->serial];
               if(k>-1) removeflag[k] = 2;
            }
            for(int s = 0; s < allerrorCount-2; s++){
               k = geno[ped[i].sibs[s]->serial];
               if(k>-1) removeflag[k] = 2;
            }
         }
      }
   }
   for(int i = 0; i < Ltrio.Length()/3; i++){
      id1 = geno[Ltrio[i*3+1]];
      id2 = geno[Ltrio[i*3+2]];
      id3 = geno[Ltrio[i*3]];   // id3 is the child
      for(int b = 0; b < shortCount; b++){
         MI_trio = GG[0][id1][b] & GG[0][id2][b] & (~(GG[1][id1][b] ^ GG[1][id2][b]))
            & (~GG[0][id3][b]) & GG[1][id3][b];
         notMissing_trio = (GG[0][id1][b] | GG[1][id1][b]) & (GG[0][id2][b] | GG[1][id2][b])
            & (GG[0][id3][b] | GG[1][id3][b]);
         MItrioCount[id1] += oneoneCount[MI_trio];
         MItrioCount[id2] += oneoneCount[MI_trio];
         MItrioCount[id3] += oneoneCount[MI_trio];
         nonmissingtrioCount[id1] += oneoneCount[notMissing_trio];
         nonmissingtrioCount[id2] += oneoneCount[notMissing_trio];
         nonmissingtrioCount[id3] += oneoneCount[notMissing_trio];
      }
   }
   if(allflags&(1<<ClusterFLAG))
      for(int f = 0; f < ped.familyCount; f++){
         int k = ped[ped.families[f]->first].pid.Find("->");
         if(ped.families[f]->famid.SubStr(0, 4) == "KING" && k > -1)
            for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
               k = ped[i].pid.Find("->"); // a family could have diff family ID
                 ped[i].famid = ped[i].pid.SubStr(0,k);
               ped[i].pid = ped[i].pid.SubStr(k+2);
               if(ped[i].fatid != "0")
                  ped[i].fatid = ped[i].fatid.SubStr(k+2);
               if(ped[i].motid != "0")
                  ped[i].motid = ped[i].motid.SubStr(k+2);
            }
      }

   String filename = prefix;
   filename.Add("bySample.txt");
   FILE *fp = fopen(filename, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)filename);
   fprintf(fp, "FID IID FA MO SEX N_SNP Missing Heterozygosity");
   if(xmarkerCount)
      fprintf(fp, " N_xSNP xHeterozygosity");
   if(ymarkerCount)
      fprintf(fp, " N_ySNP N_yHetero");
   if(mtmarkerCount)
      fprintf(fp, " N_mtSNP N_mtHetero");

   if(Lpo.Length())
      fprintf(fp, " N_pair N_MIp Err_MIp");
   if(Ltrio.Length())
      fprintf(fp, " N_trio N_MIt Err_MIt");
   if(Lpo.Length())
      fprintf(fp, " MI_Removal");
   fprintf(fp, "\n");

   for(int i = 0; i < ped.count; i++){
      int k = geno[i];
      if(k == -1) continue;
      int count = pAACount[k] + pAaCount[k] + paaCount[k];
      double mean = (pAACount[k] + pAaCount[k]*0.5) / count;
      double var = pAaCount[k];
      var = mean * (1-mean) - var / count / 2;
      fprintf(fp, "%s %s %s %s %d %d %.4lf %.4lf",
         (const char*)ped[i].famid, (const char*)ped[i].pid,
         (const char*)ped[i].fatid, (const char*)ped[i].motid,
         ped[i].sex, count,
         markerCount? (markerCount-count) * 1.0 / markerCount: 0,
         count? pAaCount[k] * 1.0 / count: 0);
      if(xmarkerCount){
         int xcount = pxAACount[k] + pxAaCount[k] + pxaaCount[k];
         double xhet = xcount? pxAaCount[k] * 1.0 / xcount: 0;
         fprintf(fp, " %d %.4lf", xcount, xhet);
      }
      if(ymarkerCount){
         int ycount = pyAACount[k] + pyAaCount[k] + pyaaCount[k];
         fprintf(fp, " %d %d", ycount, pyAaCount[k]);
      }
      if(mtmarkerCount){
         int mtcount = pmtAACount[k] + pmtAaCount[k] + pmtaaCount[k];
         fprintf(fp, " %d %d", mtcount, pmtAaCount[k]);
      }
      // MI
      if(Lpo.Length())
         fprintf(fp, " %d %d %.4lf",
            nonmissingCount[k], ibs0Count[k], ER[k]);
      if(Ltrio.Length())
         fprintf(fp, " %d %d %.4lf",
            nonmissingtrioCount[k], MItrioCount[k],
            nonmissingtrioCount[k] ? MItrioCount[k]*1.0/nonmissingtrioCount[k]: 0);
      if(Lpo.Length()) fprintf(fp, " %G", removeflag[k]==2?0.5:removeflag[k]);
      fprintf(fp, "\n");
   }
   fclose(fp);
   printf("  QC-by-sample ends at %s", currentTime());
   printf("QC statistics by samples saved in file %s\n\n", (const char*)filename);
}

void Engine::QC_By_SNP()
{
   String pedfile = prefix;
   int id1, id2, id3;
   int informative;
   printf("\nOptions in effect:\n");
   printf("\t--bysnp\n");
   if(CoreCount)
      printf("\t--cpus %d\n", CoreCount);
   if(lessmemFlag)
      printf("\t--lessmem\n");
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");

   printf("QC-by-SNP starts at %s", currentTime());
   if(!(allflags & (1<<BysampleFLAG)))
      MakeFamilyForMI();

   printf("Scanning autosomes for QC-by-SNP with %d CPU cores...\n", defaultMaxCoreCount);

   char buffer[1024];
   StringArray buffers;
   buffers.Dimension(shortCount);
   for(int i = 0; i < buffers.Length(); i++)
      buffers[i].Clear();

   const int CACHESIZE=256;
   const int CACHESIZE_SNP=CACHESIZE<<4;  // 4096
   int missingCount[CACHESIZE_SNP], AACount[CACHESIZE_SNP], AaCount[CACHESIZE_SNP];
   int ibs0Count[CACHESIZE_SNP], nonmissingCount[CACHESIZE_SNP], MItrioCount[CACHESIZE_SNP], nonmissingtrioCount[CACHESIZE_SNP];
   int ibs0MZCount[CACHESIZE_SNP], ibs1MZCount[CACHESIZE_SNP], ibs2MZCount[CACHESIZE_SNP], nonmissingMZCount[CACHESIZE_SNP];
   int HomHomCount[CACHESIZE_SNP], HetHetMZCount[CACHESIZE_SNP], HetInOffspringCount[CACHESIZE_SNP];
   int TCount[CACHESIZE_SNP], NTCount[CACHESIZE_SNP], informativeCount[CACHESIZE_SNP];;

   char revbase[65536];
   for(int i = 0; i < 16; i++)
      revbase[shortbase[i]] = i;
   char rightmost[65536];
   for(int i = 0; i < 65536; i++)
      rightmost[i] = revbase[i&(-i)];
   unsigned short int word;
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
      private(word, informative, id1, id2, id3, buffer, \
      AACount, AaCount, missingCount, \
      ibs0MZCount, ibs1MZCount, nonmissingMZCount, HetHetMZCount, \
      ibs0Count, nonmissingCount, HomHomCount, \
      MItrioCount, nonmissingtrioCount, HetInOffspringCount, \
      TCount, NTCount, informativeCount)
#endif
   for(int blockb = 0; blockb < shortCount; blockb += CACHESIZE){
      int bMax = (blockb > shortCount-CACHESIZE) ? shortCount: blockb+CACHESIZE;
      for(int j = 0; j < CACHESIZE_SNP; j++)
         AACount[j] = AaCount[j] = missingCount[j] = 0;
      for(int i = 0; i < idCount; i++)
         for(int b = blockb; b < bMax; b++){
            int bb = (b-blockb)*16;
            for(word = GG[0][i][b] & GG[1][i][b]; word; word &= (word-1))
               AACount[bb+rightmost[word]] ++;
            for(word = (~GG[0][i][b]) & GG[1][i][b]; word; word &= (word-1))
               AaCount[bb+rightmost[word]] ++;
            for(word = (~GG[0][i][b]) & (~GG[1][i][b]) & 65535; word; word &= (word-1))
               missingCount[bb+rightmost[word]] ++;
         }
      // not concordant in MZ twins or duplicates
      if(L0.Length())
         for(int j = 0; j < CACHESIZE_SNP; j++)
            ibs0MZCount[j] = ibs1MZCount[j] = nonmissingMZCount[j] = HetHetMZCount[j] = 0;
      for(int i = 0; i < L0.Length()/2; i++){
         id1 = geno[L0[i*2]];
         id2 = geno[L0[i*2+1]];
         for(int b = blockb; b < bMax; b++){
            int bb = (b-blockb)*16;
            for(word = (GG[0][id1][b] | GG[1][id1][b]) & (GG[0][id2][b] | GG[1][id2][b]);
               word; word &= (word-1)) // nonmissing
               nonmissingMZCount[bb+rightmost[word]] ++;
            for(word = (~GG[0][id1][b]) & GG[1][id1][b] & (~GG[0][id2][b]) & GG[1][id2][b];
               word; word &= (word-1)) // HetHet
               HetHetMZCount[bb+rightmost[word]] ++;
            for(word = GG[0][id1][b] & GG[0][id2][b] & (GG[1][id1][b]^GG[1][id2][b]);
               word; word &= (word-1)) // IBS0
               ibs0MZCount[bb+rightmost[word]] ++;
            for(word = ((~GG[0][id1][b]) & GG[1][id1][b] & GG[0][id2][b])
               | ((~GG[0][id2][b]) & GG[1][id2][b] & GG[0][id1][b]);
               word; word &= (word-1)) // IBS1
               ibs1MZCount[bb+rightmost[word]] ++;
         }
      }
      // MI in PO
      if(Lpo.Length())
         for(int j = 0; j < CACHESIZE_SNP; j++)
            ibs0Count[j] = nonmissingCount[j] = HomHomCount[j] = 0;
      for(int i = 0; i < Lpo.Length()/2; i++){
         id1 = geno[Lpo[i*2]];
         id2 = geno[Lpo[i*2+1]];
         for(int b = blockb; b < bMax; b++){
            int bb = (b-blockb)*16;
            for(word = GG[0][id1][b] & GG[0][id2][b] & (GG[1][id1][b]^GG[1][id2][b]);
               word; word &= (word-1)) // IBS0
               ibs0Count[bb+rightmost[word]] ++;
            for(word = (GG[0][id1][b] | GG[1][id1][b]) & (GG[0][id2][b] | GG[1][id2][b]);
               word; word &= (word-1)) // nonmissing
               nonmissingCount[bb+rightmost[word]] ++;
            for(word = GG[0][id1][b] & GG[0][id2][b];
               word; word &= (word-1)) // HomHom
               HomHomCount[bb+rightmost[word]] ++;
         }
      }
      // MI in Trios
      if(Ltrio.Length())
         for(int j = 0; j < CACHESIZE_SNP; j++)
            MItrioCount[j] = nonmissingtrioCount[j] = HetInOffspringCount[j] =
               TCount[j] = NTCount[j] = informativeCount[j] = 0;
      for(int i = 0; i < Ltrio.Length()/3; i++){
         id1 = geno[Ltrio[i*3+1]];
         id2 = geno[Ltrio[i*3+2]];
         id3 = geno[Ltrio[i*3]];   // id3 is the child
         for(int b = blockb; b < bMax; b++){
            int bb = (b-blockb)*16;
            for(word = GG[0][id1][b] & GG[0][id2][b] & (~(GG[1][id1][b] ^ GG[1][id2][b]))
            & (~GG[0][id3][b]) & GG[1][id3][b];   // AA x AA -> Aa, or aa x aa -> Aa
               word; word &= (word-1)) // MI_Trio
               MItrioCount[bb+rightmost[word]] ++;
            for(word = (GG[0][id1][b] | GG[1][id1][b]) & (GG[0][id2][b] | GG[1][id2][b])
            & (GG[0][id3][b] | GG[1][id3][b]);
               word; word &= (word-1)) // nonmissing
               nonmissingtrioCount[bb+rightmost[word]] ++;
            for(word = (~GG[0][id3][b]) & GG[1][id3][b] & (GG[0][id1][b] | GG[1][id1][b]) & (GG[0][id2][b] | GG[1][id2][b]);
               word; word &= (word-1)) // HetInOffspring
               HetInOffspringCount[bb+rightmost[word]] ++;
            // Transmission Disequilibrium
            informative = (GG[0][id3][b] | GG[1][id3][b]) & // offspring not missing
            ( ((~GG[0][id1][b]) & GG[1][id1][b] & (GG[0][id2][b] | GG[1][id2][b])) |
              ((~GG[0][id2][b]) & GG[1][id2][b] & (GG[0][id1][b] | GG[1][id1][b])) );
            for(word = informative; word; word &= (word-1))
               informativeCount[bb+rightmost[word]] ++;
            for(word = (~GG[0][id1][b]) & GG[1][id1][b] &
               (~GG[0][id2][b]) & GG[1][id2][b] &
               (~GG[0][id3][b]) & GG[1][id3][b];
               word; word &= (word-1)){
               TCount[bb+rightmost[word]] ++;
               NTCount[bb+rightmost[word]] ++;
            }
            for(word = informative &   // Aa->AA or aa->Aa
              ( (~GG[0][id1][b] & GG[1][id1][b] & GG[0][id3][b] & GG[1][id3][b]) |
                (GG[0][id1][b] & ~GG[1][id1][b] & ~GG[0][id3][b] & GG[1][id3][b]) );
               word; word &= (word-1))
               TCount[bb+rightmost[word]] ++;
            for(word = informative &   // Aa->AA or aa->Aa
              ( (~GG[0][id2][b] & GG[1][id2][b] & GG[0][id3][b] & GG[1][id3][b]) |
                (GG[0][id2][b] & ~GG[1][id2][b] & ~GG[0][id3][b] & GG[1][id3][b]) );
               word; word &= (word-1))
               TCount[bb+rightmost[word]] ++;
            for(word = informative &   // Aa->AA or aa->Aa
              ( (~GG[0][id3][b] & GG[1][id3][b] & GG[0][id1][b] & GG[1][id1][b]) |
                (GG[0][id3][b] & ~GG[1][id3][b] & ~GG[0][id1][b] & GG[1][id1][b]) );
               word; word &= (word-1))
               NTCount[bb+rightmost[word]] ++;
            for(word = informative &   // Aa->AA or aa->Aa
              ( (~GG[0][id3][b] & GG[1][id3][b] & GG[0][id2][b] & GG[1][id2][b]) |
                (GG[0][id3][b] & ~GG[1][id3][b] & ~GG[0][id2][b] & GG[1][id2][b]) );
               word; word &= (word-1))
               NTCount[bb+rightmost[word]] ++;
         }
      }  // end of trio
      for(int b = blockb; b < bMax; b++)
         for(int j = 0; j < 16; j++){
            int m = b*16+j;
            if(m >= markerCount) continue;
            int jj = (b-blockb)*16+j;
            if(snpName.Length())
               sprintf(buffer, "%s", (const char*)snpName[m]);
            else
               sprintf(buffer, "SNP%d", m+1);
            buffers[b].Add(buffer);
            if(chromosomes.Length()){
               sprintf(buffer, " %d", chromosomes[m]);
               buffers[b].Add(buffer);
            }
            if(bp.Length()){
               sprintf(buffer, " %d", bp[m]);
               buffers[b].Add(buffer);
            }
            sprintf(buffer, " %s %s",
               (const char*)alleleLabel[0][m], (const char*)alleleLabel[1][m]);
            buffers[b].Add(buffer);
            double temp = missingCount[jj] == idCount? 0:
               (AACount[jj] + AaCount[jj]*0.5) / (idCount - missingCount[jj]);
            sprintf(buffer, " %.4lf %d %d %d %d %.4lf",
               temp,
               idCount - missingCount[jj], AACount[jj], AaCount[jj],
               idCount - missingCount[jj] - AACount[jj] - AaCount[jj],
               1 - missingCount[jj]*1.0/idCount);
            buffers[b].Add(buffer);
            if(L0.Length()){
               int tempcount = nonmissingMZCount[jj] - ibs0MZCount[jj] - ibs1MZCount[jj];
               sprintf(buffer, " %d %d %d %.4lf %.4lf",
                  nonmissingMZCount[jj], ibs1MZCount[jj]+HetHetMZCount[jj], nonmissingMZCount[jj] - tempcount,
                  nonmissingMZCount[jj]? 1 - tempcount * 1.0 / nonmissingMZCount[jj]: 0,
                  ibs1MZCount[jj]+HetHetMZCount[jj]? ibs1MZCount[jj] * 1.0 / (ibs1MZCount[jj] + HetHetMZCount[jj]): 0);
               buffers[b].Add(buffer);
            }
            if(Lpo.Length()){
               sprintf(buffer, " %d %d %d %.4lf %.4lf",
                  nonmissingCount[jj], HomHomCount[jj], ibs0Count[jj],
                  nonmissingCount[jj]? ibs0Count[jj] * 1.0 / nonmissingCount[jj]: 0,
                  HomHomCount[jj]? ibs0Count[jj] * 1.0 / HomHomCount[jj]: 0);
               buffers[b].Add(buffer);
            }
            if(Ltrio.Length()){
               sprintf(buffer, " %d %d %d %.4lf %.4lf %d %d %d",
                  nonmissingtrioCount[jj], HetInOffspringCount[jj], MItrioCount[jj],
                  nonmissingtrioCount[jj]? MItrioCount[jj] * 1.0 / nonmissingtrioCount[jj]: 0,
                  HetInOffspringCount[jj]? MItrioCount[jj] * 1.0 / HetInOffspringCount[jj]: 0,
                  informativeCount[jj], TCount[jj], NTCount[jj]);
               buffers[b].Add(buffer);
            }                                             
            sprintf(buffer, "\n");
            buffers[b].Add(buffer);
         }
   }  // end of blockb

   pedfile.Add("bySNP.txt");
   FILE *fp = fopen(pedfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
   fprintf(fp, "SNP");
   if(chromosomes.Length() || xbp.Length()) fprintf(fp, " Chr");
   if(bp.Length() || xbp.Length()) fprintf(fp, " Pos");
   fprintf(fp, " Label_A Label_a Freq_A N N_AA N_Aa N_aa CallRate");
   if(L0.Length()) fprintf(fp, " N_MZ N_HetMZ N_errMZ Err_InMZ Err_InHetMZ");
   if(Lpo.Length()) fprintf(fp, " N_PO N_HomPO N_errPO Err_InPO Err_InHomPO");
   if(Ltrio.Length()) fprintf(fp, " N_trio N_HetOff N_errTrio Err_InTrio Err_InHetTrio N_TD N_T N_NT");
   fprintf(fp, "\n");
   for(int b = 0; b < buffers.Length(); b++)
      buffers[b].Write(fp);

   if(xmarkerCount) printf("Scanning chromosome X for QC-by-SNP...\n");
   double frequencies[16];
   int AA, Aa, missing;
   unsigned short int HetHet, ibs0, ibs1, notMissing, HomHom, MI_trio, HetInOffspring, T, NT;
   for(int m = 0; m < xmarkerCount; m+=16){
      int b = m/16;
      for(int i = 0; i < 16; i++){
         AACount[i] = AaCount[i] = missingCount[i] = 0;
         frequencies[i] = 0.0;
      }
      for(int i = 0; i < idCount; i++){
         AA = XG[0][i][b] & XG[1][i][b];
         Aa = (~XG[0][i][b]) & XG[1][i][b];
         missing = (~XG[0][i][b]) & (~XG[1][i][b]);
         for(int j = 0; j < 16; j++){
            if(AA & shortbase[j])
               AACount[j] ++;
            else if(Aa & shortbase[j])
               AaCount[j] ++;
            else if(missing & shortbase[j])
               missingCount[j] ++;
         }
      }
      for(int j = 0; j < 16; j++)
         if(missingCount[j] < idCount)
            frequencies[j] = (AACount[j] + AaCount[j]*0.5) / (idCount - missingCount[j]);
      // MZ
      if(L0.Length())
         for(int j = 0; j < 16; j++)
            ibs0MZCount[j] = ibs1MZCount[j] = ibs2MZCount[j]
            = nonmissingMZCount[j] = HetHetMZCount[j] = 0;
      for(int i = 0; i < L0.Length()/2; i++){
         id1 = geno[L0[i*2]];
         id2 = geno[L0[i*2+1]];
         HetHet = (~XG[0][id1][b]) & XG[1][id1][b] & (~XG[0][id2][b]) & XG[1][id2][b];
         ibs0 = XG[0][id1][b] & XG[0][id2][b] & (XG[1][id1][b]^XG[1][id2][b]);
         ibs1 = ((~XG[0][id1][b]) & XG[1][id1][b] & XG[0][id2][b])
            | ((~XG[0][id2][b]) & XG[1][id2][b] & XG[0][id1][b]);
         notMissing = (XG[0][id1][b] | XG[1][id1][b]) & (XG[0][id2][b] | XG[1][id2][b]);
         for(int j = 0; j < 16; j++)
            if(notMissing & shortbase[j]){
               nonmissingMZCount[j] ++;
               if(ibs0 & shortbase[j])
                  ibs0MZCount[j] ++;
               else if(ibs1 & shortbase[j])
                  ibs1MZCount[j] ++;
               else if(HetHet & shortbase[j])
                  HetHetMZCount[j] ++;
         }
      }
      for(int j = 0; j < 16; j++)
         ibs2MZCount[j] = nonmissingMZCount[j] - ibs0MZCount[j] - ibs1MZCount[j];
      // MI in PO
      if(Lpo.Length())
         for(int j = 0; j < 16; j++)
            ibs0Count[j] = nonmissingCount[j] = HomHomCount[j] = 0;
      for(int i = 0; i < Lpo.Length()/2; i++){
         if(ped[Lpo[i*2]].sex==1 && ped[Lpo[i*2+1]].sex==1)
            continue;   // dad does not pass genes to son
         if(ped[Lpo[i*2]].sex + ped[Lpo[i*2+1]].sex == 1)
            continue;   // sex unknown, cannot check
         id1 = geno[Lpo[i*2]];
         id2 = geno[Lpo[i*2+1]];
         ibs0 = XG[0][id1][b] & XG[0][id2][b] & (XG[1][id1][b] ^ XG[1][id2][b]);
         notMissing = (XG[0][id1][b] | XG[1][id1][b]) & (XG[0][id2][b] | XG[1][id2][b]);
         HomHom = XG[0][id1][b] & XG[0][id2][b];
         for(int j = 0; j < 16; j++){
            if(ibs0 & shortbase[j])
               ibs0Count[j] ++;
            if(notMissing & shortbase[j])
               nonmissingCount[j] ++;
            if(HomHom & shortbase[j])
               HomHomCount[j] ++;
         }
      }
      if(Ltrio.Length())
         for(int j = 0; j < 16; j++)
            MItrioCount[j] = nonmissingtrioCount[j] = HetInOffspringCount[j] =
               TCount[j] = NTCount[j] = informativeCount[j] = 0;
      for(int i = 0; i < Ltrio.Length()/3; i++){
         id1 = geno[Ltrio[i*3+1]];
         id2 = geno[Ltrio[i*3+2]];
         id3 = geno[Ltrio[i*3]];   // id3 is the child
         MI_trio = XG[0][id1][b] & XG[0][id2][b] & (~(XG[1][id1][b] ^ XG[1][id2][b]))
            & (~XG[0][id3][b]) & XG[1][id3][b];
         notMissing = (XG[0][id1][b] | XG[1][id1][b]) & (XG[0][id2][b] | XG[1][id2][b])
            & (XG[0][id3][b] | XG[1][id3][b]);
         HetInOffspring = (~XG[0][id3][b]) & XG[1][id3][b] & notMissing;
         for(int j = 0; j < 16; j++){
            if(MI_trio & shortbase[j])
               MItrioCount[j] ++;
            if(notMissing & shortbase[j])
               nonmissingtrioCount[j] ++;
            if(HetInOffspring & shortbase[j])
               HetInOffspringCount[j] ++;
         }
         // TD
         informative = (XG[0][id3][b] | XG[1][id3][b]) & // offspring not missing
            ( (~XG[0][id1][b] & XG[1][id1][b] & (XG[0][id2][b] | XG[1][id2][b])) |
              (~XG[0][id2][b] & XG[1][id2][b] & (XG[0][id1][b] | XG[1][id1][b])) );
         for(int j = 0; j < 16; j++)
            if(informative & shortbase[j])
               informativeCount[j]++;
         for(int k = 0; k < 2; k++){
            id1 = geno[Ltrio[i*3+1+k]];
            T = informative &   // Aa->AA or aa->Aa
              ( (~XG[0][id1][b] & XG[1][id1][b] & XG[0][id3][b] & XG[1][id3][b]) |
                (XG[0][id1][b] & ~XG[1][id1][b] & ~XG[0][id3][b] & XG[1][id3][b]) );
            NT = informative &   // AA->Aa or Aa->aa
              ( (~XG[0][id3][b] & XG[1][id3][b] & XG[0][id1][b] & XG[1][id1][b]) |
                (XG[0][id3][b] & ~XG[1][id3][b] & ~XG[0][id1][b] & XG[1][id1][b]) );
            if(T)
               for(int j = 0; j < 16; j++)
                  if(T & shortbase[j])
                     TCount[j]++;
            if(NT)
               for(int j = 0; j < 16; j++)
                  if(NT & shortbase[j])
                     NTCount[j]++;
         }
      }
      for(int j = 0; j < 16; j++){
         if(m+j >= xmarkerCount) continue;
         if(xsnpName.Length())
            fprintf(fp, "%s", (const char*)xsnpName[m+j]);
         else
            fprintf(fp, "SNPX%d", m+j+1);
         fprintf(fp, " X");
         if(xbp.Length())
            fprintf(fp, " %d", xbp[m+j]);
         fprintf(fp, " %s %s %.4lf %d %d %d %d %.4lf",
            (const char*)xalleleLabel[0][m+j], (const char*)xalleleLabel[1][m+j],
            frequencies[j],
            idCount - missingCount[j], AACount[j], AaCount[j],
            idCount - missingCount[j] - AACount[j] - AaCount[j],
            1 - missingCount[j]*1.0/idCount);
         if(L0.Length())
            fprintf(fp, " %d %d %d %.4lf %.4lf",
               nonmissingMZCount[j], ibs1MZCount[j]+HetHetMZCount[j], nonmissingMZCount[j] - ibs2MZCount[j],
               nonmissingMZCount[j]? 1-ibs2MZCount[j] * 1.0 / nonmissingMZCount[j]: 0,
               ibs1MZCount[j]+HetHetMZCount[j]? ibs1MZCount[j] * 1.0 / (ibs1MZCount[j] + HetHetMZCount[j]): 0);
         if(Lpo.Length())
            fprintf(fp, " %d %d %d %.4lf %.4lf",
               nonmissingCount[j], HomHomCount[j], ibs0Count[j],
               nonmissingCount[j]? ibs0Count[j] * 1.0 / nonmissingCount[j]: 0,
               HomHomCount[j]? ibs0Count[j] * 1.0 / HomHomCount[j]: 0);
         if(Ltrio.Length())
            fprintf(fp, " %d %d %d %.4lf %.4lf %d %d %d",
               nonmissingtrioCount[j], HetInOffspringCount[j], MItrioCount[j],
               nonmissingtrioCount[j]? MItrioCount[j] * 1.0 / nonmissingtrioCount[j]: 0,
               HetInOffspringCount[j]? MItrioCount[j] * 1.0 / HetInOffspringCount[j]: 0,
               informativeCount[j], TCount[j], NTCount[j]);
         fprintf(fp, "\n");
      }
   }

   if(ymarkerCount) printf("Scanning chromosome Y for QC-by-SNP...\n");
   for(int m = 0; m < ymarkerCount; m+=16){
      int b = m/16;
      for(int i = 0; i < 16; i++){
         AACount[i] = AaCount[i] = missingCount[i] = 0;
         frequencies[i] = 0.0;
      }
      for(int i = 0; i < idCount; i++){
         AA = YG[0][i][b] & YG[1][i][b];
         Aa = (~YG[0][i][b]) & YG[1][i][b];
         missing = (~YG[0][i][b]) & (~YG[1][i][b]);
         for(int j = 0; j < 16; j++){
            if(AA & shortbase[j])
               AACount[j] ++;
            else if(Aa & shortbase[j])
               AaCount[j] ++;
            else if(missing & shortbase[j])
               missingCount[j] ++;
         }
      }
      for(int j = 0; j < 16; j++)
         if(missingCount[j] < idCount)
            frequencies[j] = (AACount[j] + AaCount[j]*0.5) / (idCount - missingCount[j]);
      // MZ
      if(L0.Length())
         for(int j = 0; j < 16; j++)
            ibs0MZCount[j] = ibs1MZCount[j] = ibs2MZCount[j]
            = nonmissingMZCount[j] = HetHetMZCount[j] = 0;
      for(int i = 0; i < L0.Length()/2; i++){
         id1 = geno[L0[i*2]];
         id2 = geno[L0[i*2+1]];
         HetHet = (~YG[0][id1][b]) & YG[1][id1][b] & (~YG[0][id2][b]) & YG[1][id2][b];
         ibs0 = YG[0][id1][b] & YG[0][id2][b] & (YG[1][id1][b]^YG[1][id2][b]);
         ibs1 = ((~YG[0][id1][b]) & YG[1][id1][b] & YG[0][id2][b])
            | ((~YG[0][id2][b]) & YG[1][id2][b] & YG[0][id1][b]);
         notMissing = (YG[0][id1][b] | YG[1][id1][b]) & (YG[0][id2][b] | YG[1][id2][b]);
         for(int j = 0; j < 16; j++)
            if(notMissing & shortbase[j]){
               nonmissingMZCount[j] ++;
               if(ibs0 & shortbase[j])
                  ibs0MZCount[j] ++;
               else if(ibs1 & shortbase[j])
                  ibs1MZCount[j] ++;
               else if(HetHet & shortbase[j])
                  HetHetMZCount[j] ++;
            }
      }
      for(int j = 0; j < 16; j++)
         ibs2MZCount[j] = nonmissingMZCount[j] - ibs0MZCount[j] - ibs1MZCount[j];
      // MI in PO
      if(Lpo.Length())
         for(int j = 0; j < 16; j++)
            ibs0Count[j] = nonmissingCount[j] = HomHomCount[j] = 0;
      for(int i = 0; i < Lpo.Length()/2; i++){
         id1 = geno[Lpo[i*2]];
         id2 = geno[Lpo[i*2+1]];
         ibs0 = YG[0][id1][b] & YG[0][id2][b] & (YG[1][id1][b] ^ YG[1][id2][b]);
         notMissing = (YG[0][id1][b] | YG[1][id1][b]) & (YG[0][id2][b] | YG[1][id2][b]);
         HomHom = YG[0][id1][b] & YG[0][id2][b];
         for(int j = 0; j < 16; j++){
            if(ibs0 & shortbase[j])
               ibs0Count[j] ++;
            if(notMissing & shortbase[j])
               nonmissingCount[j] ++;
            if(HomHom & shortbase[j])
               HomHomCount[j] ++;
         }
      }
      for(int j = 0; j < 16; j++){
         if(m+j >= ymarkerCount) continue;
         if(ysnpName.Length())
            fprintf(fp, "%s", (const char*)ysnpName[m+j]);
         else
            fprintf(fp, "SNPY%d", m+j+1);
         fprintf(fp, " Y");
         if(ybp.Length())
            fprintf(fp, " %d", ybp[m+j]);
         fprintf(fp, " %s %s %.4lf %d %d %d %d %.4lf",
            (const char*)yalleleLabel[0][m+j], (const char*)yalleleLabel[1][m+j],
            frequencies[j],
            idCount - missingCount[j], AACount[j], AaCount[j],
            idCount - missingCount[j] - AACount[j] - AaCount[j],
            1 - missingCount[j]*1.0/idCount);
         if(L0.Length())
            fprintf(fp, " %d %d %d %.4lf %.4lf",
               nonmissingMZCount[j], ibs1MZCount[j]+HetHetMZCount[j], nonmissingMZCount[j] - ibs2MZCount[j],
               nonmissingMZCount[j]? 1-ibs2MZCount[j] * 1.0 / nonmissingMZCount[j]: 0,
               ibs1MZCount[j]+HetHetMZCount[j]? ibs1MZCount[j] * 1.0 / (ibs1MZCount[j] + HetHetMZCount[j]): 0);
         if(Lpo.Length())
            fprintf(fp, " %d %d %d %.4lf %.4lf",
               nonmissingCount[j], HomHomCount[j], ibs0Count[j],
               nonmissingCount[j]? ibs0Count[j] * 1.0 / nonmissingCount[j]: 0,
               HomHomCount[j]? ibs0Count[j] * 1.0 / HomHomCount[j]: 0);
         if(Ltrio.Length())
            fprintf(fp, " 0 0 0 0 0 0 0 0");
         fprintf(fp, "\n");
      }
   }
   for(int m = 0; m < mtmarkerCount; m+=16){
      int b = m/16;
      for(int i = 0; i < 16; i++){
         AACount[i] = AaCount[i] = missingCount[i] = 0;
         frequencies[i] = 0.0;
      }
      for(int i = 0; i < idCount; i++){
         AA = MG[0][i][b] & MG[1][i][b];
         Aa = (~MG[0][i][b]) & MG[1][i][b];
         missing = (~MG[0][i][b]) & (~MG[1][i][b]);
         for(int j = 0; j < 16; j++){
            if(AA & shortbase[j])
               AACount[j] ++;
            else if(Aa & shortbase[j])
               AaCount[j] ++;
            else if(missing & shortbase[j])
               missingCount[j] ++;
         }
      }

      for(int j = 0; j < 16; j++)
         if(missingCount[j] < idCount)
            frequencies[j] = (AACount[j] + AaCount[j]*0.5) / (idCount - missingCount[j]);

      // MZ
      if(L0.Length())
         for(int j = 0; j < 16; j++)
            ibs0MZCount[j] = ibs1MZCount[j] = ibs2MZCount[j]
            = nonmissingMZCount[j] = HetHetMZCount[j] = 0;
      for(int i = 0; i < L0.Length()/2; i++){
         id1 = geno[L0[i*2]];
         id2 = geno[L0[i*2+1]];
         HetHet = (~MG[0][id1][b]) & MG[1][id1][b] & (~MG[0][id2][b]) & MG[1][id2][b];
         ibs0 = MG[0][id1][b] & MG[0][id2][b] & (MG[1][id1][b]^MG[1][id2][b]);
         ibs1 = ((~MG[0][id1][b]) & MG[1][id1][b] & MG[0][id2][b])
            | ((~MG[0][id2][b]) & MG[1][id2][b] & MG[0][id1][b]);
         notMissing = (MG[0][id1][b] | MG[1][id1][b]) & (MG[0][id2][b] | MG[1][id2][b]);
         for(int j = 0; j < 16; j++)
            if(notMissing & shortbase[j]){
               nonmissingMZCount[j] ++;
               if(ibs0 & shortbase[j])
                  ibs0MZCount[j] ++;
               else if(ibs1 & shortbase[j])
                  ibs1MZCount[j] ++;
               else if(HetHet & shortbase[j])
                  HetHetMZCount[j] ++;
            }
      }
      for(int j = 0; j < 16; j++){
         if(m+j >= mtmarkerCount) continue;
         ibs2MZCount[j] = nonmissingMZCount[j] - ibs0MZCount[j] - ibs1MZCount[j];
         if(mtsnpName.Length())
            fprintf(fp, "%s", (const char*)mtsnpName[m+j]);
         else
            fprintf(fp, "SNPMT%d", m+j+1);
         fprintf(fp, " MT");
         if(mtbp.Length())
            fprintf(fp, " %d", mtbp[m+j]);
         fprintf(fp, " %s %s %.4lf %d %d %d %d %.4lf",
            (const char*)mtalleleLabel[0][m+j], (const char*)mtalleleLabel[1][m+j],
            frequencies[j],
            idCount - missingCount[j], AACount[j], AaCount[j],
            idCount - missingCount[j] - AACount[j] - AaCount[j],
            1 - missingCount[j]*1.0/idCount);
         if(L0.Length())
            fprintf(fp, " %d %d %d %.4lf %.4lf",
               nonmissingMZCount[j], ibs1MZCount[j]+HetHetMZCount[j], nonmissingMZCount[j] - ibs2MZCount[j],
               nonmissingMZCount[j]? 1-ibs2MZCount[j] * 1.0 / nonmissingMZCount[j]: 0,
               ibs1MZCount[j]+HetHetMZCount[j]? ibs1MZCount[j] * 1.0 / (ibs1MZCount[j] + HetHetMZCount[j]): 0);
         if(Lpo.Length())
            fprintf(fp, " 0 0 0 0 0");
         if(Ltrio.Length())
            fprintf(fp, " 0 0 0 0 0 0 0 0");
         fprintf(fp, "\n");
      }
   }
   fclose(fp);
   printf("QC-by-SNP ends at %s", currentTime());
   printf("QC statistics by SNPs saved in file %s\n\n", (const char*)pedfile);
}

void Engine::OutputIndividualInfo()
{
   countGenotype();
   String datfile = prefix;
   datfile.Add("af.dat");
   FILE *fp = fopen(datfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)datfile);
   fprintf(fp, "C N_SNP\n");
   fprintf(fp, "C Mean\n");
   fprintf(fp, "C Var\n");
   fprintf(fp, "C Heterozygosity\n");
   if(BetaSum.Length()){
      fprintf(fp, "C ExpectedHet\n");
//      fprintf(fp, "C ExcessHet\n");
      fprintf(fp, "C F\n");
     }
//   fprintf(fp, "C QualityScore\n");
//   if(detailFlag){
      fprintf(fp, "C CallRate\n");
      fprintf(fp, "C N_AA\n");
      fprintf(fp, "C N_Aa\n");
      fprintf(fp, "C N_aa\n");
//   }
   if(xsnpName.Length()){ // X-chromosome exists
      fprintf(fp, "C N_XSNP\n");
      fprintf(fp, "C XHeterozygosity\n");
      fprintf(fp, "C XError\n");
   }
   if(ysnpName.Length()){ // Y-chromosome exists
      fprintf(fp, "C N_YSNP\n");
      fprintf(fp, "C N_Aa\n");
   }
   if(mtsnpName.Length()){ // MT
      fprintf(fp, "C N_MTSNP\n");
      fprintf(fp, "C N_Aa\n");
   }

   fclose(fp);
   String filename = prefix;
   filename.Add("af.ped");
   fp = fopen(filename, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)filename);
   Vector pcX;
   if(BetaSum.Length())
      pcX.Dimension(BetaSum.Length());
   for(int i = 0; i < ped.count; i++){
      int k = geno[i];
      if(k == -1) {
         fprintf(fp, "%s %s %s %s %d X X X X X",
            (const char*)ped[i].famid, (const char*)ped[i].pid,
            (const char*)ped[i].fatid, (const char*)ped[i].motid,
            ped[i].sex);
         if(detailFlag)
            fprintf(fp, " X X X X");
         if(xsnpName.Length())
            fprintf(fp, " X X X");
         fprintf(fp, "\n");
         continue;
      }
      int count = pAACount[k] + pAaCount[k] + paaCount[k];
      double mean = (pAACount[k] + pAaCount[k]*0.5) / count;
      double var = pAaCount[k];
      var = mean * (1-mean) - var / count / 2;

//      double temp1 = 4.0 * pAACount[k] * paaCount[k];
//      double temp2 = pAaCount[k];
//      temp2 *= (pAaCount[k]-1.0);
//      double quality = (temp1 - temp2) / (temp1 + temp2);
      fprintf(fp, "%s %s %s %s %d %d %.4lf %.4lf %.4lf",
         (const char*)ped[i].famid, (const char*)ped[i].pid,
         (const char*)ped[i].fatid, (const char*)ped[i].motid,
         ped[i].sex, count, mean, var, pAaCount[k] * 1.0 / count);
      if(BetaSum.Length()){
         int nCov = BetaSum.Length();
         pcX[0] = 1.0;
         for(int j = 1; j < nCov; j++)
            pcX[j] = ped[i].covariates[covariatePC[j-1]];

         Vector myBetaSum = BetaSum;
         Matrix myBetaSquareSum = BetaSquareSum;
         for(int m1 = 0; m1 < shortCount; m1++){
            int m3 = ~GG[0][k][m1] & ~GG[1][k][m1] & 0xFFFF;
            if(m3==0) continue;  // missing somewhere
            for(int m2 = 0; m2 < 16; m2++)
               if(m3 & shortbase[m2]){
                  int m4 = m1*16 + m2;
                  for(int u = 0; u < nCov; u++)
                     myBetaSum[u] -= freqBeta[m4][u];
                  for(int u = 0; u < nCov; u++)
                     for(int v = 0; v < nCov; v++)
                        myBetaSquareSum[u][v] -= freqBeta[m4][u] * freqBeta[m4][v];
            }
         }
         double Evar = 0;
         for(int j = 0; j < nCov; j++)
            Evar += myBetaSum[j] * pcX[j];
         for(int u = 0; u < nCov; u++)
            for(int v = 0; v < nCov; v++)
               Evar -= myBetaSquareSum[u][v] * pcX[u] * pcX[v];
         Evar *= 2.0;
         fprintf(fp, " %.4lf %.4lf", Evar/markerCount, (Evar - pAaCount[k])/Evar);
      }
//      if(detailFlag)
         fprintf(fp, " %.4lf %d %d %d",
            count*1.0/markerCount, pAACount[k], pAaCount[k], paaCount[k]);
      if(xsnpName.Length()){
         int xcount = pxAACount[k] + pxAaCount[k] + pxaaCount[k];
         double xhet = pxAaCount[k] * 1.0 / xcount;
         fprintf(fp, " %d %.4lf %d",
            xcount, xhet, ((xhet>0.1)+1)!=ped[i].sex);
      }
      if(ysnpName.Length()){
         int ycount = pyAACount[k] + pyAaCount[k] + pyaaCount[k];
         fprintf(fp, " %d %d",
            ycount, pyAaCount[k]);
      }
      if(mtsnpName.Length()){
         int mtcount = pmtAACount[k] + pmtAaCount[k] + pmtaaCount[k];
         fprintf(fp, " %d %d",
            mtcount, pmtAaCount[k]);
      }
      fprintf(fp, "\n");
   }
   fclose(fp);
   printf("Allele frequency Statistics saved in files %s and %s\n\n",
         (const char*)datfile, (const char*)filename);
}

void Engine::countGenotype()
{
   if(pAACount==NULL) {
      pAACount = new int [idCount];
      pAaCount = new int [idCount];
      paaCount = new int [idCount];
      for(int i = 0; i < idCount; i++)
         pAACount[i] = pAaCount[i] = paaCount[i] = 0;
   }else
      return;
   if(xmarkerCount){
      if(pxAACount==0){
         pxAACount = new int [idCount];
         pxAaCount = new int [idCount];
         pxaaCount = new int [idCount];
         for(int i = 0; i < idCount; i++)
            pxAACount[i] = pxAaCount[i] = pxaaCount[i] = 0;
      }else
         return;
   }
   if(ymarkerCount){
      if(pyAACount==0){
         pyAACount = new int [idCount];
         pyAaCount = new int [idCount];
         pyaaCount = new int [idCount];
         for(int i = 0; i < idCount; i++)
            pyAACount[i] = pyAaCount[i] = pyaaCount[i] = 0;
      }else
         return;
   }
   if(mtmarkerCount){
      if(pmtAACount == NULL){
         pmtAACount = new int [idCount];
         pmtAaCount = new int [idCount];
         pmtaaCount = new int [idCount];
         for(int i = 0; i < idCount; i++)
            pmtAACount[i] = pmtAaCount[i] = pmtaaCount[i] = 0;
      }else
         return;
   }

   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   for(int k = 0; k < idCount; k++){
      if(shortFlip)
      for(int m = 0; m < shortCount; m++){
         pAACount[k] += oneoneCount[GG[0][k][m] & (GG[1][k][m]^shortFlip[m])];
         pAaCount[k] += oneoneCount[(~GG[0][k][m]) & GG[1][k][m]];
         paaCount[k] += oneoneCount[GG[0][k][m] & (~(GG[1][k][m]^shortFlip[m]))];
      }
      else
      for(int m = 0; m < shortCount; m++){
         pAACount[k] += oneoneCount[GG[0][k][m] & GG[1][k][m]];
         pAaCount[k] += oneoneCount[(~GG[0][k][m]) & GG[1][k][m]];
         paaCount[k] += oneoneCount[GG[0][k][m] & (~GG[1][k][m])];
      }
   }

   if(quality==NULL) {
      quality = new int [idCount];
      Vector score(idCount);
      for(int k = 0; k < idCount; k++){
         double temp1 = 4.0 * pAACount[k] * paaCount[k];
         double temp2 = pAaCount[k] * (pAaCount[k]-1.0);
         score[k] = (temp1 - temp2) / (temp1 + temp2);
      }
      QuickIndex idx;
      idx.Index(score);
      for(int i = 0; i < idCount; i++)
         quality[idx[i]] = i;
   }

   for(int k = 0; k < idCount; k++){
      for(int m = 0; m < xshortCount; m++){
         pxAACount[k] += oneoneCount[XG[0][k][m] & XG[1][k][m]];
         pxAaCount[k] += oneoneCount[(~XG[0][k][m]) & XG[1][k][m]];
         pxaaCount[k] += oneoneCount[XG[0][k][m] & (~XG[1][k][m])];
      }
      for(int m = 0; m < yshortCount; m++){
         pyAACount[k] += oneoneCount[YG[0][k][m] & YG[1][k][m]];
         pyAaCount[k] += oneoneCount[(~YG[0][k][m]) & YG[1][k][m]];
         pyaaCount[k] += oneoneCount[YG[0][k][m] & (~YG[1][k][m])];
      }
      for(int m = 0; m < mtshortCount; m++){
         pmtAACount[k] += oneoneCount[MG[0][k][m] & MG[1][k][m]];
         pmtAaCount[k] += oneoneCount[(~MG[0][k][m]) & MG[1][k][m]];
         pmtaaCount[k] += oneoneCount[MG[0][k][m] & (~MG[1][k][m])];
      }
   }
}


