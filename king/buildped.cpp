//////////////////////////////////////////////////////////////////////
// buildped.cpp
// (c) 2010-2018 Wei-Min Chen
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
// March 22, 2018

#include <math.h>
#include "analysis.h"
#include "Kinship.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

void Engine::internalKING(int degree)
{
   if(shortCount==0) error("No genotype data");
   printf("Autosome genotypes stored in %d", Bit64==64? longCount:shortCount);
   printf(" words for each of %d individuals.\n", idCount);

   if(bigdataFlag)
      printf("Fast algorithms for big data are used for relationship inference\n");
   else
      printf("Standard algorithms for smaller datasets are used for relationship inference\n");
/*
   int stop2 = 2304; // 256*9
   if(degree > 1) stop2 <<= 3;
   if(Bit64==64){
      if(SLG[0]==NULL)
         ConvertLGtoSLG(LG, markerCount, SLG, (stop2/4 < longCount)? (stop2/4)*64: markerCount);
   }else{
      if(SG[0]==NULL)
         ConvertGGtoSG(GG, markerCount, SG, (stop2 < shortCount)? stop2*16: markerCount);
   }
 */
    int stop1, stop2;
    if(Bit64==64){
      stop1 = 64;
      stop2 = (stop1<<3);
      if(degree == 2) stop1 <<= 3;
      if(stop1 > longCount) stop1 = longCount;
      if(stop2 > longCount) stop2 = longCount;
   }else{
      stop1 = 256;
      stop2 = (stop1<<3);
      if(degree == 2) stop1 <<= 3;
      if(stop1 > shortCount) stop1 = shortCount;
      if(stop2 > shortCount) stop2 = shortCount;
   }

   printf("Sorting autosomes...\n");
   if(Bit64==64){
      if(SLG[0]==NULL){
         if(degree==1)
            ConvertLGtoSLG(LG, markerCount, SLG, (stop2 < longCount)? (stop2<<6): markerCount);
         else if(degree==2)
            ConvertLGtoSLG(LG, markerCount, SLG, (stop1 < longCount)? (stop1<<6): markerCount);
         else
            error("Degree of relatedness not defined.");
      }
//         ConvertLGtoSLG(LG, markerCount, SLG, (stop2/4 < longCount)? (stop2/4)*64: markerCount);
   }else{
      if(SG[0]==NULL){
         if(degree==1)
            ConvertGGtoSG(GG, markerCount, SG, (stop2 < shortCount)? (stop2<<4): markerCount);
         else if(degree==2)
            ConvertGGtoSG(GG, markerCount, SG, (stop1 < shortCount)? (stop1<<4): markerCount);
         else
            error("Degree of relatedness not defined.");
      }
//      if(SG[0]==NULL)
//         ConvertGGtoSG(GG, markerCount, SG, (stop2 < shortCount)? stop2*16: markerCount);
   }

//   IntArray chrSeg;
//   double totalLength;
   bool IBDvalidFlag = false;
//   String segmessage;
   if(Bit64==64)
      IBDvalidFlag = PreSegment(/*chrSeg, totalLength, segmessage*/);
   if(!IBDvalidFlag){
      printf("%s\n", (const char*)segmessage);
      printf("  Inference will be based on kinship estimation only.\n");
   }
#ifdef _OPENMP
   printf("%d CPU cores are used to compute the pairwise kinship coefficients...\n",
      defaultMaxCoreCount);
#endif
   relativedegree = degree;
   IntArray allpairs0, allpairs1(0), allpairs;
//   unsigned long int rawrelativeCount;

   // most computationally intensive here: SCREENING RELATIVES
   if(Bit64==64)
      ScreenCloseRelativesInSubset64Bit(allpairs0);
   else
      ScreenCloseRelativesInSubset(allpairs0);
   //****************SCREENING RELATIVES end************************

   int id1, id2;
   for(int i = 0; i < allpairs0.Length()/2; i++){
      id1 = allpairs0[i*2];
      id2 = allpairs0[i*2+1];
      if(ped[phenoid[id1]].famid != ped[phenoid[id2]].famid){
         allpairs1.Push(id1);
         allpairs1.Push(id2);
      }
   }
   unsigned long int midrelativeCount = allpairs1.Length() / 2;
   if(midrelativeCount==0) return;

   double lowerbound = pow(2.0, -(relativedegree+1.5));
   Vector pM(midrelativeCount);
   for(int i = 0; i < midrelativeCount; i++){
      pM[i] = allpairs1[i*2];
      pM[i] *= idCount;
      pM[i] += allpairs1[i*2+1];
   }
   QuickIndex index;
   index.Index(pM);
   allpairs.Dimension(midrelativeCount*2);
   for(int i = 0; i < midrelativeCount; i++){
      int k = index[i];
      allpairs[i*2] = allpairs1[k*2];
      allpairs[i*2+1] = allpairs1[k*2+1];
   }

   IntArray HetHetCounts, IBS0Counts, het1Counts, het2Counts, HomHomCounts, IBSCounts;

   //****************************Compute kinship coefficient*****************************
   if(Bit64==64)
      KinshipInSubset64Bit(allpairs, HetHetCounts, IBS0Counts, het1Counts, het2Counts, HomHomCounts, IBSCounts);
   else
      KinshipInSubset(allpairs, HetHetCounts, IBS0Counts, het1Counts, het2Counts, HomHomCounts);
   //****************************Compute kinship coefficient*****************************

   if(Bit64==32){
      unsigned long int midrelativeCount = allpairs.Length() / 2;
      if(midrelativeCount==0) return;
//      double lowerbound = pow(2.0, -(relativedegree+1.5));
      double smaller;
      double kinship;
      L0.Dimension(0);
      Lfs.Dimension(0);
      Lpo.Dimension(0);
      L2.Dimension(0);
      IntArray L1(0);
      Vector IBS0L1(0);
      double ibs0;
      for(int p = 0; p < midrelativeCount; p++){
         id1 = allpairs[p*2];
         id2 = allpairs[p*2+1];
         smaller = HetHetCounts[p] + (het1Counts[p] < het2Counts[p]? het1Counts[p]: het2Counts[p]);
         kinship = 0.5 - ((het1Counts[p]+het2Counts[p])*0.25+IBS0Counts[p])/smaller;
         int deg = int(-log(kinship)/0.6931472-0.5);;
         if(deg == 0){
            L0.Push(phenoid[id1]);
            L0.Push(phenoid[id2]);
         }else if(deg == 1){
            int IBS0Count = IBS0Counts[p];
            int notMissingCount = het1Counts[p] + het2Counts[p] + HomHomCounts[p] + HetHetCounts[p];
            ibs0 = IBS0Count*1.0/notMissingCount;
            IBS0L1.Push(ibs0);
            L1.Push(phenoid[id1]);
            L1.Push(phenoid[id2]);
         }else if(deg == 2){
            L2.Push(phenoid[id1]); L2.Push(phenoid[id2]);
         }
      }
      if(errorrateCutoff == _NAN_){
      IntArray L1Count(10);
      L1Count.Zero();
      int lowT, highT;
      for(int i = 0; i < IBS0L1.Length(); i++)
         if(IBS0L1[i] < 0.01)
            L1Count[int(IBS0L1[i]*1000)] ++;
      for(lowT=0; L1Count[lowT] && lowT < 9; lowT++);
      for(highT=9; L1Count[highT] && highT >=0; highT--);
      if(lowT<=highT)
         errorrateCutoff = (lowT+highT+1)*0.0005;
      else{
         for(highT = 9; L1Count[highT] > L1Count.Min(); highT--);
         errorrateCutoff = (highT+0.5)*0.001;
      }
      printf("Cutoff value between full siblings and parent-offspring is set at %.4f\n",
         errorrateCutoff);
      }
      for(int i = 0; i < IBS0L1.Length(); i++)
         if(IBS0L1[i] > errorrateCutoff) {
            Lfs.Push(L1[i*2]);
            Lfs.Push(L1[i*2+1]);
         }else{
            Lpo.Push(L1[i*2]);
            Lpo.Push(L1[i*2+1]);
      }
      return;
   }
   double smaller, kinship;
   Vector ibdprops, maxLengths, ibd2props, maxLengths2;
   double ibdprop, ibd2prop;
   if(IBDvalidFlag){
      //******************Compute IBD2 Segments****************************
      if(Bit64==64)
         IBDSegInSubset64Bit(allpairs, ibdprops, maxLengths, ibd2props, maxLengths2);
      else
         IBD2SegInSubset(allpairs, ibd2props, maxLengths2);
      //******************Compute IBD2 Segments****************************
   }
   L0.Dimension(0);
   Lfs.Dimension(0);
   Lpo.Dimension(0);
   L2.Dimension(0);
   for(int p = 0; p < midrelativeCount; p++){
      if(!het1Counts[p] && !het2Counts[p] && !HetHetCounts[p]) continue;
      id1 = allpairs[p*2];
      id2 = allpairs[p*2+1];
      smaller = HetHetCounts[p] + (het1Counts[p] < het2Counts[p]? het1Counts[p]: het2Counts[p]);
      kinship = 0.5 - ((het1Counts[p]+het2Counts[p])*0.25+IBS0Counts[p])/smaller;
      if(IBDvalidFlag){
         ibdprop = ibdprops[p];
         ibd2prop = ibd2props[p];
      }
      if(relativedegree == 1){
         if(IBDvalidFlag){
            double pi = ibd2prop + ibdprop*0.5;
            if((kinship < 0.125) ||   // pass if phi < 0.125
               ((kinship < lowerbound) &&   //  pass if phi<0.177 AND following
               (pi < 0.3535534 || (ibd2prop<=0.08 && ibdprop+ibd2prop<0.9))) )
               continue;
         }else{
            if(kinship < lowerbound) continue;
         }
      }else{   // degree == 2
         if(IBDvalidFlag){
            if((kinship < 0.0442) ||   // pass if phi < 0.0442
               ((kinship < lowerbound) &&   //  pass if phi<0.0884 AND following
               (ibdprop+ibd2prop <= 0.3535534)) )
               continue;
         }else{
            if(kinship < lowerbound) continue;
         }
      }
      double CHet = HetHetCounts[p] * 1.0/ (HetHetCounts[p] + het1Counts[p] + het2Counts[p]);
      if(IBDvalidFlag){
         if(CHet<0.8){
            double pi = ibd2prop + ibdprop * 0.5;
            if(pi > 0.3535534){  // 1st-degree
               if(ibdprop + ibd2prop > 0.96 || (ibdprop + ibd2prop > 0.9 && ibd2prop <= 0.08)){
                  Lpo.Push(phenoid[id1]);
                  Lpo.Push(phenoid[id2]);
               }else if(ibd2prop > 0.08){
                  Lfs.Push(phenoid[id1]);
                  Lfs.Push(phenoid[id2]);
               }else{
                  L2.Push(phenoid[id1]);
                  L2.Push(phenoid[id2]);
               }
            }else if(pi > 0.1767767){  // 2nd-degree
               if(pi > 0.32 && ibd2prop > 0.15){
                  Lfs.Push(phenoid[id1]);
                  Lfs.Push(phenoid[id2]);
               }else{
                  L2.Push(phenoid[id1]);
                  L2.Push(phenoid[id2]);
               }
            }
         }else{
            L0.Push(phenoid[id1]);
            L0.Push(phenoid[id2]);
         }
      }  // end of if IBDvalidFlag
   }  // end of pairs
}

int Engine::ClusterFamily(int pedrebuildFlag, int degree)
{
   if(degree==0) degree = 1;
   if(degree>2){
      degree = 2;
      printf("Up to 2nd-degree relatedness (across families) is supported at the moment.\n");
   }
   if(pedrebuildFlag==1 || pedrebuildFlag==0 || (allflags&(1<<BysampleFLAG)) || (allflags&(1<<BysnpFLAG)) || unrelatedExtraction){ // will be expanded later
      printf("\nOptions in effect:\n");
      if(pedrebuildFlag==2){  // many applications
         if(unrelatedExtraction)
            printf("\t--unrelated\n");
         else {
            if(allflags&(1<<BysampleFLAG))
               printf("\t--bysample\n");
            if(allflags&(1<<BysnpFLAG))
               printf("\t--bySNP\n");
          }
      }
      if(pedrebuildFlag==1)
         printf("\t--build\n");
      else if(!unrelatedExtraction)
         printf("\t--cluster\n");
      if(degree > 1)
         printf("\t--degree 2\n");
      if(bigdataFlag)
         if(slower)
            printf("\t--slower %d\n", slower);
      if(CoreCount)
         printf("\t--cpus %d\n", CoreCount);
      if(SaveFormat == "MERLIN")
         printf("\t--merlin\n");
      if(SaveFormat == "PLINK")
         printf("\t--plink\n");
      if(prefix!="king")
         printf("\t--prefix %s\n", (const char*)prefix);
      printf("\n");
   }

   printf("Family clustering starts at %s", currentTime());
   if(idCount >= 100)  // fast computation
      internalKING(degree);
   else
      runKING();
   printf("\nClustering up to %d%s-degree relatives in families...\n",
      degree, degree==1?"st":"nd");

   uniqueIID = true;
   for(int i = 0; i < ped.count; i++){
      for(int j = i+1; j < ped.count; j++)
         if(ped[i].pid == ped[j].pid){
            uniqueIID = false;
            printf("  Individual IDs are not unique and family IDs will be used as well.\n");
            printf("  E.g., FAM %s IID %s and FAM %s IID %s have the same individual ID\n",
               (const char*)ped[i].famid, (const char*)ped[i].pid,
               (const char*)ped[j].famid, (const char*)ped[j].pid);
            break;
         }
      if(!uniqueIID) break;
   }
   if(uniqueIID)
      printf("Individual IDs are unique across all families.\n");

   StringArray oldID(ped.count);
   for(int i = 0; i < ped.count; i++){
      oldID[i] = ped[i].famid;
      oldID[i] += "->";
      oldID[i] += ped[i].pid;
   }

//   countGenotype();
   IntArray Lr = L0;
   Lr.Stack(Lpo);
   Lr.Stack(Lfs);
   if(degree > 1)
      Lr.Stack(L2);

   IntArray cluster[65535];
   int clusterCount = 0;
   String temp1, temp2;

   if(Lr.Length() == 0){
      printf("No families were found to be connected.\n");
      if(pedrebuildFlag==0) return 0;
   }else{ // clusters exist
      IntArray afterCount(6);
      afterCount[0] = L0.Length()/2;
      afterCount[1] = (Lpo.Length() + Lfs.Length())/2;
      afterCount[2] = L2.Length()/2;
//      afterCount[3] = L3.Length()/2;
      afterCount[4] = 0;
      afterCount[5] = Lpo.Length()/2;
      printRelationship(NULL, afterCount);

      IntArray famserial(ped.count);
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
            famserial[i] = f;

      int fid1, fid2;
      int exist[2];
      int smaller, larger;
      for(int i = 0; i < Lr.Length()/2; i++){
         exist[0] = exist[1] = -1;
         fid1 = famserial[Lr[i*2]];
         fid2 = famserial[Lr[i*2+1]];
         for(int j = 0; j < clusterCount; j++)
            if(cluster[j].Find(fid1)>-1 && cluster[j].Find(fid2)>-1){
               exist[0] = exist[1] = j;
               break;
            }else if(cluster[j].Find(fid1)>-1)
               exist[0] = j;
            else if(cluster[j].Find(fid2)>-1)
               exist[1] = j;
         if(exist[0] == -1 && exist[1] == -1){
            cluster[clusterCount].Dimension(0);
            cluster[clusterCount].Push(fid1);
            cluster[clusterCount].Push(fid2);
            clusterCount++;
         }else if(exist[0] > -1 && exist[1] > -1){
            if(exist[0]==exist[1]) continue;
            // combine exist[0] and exist[1];
            if(exist[0] < exist[1]) {
               smaller=exist[0]; larger=exist[1];
            }else{
               smaller=exist[1]; larger=exist[0];
            }
            cluster[smaller].Stack(cluster[larger]);
            if(larger < clusterCount-1)
               cluster[larger] = cluster[clusterCount-1];
            clusterCount--;
         }else if(exist[0] > -1)
            cluster[exist[0]].Push(fid2);
         else
            cluster[exist[1]].Push(fid1);
      }
      IntArray clusterID(ped.familyCount);
      clusterID.Set(-1);
      for(int i = 0; i < clusterCount; i++)
         for(int j = 0; j < cluster[i].Length(); j++)
            clusterID[cluster[i][j]] = i;
      if(clusterCount > 0 && clusterCount < 50){
         printf("The following families are found to be connected\n");
         printf("  %-10s%-50s\n", "NewFamID", "OriginalFamID");
         for(int i = 0; i < clusterCount; i++){
            printf("  KING%-6d", i+1);
            printf("%s", (const char*)ped.families[cluster[i][0]]->famid);
            for(int j = 1; j < cluster[i].Length(); j++)
               printf(",%s", (const char*)ped.families[cluster[i][j]]->famid);
            printf("\n");
         }
         printf("\n");
      }else if(clusterCount >= 50)
         printf("Families are clustered into %d new families\n", clusterCount);

      if(pedrebuildFlag==0){ // clustering only
         bool updateidFlag;
         String updateidsfile = prefix;
         updateidsfile.Add("updateids.txt");
         FILE *fp = fopen((const char*)updateidsfile, "wt");
         for(int f = 0; f < ped.familyCount; f++)
            if(clusterID[f] != -1){
               temp1 = ped.families[f]->famid;
               temp1.Add("->");
               for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
                  updateidFlag = true;
                  if(ped[i].pid.SubStr(0,4)=="KING" && ped[i].ngeno==0)
                     updateidFlag = false;
                  if(updateidFlag)
                     fprintf(fp, "%s\t%s", (const char*)ped[i].famid, (const char*)ped[i].pid);
                  ped[i].famid = "KING";
                  ped[i].famid += (clusterID[f]+1);
                  if(!uniqueIID){
                     temp2 = temp1;
                     temp2.Add(ped[i].pid);
                     ped[i].pid = temp2;
                     if(ped[i].fatid != "0"){
                        temp2 = temp1;
                        temp2.Add(ped[i].fatid);
                        ped[i].fatid = temp2;
                     }
                     if(ped[i].motid != "0"){
                        temp2 = temp1;
                        temp2.Add(ped[i].motid);
                        ped[i].motid = temp2;
                     }
                  }
                  if(updateidFlag)
                     fprintf(fp, "\t%s\t%s\n", (const char*)ped[i].famid, (const char*)ped[i].pid);
               }
            }
         fclose(fp);
         printf("Update-ID information is saved in file %s\n\n",
            (const char*)updateidsfile);
         temp1 = prefix;
         if(SaveFormat == "MERLIN")
            WriteMerlin();
         if(SaveFormat == "PLINK"){ // PLINK output by default
            temp1.Add("cluster");
            WritePlinkBinary(temp1);
         }
         printf("KING cluster analysis ends at %s", currentTime());
         return 1;
      } // end of if cluster only
   }  // end of if Lr.Length()
   String newName, tempName;

   if(FID.Length() == 0){ // read from Merlin format input
      sampleName.Dimension(0);
      FID.Dimension(0);
      PID.Dimension(0);
      FA.Dimension(0);
      MO.Dimension(0);
      SEX.Dimension(0);
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = 0; i < id[f].Length(); i++){
            tempName = ped[id[f][i]].famid;
            tempName += "->";//"_";
            tempName += ped[id[f][i]].pid;
            sampleName.Push(tempName);
            FID.Push(ped[id[f][i]].famid);
            PID.Push(ped[id[f][i]].pid);
            FA.Push(ped[id[f][i]].fatid);
            MO.Push(ped[id[f][i]].motid);
            tempName = ped[id[f][i]].sex;
            SEX.Push(tempName);
         }
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
            if(geno[i]==-1){
               FID.Push(ped[i].famid);
               PID.Push(ped[i].pid);
               FA.Push(ped[i].fatid);
               MO.Push(ped[i].motid);
               tempName = ped[i].sex;
               SEX.Push(tempName);
            }
   }

   String updateidsfile = prefix;
   updateidsfile.Add("updateids.txt");
   FILE *fp;
   if((pedrebuildFlag==1) && clusterCount)   // option --build, save update-ids file
      fp = fopen((const char*)updateidsfile, "wt");
   if(FID.Length() >= sampleName.Length())
      for(int i = 0; i < clusterCount; i++){
         newName = "KING";
         newName += (i+1);
         for(int j = 0; j < cluster[i].Length(); j++){
            temp1 = ped.families[cluster[i][j]]->famid;
            for(int k = 0; k < FID.Length(); k++){
               if(FID[k] != temp1) continue;
               if(pedrebuildFlag==1)
                  fprintf(fp, "%s\t%s", (const char*)FID[k], (const char*)PID[k]);
               FID[k] = newName;
               if((!uniqueIID) || (pedrebuildFlag==2)){
                  temp2 = PID[k];
                  PID[k] = temp1;
                  PID[k].Add("->");//("_");
                  PID[k].Add(temp2);
                  if(FA[k] != "0"){
                     temp2 = FA[k];
                     FA[k] = temp1;
                     FA[k].Add("->");//("_");
                     FA[k].Add(temp2);
                  }
                  if(MO[k] != "0"){
                     temp2 = MO[k];
                     MO[k] = temp1;
                     MO[k].Add("->");//("_");
                     MO[k].Add(temp2);
                  }
               }
               if(pedrebuildFlag==1)
                  fprintf(fp, "\t%s\t%s\n", (const char*)FID[k], (const char*)PID[k]);
               if(k < sampleName.Length()){
                  sampleName[k] = FID[k];
                  sampleName[k].Add("->");//("_");
                  sampleName[k].Add(PID[k]);
               }
            }
         }
   }
   if((pedrebuildFlag==1) && clusterCount) {
      fclose(fp);
      printf("Update-ID information is saved in file %s\n\n", (const char*)updateidsfile);
   }
   ped.familyCount = 0;
   ped.count = 0;
   MakePed();
   return 1;
}

void Engine::rebuild()
{
   printf("Pedigree reconstruction starts at %s", currentTime());
//   printf("Cutoff value to distinguish PO and FS is set at IBS0=%.4lf\n", errorrateCutoff);

   printf("Reconstructing pedigree...\n");
   IntArray nofix(0);

   cAge = ped.covariateNames.SlowFind("AGE");
   if(cAge > -1)
      printf("Covariate %s is used for pedigree reconstruction.\n",
         (const char*)ped.covariateNames[cAge]);
   else
      printf("Age information not provided.\n");

   // missingBase shouldn't be in the range of IID; default 100000
   for(missingBase=100000; ;missingBase += 100000){
      bool overlapFlag = false;
      for(int i = 0; i < ped.count; i++)
         if(int(ped[i].pid) >= missingBase && int(ped[i].pid) < missingBase+100000){
            overlapFlag = true;
            break;
         }
      if(!overlapFlag) break;
   }

//   IntArray chrSeg;
//   double totalLength;
//   String segmessage;
   if(Bit64==64)
      PreSegment(/*chrSeg, totalLength, segmessage*/);

   String updatefile = prefix;
   updatefile.Add("updateparents.txt");
   FILE *fp = fopen((const char*)updatefile, "wt");
   String message;
   for(int f = 0; f < ped.familyCount; f++){ // ready to be parallalized later
      if(!BuildOneFamily(f, chrSeg, totalLength, message)) // no pedigree reconstruction in this family
         nofix.Push(f);
      else{
         for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
            if(ped[i].pid.SubStr(0,4) == "KING" && ped[i].ngeno == 0) continue;
            fprintf(fp, "%s\t%s\t%s\t%s\n",
               (const char*)ped[i].famid, (const char*)ped[i].pid,
               (ped[i].fatid.SubStr(0, 4) == "KING" && ped[i].father && ped[i].father->ngeno==0)? "0": (const char*)ped[i].fatid,
               (ped[i].motid.SubStr(0, 4) == "KING" && ped[i].mother && ped[i].mother->ngeno==0)? "0": (const char*)ped[i].motid);
         }
         printf("%s", (const char*)message);
      }
   }
   printf("\n");
   
   fclose(fp);
   printf("\n");
   printf("Update-parent information is saved in file %s\n",
      (const char*)updatefile);

   String temp = prefix;
   if(SaveFormat == "MERLIN"){
      printf("Start writing reconstructed pedigrees in MERLIN format...\n");
      WriteMerlin();
   }
   if(SaveFormat == "PLINK"){
      printf("Start writing reconstructed pedigrees in PLINK format...\n");
      temp.Add("build");
      WritePlinkBinary(temp);
   }
   if(SaveFormat == "KING"){
      printf("Start writing reconstructed pedigrees in KING format...\n");
      temp.Add("build.king");
      WriteKingBinary(temp);
   }
   printf("Pedigree reconstruction ends at %s", currentTime());
}

void Engine::rebuild_semifamily()
{
   printf("Pedigree reconstruction starts at %s", currentTime());
   IntArray nofix(0);

   cAge = ped.covariateNames.SlowFind("AGE");
   if(cAge > -1)
      printf("Covariate %s is used for pedigree reconstruction.\n",
         (const char*)ped.covariateNames[cAge]);
   else
      printf("Age information not provided.\n");
//   missingBase = 100000;
   for(missingBase=100000; ;missingBase += 100000){
      bool overlapFlag = false;
      for(int i = 0; i < ped.count; i++)
         if(int(ped[i].pid) >= missingBase && int(ped[i].pid) < missingBase+100000){
            overlapFlag = true;
            break;
         }
      if(!overlapFlag) break;
   }

//   IntArray chrSeg;
//   double totalLength;
//   String segmessage;
   if(Bit64==64)
      PreSegment(/*chrSeg, totalLength, segmessage*/);

//   if(errorrateCutoff == _NAN_)
//      errorrateCutoff = 0.008;
//   printf("Cutoff value to distinguish PO and FS is set at IBS0=%.4lf\n", errorrateCutoff);
   String updatefile = prefix;
   updatefile.Add("updateparents.txt");
   FILE *fp = fopen((const char*)updatefile, "wt");
   String message;
   for(int f = 0; f < ped.familyCount; f++){
      if(!BuildOneFamily(f, chrSeg, totalLength, message))
         nofix.Push(f);
      else
         for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
            if(ped[i].pid.SubStr(0,4) == "KING" && ped[i].ngeno == 0) continue;
            fprintf(fp, "%s\t%s\t%s\t%s\n",
               (const char*)ped[i].famid, (const char*)ped[i].pid,
               (ped[i].fatid.SubStr(0, 4) == "KING" && ped[i].father && ped[i].father->ngeno==0)? "0": (const char*)ped[i].fatid,
               (ped[i].motid.SubStr(0, 4) == "KING" && ped[i].mother && ped[i].mother->ngeno==0)? "0": (const char*)ped[i].motid);
         }
   }
   fclose(fp);
   printf("\n");
   printf("Update-parent information is saved in file %s\n",
      (const char*)updatefile);

   String temp = prefix;
   if(SaveFormat == "MERLIN"){
      printf("Start writing reconstructed pedigrees in MERLIN format...\n");
      WriteMerlin();
   }
   if(SaveFormat == "PLINK"){
      printf("Start writing reconstructed pedigrees in PLINK format...\n");
      temp.Add("build");
      WritePlinkBinary(temp);
   }
   if(SaveFormat == "KING"){
      printf("Start writing reconstructed pedigrees in KING format...\n");
      temp.Add("build.king");
      WriteKingBinary(temp);
   }
   printf("Pedigree reconstruction ends at %s", currentTime());
}

int Engine::BuildOneFamily(int f, IntArray chrSeg, double totalLength, String & message)
{
   message.Clear();
   char buffer[256];

   Family *pf = ped.families[f];
   if(id[f].Length() < 2) return 0;

   IntArray valid(pf->count);
   valid.Set(1);
   L0.Dimension(0);
   Lfs.Dimension(0);
   Lpo.Dimension(0);
   L2.Dimension(0);
   IntArray L3(0);
   IntArray L4(0);
   Matrix relationship(pf->count, pf->count);
   relationship.Zero();
   Kinship kin;
   int id1, id2;
   double kinship;
   IntArray *poConnection, *d2Connection, *d3Connection;
   int parent, offspring;

   poConnection = new IntArray[pf->count];
   d2Connection = new IntArray[pf->count];
   d3Connection = new IntArray[pf->count];
   for(int i = 0; i < pf->count; i++){
      poConnection[i].Dimension(0);
      d2Connection[i].Dimension(0);
      d3Connection[i].Dimension(0);
   }

   IntArray pairs(0);
   kin.Setup(*ped.families[f]);
   for(int i = 0; i < id[f].Length(); i++)
      for(int j = i+1; j < id[f].Length(); j++)
         if(kin(ped[id[f][i]], ped[id[f][j]]) > 0)
            relationship[id[f][i]-pf->first][id[f][j]-pf->first] =
               relationship[id[f][j]-pf->first][id[f][i]-pf->first] =
               kin(ped[id[f][i]], ped[id[f][j]]);
         else{
            pairs.Push(geno[id[f][i]]);
            pairs.Push(geno[id[f][j]]);
         }
   int pairCount = pairs.Length()/2;
   if(pairCount==0) return 0;
   IntArray HetHetCounts, IBS0Counts, het1Counts, het2Counts, HomHomCounts, IBSCounts;
   //****************************Compute kinship coefficient*****************************
   if(Bit64==64)
      KinshipInSubset64Bit(pairs, HetHetCounts, IBS0Counts, het1Counts, het2Counts, HomHomCounts, IBSCounts);
   else
      KinshipInSubset(pairs, HetHetCounts, IBS0Counts, het1Counts, het2Counts, HomHomCounts);
   //****************************Compute kinship coefficient*****************************

   Vector ibdprops, maxLengths, ibd2props, maxLengths2;
   bool IBDvalidFlag = chrSeg.Length()>0;
   if(IBDvalidFlag){
      //****************Compute IBD2 Segments****************************
      if(Bit64==64)
         IBDSegInSubset64Bit(pairs, ibdprops, maxLengths, ibd2props, maxLengths2);
      else
         IBD2SegInSubset(pairs, ibd2props, maxLengths2);
      //****************Compute IBD2 Segments****************************
   }
   for(int p = 0; p < pairCount; p++){
      if(!het1Counts[p] && !het2Counts[p] && !HetHetCounts[p]) continue;
      id1 = pairs[p*2]; id2 = pairs[p*2+1];
      kinship = (HetHetCounts[p] - IBS0Counts[p]*2.0) / (HetHetCounts[p]*2+het1Counts[p]+het2Counts[p]);
      relationship[phenoid[id1]-pf->first][phenoid[id2]-pf->first] =
         relationship[phenoid[id2]-pf->first][phenoid[id1]-pf->first] = kinship;
      double CHet = HetHetCounts[p] * 1.0 / (HetHetCounts[p]+het1Counts[p]+het2Counts[p]);
      if(IBDvalidFlag){
         double ibdprop = ibdprops[p];
         double ibd2prop = ibd2props[p];
         if(CHet<0.8){
            double pi = ibd2prop + ibdprop * 0.5;
            if(pi > 0.3535534){  // 1st-degree
               if(ibdprop + ibd2prop > 0.96 || (ibdprop + ibd2prop > 0.9 && ibd2prop <= 0.08)){
                  Lpo.Push(phenoid[id1]-pf->first); Lpo.Push(phenoid[id2]-pf->first);
               }else if(ibd2prop > 0.15 || pi > 0.4){
                  Lfs.Push(phenoid[id1]-pf->first); Lfs.Push(phenoid[id2]-pf->first);
               }
            }else if(pi > 0.1767767){  // 2nd-degree
               if(ibd2prop < 0.08){
                  L2.Push(phenoid[id1]-pf->first); L2.Push(phenoid[id2]-pf->first);
                  d2Connection[phenoid[id1]-pf->first].Push(phenoid[id2]-pf->first);
                  d2Connection[phenoid[id2]-pf->first].Push(phenoid[id1]-pf->first);
               }
            }else if(pi > 0.08838835){ // 3rd-degree
               L3.Push(phenoid[id1]-pf->first); L3.Push(phenoid[id2]-pf->first);
               d3Connection[phenoid[id1]-pf->first].Push(phenoid[id2]-pf->first);
               d3Connection[phenoid[id2]-pf->first].Push(phenoid[id1]-pf->first);
            }
         }else{ // Duplicate
            L0.Push(phenoid[id1]-pf->first); L0.Push(phenoid[id2]-pf->first);
         }
      }
   }  // end of pairs

   if(!L0.Length() && !Lpo.Length() && !Lfs.Length()) return 0;
   sprintf(buffer, "Family %s:\n", (const char*)pf->famid);
   message += buffer;

   IntArray Lr = L0;
   Lr.Stack(Lpo);
   Lr.Stack(Lfs);
   Lr.Stack(L2);
   Lr.Stack(L3);
   int fam1, fam2;
   IntArray childCount(pf->count);
   IntArray sibship[200];
   int sibshipCount;
   int s1, s2;
   int fa, mo;
   IntArray pat, mat;
   int smaller, larger;
   String tempS;

   if(L0.Length()){
   childCount.Zero();
   for(int i = pf->first; i <= pf->last; i++){
      if(ped[i].father) childCount[ped[i].father->serial-pf->first] ++;
      if(ped[i].mother) childCount[ped[i].mother->serial-pf->first] ++;
   }
   for(int j = 0; j < L0.Length()/2; j++){
      if(!valid[L0[j*2]] || !valid[L0[j*2+1]]) continue;
      if(childCount[L0[j*2]] < childCount[L0[j*2+1]]){
         id2 = L0[j*2]; // to be removed
         id1 = L0[j*2+1];
      }else if(childCount[L0[j*2+1]] < childCount[L0[j*2]]){
         id1 = L0[j*2];
         id2 = L0[j*2+1]; // to be removed
      }else{
         if(ped[L0[j*2]].sibCount <= 1 && ped[L0[j*2+1]].sibCount <= 1){
            if(childCount[L0[j*2]] <= 1 && childCount[L0[j*2+1]] > 1){
               id2 = L0[j*2]; id1 = L0[j*2+1];
            }else if(childCount[L0[j*2]] > 1 && childCount[L0[j*2+1]] <= 1){
               id2 = L0[j*2+1]; id1 = L0[j*2];
            }else{
               id2 = L0[j*2]; id1 = L0[j*2+1];
            }
         }else if(ped[L0[j*2]].sibCount < ped[L0[j*2+1]].sibCount){
            id2 = L0[j*2]; id1 = L0[j*2+1];
         }else if(ped[L0[j*2+1]].sibCount < ped[L0[j*2]].sibCount){
            id2 = L0[j*2+1]; id1 = L0[j*2];
         }else{
            id2 = L0[j*2]; id1 = L0[j*2+1];
         }
      }
      valid[id2] = 0;
      sprintf(buffer, "  Duplicate %s (of %s) is removed.\n",
         (const char*)ped[id2+pf->first].pid, (const char*)ped[id1+pf->first].pid);
      message += buffer;
      for(int i = pf->first; i <= pf->last; i++){
         if(valid[i-pf->first] && ped[id2+pf->first].sex == 1 && ped[i].father && ped[i].father->serial==id2+pf->first){
            ped[i].father = &ped[id1+pf->first];
            ped[i].fatid = ped[id1+pf->first].pid;
         }else if(valid[i-pf->first] && ped[id2+pf->first].sex == 2 && ped[i].mother && ped[i].mother->serial==id2+pf->first){
            ped[i].mother = &ped[id1+pf->first];
            ped[i].motid = ped[id1+pf->first].pid;
         }
      }
      if(!ped[id1+pf->first].father && !ped[id1+pf->first].mother){
         ped[id1+pf->first].father = ped[id2+pf->first].father;
         ped[id1+pf->first].fatid = ped[id2+pf->first].fatid;
         ped[id1+pf->first].mother = ped[id2+pf->first].mother;
         ped[id1+pf->first].motid = ped[id2+pf->first].motid;
      }
   }
   }

   // FS
   pat.Dimension(0); mat.Dimension(0);
   sibshipCount = 0;
   for(int i = pf->first; i <= pf->last; i++)
      if(ped[i].sibCount > 1 && ped[i].sibs[0]->serial == i){
         sibship[sibshipCount].Dimension(0);
         for(int s = 0; s < ped[i].sibCount; s++)
            if(valid[ped[i].sibs[s]->serial-pf->first])
               sibship[sibshipCount].Push(ped[i].sibs[s]->serial-pf->first);
         ped[i].sibCount = sibship[sibshipCount].Length();
         if(ped[i].sibCount < 2) continue;
         pat.Push(ped[i].father->serial-pf->first);
         mat.Push(ped[i].mother->serial-pf->first);
         sibshipCount++;
      }
   int newparent = 0;

   for(int j = 0; j < Lfs.Length() / 2; j++){
      if(valid[Lfs[j*2]] == 0 || valid[Lfs[j*2+1]] == 0) continue;
      s1 = s2 = -1;
      for(int k = 0; k < sibshipCount; k++){
         if(sibship[k].Find(Lfs[j*2]) > -1) s1 = k;
         if(sibship[k].Find(Lfs[j*2+1]) > -1) s2 = k;
      }
      fa = mo = -1;
      if(ped[Lfs[j*2]+pf->first].father && valid[ped[Lfs[j*2]+pf->first].father->serial-pf->first]
         && ped[Lfs[j*2+1]+pf->first].father && valid[ped[Lfs[j*2+1]+pf->first].father->serial-pf->first]){
         if(ped[Lfs[j*2]+pf->first].father->ngeno > ped[Lfs[j*2+1]+pf->first].father->ngeno)
            fa = ped[Lfs[j*2]+pf->first].father->serial;
         else
            fa = ped[Lfs[j*2+1]+pf->first].father->serial;
      }else if(ped[Lfs[j*2]+pf->first].father && valid[ped[Lfs[j*2]+pf->first].father->serial-pf->first])
         fa = ped[Lfs[j*2]+pf->first].father->serial;
      else if(ped[Lfs[j*2+1]+pf->first].father && valid[ped[Lfs[j*2+1]+pf->first].father->serial-pf->first])
         fa = ped[Lfs[j*2+1]+pf->first].father->serial;
      if(ped[Lfs[j*2]+pf->first].mother && valid[ped[Lfs[j*2]+pf->first].mother->serial-pf->first]
         && ped[Lfs[j*2+1]+pf->first].mother && valid[ped[Lfs[j*2+1]+pf->first].mother->serial-pf->first]){
         if(ped[Lfs[j*2]+pf->first].mother->ngeno > ped[Lfs[j*2+1]+pf->first].mother->ngeno)
            mo = ped[Lfs[j*2]+pf->first].mother->serial;
         else
            mo = ped[Lfs[j*2+1]+pf->first].mother->serial;
      }else if(ped[Lfs[j*2]+pf->first].mother && valid[ped[Lfs[j*2]+pf->first].mother->serial-pf->first])
         mo = ped[Lfs[j*2]+pf->first].mother->serial;
      else if(ped[Lfs[j*2+1]+pf->first].mother && valid[ped[Lfs[j*2+1]+pf->first].mother->serial-pf->first])
         mo = ped[Lfs[j*2+1]+pf->first].mother->serial;
      if(s1 > -1 && s2 > -1){
         if(s1==s2) continue;
         // combine two sibships
         if(s1 < s2) {
            smaller=s1; larger=s2;
         }else{
            smaller=s2; larger=s1;
         }
         sprintf(buffer, "  Sibship (%s",
            (const char*)ped[sibship[s1][0]+pf->first].pid);
         message += buffer;
         for(int s = 1; s < sibship[s1].Length(); s++){
            sprintf(buffer, " %s", (const char*)ped[sibship[s1][s]+pf->first].pid);
            message += buffer;
         }
         sprintf(buffer, ") and sibship (%s",
            (const char*)ped[sibship[s2][0]+pf->first].pid);
         message += buffer;
         for(int s = 1; s < sibship[s2].Length(); s++){
            sprintf(buffer, " %s", (const char*)ped[sibship[s2][s]+pf->first].pid);
            message += buffer;
         }
         sprintf(buffer, ") are combined\n");
         message += buffer;
         sibship[smaller].Stack(sibship[larger]);
         if(fa > -1)
            pat[smaller] = fa-pf->first;
         if(mo > -1)
            mat[smaller] = mo-pf->first;
         if(larger < sibshipCount-1){
            sibship[larger] = sibship[sibshipCount-1];
            pat[larger] = pat[sibshipCount-1];
            mat[larger] = mat[sibshipCount-1];
         }
         sibshipCount--;
         pat.Delete(sibshipCount);
         mat.Delete(sibshipCount);
      }else if(s1 > -1){ // fs2 join in sibship s1
         sibship[s1].Push(Lfs[j*2+1]);
         if(fa > -1) pat[s1] = fa-pf->first;
         if(mo > -1) mat[s1] = mo-pf->first;
         sprintf(buffer, "  %s joins in sibship (%s",
            (const char*)ped[Lfs[j*2+1]+pf->first].pid,
            (const char*)ped[sibship[s1][0]+pf->first].pid);
         message += buffer;
         for(int s = 1; s < sibship[s1].Length()-1; s++){
            sprintf(buffer, " %s", (const char*)ped[sibship[s1][s]+pf->first].pid);
            message += buffer;
         }
         sprintf(buffer, ")\n");
         message += buffer;
      }else if(s2 > -1){ // fs1 join in sibship s2
         sibship[s2].Push(Lfs[j*2]);
         if(fa > -1) pat[s2] = fa-pf->first;
         if(mo > -1) mat[s2] = mo-pf->first;
         sprintf(buffer, "  %s joins in sibship (%s",
            (const char*)ped[Lfs[j*2]+pf->first].pid,
            (const char*)ped[sibship[s2][0]+pf->first].pid);
         message += buffer;
         for(int s = 1; s < sibship[s2].Length()-1; s++){
            sprintf(buffer, " %s", (const char*)ped[sibship[s2][s]+pf->first].pid);
            message += buffer;
         }
         sprintf(buffer, ")\n");
         message += buffer;
      }else{
         // create sibship
         sibship[sibshipCount].Dimension(0);
         sibship[sibshipCount].Push(Lfs[j*2]);
         sibship[sibshipCount].Push(Lfs[j*2+1]);
         sibshipCount++;

         sprintf(buffer, "  Sibship (%s",
            (const char*)ped[sibship[sibshipCount-1][0]+pf->first].pid);
         message += buffer;
         for(int s = 1; s < sibship[sibshipCount-1].Length(); s++){
            sprintf(buffer, " %s", (const char*)ped[sibship[sibshipCount-1][s]+pf->first].pid);
            message += buffer;
         }
         sprintf(buffer, ")'s parents are (");
         message += buffer;
         if(fa > -1){
            pat.Push(fa-pf->first);
            sprintf(buffer, "%s", (const char*)ped[fa].pid);
            message += buffer;
         }else{
            pat.Push(missingBase + newparent);
            inclusionList[0].Push(ped.families[f]->famid);
            tempS = missingBase + newparent;
            sprintf(buffer, "%s", (const char*)tempS);
            message += buffer;
            newparent++;
            inclusionList[1].Push(tempS);
            tempS = 1;
            inclusionList[2].Push(tempS);
         }
         if(mo > -1){
            mat.Push(mo-pf->first);
            sprintf(buffer, " %s)\n", (const char*)ped[mo].pid);
            message += buffer;
         }else{
            mat.Push(missingBase + newparent);
            inclusionList[0].Push(ped.families[f]->famid);
            tempS = missingBase + newparent;
            sprintf(buffer, " %s)\n", (const char*)tempS);
            message += buffer;
            newparent++;
            inclusionList[1].Push(tempS);
            tempS = 2;
            inclusionList[2].Push(tempS);
         }
      }
   }
   for(int j = 0; j < sibshipCount; j++){
      // KING -> newparent
      if(pat[j] < missingBase && mat[j] < missingBase
         && pat[j] > -1 && mat[j] > -1){
         if(ped[pat[j]+pf->first].pid.SubStr(0, 4)=="KING" && ped[pat[j]+pf->first].ngeno==0){
            ped[pat[j]+pf->first].pid = missingBase + newparent;
            newparent ++;
         }
         if(ped[mat[j]+pf->first].pid.SubStr(0, 4)=="KING" && ped[mat[j]+pf->first].ngeno==0){
            ped[mat[j]+pf->first].pid = missingBase + newparent;
            newparent ++;
         }
      }
      for(int k = 0; k < sibship[j].Length(); k++)
         if(pat[j] < missingBase && mat[j] < missingBase
            && pat[j] > -1 && mat[j] > -1){
            ped[sibship[j][k]+pf->first].father = &ped[pat[j]+pf->first];
            ped[sibship[j][k]+pf->first].mother = &ped[mat[j]+pf->first];
            ped[sibship[j][k]+pf->first].fatid = ped[pat[j]+pf->first].pid;
            ped[sibship[j][k]+pf->first].motid = ped[mat[j]+pf->first].pid;
         }else{
            ped[sibship[j][k]+pf->first].fatid = pat[j];
            ped[sibship[j][k]+pf->first].motid = mat[j];
         }

   }
   // parent-offspring
   for(int j = 0; j < Lpo.Length() / 2; j++){
      if(valid[Lpo[j*2]]==0 || valid[Lpo[j*2+1]] == 0) continue;
      if(ped[Lpo[j*2]+pf->first].fatid == ped[Lpo[j*2+1]+pf->first].pid ||
         ped[Lpo[j*2]+pf->first].motid == ped[Lpo[j*2+1]+pf->first].pid ||
         ped[Lpo[j*2+1]+pf->first].fatid == ped[Lpo[j*2]+pf->first].pid ||
         ped[Lpo[j*2+1]+pf->first].motid == ped[Lpo[j*2]+pf->first].pid) continue;
      sprintf(buffer, "  Reconstruct parent-offspring pair (%s, %s)...\n",
         (const char*)ped[Lpo[j*2]+pf->first].pid,
         (const char*)ped[Lpo[j*2+1]+pf->first].pid);
      message += buffer;
      s1 = s2 = -1;
      for(int k = 0; k < sibshipCount; k++){
         if(sibship[k].Find(Lpo[j*2]) > -1) s1 = k;
         if(sibship[k].Find(Lpo[j*2+1]) > -1) s2 = k;
      }
      parent = -1;
      offspring = -1;
      int sibshiplist = -1;
      if(s1 == -1 && s2 == -1){ // singletons
         if(ped[Lpo[j*2]+pf->first].sex == 1 && ped[Lpo[j*2+1]+pf->first].father
            && valid[ped[Lpo[j*2+1]+pf->first].father->serial-pf->first]
            && ped[Lpo[j*2+1]+pf->first].father->ngeno >= MINSNPCOUNT){
            // j*2 is a male and j*2+1 has a father then j*2+1 is parent of j*2
               parent = Lpo[j*2+1];
               offspring = Lpo[j*2];
         }else if(ped[Lpo[j*2]+pf->first].sex == 1 && ped[Lpo[j*2+1]+pf->first].mother
            && valid[ped[Lpo[j*2+1]+pf->first].mother->serial-pf->first]
            && ped[Lpo[j*2+1]+pf->first].mother->ngeno >= MINSNPCOUNT ){
            // check kinship between Lpo[j*2] and mother
               for(int k = 0; k < Lr.Length()/2; k++)
                  if( (Lpo[j*2] == Lr[k*2]
                  && ped[Lpo[j*2+1]+pf->first].mother->serial == Lr[k*2+1]+pf->first) ||
                  (Lpo[j*2] == Lr[k*2+1]
                  && ped[Lpo[j*2+1]+pf->first].mother->serial == Lr[k*2]+pf->first) ){
                     parent = Lpo[j*2+1];
                     offspring = Lpo[j*2];
                     break;
                  }
               if(parent == -1){ // Lpo[j*2] and mother are unrelated
                  parent = Lpo[j*2];
                  offspring = Lpo[j*2+1];
               }
         }else if(ped[Lpo[j*2]+pf->first].sex == 2 && ped[Lpo[j*2+1]+pf->first].mother
            && valid[ped[Lpo[j*2+1]+pf->first].mother->serial-pf->first]
            && ped[Lpo[j*2+1]+pf->first].mother->ngeno >= MINSNPCOUNT){
            // j*2 is a female and j*2+1 has a mother then j*2+1 is parent of j*2
               parent = Lpo[j*2+1];
               offspring = Lpo[j*2];
         }else if(ped[Lpo[j*2]+pf->first].sex == 2 && ped[Lpo[j*2+1]+pf->first].father
            && valid[ped[Lpo[j*2+1]+pf->first].father->serial-pf->first]
            && ped[Lpo[j*2+1]+pf->first].father->ngeno >= MINSNPCOUNT ){
            // check kinship between Lpo[j*2] and father
               for(int k = 0; k < Lr.Length()/2; k++)
                  if( (Lpo[j*2] == Lr[k*2] &&
                     ped[Lpo[j*2+1]+pf->first].father->serial == Lr[k*2+1]+pf->first) ||
                     (Lpo[j*2] == Lr[k*2+1]
                     && ped[Lpo[j*2+1]+pf->first].father->serial == Lr[k*2]+pf->first) ){
                        parent = Lpo[j*2+1];
                        offspring = Lpo[j*2];
                        break;
                     }
               if(parent == -1){ // Lpo[j*2] and father are unrelated
                  parent = Lpo[j*2];
                  offspring = Lpo[j*2+1];
               }
         }else if(ped[Lpo[j*2+1]+pf->first].sex == 1
            && ped[Lpo[j*2]+pf->first].father
            && valid[ped[Lpo[j*2]+pf->first].father->serial-pf->first]
            && ped[Lpo[j*2]+pf->first].father->ngeno >= MINSNPCOUNT){
            // j*2+1 is a male and j*2 has a father then j*2 is parent of j*2+1
               parent = Lpo[j*2];
               offspring = Lpo[j*2+1];
         }else if(ped[Lpo[j*2+1]+pf->first].sex == 1
            && ped[Lpo[j*2]+pf->first].father
            && ped[Lpo[j*2]+pf->first].mother && valid[ped[Lpo[j*2]+pf->first].mother->serial-pf->first]
            && ped[Lpo[j*2]+pf->first].mother->ngeno >= MINSNPCOUNT ){
            // check kinship between Lpo[j*2+1] and mother
               for(int k = 0; k < Lr.Length()/2; k++)
                  if( (Lpo[j*2+1] == Lr[k*2]
                     && ped[Lpo[j*2]+pf->first].mother->serial == Lr[k*2+1]+pf->first) ||
                     (Lpo[j*2+1] == Lr[k*2+1]
                     && ped[Lpo[j*2]+pf->first].mother->serial == Lr[k*2]+pf->first) ){
                        parent = Lpo[j*2];
                        offspring = Lpo[j*2+1];
                        break;
                     }
               if(parent == -1){ // Lpo[j*2+1] and mother are unrelated
                  parent = Lpo[j*2+1];
                  offspring = Lpo[j*2];
               }
         }else if(ped[Lpo[j*2+1]+pf->first].sex == 2 && ped[Lpo[j*2]+pf->first].mother
            && valid[ped[Lpo[j*2]+pf->first].mother->serial-pf->first]
            && ped[Lpo[j*2]+pf->first].mother->ngeno >= MINSNPCOUNT){
            // j*2+1 is a female and j*2 has a mother then j*2 is parent of j*2+1
               parent = Lpo[j*2];
               offspring = Lpo[j*2+1];
         }else if(ped[Lpo[j*2+1]+pf->first].sex == 2 && ped[Lpo[j*2]+pf->first].mother
            && ped[Lpo[j*2]+pf->first].father && valid[ped[Lpo[j*2]+pf->first].father->serial-pf->first]
            && ped[Lpo[j*2]+pf->first].father->ngeno >= MINSNPCOUNT ){
            // check kinship between Lpo[j*2+1] and father
               for(int k = 0; k < Lr.Length()/2; k++)
                  if( (Lpo[j*2+1] == Lr[k*2]
                  && ped[Lpo[j*2]+pf->first].father->serial == Lr[k*2+1]+pf->first) ||
                  (Lpo[j*2+1] == Lr[k*2+1]
                  && ped[Lpo[j*2]+pf->first].father->serial == Lr[k*2]+pf->first) ){
                     parent = Lpo[j*2];
                     offspring = Lpo[j*2+1];
                     break;
                  }
               if(parent == -1){ // Lpo[j*2+1] and father are unrelated
                  parent = Lpo[j*2+1];
                  offspring = Lpo[j*2];
               }
         }else{   // none of the PO pair has known parents
            for(int k = 0; k < d2Connection[Lpo[j*2]].Length(); k++)
               if(relationship[Lpo[j*2+1]][d2Connection[Lpo[j*2]][k]] < 0.0375){
                  // Lpo[j*2+1] and Lpo[j*2]'s relative are unrelated
                  parent = Lpo[j*2+1];  // then Lpo[j*2+1] is the parent
                  offspring = Lpo[j*2];
               }
            for(int k = 0; k < d2Connection[Lpo[j*2+1]].Length(); k++)
               if(relationship[Lpo[j*2]][d2Connection[Lpo[j*2+1]][k]] < 0.0375){
                  // Lpo[j*2] and Lpo[j*2+1]'s relative are unrelated
                  parent = Lpo[j*2];  // then Lpo[j*2] is the parent
                  offspring = Lpo[j*2+1];
               }
            if(parent==-1){// Age information should be used here
               // Whoever older is the parent
               // otherwise if age is unknown
               if(cAge == -1 || ped[Lpo[j*2]+pf->first].covariates[cAge] == _NAN_ ||
                  ped[Lpo[j*2+1]+pf->first].covariates[cAge] == _NAN_){
                  poConnection[Lpo[j*2]].Push(Lpo[j*2+1]);
                  poConnection[Lpo[j*2+1]].Push(Lpo[j*2]);
               }else{// Age helps here
                  if(ped[Lpo[j*2]+pf->first].covariates[cAge] > ped[Lpo[j*2+1]+pf->first].covariates[cAge]+10){
                     parent = Lpo[j*2];
                     offspring = Lpo[j*2+1];
                     sprintf(buffer, "  Age information is used\n");
                     message += buffer;
                  }else if(ped[Lpo[j*2+1]+pf->first].covariates[cAge] > ped[Lpo[j*2]+pf->first].covariates[cAge]+10){
                     parent = Lpo[j*2+1];
                     offspring = Lpo[j*2];
                     sprintf(buffer, "  Age information is used\n");
                     message += buffer;
                  }else{
                     poConnection[Lpo[j*2]].Push(Lpo[j*2+1]);
                     poConnection[Lpo[j*2+1]].Push(Lpo[j*2]);
                  }
               }
            }
         }
      }else if(s1 > -1){   // check if Lpo[j*2+1] is P or O
         int r = sibship[s1].Find(Lpo[j*2]);
         int t = (r==0? 1: 0);
         for(int k = 0; k < Lpo.Length()/2; k++)
            if( (Lpo[k*2+1]==Lpo[j*2+1] && Lpo[k*2]==sibship[s1][t]) ||
               (Lpo[k*2]==Lpo[j*2+1] && Lpo[k*2+1]==sibship[s1][t]) )
               parent = Lpo[j*2+1];
         if(parent == -1){// Lpo[j*2+1] is offspring of Lpo[j*2]
            parent = Lpo[j*2];
            offspring = Lpo[j*2+1];
         }else // Lpo[j*2+1] is the parent
            sibshiplist = s1;
         sprintf(buffer, "  %s's sibship is used to determine the parent/offspring\n",
            (const char*)ped[Lpo[j*2]+pf->first].pid);
         message += buffer;
      }else if(s2 > -1){   // check if Lpo[j*2] is P or O
         int r = sibship[s2].Find(Lpo[j*2+1]);
         int t = (r==0? 1: 0);
         for(int k = 0; k < Lpo.Length()/2; k++)
            if( (Lpo[k*2+1]==Lpo[j*2] && Lpo[k*2]==sibship[s2][t]) ||
            (Lpo[k*2]==Lpo[j*2] && Lpo[k*2+1]==sibship[s2][t]) )
               parent = Lpo[j*2];
         if(parent == -1){// Lpo[j*2] is offspring of Lpo[j*2+1]
            parent = Lpo[j*2+1];
            offspring = Lpo[j*2];
         }else
            sibshiplist = s2;
         sprintf(buffer, "  %s's sibship is used to determine the parent/offspring\n",
            (const char*)ped[Lpo[j*2+1]+pf->first].pid);
         message += buffer;
      }
      if(parent == -1) {
         if(cAge == -1 || ped[Lpo[j*2]+pf->first].covariates[cAge] == _NAN_ ||
            ped[Lpo[j*2+1]+pf->first].covariates[cAge] == _NAN_)
            continue;
         // Age helps here
         if(ped[Lpo[j*2]+pf->first].covariates[cAge] > ped[Lpo[j*2+1]+pf->first].covariates[cAge]+10){
            parent = Lpo[j*2];
            offspring = Lpo[j*2+1];
            sprintf(buffer, "  Information of covariate %s is used\n",
               (const char*)ped.covariateNames[cAge]);
            message += buffer;
         }else if(ped[Lpo[j*2+1]+pf->first].covariates[cAge] > ped[Lpo[j*2]+pf->first].covariates[cAge]+10){
            parent = Lpo[j*2+1];
            offspring = Lpo[j*2];
            sprintf(buffer, "  Information of covariate %s is used\n",
                (const char*)ped.covariateNames[cAge]);
            message += buffer;
         }else
            continue;
      }
      if(ped[parent+pf->first].sex == 1){
         if(offspring > -1){
            tempS = ped[offspring+pf->first].fatid;
            if(tempS == "0")
               sprintf(buffer, "  %s is now father of %s\n",
               (const char*)ped[parent+pf->first].pid, (const char*)ped[offspring+pf->first].pid);
            else
               sprintf(buffer, "  %s is now father of %s, replacing %s (%d genotypes)\n",
                  (const char*)ped[parent+pf->first].pid,
                  (const char*)ped[offspring+pf->first].pid,
                  (const char*)tempS,
                  ped[offspring+pf->first].father ?
                  ped[ped[offspring+pf->first].father->serial].ngeno: 0);
            message += buffer;
            ped[offspring+pf->first].fatid = ped[parent+pf->first].pid;
            ped[offspring+pf->first].father = &ped[parent+pf->first];
            if(ped[offspring+pf->first].motid=="0"){
               inclusionList[0].Push(ped.families[f]->famid);
               tempS = missingBase + newparent;
               sprintf(buffer, "    %s is created as %s's mother.\n",
                  (const char*)tempS, (const char*)ped[offspring+pf->first].pid);
               message += buffer;
               ped[offspring+pf->first].motid = tempS;
               newparent++;
               inclusionList[1].Push(tempS);
               tempS = 2;
               inclusionList[2].Push(tempS);
            }
         }else{
            tempS = ped[sibship[sibshiplist][0]+pf->first].fatid;
            if(tempS == "0")
               sprintf(buffer, "  %s is now father of %s's sibship\n",
                  (const char*)ped[parent+pf->first].pid,
                  (const char*)ped[sibship[sibshiplist][0]+pf->first].pid);
            else
               sprintf(buffer, "  %s is now father of %s's sibship, replacing %s (%d genotypes)\n",
                  (const char*)ped[parent+pf->first].pid,
                  (const char*)ped[sibship[sibshiplist][0]+pf->first].pid,
                  (const char*)tempS,
                  ped[sibship[sibshiplist][0]+pf->first].father?
                  ped[ped[sibship[sibshiplist][0]+pf->first].father->serial].ngeno: 0);
            message += buffer;
            for(int k = 0; k < sibship[sibshiplist].Length(); k++){
               ped[sibship[sibshiplist][k]+pf->first].fatid = ped[parent+pf->first].pid;
               ped[sibship[sibshiplist][k]+pf->first].father = &ped[parent+pf->first];
            }
         }
      }else{
         if(offspring > -1){
            tempS = ped[offspring+pf->first].motid;
            if(tempS == "0")
               sprintf(buffer, "  %s is now mother of %s\n",
               (const char*)ped[parent+pf->first].pid, (const char*)ped[offspring+pf->first].pid);
            else
               sprintf(buffer, "  %s is now mother of %s, replacing %s (%d genotypes)\n",
                  (const char*)ped[parent+pf->first].pid,
                  (const char*)ped[offspring+pf->first].pid,
                  (const char*)tempS,
                  ped[offspring+pf->first].mother?
                  ped[ped[offspring+pf->first].mother->serial].ngeno:0);
            message += buffer;
            ped[offspring+pf->first].motid = ped[parent+pf->first].pid;
            ped[offspring+pf->first].mother = &ped[parent+pf->first];
            if(ped[offspring+pf->first].fatid=="0"){
               inclusionList[0].Push(ped.families[f]->famid);
               tempS = missingBase + newparent;
               sprintf(buffer, "    %s is created as %s's father.\n",
                  (const char*)tempS, (const char*)ped[offspring+pf->first].pid);
               message += buffer;
               ped[offspring+pf->first].fatid = tempS;
               newparent++;
               inclusionList[1].Push(tempS);
               tempS = 1;
               inclusionList[2].Push(tempS);
            }
         }else{
            tempS = ped[sibship[sibshiplist][0]+pf->first].motid;
            if(tempS == "0")
               sprintf(buffer, "  %s is now mother of %s's sibship\n",
                  (const char*)ped[parent+pf->first].pid,
                  (const char*)ped[sibship[sibshiplist][0]+pf->first].pid);
            else
               sprintf(buffer, "  %s is now mother of %s's sibship, replacing %s (%d genotypes)\n",
                  (const char*)ped[parent+pf->first].pid,
                  (const char*)ped[sibship[sibshiplist][0]+pf->first].pid,
                  (const char*)tempS,
                  ped[sibship[sibshiplist][0]+pf->first].mother?
                  ped[ped[sibship[sibshiplist][0]+pf->first].mother->serial].ngeno:0);
            message += buffer;
            for(int k = 0; k < sibship[sibshiplist].Length(); k++){
               ped[sibship[sibshiplist][k]+pf->first].motid = ped[parent+pf->first].pid;
               ped[sibship[sibshiplist][k]+pf->first].mother = &ped[parent+pf->first];
            }
         }
      }
   }
   for(int i = 0; i < pf->count; i++){
      int couple1 = -1; int couple2 = -1;
      for(int j = 0; j < poConnection[i].Length(); j++)
         for(int k = j+1; k < poConnection[i].Length(); k++)
            if(relationship[poConnection[i][j]][poConnection[i][k]] < 0.0625){
               // two PO are unrelated
               if(ped[poConnection[i][j]+pf->first].sex==1){
                  couple1 = poConnection[i][j];
                  couple2 = poConnection[i][k];
               }else{
                  couple1 = poConnection[i][k];
                  couple2 = poConnection[i][j];
               }
            }
      if(couple1 > -1){ // two parents identified
            ped[i+pf->first].fatid = ped[couple1+pf->first].pid;
            ped[i+pf->first].motid = ped[couple2+pf->first].pid;
            sprintf(buffer, "  %s and %s are %s's parents\n",
               (const char*)ped[i+pf->first].fatid,
               (const char*)ped[i+pf->first].motid,
               (const char*)ped[i+pf->first].pid);
            message += buffer;
            for(int k = 0; k < poConnection[i].Length(); k++){
               if(poConnection[i][k] == couple1 || poConnection[i][k] == couple2)
                  continue;
               if(ped[i+pf->first].sex == 1)
                  ped[poConnection[i][k] + pf->first].fatid = i;
               else
                  ped[poConnection[i][k] + pf->first].motid = i;
            }
      }else if(d2Connection[i].Length()){
         for(int k = 0; k < poConnection[i].Length(); k++)
            if(relationship[d2Connection[i][0]][poConnection[i][k]] > 0.0375 &&
               relationship[d2Connection[i][0]][poConnection[i][k]] < 0.0884){
                  parent = i;
                  offspring = poConnection[i][k];
                  if(ped[i+pf->first].sex == 1){
                     ped[offspring + pf->first].fatid = ped[parent+pf->first].pid;
                     sprintf(buffer, "  %s is now father of %s\n",
                        (const char*)ped[parent+pf->first].pid,
                        (const char*)ped[offspring+pf->first].pid);
                     message += buffer;
                     if(ped[poConnection[i][k] + pf->first].motid=="0"){
                        inclusionList[0].Push(ped.families[f]->famid);
                        tempS = missingBase + newparent;
                        sprintf(buffer, "    %s is created as %s's mother.\n",
                        (const char*)tempS, (const char*)ped[offspring+pf->first].pid);
                        message += buffer;
                        ped[offspring+pf->first].motid = tempS;
                        newparent++;
                        inclusionList[1].Push(tempS);
                        tempS = 2;
                        inclusionList[2].Push(tempS);
                     }
                  }else{
                     ped[offspring + pf->first].motid = ped[parent+pf->first].pid;
                     sprintf(buffer, "  %s is now mother of %s\n",
                        (const char*)ped[parent+pf->first].pid,
                        (const char*)ped[offspring+pf->first].pid);
                     message += buffer;
                     if(ped[poConnection[i][k] + pf->first].fatid=="0"){
                        inclusionList[0].Push(ped.families[f]->famid);
                        tempS = missingBase + newparent;
                        sprintf(buffer, "    %s is created as %s's father.\n",
                        (const char*)tempS, (const char*)ped[offspring+pf->first].pid);
                        message += buffer;
                        ped[offspring+pf->first].fatid = tempS;
                        newparent++;
                        inclusionList[1].Push(tempS);
                        tempS = 1;
                        inclusionList[2].Push(tempS);
                     }
                  }  // end of if sex
         }
      }  // end of if there are two parents
   }  // end of loop over each person
   for(int i = pf->first; i <= pf->last; i++)
      if(!valid[i-pf->first]){
         tempS = ped[i].famid;
         tempS += "->";// '->'
         tempS += ped[i].pid;
         exclusionList.Push(tempS);
      }

   delete []poConnection;
   delete []d2Connection;
   delete []d3Connection;
   missingBase += newparent;
   return 1;
}


