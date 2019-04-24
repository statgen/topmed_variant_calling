//////////////////////////////////////////////////////////////////////
// autoQC.cpp
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
// March 6, 2019

#include "analysis.h"
#include "Kinship.h"
#include "QuickIndex.h"
#include "rplot.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

void Engine::autoQC(double samplecallrate, double snpcallrate)
{
   printf("Autosome genotypes stored in %d words for each of %d individuals.\n",
      shortCount, idCount);
   printf("\nOptions in effect:\n");
   printf("\t--autoQC\n");
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\nAuto-QC starts at %s\n", currentTime());

   if(idMask == NULL)
      idMask = new char[idCount];
   char *tempidMask = new char[idCount];
   for(int i = 0; i < idCount; i++)
      idMask[i] = 0;
   if(gMask == NULL && shortCount) {
      gMask = new unsigned short int[shortCount];
      for(int m = 0; m < shortCount; m++) gMask[m] = 0;
      for(int j = 0; j < 16-markerCount%16; j++)
         gMask[shortCount-1] |= shortbase[15-j];
   }
   if(xMask == NULL && xshortCount) {
      xMask = new unsigned short int[xshortCount];
      for(int m = 0; m < xshortCount; m++) xMask[m] = 0;
      for(int j = 0; j < 16-xmarkerCount%16; j++)
         xMask[xshortCount-1] |= shortbase[15-j];
   }
   if(yMask == NULL && yshortCount) {
      yMask = new unsigned short int[yshortCount];
      for(int m = 0; m < yshortCount; m++) yMask[m] = 0;
      for(int j = 0; j < 16-ymarkerCount%16; j++)
         yMask[yshortCount-1] |= shortbase[15-j];
   }
   unsigned char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   IntArray snptoberemovedTMP, mlog(0);
   IntArray sampletoberemoved(0), sampletoberemovedTMP, nlog(0);
   int ncount = idCount;
   int mcount = markerCount;
   int xmcount = xmarkerCount;
   int ymcount = ymarkerCount;
   IntArray updatesex(0);
   mlog.Push(markerCount + xmarkerCount + ymarkerCount + mtmarkerCount);
   nlog.Push(idCount);

   String snpremovalfile = prefix;
   snpremovalfile.Add("_autoQC_snptoberemoved.txt");
   FILE *fp1 = fopen(snpremovalfile, "wt");
   if(fp1 == NULL) error("Cannot open %s to write.", (const char*)snpremovalfile);
   fprintf(fp1, "SNP\tREASON\n");
   String sampleremovalfile = prefix;
   sampleremovalfile.Add("_autoQC_sampletoberemoved.txt");
   FILE *fp2 = fopen(sampleremovalfile, "wt");
   if(fp2 == NULL) error("Cannot open %s to write.", (const char*)sampleremovalfile);
   fprintf(fp2, "FID\tIID\tREASON\n");

   double prefilter=0.8;
   if(snpcallrate < 0.9)
      prefilter=int(snpcallrate*10)*0.1-0.1; // eg., 85% -> 70%
   printf("Auto-QC step 1: Apply SNP call rate filter %.1lf%% on %d SNPs (in %d samples)\n",
      prefilter*100, mcount, ncount);
   CallRate_SNP(prefilter, snptoberemovedTMP);
   mcount -= snptoberemovedTMP.Length();
   int count = snptoberemovedTMP.Length();
   for(int m = 0; m < snptoberemovedTMP.Length(); m++)
      fprintf(fp1, "%s\tCallRateLessThan%d\n",
//         (const char*)snpName[bigdataIdx[markerCount-1-snptoberemovedTMP[m]]],
         (const char*)snpName[snptoberemovedTMP[m]],
         int(prefilter*100+0.5));
   printf("  %d autosome SNPs have call rate < %.1lf%%\n",
      snptoberemovedTMP.Length(), prefilter*100);

   CallRate_xSNP(prefilter, snptoberemovedTMP);
   xmcount -= snptoberemovedTMP.Length();
   count += snptoberemovedTMP.Length();
   for(int m = 0; m < snptoberemovedTMP.Length(); m++)
      fprintf(fp1, "%s\tCallRateLessThan%d\n",
         (const char*)xsnpName[snptoberemovedTMP[m]],
         int(prefilter*100+0.5));
   printf("  %d X-chr SNPs have call rate < %.1lf%%\n",
      snptoberemovedTMP.Length(), prefilter*100);
   mlog.Push(count);

   Monomorphic_SNP(snptoberemovedTMP);
   for(int m = 0; m < snptoberemovedTMP.Length(); m++)
      if(snptoberemovedTMP[m] < markerCount){
         fprintf(fp1, "%s\tMonomorphic\n",
//         (const char*)snpName[bigdataIdx[markerCount-1-snptoberemovedTMP[m]]]);
         (const char*)snpName[snptoberemovedTMP[m]]);
         mcount --;
      }else if(snptoberemovedTMP[m] < shortCount*16+xmarkerCount){
         fprintf(fp1, "%s\tMonomorphic\n",
         (const char*)xsnpName[snptoberemovedTMP[m]-shortCount*16]);
         xmcount --;
      }else{
         fprintf(fp1, "%s\tMonomorphic\n",
         (const char*)ysnpName[snptoberemovedTMP[m]-(shortCount+xshortCount)*16]);
         ymcount --;
      }
   printf("  %d SNPs are monomorphic\n", snptoberemovedTMP.Length());
   mlog.Push(snptoberemovedTMP.Length());

   printf("\nAuto-QC step 2: Apply sample call rate filter %.1lf%% on %d samples (with %d SNPs)\n",
      samplecallrate*100, ncount, mcount);
   CallRate_Sample(samplecallrate, sampletoberemovedTMP);
   for(int i = 0; i < sampletoberemovedTMP.Length(); i++)
      idMask[sampletoberemovedTMP[i]] = 1;   // sampletoberemoved[i] is now removed
   printf("  %d samples have call rate < %.1lf%%\n",
      sampletoberemovedTMP.Length(), samplecallrate*100);
   for(int i = 0; i < sampletoberemovedTMP.Length(); i++)
      fprintf(fp2, "%s\t%s\tMissingMoreThan%d\n",
         (const char*)ped[phenoid[sampletoberemovedTMP[i]]].famid,
         (const char*)ped[phenoid[sampletoberemovedTMP[i]]].pid,
         int((1-samplecallrate)*100+0.5));
   ncount -= sampletoberemovedTMP.Length();
   nlog.Push(sampletoberemovedTMP.Length());

   printf("\nAuto-QC step 3: Apply SNP call rate filter %.1lf%% on %d SNPs (in %d samples)\n",
      snpcallrate*100, mcount, ncount);
   CallRate_SNP(snpcallrate, snptoberemovedTMP);
   printf("  %d SNPs have call rate < %.1lf%%\n",
      snptoberemovedTMP.Length(), snpcallrate*100);
   mcount -= snptoberemovedTMP.Length();
   count = snptoberemovedTMP.Length();
   for(int m = 0; m < snptoberemovedTMP.Length(); m++)
      fprintf(fp1, "%s\tCallRateLessThan%d\n",
//         (const char*)snpName[bigdataIdx[markerCount-1-snptoberemovedTMP[m]]],
         (const char*)snpName[snptoberemovedTMP[m]],
         int(snpcallrate*100+0.5));
   CallRate_xSNP(snpcallrate, snptoberemovedTMP);
   printf("  %d chr-X SNPs have call rate < %.1lf%%\n",
      snptoberemovedTMP.Length(), snpcallrate*100);
   xmcount -= snptoberemovedTMP.Length();
   count += snptoberemovedTMP.Length();
   for(int m = 0; m < snptoberemovedTMP.Length(); m++)
      fprintf(fp1, "%s\tCallRateLessThan%d\n",
         (const char*)xsnpName[snptoberemovedTMP[m]],
         int(snpcallrate*100+0.5));
   mlog.Push(count);

   if(xmarkerCount && ymarkerCount){
      printf("\nAuto-QC step 4: Apply call rate filters on %d Y-chr SNPs\n", ymarkerCount);
      printf("\n  Step 4a: Apply Y-chr call rate filter %.1lf%% in males\n", prefilter*100);
      for(int i = 0; i < idCount; i++)
         tempidMask[i] = idMask[i];
      for(int i = 0; i < idCount; i++)
         if(ped[phenoid[i]].sex!=1)
            idMask[i] = 1;
      CallRate_ySNP(prefilter, snptoberemovedTMP);
      for(int m = 0; m < snptoberemovedTMP.Length(); m++)
         fprintf(fp1, "%s\tCallRateLessThan%d\n",
            (const char*)ysnpName[snptoberemovedTMP[m]], int(snpcallrate*100+0.5));
      printf("  %d chr-Y SNPs have call rate < %.1lf%% in males\n",
         snptoberemovedTMP.Length(), prefilter*100);
      ymcount -= snptoberemovedTMP.Length();
      mlog.Push(snptoberemovedTMP.Length());

      printf("\n  Step 4b: Apply X-chr heterozygosity filter 5%% in males\n");
      xHeterozygosity_SNP(snptoberemovedTMP, 0.05);
      for(int m = 0; m < snptoberemovedTMP.Length(); m++)
         fprintf(fp1, "%s\txHeterozygosityInMale\n",
            (const char*)xsnpName[snptoberemovedTMP[m]]);
      printf("  %d X-chr SNPs have heterozygosity > 5%% in males\n",
         snptoberemovedTMP.Length());
      xmcount -= snptoberemovedTMP.Length();
      mlog.Push(snptoberemovedTMP.Length());
      for(int i = 0; i < idCount; i++)
         idMask[i] = tempidMask[i];

      printf("\n  Step 4c: Apply Y-chr call rate filter 10%% in females\n");
      for(int i = 0; i < idCount; i++)
         if(ped[phenoid[i]].sex!=2)
            idMask[i] = 1;
      CallRate_ySNP(0.1, snptoberemovedTMP, false);
      for(int m = 0; m < snptoberemovedTMP.Length(); m++)
         fprintf(fp1, "%s\tYSNPInFemales\n", (const char*)ysnpName[snptoberemovedTMP[m]]);
      printf("  %d chr-Y SNPs have call rate > 10%% in females\n",
         snptoberemovedTMP.Length());
      ymcount -= snptoberemovedTMP.Length();
      mlog.Push(snptoberemovedTMP.Length());
      for(int i = 0; i < idCount; i++)
         idMask[i] = tempidMask[i];

      printf("\nAuto-QC step 5: Gender QC on %d samples\n",
         idCount - sampletoberemovedTMP.Length());

      printf("\n  Step 5a: Determine thresholds in %d Y-chr SNPs for gender checking\n", ymcount);
      int ysnps = ymcount;
      int xsnps = xmcount;
      double cutoff_ysnps = ysnps * 0.5;
      double cutoff_left = ysnps * 0.3333;
      double cutoff_right = ysnps * 0.6667;
      printf("  %.1lf is used as the cutoff value for the chr-Y SNP count between males vs females\n",
         cutoff_ysnps);

      printf("\n  Step 5b: Gender checking\n");
      IntArray maleerror(0);
      IntArray femaleerror(0);
      IntArray gendererror(0);
      int xmissingCount, ymissingCount, xhetCount;
      IntArray plotx_ysnpcount(0);
      Vector ploty_xhet(0);
      IntArray plotz_sex(0);
      Vector tempArray(0);
      for(int id = 0; id < idCount; id++){
         if(idMask[id]) continue;
         xmissingCount = ymissingCount = xhetCount = 0;
         for(int m = 0; m < yshortCount; m++)
            ymissingCount += oneoneCount[(~YG[0][id][m])&(~YG[1][id][m])&(~yMask[m])&65535];
         for(int m = 0; m < xshortCount; m++){
            xhetCount += oneoneCount[(~XG[0][id][m])&XG[1][id][m]&(~xMask[m])&65535];
            xmissingCount += oneoneCount[(~XG[0][id][m])&(~XG[1][id][m])&(~xMask[m])&65535];
         }
         plotx_ysnpcount.Push(ysnps - ymissingCount);
         double temp = xhetCount*1.0/(xsnps-xmissingCount);
         ploty_xhet.Push(temp);
         if(ped[phenoid[id]].sex == 2 && temp < 0.1 && temp > 0.05) tempArray.Push(temp);
      }
      double xHeterozygosity = 0.1;
      for(int i = 0; i < tempArray.Length(); i++)
         if(tempArray[i] < xHeterozygosity) xHeterozygosity = tempArray[i];
      xHeterozygosity = 0.01 * int(xHeterozygosity*100);
      printf("  X-chr heterozygosity is set as %.2lf\n", xHeterozygosity);

      int index = 0;
      for(int id = 0; id < idCount; id++){
         if(idMask[id]) continue;
         plotz_sex.Push(ped[phenoid[id]].sex);
         if(ped[phenoid[id]].sex == 2){   // labeled as a female
            if(plotx_ysnpcount[index] > cutoff_ysnps)
               femaleerror.Push(id);
            else if(plotx_ysnpcount[index] > cutoff_left || ploty_xhet[index] < xHeterozygosity)
               gendererror.Push(id);
         }else if(ped[phenoid[id]].sex == 1){ // labeled as a male
            if(plotx_ysnpcount[index] < cutoff_ysnps)
               maleerror.Push(id);
            else if(plotx_ysnpcount[index] < cutoff_right || ploty_xhet[index] > xHeterozygosity)
               gendererror.Push(id);
         }else if(ped[phenoid[id]].sex == 0){
            updatesex.Push(id);
            updatesex.Push((plotx_ysnpcount[index] < cutoff_ysnps)+1);
         }
         index ++;
      }
      printf("  %d (reported) males have gender errors (less than %.1lf Y-chr SNPs)\n",
         maleerror.Length(), cutoff_ysnps);
      printf("  %d (reported) females have gender errors (more than %.1lf Y-chr SNPs)\n",
         femaleerror.Length(), cutoff_ysnps);
      printf("  %d samples have additional gender errors (according to X-Chr heterozygosity filter %.3lf)\n",
         gendererror.Length(), xHeterozygosity);
      printf("  Genders of %d samples are not reported but will be inferred now\n", updatesex.Length()/2);
      if(rplotFlag)
         plotGenderError(prefix, plotx_ysnpcount, ploty_xhet, plotz_sex, xHeterozygosity, maleerror.Length()+femaleerror.Length());

      for(int i = 0; i < maleerror.Length(); i++)
         fprintf(fp2, "%s\t%s\tMislabeledAsMale\n",
            (const char*)ped[phenoid[maleerror[i]]].famid,
            (const char*)ped[phenoid[maleerror[i]]].pid);
      for(int i = 0; i < femaleerror.Length(); i++)
         fprintf(fp2, "%s\t%s\tMislabeledAsFemale\n",
            (const char*)ped[phenoid[femaleerror[i]]].famid,
            (const char*)ped[phenoid[femaleerror[i]]].pid);
      if(xHeterozygosity > 0 && xHeterozygosity < 1)
         for(int i = 0; i < gendererror.Length(); i++)
            fprintf(fp2, "%s\t%s\tGenderQC\n",
            (const char*)ped[phenoid[gendererror[i]]].famid,
            (const char*)ped[phenoid[gendererror[i]]].pid);

      IntArray allgendererror(0);
      for(int i = 0; i < maleerror.Length(); i++) allgendererror.Push(maleerror[i]);
      for(int i = 0; i < femaleerror.Length(); i++) allgendererror.Push(femaleerror[i]);
      if(xHeterozygosity > 0 && xHeterozygosity < 1)
         for(int i = 0; i < gendererror.Length(); i++) allgendererror.Push(gendererror[i]);

      for(int i = 0; i < allgendererror.Length(); i++)
         idMask[allgendererror[i]] = 1;
      ncount -= allgendererror.Length();
      nlog.Push(maleerror.Length());
      nlog.Push(femaleerror.Length());
      nlog.Push(gendererror.Length());

      printf("\nAuto-QC step 6: Apply call rate filters on %d Y-chr SNPs\n", ymcount);
      printf("\n  Step 6a: Apply Y-chr call rate filter %.1lf%% in males\n", snpcallrate*100);
      for(int i = 0; i < idCount; i++)
         tempidMask[i] = idMask[i];
      for(int i = 0; i < idCount; i++)
         if(ped[phenoid[i]].sex!=1)
            idMask[i] = 1;
      CallRate_ySNP(snpcallrate, snptoberemovedTMP);
      ymcount -= snptoberemovedTMP.Length();
      mlog.Push(snptoberemovedTMP.Length());
      for(int m = 0; m < snptoberemovedTMP.Length(); m++)
         fprintf(fp1, "%s\tCallRateLessThan%d\n",
            (const char*)ysnpName[snptoberemovedTMP[m]],
            int(snpcallrate*100+0.5));
      printf("  %d chr-Y SNPs have call rate < %.1lf%% in males\n",
         snptoberemovedTMP.Length(), snpcallrate*100);

      printf("\n  Step 6b: Apply X-chr heterozygosity filter 1%% in males\n");
      xHeterozygosity_SNP(snptoberemovedTMP, 0.01);
      for(int m = 0; m < snptoberemovedTMP.Length(); m++)
         fprintf(fp1, "%s\txHeterozygosityInMale\n",
            (const char*)xsnpName[snptoberemovedTMP[m]]);
      printf("  %d X-chr SNPs have heterozygosity > 1%% in males\n",
         snptoberemovedTMP.Length());
      xmcount -= snptoberemovedTMP.Length();
      mlog.Push(snptoberemovedTMP.Length());
      for(int i = 0; i < idCount; i++)
         idMask[i] = tempidMask[i];

      printf("\n  Step 6c: Apply Y-chr call rate filter 2%% in females\n");
      for(int i = 0; i < idCount; i++)
         if(ped[phenoid[i]].sex!=2)
            idMask[i] = 1;
      CallRate_ySNP(0.02, snptoberemovedTMP, false);
      ymcount -= snptoberemovedTMP.Length();
      mlog.Push(snptoberemovedTMP.Length());
      for(int m = 0; m < snptoberemovedTMP.Length(); m++)
         fprintf(fp1, "%s\tYSNPInFemales\n",
            (const char*)ysnpName[snptoberemovedTMP[m]]);
      printf("  %d chr-Y SNPs have call rate > 2%% in females\n",
         snptoberemovedTMP.Length());
      for(int i = 0; i < idCount; i++)
         idMask[i] = tempidMask[i];
   }else{
      if(xmarkerCount==0)
         printf("X-Chr SNPs are not available. Gender QC is skipped.\n");
      if(ymarkerCount==0)
         printf("Y-Chr SNPs are not available. Gender QC is skipped.\n");
   }

   printf("\nAuto-QC step 7: Final check\n");
   ncount = idCount;
   for(int i = 0; i < idCount; i++)
      if(idMask[i]) ncount--;
   printf("  %d samples", ncount);
   mcount = shortCount*16;
   for(int m = 0; m < shortCount; m++)
      if(gMask[m]) mcount -= oneoneCount[gMask[m]];
   printf(", %d autosome SNPs", mcount);
   xmcount = xshortCount*16;
   for(int m = 0; m < xshortCount; m++)
      if(xMask[m]) xmcount -= oneoneCount[xMask[m]];
   if(xmcount) printf(", %d X-chr SNPs", xmcount);
   ymcount = yshortCount*16;
   for(int m = 0; m < yshortCount; m++)
      if(yMask[m]) ymcount -= oneoneCount[yMask[m]];
   if(ymcount) printf(", %d Y-chr SNPs", ymcount);
   printf("\n");
   fclose(fp1);
   fclose(fp2);

   printf("\nAuto-QC step 8: QC Summary Report\n\n");
   char buffer[1024];
   StringArray buffers(0);
   sprintf(buffer, "%-5s%-55s%-10s%-10s", "Step", "Description", "Subjects", "SNPs");
   buffers.Add(buffer);
   sprintf(buffer, "%-5d%-55s%-10d%-10d", 1, "Raw data counts", nlog[0], mlog[0]);
   buffers.Add(buffer);
   sprintf(buffer, "%-5.1f%-55s%10s(%-d)", 1.1, "SNPs with very low call rate < 80% (removed)", "", mlog[1]);
   buffers.Add(buffer);
   sprintf(buffer, "%-5.1f%-55s%10s(%-d)", 1.2, "Monomorphic SNPs (removed)", "", mlog[2]);
   buffers.Add(buffer);
   sprintf(buffer, "%-5.1f%-55s(%d)", 1.3, "Sample call rate < 95\% (removed)", nlog[1]);
   buffers.Add(buffer);
   sprintf(buffer, "%-5.1f%-55s%10s(%-d)", 1.4, "SNPs with call rate < 95% (removed)", "", mlog[3]);
   buffers.Add(buffer);
   if(mlog.Length()>4){
      sprintf(buffer, "%-5d%-55s%-10d%-10d", 2, "data counts for gender error checking",
         nlog[0]-nlog[1], mlog[0]-mlog[1]-mlog[2]-mlog[3]);
      buffers.Add(buffer);
      sprintf(buffer, "%-5.1f%-55s%10s(%-d)", 2.1, "Y-chr SNPs with call rate < 80% in men (removed)", "", mlog[4]);
      buffers.Add(buffer);
      sprintf(buffer, "%-5.1f%-55s%10s(%-d)", 2.2, "X-chr SNPs with heterozygosity > 5% in men (removed)", "", mlog[5]);
      buffers.Add(buffer);
      sprintf(buffer, "%-5.1f%-55s%10s(%-d)", 2.3, "Y-chr SNPs with genotypes in >10% women (removed)", "", mlog[6]);
      buffers.Add(buffer);
      sprintf(buffer, "%-5.1f%-55s(%d)", 2.4,"Mislabeled as male (removed)", nlog[2]);
      buffers.Add(buffer);
      sprintf(buffer, "%-5.1f%-55s(%d)", 2.5,"Mislabeled as female (removed)", nlog[3]);
      buffers.Add(buffer);
      sprintf(buffer, "%-5.1f%-55s(%d)", 2.6,"Suspicious gender error (removed)", nlog[4]);
      buffers.Add(buffer);
      sprintf(buffer, "%-5.1f%-55s%10s(%-d)", 2.7, "Y-chr SNPs with call rate < 95% in men (removed)", "", mlog[7]);
      buffers.Add(buffer);
      sprintf(buffer, "%-5.1f%-55s%10s(%-d)", 2.8, "X-chr SNPs with heterozygosity > 1% in men (removed)", "", mlog[8]);
      buffers.Add(buffer);
      sprintf(buffer, "%-5.1f%-55s%10s(%-d)", 2.9, "Y-chr SNPs with genotypes in >2% women (removed)", "", mlog[9]);
      buffers.Add(buffer);
   }
   sprintf(buffer, "%-5d%-55s", 3, "Generate Final Study Files");
   buffers.Add(buffer);
   sprintf(buffer, "%5s%-55s%-10d%-10d", "", "Final QC'ed data",
      ncount, mcount+xmcount+ymcount+mtmarkerCount);
   buffers.Add(buffer);

   buffers.Write(stdout);
   String QCsummaryfile = prefix;
   QCsummaryfile.Add("_autoQC_Summary.txt");
   FILE *fp = fopen(QCsummaryfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)QCsummaryfile);
   buffers.Write(fp);
   fclose(fp);

   if(gMask) delete []gMask;
   if(xMask) delete []xMask;
   if(yMask) delete []yMask;
   if(idMask) delete []idMask;

   printf("\nQC summary report saved in %s\n", (const char*)QCsummaryfile);
   printf("SNP-removal QC file saved in %s\n", (const char*)snpremovalfile);
   printf("Sample-removal QC file saved in %s\n", (const char*)sampleremovalfile);
   if(updatesex.Length()){
      String updatesexfile=prefix;
      updatesexfile.Add("_autoQC_updatesex.txt");
      FILE *fp = fopen(updatesexfile, "wt");
      if(fp == NULL) error("Cannot open %s to write.", (const char*)updatesexfile);
      for(int i = 0; i < updatesex.Length()/2; i++){
      fprintf(fp, "%s\t%s\t%d\n",
         (const char*)ped[phenoid[updatesex[i*2]]].famid,
         (const char*)ped[phenoid[updatesex[i*2]]].pid,
         updatesex[i*2+1]);
      }
      fclose(fp);
      printf("Update-sex QC file saved in %s\n", (const char*)updatesexfile);
   }

   printf("\nAuto-QC ends at %s", currentTime());
}

void Engine::xHeterozygosity_SNP(IntArray & removeList, double xHeterozygosity)
{
   removeList.Dimension(0);
   unsigned char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   int het, missing, fixedMask;
   int hetCount[16], missingCount[16];
   int id;
   IntArray idIndex(0);
   for(int i = 0; i < idCount; i++)
      if(idMask[i]==0) idIndex.Push(i);
   int indexCount = idIndex.Length();
   IntArray *wordremoval = new IntArray[xshortCount];
   for(int m = 0; m < xshortCount; m++)
      wordremoval[m].Dimension(0);
   for(int m = 0; m < xshortCount; m++){
      if(xMask[m]==65535) continue;  // all SNPs masked
      fixedMask = (~xMask[m]) & 65535;
      for(int j = 0; j < 16; j++)
         hetCount[j] = missingCount[j] = 0;
      for(int i = 0; i < indexCount; i++){
         id = idIndex[i];
         het = (~XG[0][id][m]) & XG[1][id][m] & fixedMask;
         missing = (~XG[0][id][m]) & (~XG[1][id][m]) & fixedMask;
         for(int j = 0; j < 16; j++)
            if(het & shortbase[j])
               hetCount[j] ++;
            else if(missing & shortbase[j])
               missingCount[j] ++;
      }
      for(int j = 0; j < 16; j++)
         if(hetCount[j] > (indexCount-missingCount[j])*xHeterozygosity){
            if((~xMask[m]) & shortbase[j])
               wordremoval[m].Push(j);
            xMask[m] |= shortbase[j];
         }
   }
   for(int m = 0; m < xshortCount; m++)
      if(xMask[m])
         for(int j = 0; j < wordremoval[m].Length(); j++)
            removeList.Push(m*16+wordremoval[m][j]);
   delete []wordremoval;
}

void Engine::Monomorphic_SNP(IntArray & removeList)
{
   removeList.Dimension(0);
   unsigned char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   int het, AA, aa, fixedMask;
   int hetCount[16], AACount[16], aaCount[16];
   int id;
   IntArray idIndex(0);
   for(int i = 0; i < idCount; i++)
      if(idMask[i]==0) idIndex.Push(i);
   int indexCount = idIndex.Length();
   IntArray *wordremoval = new IntArray[shortCount+xshortCount+yshortCount];
   for(int m = 0; m < shortCount+xshortCount+yshortCount; m++)
      wordremoval[m].Dimension(0);
   for(int i = 0; i < indexCount; i++)
      id = idIndex[i];
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultEfficientCoreCount) \
      private(id, het, AA, aa, hetCount, AACount, aaCount, fixedMask)
#endif
   for(int m = 0; m < shortCount; m++){
      if(gMask[m]==65535) continue;  // all SNPs masked
      fixedMask = (~gMask[m]) & 65535;
      for(int j = 0; j < 16; j++)
         hetCount[j] = AACount[j] = aaCount[j] = 0;
      for(int i = 0; i < indexCount; i++){
         id = idIndex[i];
         het = (~GG[0][id][m]) & GG[1][id][m] & fixedMask;
         AA = GG[0][id][m] & GG[1][id][m] & fixedMask;
         aa = GG[0][id][m] & (~GG[1][id][m]) & fixedMask;
         for(int j = 0; j < 16; j++)
            if(AA & shortbase[j])
               AACount[j] ++;
            else if(aa & shortbase[j])
               aaCount[j] ++;
            else if(het & shortbase[j])
               hetCount[j] ++;
      }
      for(int j = 0; j < 16; j++)
         if(hetCount[j] + AACount[j] == 0 || hetCount[j] + aaCount[j] == 0){
            if((~gMask[m]) & shortbase[j])
               wordremoval[m].Push(j);
            gMask[m] |= shortbase[j];
         }
   }
   for(int m = 0; m < shortCount; m++)
      if(gMask[m])
         for(int j = 0; j < wordremoval[m].Length(); j++)
            removeList.Push(m*16+wordremoval[m][j]);

      for(int m = 0; m < xshortCount; m++){
         if(xMask[m]==65535) continue;  // all SNPs masked
         fixedMask = (~xMask[m]) & 65535;
         for(int j = 0; j < 16; j++)
            hetCount[j] = AACount[j] = aaCount[j] = 0;
         for(int i = 0; i < indexCount; i++){
            id = idIndex[i];
            het = (~XG[0][id][m]) & XG[1][id][m] & fixedMask;
            AA = XG[0][id][m] & XG[1][id][m] & fixedMask;
            aa = XG[0][id][m] & (~XG[1][id][m]) & fixedMask;
            for(int j = 0; j < 16; j++)
               if(AA & shortbase[j])
                  AACount[j] ++;
               else if(aa & shortbase[j])
                  aaCount[j] ++;
               else if(het & shortbase[j])
                  hetCount[j] ++;
         }
         for(int j = 0; j < 16; j++)
            if(hetCount[j] + AACount[j] == 0 || hetCount[j] + aaCount[j] == 0){
               if((~xMask[m]) & shortbase[j])
                  wordremoval[shortCount+m].Push(j);
               xMask[m] |= shortbase[j];
            }
      }
      for(int m = 0; m < xshortCount; m++)
         if(xMask[m])
            for(int j = 0; j < wordremoval[shortCount+m].Length(); j++)
               removeList.Push((shortCount+m)*16+wordremoval[shortCount+m][j]);

      for(int m = 0; m < yshortCount; m++){
         if(yMask[m]==65535) continue;  // all SNPs masked
         fixedMask = (~yMask[m]) & 65535;
         for(int j = 0; j < 16; j++)
            hetCount[j] = AACount[j] = aaCount[j] = 0;
         for(int i = 0; i < indexCount; i++){
            id = idIndex[i];
            het = (~YG[0][id][m]) & YG[1][id][m] & fixedMask;
            AA = YG[0][id][m] & YG[1][id][m] & fixedMask;
            aa = YG[0][id][m] & (~YG[1][id][m]) & fixedMask;
            for(int j = 0; j < 16; j++)
               if(AA & shortbase[j])
                  AACount[j] ++;
               else if(aa & shortbase[j])
                  aaCount[j] ++;
               else if(het & shortbase[j])
                  hetCount[j] ++;
         }
         for(int j = 0; j < 16; j++)
            if(hetCount[j] + AACount[j] == 0 || hetCount[j] + aaCount[j] == 0){
               if((~yMask[m]) & shortbase[j])   // not masked before
                  wordremoval[shortCount+xshortCount+m].Push(j);
               yMask[m] |= shortbase[j];
            }
      }
      for(int m = 0; m < yshortCount; m++)
         if(yMask[m])
            for(int j = 0; j < wordremoval[shortCount+xshortCount+m].Length(); j++)
               removeList.Push((shortCount+xshortCount+m)*16+wordremoval[shortCount+xshortCount+m][j]);
   delete []wordremoval;
}

void Engine::CallRate_SNP(double rateFilter, IntArray & removeList)
{
   removeList.Dimension(0);
   unsigned char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   int missing, fixedMask;
   int missingCount[16];
   int id;
   IntArray idIndex(0);
   for(int i = 0; i < idCount; i++)
      if(idMask[i]==0) idIndex.Push(i);
   int indexCount = idIndex.Length();
   double cutoff = indexCount * (1-rateFilter);
   IntArray *wordremoval = new IntArray[shortCount];
   for(int m = 0; m < shortCount; m++)
      wordremoval[m].Dimension(0);
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultEfficientCoreCount) \
      private(id, missing, missingCount, fixedMask)
#endif
   for(int m = 0; m < shortCount; m++){
      if(gMask[m]==65535) continue;  // all SNPs masked
      // pool all genotypes from 16 SNPs together
      fixedMask = (~gMask[m]) & 65535;
      missing = 0;
      for(int i = 0; (i < indexCount) && (missing <= cutoff); i++){
         id = idIndex[i];
         missing += oneoneCount[(~GG[0][id][m]) & (~GG[1][id][m]) & fixedMask];
      }
      if(missing <= cutoff) continue;  // all missingness together still small
      for(int j = 0; j < 16; j++)
         missingCount[j] = 0;
      for(int i = 0; i < indexCount; i++){
         id = idIndex[i];
         missing = (~GG[0][id][m]) & (~GG[1][id][m]) & fixedMask;
         if(missing)  // if not every SNP is called for this id
            for(int j = 0; j < 16; j++)
               if(missing & shortbase[j])
                  missingCount[j] ++;
      }
      for(int j = 0; j < 16; j++)
         if(missingCount[j] > cutoff){
            if((~gMask[m]) & shortbase[j])
               wordremoval[m].Push(j);
            gMask[m] |= shortbase[j];
         }
   }
   for(int m = 0; m < shortCount; m++)
      if(gMask[m])
         for(int j = 0; j < wordremoval[m].Length(); j++)
            removeList.Push(m*16+wordremoval[m][j]);
   delete []wordremoval;
}

void Engine::CallRate_xSNP(double rateFilter, IntArray & removeList)
{
   removeList.Dimension(0);
   unsigned char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   int missing, fixedMask;
   int missingCount[16];
   int id;
   IntArray idIndex(0);
   for(int i = 0; i < idCount; i++)
      if(idMask[i]==0) idIndex.Push(i);
   int indexCount = idIndex.Length();
   double cutoff = indexCount * (1-rateFilter);
   IntArray *wordremoval = new IntArray[xshortCount];
   for(int m = 0; m < xshortCount; m++)
      wordremoval[m].Dimension(0);
   for(int m = 0; m < xshortCount; m++){
      if(xMask[m]==65535) continue;  // all SNPs masked
      fixedMask = (~xMask[m]) & 65535;
      // pool all genotypes from 16 SNPs together
      missing = 0;
      for(int i = 0; (i < indexCount) && (missing <= cutoff); i++){
         id = idIndex[i];
         missing += oneoneCount[(~XG[0][id][m]) & (~XG[1][id][m]) & fixedMask];
      }
      if(missing <= cutoff) continue;  // all missingness together still small
      for(int j = 0; j < 16; j++)
         missingCount[j] = 0;
      for(int i = 0; i < indexCount; i++){
         id = idIndex[i];
         missing = (~XG[0][id][m]) & (~XG[1][id][m]) & fixedMask;
         if(missing)  // if not every SNP is called for this id
            for(int j = 0; j < 16; j++)
               if(missing & shortbase[j])
                  missingCount[j] ++;
      }
      for(int j = 0; j < 16; j++)
         if(missingCount[j] > cutoff){
            if((~xMask[m]) & shortbase[j])
               wordremoval[m].Push(j);
            xMask[m] |= shortbase[j];
         }
   }
   for(int m = 0; m < xshortCount; m++)
      if(xMask[m])
         for(int j = 0; j < wordremoval[m].Length(); j++)
            removeList.Push(m*16+wordremoval[m][j]);
   delete []wordremoval;
}

void Engine::CallRate_ySNP(double rateFilter, IntArray & removeList, bool lessthanFlag)
{
   removeList.Dimension(0);
   unsigned char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   int missing, fixedMask;
   int missingCount[16];
   int id;
   IntArray idIndex(0);
   for(int i = 0; i < idCount; i++)
      if(idMask[i]==0) idIndex.Push(i);
   int indexCount = idIndex.Length();
   double cutoff = indexCount * (1-rateFilter);
   IntArray *wordremoval = new IntArray[yshortCount];
   for(int m = 0; m < yshortCount; m++)
      wordremoval[m].Dimension(0);
   for(int m = 0; m < yshortCount; m++){
      if(yMask[m]==65535) continue;  // all SNPs masked
      fixedMask = (~yMask[m]) & 65535;
      for(int j = 0; j < 16; j++)
         missingCount[j] = 0;
      for(int i = 0; i < indexCount; i++){
         id = idIndex[i];
         missing = (~YG[0][id][m]) & (~YG[1][id][m]) & fixedMask;
         if(missing)  // if not every SNP is called for this id
            for(int j = 0; j < 16; j++)
               if(missing & shortbase[j])
                  missingCount[j] ++;
      }
      for(int j = 0; j < 16; j++)
         if( (lessthanFlag && (missingCount[j] > cutoff) ) ||
            ((!lessthanFlag) && (missingCount[j] < cutoff) ) ){
            if((~yMask[m]) & shortbase[j])
               wordremoval[m].Push(j);
            yMask[m] |= shortbase[j];
         }
   }
   for(int m = 0; m < yshortCount; m++)
      if(yMask[m])
         for(int j = 0; j < wordremoval[m].Length(); j++)
            removeList.Push(m*16+wordremoval[m][j]);
   delete []wordremoval;
}

void Engine::CallRate_Sample(double rateFilter, IntArray & removeList)
{
   removeList.Dimension(0);
   IntArray idIndex(0);
   for(int i = 0; i < idCount; i++)
      if(idMask[i]==0) idIndex.Push(i);
   int indexCount = idIndex.Length();
   unsigned char oneoneCount[65536];
   for(int i = 0; i < 65536; i++){
      oneoneCount[i] = 0;
      for(int j = 0; j < 16; j++)
         if(i & shortbase[j]) oneoneCount[i]++;
   }
   int validCount = 0;
   for(int m = 0; m < shortCount; m++)
      validCount += oneoneCount[gMask[m]];   // invalidCount
   validCount = shortCount*16 - validCount;
   double cutoff = validCount * (1-rateFilter);

   int missingCount;
   int id;
   for(int i = 0; i < indexCount; i++){
      id = idIndex[i];
      missingCount = 0;
      for(int m = 0; (m < shortCount) && (missingCount <= cutoff); m++)
         missingCount += oneoneCount[(~GG[0][id][m])&(~GG[1][id][m])&(~gMask[m])&65535];
      if(missingCount > cutoff) removeList.Push(id);
   }
}


