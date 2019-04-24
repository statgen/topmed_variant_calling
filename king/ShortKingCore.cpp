//////////////////////////////////////////////////////////////////////
// ShortKingCore.cpp
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

#include <math.h>
#include <time.h>
#include "KingCore.h"
#include "Kinship.h"
#include "MathStats.h"
#include "MathSVD.h"
#include "QuickIndex.h"
#include "MathCholesky.h"

bool AlleleClassDifferent(char a, char b)
{
   if( ((a=='A' || a=='T') && (b=='C' || b=='G')) ||
       ((b=='A' || b=='T') && (a=='C' || a=='G')) ||
       ((a=='a' || a=='t') && (b=='c' || b=='g')) ||
       ((b=='a' || b=='t') && (a=='c' || a=='g')) )
      return true;
   else
      return false;
}

void KingEngine::autoflip()
{
   StringArray label;
   int ambiguityCount = 0;
   int ambiguity2Count = 0;
   int flipCount = 0;
   int newcode[4];
   String flipfile = prefix;
   flipfile.Add("flip.txt");
   FILE *fp = fopen(flipfile, "wt");
   fprintf(fp, "SNP\tFlipFlag\tAlleleCount\n");
   for(int m = 0; m < ped.markerCount; m++){
      if(ped.GetMarkerInfo(m)->alleleLabels.Length() == 3){
         label = ped.GetMarkerInfo(m)->alleleLabels;
         if(!AlleleClassDifferent(label[1][0], label[2][0])) { // same class
            ped.GetMarkerInfo(m)->alleleLabels.Dimension(1);
            ambiguityCount ++;
            fprintf(fp, "%s\t-9\t2\n", (const char*)ped.markerNames[m]);
         }else // different class: no allele flip needed
            fprintf(fp, "%s\t0\t2\n", (const char*)ped.markerNames[m]);
      }else if(ped.GetMarkerInfo(m)->alleleLabels.Length() == 4){ // allele flip
         ped.GetMarkerInfo(m)->alleleLabels.Dimension(1);
         ambiguity2Count++;
         fprintf(fp, "%s\t1\t3\n", (const char*)ped.markerNames[m]);
      }else if(ped.GetMarkerInfo(m)->alleleLabels.Length() == 5){ // allele flip
         label = ped.GetMarkerInfo(m)->alleleLabels;
         if( (!AlleleClassDifferent(label[1][0], label[2][0])) ||
             (!AlleleClassDifferent(label[3][0], label[4][0])) ){   // error
             ped.GetMarkerInfo(m)->alleleLabels.Dimension(1);
             ambiguity2Count ++;
             fprintf(fp, "%s\t-1\t4\n", (const char*)ped.markerNames[m]);
             continue;
         }
         fprintf(fp, "%s\t1\t4\n", (const char*)ped.markerNames[m]);
         if(AlleleClassDifferent(label[3][0], label[1][0])){
            newcode[0] = 2;
            newcode[1] = 1;
         }else{
            newcode[0] = 1;
            newcode[1] = 2;
         }
         ped.GetMarkerInfo(m)->alleleLabels.Dimension(3);
         for(int i = 0; i < ped.count; i++)
            if(ped[i].markers[m].Lo() > 2 || ped[i].markers[m].Hi() > 2)
               ped[i].markers[m].AssignGenotype(newcode[ped[i].markers[m].Lo()-3],
                  newcode[ped[i].markers[m].Hi()-3]);
         flipCount ++;
      }else // error
         fprintf(fp, "%s\t-1\t0\n", (const char*)ped.markerNames[m]);
   }
   fclose(fp);
   printf("Allele flip status is saved in file %s\n", (const char*)flipfile);
   printf("%d SNPs with genotype AT or CG are removed.\n", ambiguityCount);
   printf("%d SNPs are flipped.\n", flipCount);
   if(ambiguity2Count)
      printf("%d SNPs are flipped but cannot be fixed.\n", ambiguity2Count);
}


void KingEngine::MakePed()
{
   StringArray tokens;
   String pedfile = prefix;
   pedfile.Add("$TMP$.ped");
   String datfile = prefix;
   datfile.Add("$TMP$.dat");
   FILE *fp = fopen((const char*)pedfile, "wt");
   if(fp==NULL) error("Cannot open %s to write.", (const char*)pedfile);
   if(FID.Length() >= sampleName.Length())
      for(int i = 0; i < FID.Length(); i++){
         fprintf(fp, "%s %s %s %s %s",
         (const char*)FID[i], (const char*)PID[i], (const char*)FA[i],
         (const char*)MO[i], (const char*)SEX[i]);
         for(int t = 0; t < affectionNames.Length(); t++)
            fprintf(fp, " %d", (int)affections[t][i]);
         for(int t = 0; t < traitNames.Length(); t++)
            if(traits[t][i] == _NAN_)
               fprintf(fp, " X");
            else
               fprintf(fp, " %lf", traits[t][i]);
         for(int t = 0; t < covariateNames.Length(); t++)
            if(covariates[t][i] == _NAN_)
               fprintf(fp, " X");
            else
               fprintf(fp, " %lf", covariates[t][i]);
         if(ZG.Length()) fprintf(fp, " %s", (const char*)ZG[i]);
         fprintf(fp, "\n");
      }
   else
      for(int i = 0; i < sampleName.Length(); i++){
         int m = sampleName[i].Find("->");
         if(m>-1)
            fprintf(fp, "%s %s 0 0 0\n",
            (const char*)sampleName[i].SubStr(0, m), (const char*)sampleName[i].SubStr(m+2));
         else{
            fprintf(fp, "%s 1 0 0 0\n", (const char*)sampleName[i]);
             printf("Warning: sample %s has no ->\n", (const char*)sampleName[i]);
         }
      }
   fclose(fp);
   fp = fopen((const char*)datfile, "wt");
   if(fp==NULL) error("Cannot open %s to write.", (const char*)datfile);
   for(int t = 0; t < affectionNames.Length(); t++)
      fprintf(fp, "A %s\n", (const char*)affectionNames[t]);
   for(int t = 0; t < traitNames.Length(); t++)
      fprintf(fp, "T %s\n", (const char*)traitNames[t]);
   for(int t = 0; t < covariateNames.Length(); t++)
      fprintf(fp, "C %s\n", (const char*)covariateNames[t]);
   if(ZG.Length()) fprintf(fp, "Z twin\n");
   fclose(fp);

   ped.Prepare(datfile);
   ped.Load(pedfile);
   remove(datfile);
   remove(pedfile);

   String tempName;
   if(id) delete []id;
   id = new IntArray [ped.familyCount];
   geno.Dimension(ped.count);
   geno.Set(-1);
   phenoid.Dimension(sampleName.Length());
   phenoid.Set(-1);
   int k;
   int indexSampleName = 0;
   StringIntHash sampleHash;
   for(int i = 0; i < sampleName.Length(); i++)
      sampleHash.SetInteger(sampleName[i], i);
   for(int f = 0; f < ped.familyCount; f++){
      id[f].Dimension(0);
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
         tempName = ped.families[f]->famid;
         tempName.Add("->"); // '->'
         tempName.Add(ped[i].pid);
         k = sampleHash.Integer(tempName);
         if(k > -1){
            id[f].Push(i);
            geno[i] = k;
            phenoid[k] = i;
            ped[i].ngeno = markerCount+xmarkerCount;
         }else
            ped[i].ngeno = 0;
      }
   }
   if(phenoid.Find(-1)>-1)
      error("Sample %s's genotypes do not correspond to an internal pedigree ID",
         (const char*)sampleName[phenoid.Find(-1)]);
//   printf("Debug: ID matching ends at %s", currentTime());
}

void KingEngine::ReadKingBinaryData(const char *datafile)
{
   String line;
   StringArray tokens;
   FILE *input = fopen(datafile, "rb");
   if(input==NULL)
      error("Cannot read from file %s\n", datafile);
   printf("Loading data in KING binary format...\n");
   tokens.Clear();
   line.ReadLine(input);
   tokens.AddTokens(line);
   if(tokens.Length() < 3 || tokens.Length() > 6)
      error("Line 1: genotype file %s not in KING binary format", datafile);
   markerCount = xmarkerCount = ymarkerCount = mtmarkerCount = 0;
   int headerCount = tokens[0].AsInteger();
   printf("Loading %d lines of headers...\n", headerCount);
   idCount = tokens[1].AsInteger();
   markerCount = tokens[2].AsInteger();
   if(tokens.Length() >= 4) xmarkerCount = tokens[3].AsInteger();
   if(tokens.Length() >= 5) ymarkerCount = tokens[4].AsInteger();
   if(tokens.Length() >= 6) mtmarkerCount = tokens[5].AsInteger();

   int totalmarkerCount = markerCount + xmarkerCount + ymarkerCount + mtmarkerCount;
   if(headerCount < 1 || idCount < 1 || totalmarkerCount < 1)
      error("Line 1: genotype file %s not in KING binary format", datafile);
   printf(" %d samples, %d autosome SNPs", idCount, markerCount);
   if(xmarkerCount) printf(", %d X-chromosome SNPs", xmarkerCount);
   if(ymarkerCount) printf(", %d Y-chromosome SNPs", ymarkerCount);
   if(mtmarkerCount) printf(", %d mitochondrial SNPs", mtmarkerCount);
   printf("\n");
   shortCount = (markerCount+15)/16;
   xshortCount = (xmarkerCount+15)/16;
   yshortCount = (ymarkerCount+15)/16;
   mtshortCount = (mtmarkerCount+15)/16;

   if(!genotypeOnly){
      for(int i = 0; i < 2; i++){
         alleleLabel[i].Dimension(markerCount);
         for(int m = 0; m < markerCount; m++)
            alleleLabel[i][m] = i+1;
      }
      if(!autosomeOnly){
         for(int i = 0; i < 2; i++){
            xalleleLabel[i].Dimension(xmarkerCount);
            for(int m = 0; m < xmarkerCount; m++)
               xalleleLabel[i][m] = i+1;
         }
         for(int i = 0; i < 2; i++){
            yalleleLabel[i].Dimension(ymarkerCount);
            for(int m = 0; m < ymarkerCount; m++)
               yalleleLabel[i][m] = i+1;
         }
         for(int i = 0; i < 2; i++){
            mtalleleLabel[i].Dimension(mtmarkerCount);
            for(int m = 0; m < mtmarkerCount; m++)
               mtalleleLabel[i][m] = i+1;
         }
      }
   }
   FID.Dimension(0);
   PID.Dimension(0);
   FA.Dimension(0);
   MO.Dimension(0);
   SEX.Dimension(0);
   ZG.Dimension(0);
   chromosomes.Dimension(0);
   bp.Dimension(0);  //positions.Dimension(0);
   snpName.Dimension(0);
   xbp.Dimension(0); //xpositions.Dimension(0);
   xsnpName.Dimension(0);
   ybp.Dimension(0); //ypositions.Dimension(0);
   ysnpName.Dimension(0);
   mtbp.Dimension(0);   //mtpositions.Dimension(0);
   mtsnpName.Dimension(0);
   sampleName.Dimension(0);
   affectionNames.Dimension(0);
   traitNames.Dimension(0);
   covariateNames.Dimension(0);
   StringArray bufferT(0), bufferA(0), bufferC(0);

   int buffersize = 60000;
   if(totalmarkerCount > 2000)
      buffersize = totalmarkerCount * 30;
   if(idCount > totalmarkerCount)
      buffersize = idCount * 30;
   char *buffer = new char[buffersize];
   for(int k = 2; k <= headerCount; k++){
      fgets(buffer, buffersize, input);
      if(buffer[0]=='\0'){
         printf("  Line %d is skipped.\n", k);
         continue;
      }
      if(genotypeOnly && (buffer[0] == 'M' || buffer[0] == 'T' || buffer[0] == 'A' || buffer[0] == 'C'))
         continue;
      tokens.Clear();
      tokens.AddTokens(buffer);
      if(tokens.Length() < 3 || tokens[0].Length() != 1 ||
         (tokens[0][0]!= 'M' && tokens[0][0]!='S' && tokens[0][0]!='T'
            && tokens[0][0]!='A' && tokens[0][0]!='C' && tokens[0][0]!='Z')){
         printf("  Line %d (%s %s) is skipped for lack of detailed information.\n",
            k, (const char*)tokens[0], (const char*)tokens[1]);
         continue;
      }
      switch(tokens[0][0]){
         case 'M':   // Marker
            if(tokens[1] == "KING_ALLELE"){
               if(tokens.Length()!=2*(markerCount+xmarkerCount+ymarkerCount+mtmarkerCount)+2)
                  error("Line %d: total number of allele labels %d is different from twice of SNP count %d",
                     k, tokens.Length(), markerCount+xmarkerCount+ymarkerCount+mtmarkerCount);
            }else if(tokens[1] == "KING_CHROMOSOME"){
               if(tokens.Length() < markerCount+2)
                  error("Line %d: total number of SNPs %d is less than %d",
                     k, tokens.Length()-2, markerCount);
            }else{
               if(tokens.Length() != markerCount+xmarkerCount+ymarkerCount+mtmarkerCount+2)
                  error("Line %d: total number of SNPs %d is different from %d",
                     k, tokens.Length()-2, markerCount+xmarkerCount+ymarkerCount+mtmarkerCount);
            }
            if(tokens[1] == "KING_SNP"){ // snpName
               for(int i = 0; i < markerCount; i++)
                  snpName.Push(tokens[i+2]);
               for(int i = 0; i < xmarkerCount; i++)
                  xsnpName.Push(tokens[markerCount+2+i]);
               for(int i = 0; i < ymarkerCount; i++)
                  ysnpName.Push(tokens[markerCount+xmarkerCount+2+i]);
               for(int i = 0; i < mtmarkerCount; i++)
                  mtsnpName.Push(tokens[markerCount+xmarkerCount+ymarkerCount+2+i]);
            }else if(tokens[1] == "KING_ALLELE"){ // allele label
               for(int m = 0; m < markerCount; m++){
                  alleleLabel[0][m] = tokens[m*2+2];
                  alleleLabel[1][m] = tokens[m*2+3];
               }
               for(int m = 0; m < xmarkerCount; m++){
                  xalleleLabel[0][m] = tokens[markerCount*2+m*2+2];
                  xalleleLabel[1][m] = tokens[markerCount*2+m*2+3];
               }
               for(int m = 0; m < ymarkerCount; m++){
                  yalleleLabel[0][m] = tokens[(markerCount+xmarkerCount)*2+m*2+2];
                  yalleleLabel[1][m] = tokens[(markerCount+xmarkerCount)*2+m*2+3];
               }
               for(int m = 0; m < mtmarkerCount; m++){
                  mtalleleLabel[0][m] = tokens[(markerCount+xmarkerCount+ymarkerCount)*2+m*2+2];
                  mtalleleLabel[1][m] = tokens[(markerCount+xmarkerCount+ymarkerCount)*2+m*2+3];
               }
            }else if(tokens[1] == "KING_CHROMOSOME")
               for(int i = 2; i < markerCount+2; i++)
                  chromosomes.Push(tokens[i].AsInteger());
            else if(tokens[1] == "KING_POSITION"){
               for(int i = 2; i < markerCount+2; i++){
                  bp.Push(int(tokens[i].AsDouble()*1000000+0.5)); //positions.Push(tokens[i].AsDouble());
               }
               for(int i = markerCount+2; i < markerCount+xmarkerCount+2; i++){
                  xbp.Push(int(tokens[i].AsDouble()*1000000+0.5));   //xpositions.Push(tokens[i].AsDouble());
               }
               for(int i = markerCount+xmarkerCount+2; i < markerCount+xmarkerCount+ymarkerCount+2; i++){
                  ybp.Push(int(tokens[i].AsDouble()*1000000+0.5));   //ypositions.Push(tokens[i].AsDouble());
               }
               for(int i = markerCount+xmarkerCount+ymarkerCount+2; i < markerCount+xmarkerCount+ymarkerCount+mtmarkerCount+2; i++){
                  mtbp.Push(int(tokens[i].AsDouble()*1000000+0.5));  //mtpositions.Push(tokens[i].AsDouble()); 
               }
            }
            break;
         case 'S':
            if(tokens[1] == "KING_SAMPLE"){   // sample name
               if(tokens.Length() != idCount+2)
                  error("Line 4: total number of samples %d is different from %d",
                     tokens.Length()-2, idCount);
               for(int i = 2; i < tokens.Length(); i++)
                  sampleName.Push(tokens[i]);
            }else if(tokens[1] == "KING_FID")
               for(int i = 2; i < tokens.Length(); i++)
                  FID.Push(tokens[i]);
            else if(tokens[1] == "KING_PID")
               for(int i = 2; i < tokens.Length(); i++)
                  PID.Push(tokens[i]);
            else if(tokens[1] == "KING_FATHER")
               for(int i = 2; i < tokens.Length(); i++)
                  FA.Push(tokens[i]);
            else if(tokens[1] == "KING_MOTHER")
               for(int i = 2; i < tokens.Length(); i++)
                  MO.Push(tokens[i]);
            else if(tokens[1] == "KING_SEX")
               for(int i = 2; i < tokens.Length(); i++)
                  SEX.Push(tokens[i]);
            else
               printf("  Line %d skipped.\n", k);
            break;
         case 'T':
            bufferT.Push(line);
            break;
         case 'A':
            bufferA.Push(line);
            break;
         case 'C':
            bufferC.Push(line);
            break;
         case 'Z': // zygosity
            for(int i = 2; i < tokens.Length(); i++)
               ZG.Push(tokens[i]);
            break;
         default: printf("%c not defined", tokens[0][0]);
         break;
      }
   }
   delete []buffer;
   if(bufferT.Length()){
      tokens.Clear();
      tokens.AddTokens(bufferT[0]);
      traits.Dimension(bufferT.Length(), tokens.Length()-2);
      for(int t = 0; t < bufferT.Length(); t++){
         tokens.Clear();
         tokens.AddTokens(bufferT[t]);
         traitNames.Push(tokens[1]);
         for(int i = 2; i < tokens.Length(); i++)
            if(tokens[i] != "X" && tokens[i] != "x")
               traits[t][i-2] = tokens[i].AsDouble();
            else
               traits[t][i-2] = _NAN_;
      }
   }
   if(bufferC.Length()){
      tokens.Clear();
      tokens.AddTokens(bufferC[0]);
      covariates.Dimension(bufferC.Length(), tokens.Length()-2);
      for(int t = 0; t < bufferC.Length(); t++){
         tokens.Clear();
         tokens.AddTokens(bufferC[t]);
         covariateNames.Push(tokens[1]);
         for(int i = 2; i < tokens.Length(); i++)
            if(tokens[i] != "X" && tokens[i] != "x")
               covariates[t][i-2] = tokens[i].AsDouble();
            else
               covariates[t][i-2] = _NAN_;
      }
   }
   if(bufferA.Length()){
      tokens.Clear();               
      tokens.AddTokens(bufferA[0]);
      affections.Dimension(bufferA.Length(), tokens.Length()-2);
      for(int t = 0; t < bufferA.Length(); t++){
         tokens.Clear();
         tokens.AddTokens(bufferA[t]);
         affectionNames.Push(tokens[1]);
         for(int i = 2; i < tokens.Length(); i++)
            if(tokens[i] != "X" && tokens[i] != "x")
               affections[t][i-2] = tokens[i].AsInteger();
            else
               affections[t][i-2] = 0;
      }
   }
   if(sampleName.Length()==0){
      String KING("KING");
      KING.Add("->");  // '->'
      String temp1, temp2;
      for(int i = 0; i < idCount; i++){
         temp2 = i;
         temp1 = KING;
         sampleName.Push(temp1.Add(temp2));
      }
   }

   // line: genotypes per sample
   String king;
   king.Clear();
   for(int i = 0; i < 4; i++) king.Add(fgetc(input));
   if(king != "KING")
      error("Line %d: %s", headerCount+1, (const char*)king);
   printf("Loading genotypes...\n");
   fgetc(input); fgetc(input);   // skip KING format version
   for(int j = 0; j < 2; j++){
      GG[j] = new unsigned short int * [idCount];
      for(int i = 0; i < idCount; i++)
         GG[j][i] = new unsigned short int [shortCount];
   }
   if(autosomeOnly){
      if(xshortCount+yshortCount+mtshortCount>0)
         for(int i = 0; i < idCount; i++){
            fread(GG[0][i], sizeof(unsigned short int), shortCount, input);
            fread(GG[1][i], sizeof(unsigned short int), shortCount, input);
            fseek(input, (xshortCount+yshortCount+mtshortCount)*2*sizeof(unsigned short int), SEEK_CUR);
         }
      else
         for(int i = 0; i < idCount; i++){
            fread(GG[0][i], sizeof(unsigned short int), shortCount, input);
            fread(GG[1][i], sizeof(unsigned short int), shortCount, input);
         }
   }else{
      if(xshortCount)
         for(int j = 0; j < 2; j++){
            XG[j] = new unsigned short int * [idCount];
            for(int i = 0; i < idCount; i++)
               XG[j][i] = new unsigned short int [xshortCount];
         }
      if(yshortCount)
         for(int j = 0; j < 2; j++){
            YG[j] = new unsigned short int * [idCount];
            for(int i = 0; i < idCount; i++)
               YG[j][i] = new unsigned short int [yshortCount];
         }
      if(mtshortCount)
         for(int j = 0; j < 2; j++){
            MG[j] = new unsigned short int * [idCount];
            for(int i = 0; i < idCount; i++)
               MG[j][i] = new unsigned short int [mtshortCount];
         }
      for(int i = 0; i < idCount; i++){
         if(shortCount){
            fread(GG[0][i], sizeof(unsigned short int), shortCount, input);
            fread(GG[1][i], sizeof(unsigned short int), shortCount, input);
         }
         if(xshortCount){
            fread(XG[0][i], sizeof(unsigned short int), xshortCount, input);
            fread(XG[1][i], sizeof(unsigned short int), xshortCount, input);
         }
         if(yshortCount){
            fread(YG[0][i], sizeof(unsigned short int), yshortCount, input);
            fread(YG[1][i], sizeof(unsigned short int), yshortCount, input);
         }
         if(mtshortCount){
            fread(MG[0][i], sizeof(unsigned short int), mtshortCount, input);
            fread(MG[1][i], sizeof(unsigned short int), mtshortCount, input);
         }
      }
   }
   fclose(input);
   printf("Genotype data loaded.\n");
   MakePed();
}

void KingEngine::AfterBinaryLoaded(void)
{
   unsigned short int AA, Aa, aa, missing, wipe;
   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++){
      oneoneCount[i] = 0;
      for(int j = 0; j < 16; j++)
         if(i & shortbase[j]) oneoneCount[i]++;
   }
   int count = 0;
   for(int i = 0; i < ped.count; i++){
      ped[i].ngeno = 0;
      if(geno[i]!=-1)
         for(int b = 0; b < shortCount; b++)
            ped[i].ngeno += oneoneCount[GG[0][geno[i]][b] | GG[1][geno[i]][b]];
      if(ped[i].ngeno) count++;
   }
   if(count < idCount)
      printf("Among %d samples with genotypes, %d have 100%% missing genotypes on autosomes.\n",
         idCount, idCount - count);
   if(genoFilterFlag){
      int wipeCount = 0;
      for(int b = 0; b < shortCount; b++){
         AA = Aa = aa = 0;
         for(int i = 0; i < idCount; i++){
            AA |= GG[0][i][b] & GG[1][i][b];
            Aa |= (~GG[0][i][b]) & GG[1][i][b];
            aa |= GG[0][i][b] & (~GG[1][i][b]);
         }
         wipe = ~(AA & Aa & aa);
         if(wipe)
            for(int i = 0; i < idCount; i++){
               GG[0][i][b] &= (~wipe);
               GG[1][i][b] &= (~wipe);
            }
         wipeCount += oneoneCount[wipe];
      }
      if(wipeCount)
         printf("  %d SNPs with less than 3 genotypes are wiped out.\n", wipeCount);
   }
   if(individualInfo || detailFlag || diagFlag){
      int flipCount[16]; // AACount - aaCount;
      shortFlip = new unsigned short int[shortCount];
      if(shortFlip == NULL) error("Cannot allocate memory for flip signs.");
      for(int b = 0; b < shortCount; b++){
         shortFlip[b] = 0;
         for(int i = 0; i < 16; i++) flipCount[i] = 0;
         for(int i = 0; i < idCount; i++){
            AA = GG[0][i][b] & GG[1][i][b];
            aa = GG[0][i][b] & (~GG[1][i][b]);
            for(int i = 0; i < 16; i++)
               if(AA & shortbase[i])
                  flipCount[i] ++;
               else if(aa & shortbase[i])
                  flipCount[i] --;
         }
         for(int i = 0; i < 16; i++)
            if(flipCount[i] > 0)
               shortFlip[b] |= shortbase[i];
      }
   }
   if(mafFilterFlag){
      Vector frequencies;
      ComputeAlleleFrequency(frequencies);
      int wipeCount = 0;
      int pos;
      for(int b = 0; b < shortCount; b++){
         if(mafFilterFlag && minMAF > 1E-10){
            wipe = 0;
            for(int j = 0; j < 16; j++)
               if(frequencies[pos+j] > 1-minMAF || frequencies[pos+j] < minMAF){
                  wipe |= shortbase[j];
                  frequencies[pos+j] = 0;
               }
            if(wipe){
               for(int i = 0; i < idCount; i++){
                  GG[0][i][b] &= (~wipe);
                  GG[1][i][b] &= (~wipe);
               }
               wipeCount += oneoneCount[wipe];
            }
         }
      }
      if(wipeCount)
         printf("  %d SNPs with MAF < %lf are wiped out.\n", wipeCount, minMAF);
      printf("Allele frequencies by SNP are calculated.\n");
   }
   if(callrate != _NAN_){
      int wipeCount = 0;
      int missingCount[16];
      for(int b = 0; b < shortCount; b++){
         for(int i = 0; i < 16; i++) missingCount[i] = 0;
         for(int i = 0; i < idCount; i++){
            missing = (~(GG[0][i][b] | GG[1][i][b])) & 65535;
            if(missing)
               for(int i = 0; i < 16; i++)
                  if(missing & shortbase[i])
                     missingCount[i] ++;
         }
         wipe = 0;
         for(int i = 0; i < 16; i++)
            if(missingCount[i] > idCount*(1-callrate))
               wipe |= shortbase[i];
         if(wipe){
            for(int i = 0; i < idCount; i++){
               GG[0][i][b] &= (~wipe);
               GG[1][i][b] &= (~wipe);
            }
            wipeCount += oneoneCount[wipe];
         }
      }
      if(wipeCount)
         printf("  %d autosome SNPs with call rate < %.1lf are wiped out.\n", wipeCount, callrate*100);
   }

   if(xshortCount){
      Vector xfrequencies;
      ComputeAlleleFrequencyInX(xfrequencies);
      int pos;
      int wipeCount = 0;
      for(int b = 0; b < xshortCount; b++){
         if(mafFilterFlag && minMAF > 1E-10){
            wipe = 0;
            for(int j = 0; j < 16; j++)
               if(xfrequencies[pos+j] > 1-minMAF || xfrequencies[pos+j] < minMAF){
                  wipe |= shortbase[j];
                  xfrequencies[pos+j] = 0;
               }
            if(wipe){
               for(int i = 0; i < idCount; i++){
                  XG[0][i][b] &= (~wipe);
                  XG[1][i][b] &= (~wipe);
               }
               wipeCount += oneoneCount[wipe];
            }
         }
      }
      if(wipeCount)
         printf("  %d X-chromosome SNPs with MAF < %lf are wiped out.\n",
            wipeCount, minMAF);
   }

   if(!pedstatFlag) return;
   familyFlag = false;
   for(int f = 0; f < ped.familyCount; f++){
      int count = 0;
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(ped[i].ngeno) count++;
      if(count > 1) {
         familyFlag = true;
         break;
      }
   }
   if(familyFlag){
      IntArray pedsize(11);
      pedsize.Zero();
      int count;
      int maxcount = 0;
      int maxf = -1;
      for(int f = 0; f < ped.familyCount; f++){
         count = 0;
         for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
            if(ped[i].ngeno) count++;
         if(count == 0) continue;
         if(count < 11) pedsize[count-1]++;
         else pedsize[10]++;
         if(count > maxcount) {
            maxcount = count;
            maxf = f;
         }
      }
      printf("\nPedigree information (%d families, before relationship inference)\n",
         pedsize.Sum());
      printf("%10s", "PedSize");
      for(int i = 1; i <= 10 && i <= maxcount; i++)
         printf("%5d", i);
      if(maxcount > 10) printf("%5s", "11+");
      printf("\n");
      printf(" =========");
      for(int i = 1; i <= 11 && i <= maxcount; i++)
         printf("=====");
      printf("\n");
      printf("%10s", "PedCount");
      for(int i = 0; i < 11 && i < maxcount; i++)
         printf("%5d", pedsize[i]);
      printf("\n\n");
      if(maxcount > 10)
         printf("The largest pedigree is family %s, with %d individuals\n\n",
         (const char*)ped.families[maxf]->famid , maxcount);

      IntArray nuclear[3];
      for(int i = 0; i < 3; i++){
         nuclear[i].Dimension(6);
         nuclear[i].Zero();
      }
      int pcount, affcount;

      String outfile;
      outfile = prefix;
      outfile.Add("faminfo.txt");
      FILE *fp = fopen((const char*)outfile, "wt");
      if(fp==NULL) error("Cannot open %s to write", (const char*)outfile);
      fprintf(fp, "FAMID\tPedSize\tGenotyped\tParentCount\tOffspringCount");
      if(ped.affectionCount) fprintf(fp, "\tAffectedOffspringCount");
      fprintf(fp, "\n");
      for(int f=0; f < ped.familyCount; f++)
         for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
            if(ped[i].sibCount && ped[i].sibs[0]->serial == i){
               pcount = 0;
               if(ped[i].father->ngeno) pcount++;
               if(ped[i].mother->ngeno) pcount++;
               count = ped[i].sibCount;

               affcount = 0;
               if(ped[i].affectionCount > 0)
               for(int k = 0; k < ped[i].sibCount; k++)
                  if(ped[i].sibs[k]->affections[0]==2)
                     affcount++;
               if(count > 1 || pcount){
                  nuclear[pcount][(count>5?6:count)-1]++;
                  fprintf(fp, "%s\t%d\t%d\t%d\t%d\t%d\n",
                  (const char*)ped.families[f]->famid, ped.families[f]->count,
                  count+pcount, pcount, count, affcount);
               }
            }
      fclose(fp);
      printf("Nuclear family information (#parent vs. #offspring)\n");
      printf("%10s%5s%5s%5s%5s%5s%5s%10s\n",
         "#parents", "1","2","3","4","5","6+", "Total");
      printf("  =================================================\n");
      int total;
      for(int i = 0; i < 3; i++){
         printf("%10d", i);
         total = 0;
         for(int j = 0; j < 6; j++){
            printf("%5d", nuclear[i][j]);
            total += nuclear[i][j];
         }
         printf("%10d", total);
         printf("\n");
      }
      printf("Nuclear family information saved in file %s\n", (const char*)outfile);
      printf("\n");

      if(ped.affectionCount==0) return;
      // relatives
      IntArray affrelCount(ped.count);
      affrelCount.Zero();
      IntArray unaffrelCount(ped.count);
      unaffrelCount.Zero();
      for(int f = 0; f < ped.familyCount; f++){
         for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
            if(ped[i].sibCount && ped[i].sibs[0]->serial == i){
               int N_affParent = 0;
               int N_affSibling = 0;
               int N_unaffParent = 0;
               int N_unaffSibling = 0;
               if(ped[i].father){
                  if(ped[i].father->affections[0]==2){
                     N_affParent++;
                  }else if(ped[i].father->affections[0]==1){
                     N_unaffParent++;
                  }
               }
               if(ped[i].mother){
                  if(ped[i].mother->affections[0]==2)
                     N_affParent++;
                  else if(ped[i].mother->affections[0]==1)
                     N_unaffParent++;
               }
               for(int s = 0; s < ped[i].sibCount; s++){
                  if(ped[i].sibs[s]->affections[0]==2) {
                     N_affSibling++;
                     if(ped[i].sibs[s]->father)
                        affrelCount[ped[i].sibs[s]->father->serial]++;
                     if(ped[i].sibs[s]->mother)
                        affrelCount[ped[i].sibs[s]->mother->serial]++;
                  }else if(ped[i].sibs[s]->affections[0]==1){
                     N_unaffSibling++;
                     if(ped[i].sibs[s]->father)
                        unaffrelCount[ped[i].sibs[s]->father->serial]++;
                     if(ped[i].sibs[s]->mother)
                        unaffrelCount[ped[i].sibs[s]->mother->serial]++;
                  }
               }
               for(int s = 0; s < ped[i].sibCount; s++){
                  if(ped[i].sibs[s]->affections[0]==2)
                     affrelCount[ped[i].sibs[s]->serial] += N_affParent + N_affSibling - 1;
                  else
                     affrelCount[ped[i].sibs[s]->serial] += N_affParent + N_affSibling;
                  if(ped[i].sibs[s]->affections[0]==1)
                     unaffrelCount[ped[i].sibs[s]->serial] += N_unaffParent + N_unaffSibling - 1;
                  else
                     unaffrelCount[ped[i].sibs[s]->serial] += N_unaffParent + N_unaffSibling;
               }
            }
      }
      outfile = prefix;
      outfile.Add("relinfo.txt");
      fp = fopen((const char*)outfile, "wt");
      if(fp==NULL) error("Cannot open %s to write", (const char*)outfile);
      fprintf(fp, "FID\tIID\tFA\tMO\tSEX\tN_aff\tN_unaff\n");
      for(int i = 0; i < ped.count; i++)
         fprintf(fp, "%s\t%s\t%s\t%s\t%d\t%d\t%d\n",
            (const char*)ped[i].famid, (const char*)ped[i].pid,
            (const char*)ped[i].fatid, (const char*)ped[i].motid,
            ped[i].sex, affrelCount[i], unaffrelCount[i]);
      fclose(fp);
      printf("# of 1st-degree relative information saved in file %s.\n\n", (const char*)outfile);
   }else
      printf("The data consist of %d unrelated individuals\n", ped.familyCount);
}


void KingEngine::BuildShortBinary()
{
   IntArray markers, xmarkers;
   
   id = new IntArray [ped.familyCount];
   geno.Dimension(ped.count);
   geno.Set(-1);
   Vector frequencies(markerCount);
   frequencies.Zero();
   Vector xfrequencies(xmarkerCount);
   xfrequencies.Zero();

   IntArray tooShortList(0);
   idCount = 0;
   for(int f = 0; f < ped.familyCount; f++){
      id[f].Dimension(0);
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
//         if(ped[i].ngeno >= MINSNPCOUNT || ((ibsFlag||QCbySample||QCbySNP) && ped[i].ngeno > 0)){
         if(ped[i].ngeno >= MINSNPCOUNT ||
         ((!(allflags&(1<<KinshipFLAG)) || (allflags&(1<<BysampleFLAG)) || (allflags&(1<<BysnpFLAG))) && ped[i].ngeno > 0)){
            id[f].Push(i);
            geno[i] = idCount;
            idCount++;
         }else if(ped[i].ngeno > 0)
            tooShortList.Push(i);
   }
   if(tooShortList.Length()){
      printf("The following %d samples are removed for having less than %d SNPs.\n",
         tooShortList.Length(), MINSNPCOUNT);
      for(int i = 0; i < tooShortList.Length(); i++)
         printf("%s->%s\t",
            (const char*)ped[tooShortList[i]].famid,
            (const char*)ped[tooShortList[i]].pid);
      printf("\n");
   }
   if(idCount == 0)
      error("There are no samples with sufficient amount of SNPs for your analysis.\nThe only analysis that allows low throughput SNP data is --ibs.");
   else if(idCount == 1/* && (analysisFlag || !notrimFlag)*/)
      warning("There is only 1 sample in the data. It's recommended to have at least 2 samples.");

   markerCount = xmarkerCount = ymarkerCount = mtmarkerCount = 0;
//  frequencies = xfrequencies = NULL;
   markers.Dimension(0);
   xmarkers.Dimension(0);
   IntArray ymarkers(0);
   IntArray mtmarkers(0);

   if(autoflipFlag) autoflip();

   int nonsnp0Count = 0; int nonsnp3Count = 0; int nonAutosomeCount = 0;
   int notspecifiedCount = 0;
   int chr;
   if( mafFilterFlag || (individualInfo && !freqLoaded) || callrate != _NAN_ || diagFlag){
      Vector freqs(0);
      Vector xfreqs(0);
      double freq; int freqCount; int mafCount = 0; int lowcallCount = 0;
      for(int m = 0; m < ped.markerCount; m++){
         if(ped.GetMarkerInfo(m)->alleleLabels.Length() < 3) {
            nonsnp0Count++;
            continue;
         }else if(ped.GetMarkerInfo(m)->alleleLabels.Length() > 3) {
            nonsnp3Count++;
            continue;
         }
         chr = ped.GetMarkerInfo(m)->chromosome;
         double pos = ped.GetMarkerInfo(m)->position*100;
         if(chr == 999)
            ped.GetMarkerInfo(m)->chromosome = chr = SEXCHR;
         if(chr > SEXCHR+3 || chr < -1) {
            nonAutosomeCount++;
            continue;// only autosomes
         }
         if(chrList.Length() && (chrList.Find(chr)==-1
            || (start != _NAN_ && pos < start)
            || (stop != _NAN_ && pos > stop) )){
            notspecifiedCount = 0;
            continue;
         }
         freq = 0;
         freqCount = 0;
         for(int f = 0; f < ped.familyCount; f++)
            for(int i = 0; i < id[f].Length(); i++)
               if(ped[id[f][i]].markers[m].isKnown()){
                  freq += ped[id[f][i]].markers[m].countAlleles(1);
                  freqCount += 2;
               }
         if(freqCount) freq /= freqCount;
         if(mafFilterFlag)
            if(freq < minMAF || freq > 1 - minMAF) {mafCount++; continue;}
         if(callrate != _NAN_ && chr != 24)
            if(freqCount/2 < idCount*callrate) {lowcallCount++; continue;}
         if(chr == SEXCHR){
            xmarkers.Push(m);
            xfreqs.Push(freq);
            xmarkerCount++;
         }else if(chr == SEXCHR+1){ // Y-chromosome
            ymarkers.Push(m);
            ymarkerCount++;
         }else if(chr == SEXCHR+3){ // MT
            mtmarkers.Push(m);
            mtmarkerCount++;
         }else{   // autosome
            markers.Push(m);
            freqs.Push(freq);
            markerCount++;
         }
      }
      if(mafFilterFlag || individualInfo || diagFlag){ // not for callrate
         for(int m = 0; m < markerCount; m++)
            frequencies[m] = freqs[m];
         for(int m = 0; m < xmarkerCount; m++)
            xfrequencies[m] = xfreqs[m];
      }
      if(mafCount)
         printf("  %d SNPs with MAF < %lf are removed.\n", mafCount, minMAF);
      if(lowcallCount)
         printf("  %d SNPs with call rate < %lf are removed.\n", lowcallCount, callrate);
   }else
      for(int m = 0; m < ped.markerCount; m++){
         if(ped.GetMarkerInfo(m)->alleleLabels.Length() < 3) {
//            if(analysisFlag || !notrimFlag){
               nonsnp0Count++;
               continue;
//            }
         }else if(ped.GetMarkerInfo(m)->alleleLabels.Length() > 3) {
            nonsnp3Count++;
            continue;
         }
         chr = ped.GetMarkerInfo(m)->chromosome;
         double pos = ped.GetMarkerInfo(m)->position*100;
         if(chr == 999)
            ped.GetMarkerInfo(m)->chromosome = chr = SEXCHR;
         if(chr > SEXCHR+3 || chr < -1){
            nonAutosomeCount ++;
            continue;// only autosomes
         }
         if(chrList.Length() && (chrList.Find(chr)==-1
            || (start != _NAN_ && pos < start)
            || (stop != _NAN_ && pos > stop) )){
            notspecifiedCount = 0;
            continue;
         }
         if(chr == SEXCHR){
            xmarkers.Push(m);
            xmarkerCount++;
         }else if(chr == SEXCHR+1){ // Y-chromosome
            ymarkers.Push(m);
            ymarkerCount++;
         }else if(chr == SEXCHR+3){ // MT
            mtmarkers.Push(m);
            mtmarkerCount++;
         }else{
            markers.Push(m);
            markerCount++;
         }
      }

   for(int i = 0; i < 2; i++){
      if(markerCount) alleleLabel[i].Dimension(markerCount);
      if(xmarkerCount) xalleleLabel[i].Dimension(xmarkerCount);
      if(ymarkerCount) yalleleLabel[i].Dimension(ymarkerCount);
      if(mtmarkerCount) mtalleleLabel[i].Dimension(mtmarkerCount);
   }
   bp.Dimension(0);  //positions.Dimension(0);
   xbp.Dimension(0); //xpositions.Dimension(0);
   ybp.Dimension(0); //ypositions.Dimension(0);
   mtbp.Dimension(0);   //mtpositions.Dimension(0); 
   snpName.Dimension(markerCount);
   if(xmarkerCount) xsnpName.Dimension(xmarkerCount);
   if(ymarkerCount) ysnpName.Dimension(ymarkerCount);
   if(mtmarkerCount) mtsnpName.Dimension(mtmarkerCount);
   for(int m = 0; m < markerCount; m++){
      chr = ped.GetMarkerInfo(markers[m])->chromosome;
      if(chr < 0)
         chromosomes.Push(0);
      else
         chromosomes.Push(chr);
      //positions.Push(ped.GetMarkerInfo(markers[m])->position*100);
      bp.Push(int(ped.GetMarkerInfo(markers[m])->position*100000000+0.5));
      snpName[m] = ped.markerNames[markers[m]];
   }
/*
   if(analysisFlag || !notrimFlag)
      for(int m = 0; m < markerCount; m++){
         alleleLabel[0][m] = ped.GetMarkerInfo(markers[m])->alleleLabels[1];
         alleleLabel[1][m] = ped.GetMarkerInfo(markers[m])->alleleLabels[2];
      }
   else{ */
      for(int m = 0; m < markerCount; m++)
         if(ped.GetMarkerInfo(markers[m])->alleleLabels.Length() > 2){
            alleleLabel[0][m] = ped.GetMarkerInfo(markers[m])->alleleLabels[1];
            alleleLabel[1][m] = ped.GetMarkerInfo(markers[m])->alleleLabels[2];
         }else if(ped.GetMarkerInfo(markers[m])->alleleLabels.Length()==2){
            alleleLabel[0][m] = ped.GetMarkerInfo(markers[m])->alleleLabels[1];
            alleleLabel[1][m] = "0";
         }else{
            alleleLabel[0][m] = "0";
            alleleLabel[1][m] = "0";
         }
//   }
   for(int m = 0; m < xmarkerCount; m++){
      //xpositions.Push(ped.GetMarkerInfo(xmarkers[m])->position*100);
      xbp.Push(int(ped.GetMarkerInfo(xmarkers[m])->position*100000000+0.5));
      xsnpName[m] = ped.markerNames[xmarkers[m]];
   }
/*   if(analysisFlag || !notrimFlag)
      for(int m = 0; m < xmarkerCount; m++){
         xalleleLabel[0][m] = ped.GetMarkerInfo(xmarkers[m])->alleleLabels[1];
         xalleleLabel[1][m] = ped.GetMarkerInfo(xmarkers[m])->alleleLabels[2];
      }
   else{ */
      for(int m = 0; m < xmarkerCount; m++)
         if(ped.GetMarkerInfo(xmarkers[m])->alleleLabels.Length() > 2){
            xalleleLabel[0][m] = ped.GetMarkerInfo(xmarkers[m])->alleleLabels[1];
            xalleleLabel[1][m] = ped.GetMarkerInfo(xmarkers[m])->alleleLabels[2];
         }else if(ped.GetMarkerInfo(xmarkers[m])->alleleLabels.Length()==2){
            xalleleLabel[0][m] = "0";
            xalleleLabel[1][m] = ped.GetMarkerInfo(xmarkers[m])->alleleLabels[1];
         }else{
            xalleleLabel[0][m] = "0";
            xalleleLabel[1][m] = "0";
         }
//   }
   for(int m = 0; m < ymarkerCount; m++){
      //ypositions.Push(ped.GetMarkerInfo(ymarkers[m])->position*100);
      ybp.Push(int(ped.GetMarkerInfo(ymarkers[m])->position*100000000+0.5));
      ysnpName[m] = ped.markerNames[ymarkers[m]];
   }
/*   if(analysisFlag || !notrimFlag)
      for(int m = 0; m < ymarkerCount; m++){
         yalleleLabel[0][m] = ped.GetMarkerInfo(ymarkers[m])->alleleLabels[1];
         yalleleLabel[1][m] = ped.GetMarkerInfo(ymarkers[m])->alleleLabels[2];
      }
   else{*/
      for(int m = 0; m < ymarkerCount; m++)
         if(ped.GetMarkerInfo(ymarkers[m])->alleleLabels.Length() > 2){
            yalleleLabel[0][m] = ped.GetMarkerInfo(ymarkers[m])->alleleLabels[1];
            yalleleLabel[1][m] = ped.GetMarkerInfo(ymarkers[m])->alleleLabels[2];
         }else if(ped.GetMarkerInfo(ymarkers[m])->alleleLabels.Length()==2){
            yalleleLabel[0][m] = "0";
            yalleleLabel[1][m] = ped.GetMarkerInfo(xmarkers[m])->alleleLabels[1];
         }else{
            yalleleLabel[0][m] = "0";
            yalleleLabel[1][m] = "0";
         }
//   }

   for(int m = 0; m < mtmarkerCount; m++){
      //mtpositions.Push(ped.GetMarkerInfo(mtmarkers[m])->position*100);
      mtbp.Push(int(ped.GetMarkerInfo(mtmarkers[m])->position*100000000+0.5));
      mtsnpName[m] = ped.markerNames[mtmarkers[m]];
   }
/*   if(analysisFlag || !notrimFlag)
      for(int m = 0; m < mtmarkerCount; m++){
         mtalleleLabel[0][m] = ped.GetMarkerInfo(mtmarkers[m])->alleleLabels[1];
         mtalleleLabel[1][m] = ped.GetMarkerInfo(mtmarkers[m])->alleleLabels[2];
      }
   else{*/
      for(int m = 0; m < mtmarkerCount; m++)
         if(ped.GetMarkerInfo(mtmarkers[m])->alleleLabels.Length() > 2){
            mtalleleLabel[0][m] = ped.GetMarkerInfo(mtmarkers[m])->alleleLabels[1];
            mtalleleLabel[1][m] = ped.GetMarkerInfo(mtmarkers[m])->alleleLabels[2];
         }else if(ped.GetMarkerInfo(mtmarkers[m])->alleleLabels.Length()==2){
            mtalleleLabel[0][m] = "0";
            mtalleleLabel[1][m] = ped.GetMarkerInfo(mtmarkers[m])->alleleLabels[1];
         }else{
            mtalleleLabel[0][m] = "0";
            mtalleleLabel[1][m] = "0";
         }
//   }

   if(nonsnp0Count/* && (analysisFlag || !notrimFlag)*/)
      printf("  %d markers with 1 allele are removed.\n", nonsnp0Count);
   if(nonsnp3Count)
      printf("  %d markers with >2 alleles are removed.\n", nonsnp3Count);
   if(nonAutosomeCount)
      printf("  %d SNPs on unknown chromosomes are removed.\n",
         nonAutosomeCount);
   if(chrList.Length() && notspecifiedCount)
      printf("  %d SNPs not specified through --chr --start --stop are removed.\n",
         notspecifiedCount);
   if(markerCount+xmarkerCount+ymarkerCount+mtmarkerCount==0)
      error("There are no SNPs (with 2 alleles) in the data.");
   else{
      if(markerCount) printf("  %d autosome", markerCount);
      if(xmarkerCount && markerCount) printf(" and");
      if(xmarkerCount) printf(" %d X-chromsome", xmarkerCount);
      if(ymarkerCount) printf(" %d Y-chromsome", ymarkerCount);
      if(mtmarkerCount) printf(" %d Mitochondrial", mtmarkerCount);
      printf(" SNPs are loaded\n");
   }
   if(markerCount){
      shortCount = (markerCount-1)/16+1;
//      byteCount = (markerCount-1)/8+1;
   }else
      shortCount = 0;
   if(xmarkerCount){
      xshortCount = (xmarkerCount-1)/16+1;
//      xbyteCount = (xmarkerCount-1)/8+1;
      for(int j = 0; j < 2; j++)
         XG[j] = new unsigned short int * [idCount];
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = 0; i < id[f].Length(); i++)
            for(int j = 0; j < 2; j++){
               XG[j][geno[id[f][i]]] = new unsigned short int [xshortCount];
               for(int m = 0; m < xshortCount; m++)
                  XG[j][geno[id[f][i]]][m] = 0;
            }
   }else
      xshortCount = 0;

   if(ymarkerCount){
      yshortCount = (ymarkerCount-1)/16+1;
      for(int j = 0; j < 2; j++)
         YG[j] = new unsigned short int * [idCount];
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = 0; i < id[f].Length(); i++)
            for(int j = 0; j < 2; j++){
               YG[j][geno[id[f][i]]] = new unsigned short int [yshortCount];
               for(int m = 0; m < yshortCount; m++)
                  YG[j][geno[id[f][i]]][m] = 0;
            }
   }else
      yshortCount = 0;

   if(mtmarkerCount){
      mtshortCount = (mtmarkerCount-1)/16+1;
      for(int j = 0; j < 2; j++)
         MG[j] = new unsigned short int * [idCount];
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = 0; i < id[f].Length(); i++)
            for(int j = 0; j < 2; j++){
               MG[j][geno[id[f][i]]] = new unsigned short int [mtshortCount];
               for(int m = 0; m < mtshortCount; m++)
                  MG[j][geno[id[f][i]]][m] = 0;
            }
   }else
      mtshortCount = 0;

   for(int j = 0; j < 2; j++)
      GG[j] = new unsigned short int * [idCount];
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = 0; i < id[f].Length(); i++)
         for(int j = 0; j < 2; j++){
            GG[j][geno[id[f][i]]] = new unsigned short int [shortCount];
            for(int m = 0; m < shortCount; m++)
               GG[j][geno[id[f][i]]][m] = 0;
         }
   int minorAllele, pos, offset;
   for(int m = 0; m < markerCount; m++){
      minorAllele = ((individualInfo||mafFilterFlag||diagFlag) && (frequencies[m] > 0.5))? 2:1;
      pos = m / 16;
      offset = m % 16;
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = 0; i < id[f].Length(); i++)
            if(ped[id[f][i]].markers[markers[m]].isKnown()){
               if(ped[id[f][i]].markers[markers[m]].isHeterozygous())  // Aa: 01
                  GG[1][geno[id[f][i]]][pos] |= shortbase[offset];
               else{
                  GG[0][geno[id[f][i]]][pos] |= shortbase[offset];  // aa: 10
                  if(ped[id[f][i]].markers[markers[m]].countAlleles(minorAllele) == 2) //AA: 11
                     GG[1][geno[id[f][i]]][pos] |= shortbase[offset];
               }
            }
   }
   for(int m = 0; m < xmarkerCount; m++){
      minorAllele = ((individualInfo||mafFilterFlag||diagFlag) && (xfrequencies[m] > 0.5))? 2:1;
      pos = m / 16;
      offset = m % 16;
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = 0; i < id[f].Length(); i++)
            if(ped[id[f][i]].markers[xmarkers[m]].isKnown()){
               if(ped[id[f][i]].markers[xmarkers[m]].isHeterozygous())  // Aa: 01
                  XG[1][geno[id[f][i]]][pos] |= shortbase[offset];
               else{
                  XG[0][geno[id[f][i]]][pos] |= shortbase[offset];  // aa: 10
                  if(ped[id[f][i]].markers[xmarkers[m]].countAlleles(minorAllele) == 2) //AA: 11
                     XG[1][geno[id[f][i]]][pos] |= shortbase[offset];
               }
            }
   }

   for(int m = 0; m < ymarkerCount; m++){
      minorAllele = 1;
      pos = m / 16;
      offset = m % 16;
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = 0; i < id[f].Length(); i++)
            if(ped[id[f][i]].markers[ymarkers[m]].isKnown()){
               if(ped[id[f][i]].markers[ymarkers[m]].isHeterozygous())  // Aa: 01
                  YG[1][geno[id[f][i]]][pos] |= shortbase[offset];
               else{
                  YG[0][geno[id[f][i]]][pos] |= shortbase[offset];  // aa: 10
                  if(ped[id[f][i]].markers[ymarkers[m]].countAlleles(minorAllele) == 2) //AA: 11
                     YG[1][geno[id[f][i]]][pos] |= shortbase[offset];
               }
            }
   }

   for(int m = 0; m < mtmarkerCount; m++){
      minorAllele = 1;
      pos = m / 16;
      offset = m % 16;
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = 0; i < id[f].Length(); i++)
            if(ped[id[f][i]].markers[mtmarkers[m]].isKnown()){
               if(ped[id[f][i]].markers[mtmarkers[m]].isHeterozygous())  // Aa: 01
                  MG[1][geno[id[f][i]]][pos] |= shortbase[offset];
               else{
                  MG[0][geno[id[f][i]]][pos] |= shortbase[offset];  // aa: 10
                  if(ped[id[f][i]].markers[mtmarkers[m]].countAlleles(minorAllele) == 2) //AA: 11
                     MG[1][geno[id[f][i]]][pos] |= shortbase[offset];
               }
            }
   }
}

char *KingEngine::currentTime()
{
   time_t timer;
   struct tm *tblock;
   timer = time(NULL);
   tblock = localtime(&timer);
   return(asctime(tblock));
}
void KingEngine::BuildBinary()
{
   IntArray markers;
   id = new IntArray [ped.familyCount];
   geno.Dimension(ped.count);
   geno.Set(-1);
   int k = 0;
   for(int f = 0; f < ped.familyCount; f++){
      id[f].Dimension(0);
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(ped[i].ngeno >= MINSNPCOUNT ){
            id[f].Push(i);
            geno[i] = k;
            k++;
         }
   }

   Vector frequencies(markerCount);
   markerCount = 0;
   markers.Dimension(0);
   int nonsnpCount = 0; int nonAutosomeCount = 0;
   if( mafFilterFlag || (individualInfo && !freqLoaded)){
      double freq; int freqCount; int mafCount = 0;
      for(int m = 0; m < ped.markerCount; m++){
         if(ped.GetMarkerInfo(m)->alleleLabels.Length() != 3) {
            nonsnpCount++;
            continue;
         }
         if(ped.GetMarkerInfo(m)->chromosome > SEXCHR-1 ||
            ped.GetMarkerInfo(m)->chromosome == 0) {
            nonAutosomeCount++;
            continue;// only autosomes
         }
         freq = 0;
         freqCount = 0;
         for(int f = 0; f < ped.familyCount; f++)
            for(int i = 0; i < id[f].Length(); i++)
               if(ped[id[f][i]].markers[m].isKnown()){
                  freq += ped[id[f][i]].markers[m].countAlleles(1);
                  freqCount += 2;
               }
         if(freqCount) freq /= freqCount;
         if(mafFilterFlag){
            if(freq < minMAF || freq > 1 - minMAF) {mafCount++; continue;}
         }
         markers.Push(m);
         frequencies[markerCount] = freq;
         markerCount++;
      }
      if(mafCount)
         printf("  %d SNPs with MAF < %lf are removed.\n", mafCount, minMAF);
   }else
      for(int m = 0; m < ped.markerCount; m++){
         if(ped.GetMarkerInfo(m)->alleleLabels.Length() != 3) {
            nonsnpCount++;
            continue;
         }
         if(ped.GetMarkerInfo(m)->chromosome > SEXCHR-1 ||
            ped.GetMarkerInfo(m)->chromosome == 0) {
            nonAutosomeCount++;
            continue;// only autosomes
         }
         markers.Push(m);
         markerCount++;
      }
   if(nonsnpCount)
      printf("  %d markers with 1 or > 2 alleles are removed.\n", nonsnpCount);

   if(nonAutosomeCount)
      printf("  %d SNPs that are not on autosomes are removed.\n", nonAutosomeCount);

   idCount = 0;
   for(int f = 0; f < ped.familyCount; f++)
      idCount += id[f].Length();
}
void KingEngine::WriteKingBinary(const char *pedfile)
{
   if(geno.Length()==0)  // binary genotype does not exist
      BuildShortBinary();
   DumpBinary(pedfile);
}

void KingEngine::DumpBinary(const char *pedfile)
{
   bool markerinfo = (bp.Length()||xbp.Length())? true: false;  // update later

   printf("Writing to file %s starts at %s",
      (const char*)pedfile, currentTime());
   IntArray valid(ped.count);
   valid.Set(1);
   String temp;
   int invalidCount = 0;
   int sampleinvalidCount = 0;
   if(exclusionList.Length()){
      for(int i = 0; i < ped.count; i++){
         temp = ped[i].famid;
         temp += "->"; // '->'
         temp += ped[i].pid;
         if(exclusionList.Find(temp)>-1 && valid[i]) {
            valid[i] = 0;
            invalidCount++;
            if(geno[i] > -1) sampleinvalidCount++;
         }
      }
      printf("%d samples are excluded.\n", invalidCount);
   }
   FILE *fp = fopen(pedfile, "wb");
   if(fp == NULL) error("Cannot open %s to write.", pedfile);
   // line 1: headerCount, idCount, markerCount
   int headerCount = 7 + markerinfo*2;
   if(snpName.Length() >= markerCount) headerCount += 2;
   headerCount += ped.affectionCount + ped.traitCount + ped.covariateCount + (ped.haveTwins>0);
   printf("  # Header lines: %d\n", headerCount);
   fprintf(fp, "%d %d %d", headerCount, idCount-sampleinvalidCount, markerCount);
   if(xmarkerCount || ymarkerCount || mtmarkerCount) {
      fprintf(fp, " %d", xmarkerCount);
      if(ymarkerCount || mtmarkerCount)
         fprintf(fp, " %d %d", ymarkerCount, mtmarkerCount);
   }
   fprintf(fp, "\n");
   int maxcount = snpName.Length()-1;
   if(snpName.Length() >= markerCount){
      // line 2: snpName
      //bigdataIdx[oldmarkerCount-1-m]
      fprintf(fp, "M KING_SNP");
      for(int m = 0; m < markerCount; m++)
//         fprintf(fp, " %s", (const char*)snpName[bigdataIdx[maxcount-m]]);
         fprintf(fp, " %s", (const char*)snpName[m]);
      for(int m = 0; m < xmarkerCount; m++)
         fprintf(fp, " %s", (const char*)xsnpName[m]);
      for(int m = 0; m < ymarkerCount; m++)
         fprintf(fp, " %s", (const char*)ysnpName[m]);
      for(int m = 0; m < mtmarkerCount; m++)
         fprintf(fp, " %s", (const char*)mtsnpName[m]);
      fprintf(fp, "\n");
      // line 3: alleleLabel
      fprintf(fp, "M KING_ALLELE ");
      for(int m = 0; m < markerCount; m++)
         fprintf(fp, " %s %s",
//            (const char*)alleleLabel[0][bigdataIdx[maxcount-m]], (const char*)alleleLabel[1][bigdataIdx[maxcount-m]]);
            (const char*)alleleLabel[0][m], (const char*)alleleLabel[1][m]);
      for(int m = 0; m < xmarkerCount; m++)
         fprintf(fp, " %s %s",
            (const char*)xalleleLabel[0][m], (const char*)xalleleLabel[1][m]);
      for(int m = 0; m < ymarkerCount; m++)
         fprintf(fp, " %s %s",
            (const char*)yalleleLabel[0][m], (const char*)yalleleLabel[1][m]);
      for(int m = 0; m < mtmarkerCount; m++)
         fprintf(fp, " %s %s",
            (const char*)mtalleleLabel[0][m], (const char*)mtalleleLabel[1][m]);
      fprintf(fp, "\n");
   }
   // line 4: sampleName
   fprintf(fp, "S KING_SAMPLE");
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = 0; i < id[f].Length(); i++)
         if(valid[id[f][i]])
            fprintf(fp, " %s->%s",
            (const char*)ped[id[f][i]].famid, (const char*)ped[id[f][i]].pid);
   fprintf(fp, "\n");
   if(markerinfo){
      // line: KING_CHROMOSOME
      fprintf(fp, "M KING_CHROMOSOME");
      for(int m = 0; m < markerCount; m++)
         fprintf(fp, " %d", chromosomes[m]);
      fprintf(fp, "\n");
      // line: KING_POSITION
      fprintf(fp, "M KING_POSITION");
      for(int m = 0; m < markerCount; m++)
         fprintf(fp, " %.6lf", bp[m]*0.000001);
        for(int m = 0; m < xmarkerCount; m++)
         fprintf(fp, " %.6lf", xbp[m]*0.000001);
      for(int m = 0; m < ymarkerCount; m++)
         fprintf(fp, " %.6lf", ybp[m]*0.000001);
      for(int m = 0; m < mtmarkerCount; m++)
         fprintf(fp, " %.6lf", mtbp[m]*0.000001);
      fprintf(fp, "\n");
   }
   // line: FID
   fprintf(fp, "S KING_FID");
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = 0; i < id[f].Length(); i++)
         if(valid[id[f][i]])
            fprintf(fp, " %s", (const char*)ped[id[f][i]].famid);
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(geno[i]==-1 && valid[i])
            fprintf(fp, " %s", (const char*)ped[i].famid);
   for(int i = 0; i < inclusionList[0].Length(); i++)
      fprintf(fp, " %s", (const char*)inclusionList[0][i]);
   fprintf(fp, "\n");
   // line: PID
   fprintf(fp, "S KING_PID");
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = 0; i < id[f].Length(); i++)
         if(valid[id[f][i]])
            fprintf(fp, " %s", (const char*)ped[id[f][i]].pid);
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(geno[i]==-1 && valid[i])
            fprintf(fp, " %s", (const char*)ped[i].pid);
   for(int i = 0; i < inclusionList[0].Length(); i++)
      fprintf(fp, " %s", (const char*)inclusionList[1][i]);
   fprintf(fp, "\n");
   // line: FATHER
   fprintf(fp, "S KING_FATHER");
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = 0; i < id[f].Length(); i++)
         if(valid[id[f][i]])
            fprintf(fp, " %s", (const char*)ped[id[f][i]].fatid);
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(geno[i]==-1 && valid[i])
            fprintf(fp, " %s", (const char*)ped[i].fatid);
   for(int i = 0; i < inclusionList[0].Length(); i++)
      fprintf(fp, " 0");
   fprintf(fp, "\n");
   // line: MOTHER
   fprintf(fp, "S KING_MOTHER");
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = 0; i < id[f].Length(); i++)
         if(valid[id[f][i]])
            fprintf(fp, " %s", (const char*)ped[id[f][i]].motid);
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(geno[i]==-1 && valid[i])
            fprintf(fp, " %s", (const char*)ped[i].motid);
   for(int i = 0; i < inclusionList[0].Length(); i++)
      fprintf(fp, " 0");
   fprintf(fp, "\n");
   // line: SEX
   fprintf(fp, "S KING_SEX");
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = 0; i < id[f].Length(); i++)
         if(valid[id[f][i]])
            fprintf(fp, " %d", ped[id[f][i]].sex);
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(geno[i]==-1 && valid[i])
            fprintf(fp, " %d", ped[i].sex);
   for(int i = 0; i < inclusionList[0].Length(); i++)
      fprintf(fp, " %s", (const char*)inclusionList[2][i]);
   fprintf(fp, "\n");
   // lines: affections
   for(int t = 0; t < ped.affectionCount; t++){
      fprintf(fp, "A %s", (const char*)ped.affectionNames[t]);
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = 0; i < id[f].Length(); i++)
            if(valid[id[f][i]]){
               if(ped[id[f][i]].isDiagnosed(t))
                  fprintf(fp, " %d", ped[id[f][i]].affections[t]);
               else
                  fprintf(fp, " 0");
            }
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
            if(geno[i]==-1 && valid[i])
               fprintf(fp, " %d", ped[i].affections[t]);
      for(int i = 0; i < inclusionList[0].Length(); i++)
         fprintf(fp, " 0");
      fprintf(fp, "\n");
   }
   // lines: traits
   for(int t = 0; t < ped.traitCount; t++){
      fprintf(fp, "T %s", (const char*)ped.traitNames[t]);
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = 0; i < id[f].Length(); i++)
            if(valid[id[f][i]]){
               if(ped[id[f][i]].isPhenotyped(t))
                  fprintf(fp, " %lf", ped[id[f][i]].traits[t]);
               else
                  fprintf(fp, " X");
            }
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
            if(valid[i]){
               if(geno[i]==-1){
                  if(ped[i].isPhenotyped(t))
                     fprintf(fp, " %lf", ped[i].traits[t]);
                  else
                     fprintf(fp, " X");
               }
            }
      for(int i = 0; i < inclusionList[0].Length(); i++)
         fprintf(fp, " X");
      fprintf(fp, "\n");
   }
   // lines: covariates
   for(int t = 0; t < ped.covariateCount; t++){
      fprintf(fp, "C %s", (const char*)ped.covariateNames[t]);
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = 0; i < id[f].Length(); i++)
            if(valid[id[f][i]]){
               if(ped[id[f][i]].isControlled(t))
                  fprintf(fp, " %lf", ped[id[f][i]].covariates[t]);
               else
                  fprintf(fp, " X");
            }
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
            if(valid[i]){
               if(geno[i]==-1){
                  if(ped[i].isControlled(t))
                     fprintf(fp, " %lf", ped[i].covariates[t]);
                  else
                     fprintf(fp, " X");
               }
            }
      for(int i = 0; i < inclusionList[0].Length(); i++)
         fprintf(fp, " X");
      fprintf(fp, "\n");
   }
   // line: twin
   if(ped.haveTwins){
      StringArray twinCodes(0);
      twinCodes.Push("0"); twinCodes.Push("MZ"); twinCodes.Push("DZ");
      fprintf(fp, "Z twin");
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = 0; i < id[f].Length(); i++)
            if(valid[id[f][i]]){
               if(ped[id[f][i]].zygosity <= 2)
                  fprintf(fp, " %s", (const char*)twinCodes[ped[id[f][i]].zygosity]);
               else
                  fprintf(fp, " %d", ped[id[f][i]].zygosity);
            }
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
            if(geno[i]==-1 && valid[i]){
               if(ped[i].zygosity <= 2)
                  fprintf(fp, " %s", (const char*)twinCodes[ped[i].zygosity]);
               else
                  fprintf(fp, " %d", ped[i].zygosity);
            }
      for(int i = 0; i < inclusionList[0].Length(); i++)
         fprintf(fp, " 0");
      fprintf(fp, "\n");
   }

   // line: genotype data
   fprintf(fp, "KING01");
   if(GG[0] || XG[0] || YG[0] || MG[0])
      for(int f = 0; f < ped.familyCount; f++)
         for(int i = 0; i < id[f].Length(); i++)
            if(valid[id[f][i]]){
               fwrite(GG[0][geno[id[f][i]]], sizeof(unsigned short int), shortCount, fp);
               fwrite(GG[1][geno[id[f][i]]], sizeof(unsigned short int), shortCount, fp);
               if(xshortCount){
                  fwrite(XG[0][geno[id[f][i]]], sizeof(unsigned short int), xshortCount, fp);
                  fwrite(XG[1][geno[id[f][i]]], sizeof(unsigned short int), xshortCount, fp);
               }
               if(yshortCount){
                  fwrite(YG[0][geno[id[f][i]]], sizeof(unsigned short int), yshortCount, fp);
                  fwrite(YG[1][geno[id[f][i]]], sizeof(unsigned short int), yshortCount, fp);
               }
               if(mtshortCount){
                  fwrite(MG[0][geno[id[f][i]]], sizeof(unsigned short int), mtshortCount, fp);
                  fwrite(MG[1][geno[id[f][i]]], sizeof(unsigned short int), mtshortCount, fp);
               }
            }
   fclose(fp);
   printf("Genotype data saved as KING binary format in file %s\n",
      (const char*)pedfile);
}


