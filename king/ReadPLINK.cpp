//////////////////////////////////////////////////////////////////////
// ReadPLINK.cpp
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

#include <string.h>
#include "analysis.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

void KingEngine::ReadMultiplePlinkBinaryBigData(const char *filename, const char *famfile, const char *bimfile, const char *covfile, const char *phefile)
{
   StringArray filenames;
   filenames.AddTokens(filename, ',');
   for(int i = 0; i < filenames.Length(); i++){
      int m = filenames[i].Find(".bed");
      if(m > -1) filenames[i]=filenames[i].SubStr(0,m);
   }
   int cohortCount = filenames.Length();
   IntArray idCounts(cohortCount);
   idCounts.Zero();
   String line, tempS;
   StringArray tokens;
   FID.Dimension(0);
   PID.Dimension(0);
   FA.Dimension(0);
   MO.Dimension(0);
   SEX.Dimension(0);
   sampleName.Dimension(0);
   affectionNames.Dimension(1);
   affectionNames[0] = "DISEASEKING";
   IntArray disease(0);
   String pedfile;
   IFILE input;

   StringArray myfilenames=filenames;
   int mycohortCount = cohortCount;
   if(famfile!=NULL && famfile[0]!='\0'){
      myfilenames.Dimension(0);
      myfilenames.AddTokens(famfile, ',');
      for(int i = 0; i < myfilenames.Length(); i++){
         int m = myfilenames[i].Find(".fam");
         if(m > -1) myfilenames[i]=myfilenames[i].SubStr(0,m);
      }
      mycohortCount = myfilenames.Length();
   }
   printf("Read in PLINK fam files\n");
   idCount = 0;
   for(int c = 0; c < mycohortCount; c++){
      pedfile = myfilenames[c];
      pedfile.Add(".fam");
      input = ifopen(pedfile, "rt");
      if(input == NULL){
         String pedfile2 = pedfile;
         pedfile2.Add(".gz");
         input = ifopen(pedfile2, "rt");
         if(input == NULL)
            error("Pedigree file %s cannot be opened", (const char*)pedfile);
         printf("\t%s...\n", (const char*)pedfile2);
      }else
         printf("\t%s...\n", (const char*)pedfile);
      while(!ifeof(input)){
         tokens.Clear();
         line.ReadLine(input);
         tokens.AddTokens(line);
         if(tokens.Length() < 6) continue;
         FID.Push(tokens[0]);
         PID.Push(tokens[1]);
         tempS = tokens[0];
         tempS.Add("->"); // '->'
         tempS.Add(tokens[1]);
         sampleName.Push(tempS);
         FA.Push(tokens[2]);
         MO.Push(tokens[3]);
         if(tokens[4].AsInteger() == -9)
            SEX.Push("0");
         else
            SEX.Push(tokens[4]);
         if(tokens[5].AsInteger() == -9)
            disease.Push(0);
         else
            disease.Push(tokens[5].AsInteger());
         idCounts[c] ++;
      }
      idCount += idCounts[c];
      ifclose(input);
   }
   StringIntHash sampleHash;
   sampleHash.Clear();
   for(int i = 0; i < idCounts[0]; i++)
      sampleHash.SetInteger(sampleName[i], i);
   int overlapCount = 0;
   int pos = idCounts[0];
   for(int c = 1; c < mycohortCount; c++){
      for(int i = 0; i < idCounts[c]; i++)
         if(sampleHash.Integer(sampleName[pos+i])>-1) overlapCount++;
      pos += idCounts[c];
   }
   if(overlapCount > (idCount-idCounts[0])*0.05){
      printf("  %d samples have identical IDs between %s.fam and %s\n",
         overlapCount, (const char*)myfilenames[0],
         (mycohortCount==2) ? (const char*)myfilenames[1]: "other families");
      String REF="REF_";
      printf("  REF_ is added to all IDs in %s.fam.\n", (const char*)myfilenames[0]);
      for(int i = 0; i < idCounts[0]; i++) PID[i] = REF+PID[i];
      String QUERY="QRY_";
      printf("  QRY_ is added to all IDs in %s.fam",(const char*)myfilenames[1]);
      pos = idCounts[0]+idCounts[1];
      for(int i = idCounts[0]; i < pos; i++)
         PID[i] = QUERY+PID[i];
      for(int c = 2; c < mycohortCount; c++){
         printf(" %s.fam", (const char*)myfilenames[c]);
         for(int i = pos; i < pos+idCounts[c]; i++)
            PID[i] = QUERY+PID[i];
         pos += idCounts[c];
      }
      printf("\n");
      for(int i = 0; i < idCount; i++){
         tempS = FID[i];
         tempS.Add("->");
         tempS.Add(PID[i]);
         sampleName[i] = tempS;
      }
   }else if(overlapCount > 0)
      error("KING cannot handle %d samples with identical IDs between %s.fam and %s\n",
         overlapCount, (const char*)myfilenames[0],
         (mycohortCount==2) ? (const char*)myfilenames[1]: "other families");

   // check if pedigrees are complete
   IntArray inclusionFA(0);
   IntArray inclusionMO(0);
   StringArray newName = sampleName;
   int newid = 1;
   for(int i = 0; i < idCount; i++){
      if(FA[i]!="0"){
         tempS = FID[i];
         tempS.Add("->");    // '->'
         tempS.Add(FA[i]);
         if(newName.Find(tempS)==-1){
            inclusionFA.Push(i);
            newName.Push(tempS);
         }
      }
      if(MO[i]!="0"){
         tempS = FID[i];
         tempS.Add("->"); // '->'
         tempS.Add(MO[i]);
         if(newName.Find(tempS)==-1){
            inclusionMO.Push(i);
            newName.Push(tempS);
         }
      }
      if(FA[i]!="0" && MO[i]=="0"){
         MO[i] = "KING";
         MO[i] += newid;
         tempS = FID[i];
         tempS.Add("->"); // '->'
         tempS.Add(MO[i]);
         newName.Push(tempS);
         inclusionMO.Push(i);
         newid++;
      }else if(FA[i]=="0" && MO[i]!="0"){
         FA[i] = "KING";
         FA[i] += newid;
         tempS = FID[i];
         tempS.Add("->"); // '->'
         tempS.Add(FA[i]);
         newName.Push(tempS);
         inclusionFA.Push(i);
         newid++;
      }
   }
   if(inclusionFA.Length() || inclusionMO.Length()) {
      for(int i = 0; i < inclusionFA.Length(); i++){
         FID.Push(FID[inclusionFA[i]]);
         PID.Push(FA[inclusionFA[i]]);
         FA.Push("0");
         MO.Push("0");
         SEX.Push("1");
         disease.Push(0);
      }
      for(int i = 0; i < inclusionMO.Length(); i++){
         FID.Push(FID[inclusionMO[i]]);
         PID.Push(MO[inclusionMO[i]]);
         FA.Push("0");
         MO.Push("0");
         SEX.Push("2");
         disease.Push(0);
      }
   }
   affections.Dimension(1, FID.Length());
   for(int i = 0; i < FID.Length(); i++)
      affections[0][i] = disease[i];
   markerCount = 10000;
   MakePed();
   markerCount = 0;
   printf("  PLINK pedigrees loaded: %d samples\n", idCount);
   bool diseaseMissing = true;
   for(int i = 0; i < ped.count; i++)
      if(ped[i].affections[0]){
         diseaseMissing = false;
         break;
      }
   if(diseaseMissing)
      ped.affectionCount = 0;

   myfilenames=filenames;
   mycohortCount = cohortCount;
   if(bimfile!=NULL && bimfile[0]!='\0'){
      myfilenames.Dimension(0);
      myfilenames.AddTokens(bimfile, ',');
      for(int i = 0; i < myfilenames.Length(); i++){
         int m = myfilenames[i].Find(".bim");
         if(m > -1) myfilenames[i]=myfilenames[i].SubStr(0,m);
      }
      mycohortCount = myfilenames.Length();
   }
   String mapfile;
   printf("Read in PLINK bim files\n");
   int chr;
   int xymarkerCount = 0;
   StringArray tempLabels[2];
   Vector tempPositions(0);   IntArray tempBP(0);
   StringArray tempNames(0);
   IntArray tempChromosomes(0);
   tempLabels[0].Dimension(0); tempLabels[1].Dimension(0);
   tempPositions.Dimension(0); tempBP.Dimension(0);
   tempNames.Dimension(0);
   tempChromosomes.Dimension(0);
   StringIntHash markerHash;
   markerHash.Clear();
   IntArray *valids = new IntArray[cohortCount];
   IntArray *flips = new IntArray[cohortCount];
   for(int c = 0; c < cohortCount; c++){
      valids[c].Dimension(0);
      flips[c].Dimension(0);
   }
   char AT[256];
   for(int i = 0; i < 256; i++)
      AT[i] = '0';
   AT['A'] = AT['T'] = AT['a'] = AT['t'] = 'A';
   AT['C'] = AT['G'] = AT['c'] = AT['g'] = 'C';
   mapfile = myfilenames[0];
   mapfile.Add(".bim");
   input = ifopen(mapfile, "rt");
   if(input == NULL){
      String mapfile2=mapfile;
      mapfile2.Add(".gz");
      input = ifopen(mapfile2, "rt");
      if(input == NULL)
         error("Map file %s cannot be opened", (const char*)mapfile);
      printf("\t%s...\n", (const char*)mapfile2);
   }else
      printf("\t%s...\n", (const char*)mapfile);
   while(!ifeof(input)){
      tokens.Clear();
      line.ReadLine(input);
      tokens.AddTokens(line);
      if(tokens.Length() < 6) continue;
      chr = tokens[0].AsInteger();
      tempChromosomes.Push(chr);
      tempNames.Push(tokens[1]);
      tempPositions.Push(tokens[3].AsDouble()*0.000001); tempBP.Push(tokens[3].AsInteger());
      tempLabels[0].Push(tokens[4]);
      tempLabels[1].Push(tokens[5]);
      flips[0].Push(0);
      if(((chr < 1) && (tokens[0]!="0")) || (chr > SEXCHR+3) // chromosome not 0-26
         || (chrList.Length() && (chrList.Find(chr)==-1 // chr not specified
         || (start != _NAN_ && tokens[3].AsDouble()*0.000001 < start) // pos not specified
         || (stop != _NAN_ && tokens[3].AsDouble()*0.000001 > stop) ) )
         || tokens[4][1]!='\0' || tokens[5][1]!='\0'  // not single letter
         || (AT[tokens[4][0]]!='A' && AT[tokens[4][0]]!='C')   // not ATCG
         || (AT[tokens[5][0]]!='A' && AT[tokens[5][0]]!='C')   // not ATCG
         || (AT[tokens[4][0]]==AT[tokens[5][0]])   // AT or CG
         ){ // skip
         valids[0].Push(0);
      }else{
         valids[0].Push(1);
         markerHash.SetInteger(tokens[1], tempNames.Length()-1);
      }
   }
   ifclose(input);
   for(int c = 1; c < mycohortCount; c++){
      mapfile = myfilenames[c];
      mapfile.Add(".bim");
      input = ifopen(mapfile, "rt");
      if(input == NULL){
         String mapfile2 = mapfile;
         mapfile2.Add(".gz");
         input = ifopen(mapfile2, "rt");
         if(input == NULL)
            error("Map file %s cannot be opened", (const char*)mapfile);
         printf("\t%s...\n", (const char*)mapfile2);
      }else
         printf("\t%s...\n", (const char*)mapfile);
      while(!ifeof(input)){
         tokens.Clear();
         line.ReadLine(input);
         tokens.AddTokens(line);
         if(tokens.Length() < 6) continue;
         int k = markerHash.Integer(tokens[1]);
         if(k==-1){  // SNP not available in cohort 1
            valids[c].Push(-1);
            flips[c].Push(0);
         }else{
            if(tokens[4][1]=='\0' && tokens[5][1]=='\0'
               && AT[int(tokens[4][0])] == AT[int(tempLabels[1][k][0])]
               && AT[int(tokens[5][0])] == AT[int(tempLabels[0][k][0])]){
               flips[c].Push(1);
               valids[c].Push(k);
               valids[0][k] |= (1<<c);
            }else if(tokens[4][1]=='\0' && tokens[5][1]=='\0'
               && AT[int(tokens[4][0])] == AT[int(tempLabels[0][k][0])]
               && AT[int(tokens[5][0])] == AT[int(tempLabels[1][k][0])]){
               flips[c].Push(0);
               valids[c].Push(k);
               valids[0][k] |= (1<<c);
            }else{   // allele labels not matching
               flips[c].Push(0);
               valids[c].Push(-1);
            }
         }
      }  // end of input
      ifclose(input);
   }   // end of cohort
   chromosomes.Dimension(0);
   bp.Dimension(0);  //positions.Dimension(0);
   snpName.Dimension(0);
   StringArray labels[2], xlabels[2], ylabels[2], mtlabels[2];
   labels[0].Dimension(0); labels[1].Dimension(0);
   xbp.Dimension(0); //xpositions.Dimension(0);
   xsnpName.Dimension(0);
   xlabels[0].Dimension(0); xlabels[1].Dimension(0);
   xbp.Dimension(0); //ypositions.Dimension(0);
   ysnpName.Dimension(0);
   ylabels[0].Dimension(0); ylabels[1].Dimension(0);
   mtbp.Dimension(0);   //mtpositions.Dimension(0);
   mtsnpName.Dimension(0);
   mtlabels[0].Dimension(0); mtlabels[1].Dimension(0);
   IntArray *validMarkers = new IntArray[cohortCount];
   IntArray *xvalidMarkers = new IntArray[cohortCount];
   IntArray *yvalidMarkers = new IntArray[cohortCount];
   IntArray *mtvalidMarkers = new IntArray[cohortCount];
   for(int c = 0; c < cohortCount; c++){
      int count = valids[c].Length();
      validMarkers[c].Dimension(count);
      xvalidMarkers[c].Dimension(count);
      yvalidMarkers[c].Dimension(count);
      mtvalidMarkers[c].Dimension(count);
      validMarkers[c].Set(-1);
      xvalidMarkers[c].Set(-1);
      yvalidMarkers[c].Set(-1);
      mtvalidMarkers[c].Set(-1);
   }
   int ALLVALID = (1<<cohortCount)-1;
   markerCount = xmarkerCount = ymarkerCount = mtmarkerCount = 0;
   for(int m = 0; m < valids[0].Length(); m++){
      if(valids[0][m]==ALLVALID){
         if(tempChromosomes[m]==SEXCHR){
            xbp.Push(tempBP[m]); //xpositions.Push(tempPositions[m]); 
            xsnpName.Push(tempNames[m]);
            xlabels[0].Push(tempLabels[0][m]);
            xlabels[1].Push(tempLabels[1][m]);
            xvalidMarkers[0][m] = xmarkerCount++;
         }else if(tempChromosomes[m]==SEXCHR+1){
            ybp.Push(tempBP[m]); //ypositions.Push(tempPositions[m]); 
            ysnpName.Push(tempNames[m]);
            ylabels[0].Push(tempLabels[0][m]);
            ylabels[1].Push(tempLabels[1][m]);
            yvalidMarkers[0][m] = ymarkerCount++;
         }else if(tempChromosomes[m]==SEXCHR+3){
            mtbp.Push(tempBP[m]);   //mtpositions.Push(tempPositions[m]); 
            mtsnpName.Push(tempNames[m]);
            mtlabels[0].Push(tempLabels[0][m]);
            mtlabels[1].Push(tempLabels[1][m]);
            mtvalidMarkers[0][m] = mtmarkerCount++;
         }else{
            bp.Push(tempBP[m]);  //positions.Push(tempPositions[m]); 
            chromosomes.Push(tempChromosomes[m]);
            snpName.Push(tempNames[m]);
            labels[0].Push(tempLabels[0][m]);
            labels[1].Push(tempLabels[1][m]);
            validMarkers[0][m] = markerCount++;
         }
      }
   }
   for(int c = 1; c < cohortCount; c++)
      for(int m = 0; m < valids[c].Length(); m++){
         int m1 = valids[c][m];
         if(m1>-1){
            if(tempChromosomes[m1]==SEXCHR)
               xvalidMarkers[c][m] = xvalidMarkers[0][m1];
            else if(tempChromosomes[m1]==SEXCHR+1)
               yvalidMarkers[c][m] = yvalidMarkers[0][m1];
            else if(tempChromosomes[m1]==SEXCHR+3)
               mtvalidMarkers[c][m] = mtvalidMarkers[0][m1];
            else
               validMarkers[c][m] = validMarkers[0][m1];
         }
      }
   if(markerCount)
      printf("  Genotype data consist of %d autosome SNPs", markerCount);
   else error("No autosome genotypes available for KING inferences.");
   if(xmarkerCount) printf(", %d X-chromosome SNPs", xmarkerCount);
   if(ymarkerCount) printf(", %d Y-chromosome SNPs", ymarkerCount);
   if(mtmarkerCount) printf(", %d mitochondrial SNPs", mtmarkerCount);
   printf("\n");
   printf("  PLINK maps loaded: %d SNPs\n", markerCount + xmarkerCount + ymarkerCount + mtmarkerCount);

   for(int i = 0; i < ped.count; i++)
      if(ped[i].ngeno == 10000)
         ped[i].ngeno = markerCount+xmarkerCount+ymarkerCount+mtmarkerCount;
   if(markerCount){
      if(Bit64==64)
         longCount = (markerCount-1)/64+1;
      shortCount = (markerCount-1)/16+1;
   }else
      longCount = shortCount = 0;
   for(int i = 0; i < 2; i++)
      alleleLabel[i] = labels[i];
   if(Bit64 == 64)
      for(int j = 0; j < 2; j++){
         LG[j] = new unsigned long long int * [idCount];
         for(int i = 0; i < idCount; i++){
            LG[j][i] = new unsigned long long int [longCount];
         }
      }
   else
      for(int j = 0; j < 2; j++){
         GG[j] = new unsigned short int * [idCount];
         for(int i = 0; i < idCount; i++){
            GG[j][i] = new unsigned short int [shortCount];
         }
      }
   if(!autosomeOnly){
      if(xmarkerCount){
         xshortCount = (xmarkerCount-1)/16+1;
         for(int j = 0; j < 2; j++)
            XG[j] = new unsigned short int * [idCount];
         for(int i = 0; i < idCount; i++)
            for(int j = 0; j < 2; j++){
               XG[j][i] = new unsigned short int [xshortCount];
               for(int m = 0; m < xshortCount; m++)
                  XG[j][i][m] = 0;
            }
         for(int i = 0; i < 2; i++)
            xalleleLabel[i] = xlabels[i];
      }else
         xshortCount = 0;
      if(ymarkerCount){
         yshortCount = (ymarkerCount-1)/16+1;
         for(int j = 0; j < 2; j++)
            YG[j] = new unsigned short int * [idCount];
         for(int i = 0; i < idCount; i++)
            for(int j = 0; j < 2; j++){
               YG[j][i] = new unsigned short int [yshortCount];
               for(int m = 0; m < yshortCount; m++)
                  YG[j][i][m] = 0;
            }
         for(int i = 0; i < 2; i++)
            yalleleLabel[i] = ylabels[i];
      }else
         yshortCount = 0;
      if(mtmarkerCount){
         mtshortCount = (mtmarkerCount-1)/16+1;
         for(int j = 0; j < 2; j++)
            MG[j] = new unsigned short int * [idCount];
         for(int i = 0; i < idCount; i++)
            for(int j = 0; j < 2; j++){
               MG[j][i] = new unsigned short int [mtshortCount];
               for(int m = 0; m < mtshortCount; m++)
                  MG[j][i][m] = 0;
            }
         for(int i = 0; i < 2; i++)
            mtalleleLabel[i] = mtlabels[i];
      }else
         mtshortCount = 0;
   }  // end of if(!autosomeOnly)
   unsigned char flipByte[256];
   for(int b = 0; b < 256; b++){
      flipByte[b] = 0;
      for(int i = 0; i < 4; i++){
         int g = (b>>(i*2))&3;
         if(g==3);      // 3->0
         else if(g==0)  // 0->3
            flipByte[b] |= (3<<(i*2));
         else           // 1->1, 2->2
            flipByte[b] |= (g<<(i*2));
      }
   }
   String bedfile;
   unsigned char bit1, bit2;
   unsigned short int intbit1, intbit2;
   int newbyte, newbit;
   long int oldbyte;
   int SNP_space;
   IntArray SNP_spaces(cohortCount);
   for(int c = 0; c < cohortCount; c++){
      SNP_spaces[c] = (idCounts[c]+3)/4;
      SNP_space += SNP_spaces[c];
   }
   unsigned char **G_SNP = new unsigned char *[Bit64==64? longCount*64: shortCount*16];
   FILE *fp;
   printf("Read in PLINK bed files\n");
   int baseID = 0;
   for(int c = 0; c < cohortCount; c++){  // Read in bed files by cohorts
      bedfile = filenames[c];
      bedfile.Add(".bed");
      fp = fopen(bedfile, "rb");
      if(fp == NULL)
         error("Binary genotype file %s cannot be opened", (const char*)bedfile);
      printf("\t%s...\n", (const char*)bedfile); fflush(stdout);
      unsigned char byte;
      byte = fgetc(fp);
      if(byte != 108)
         error("The PLINK format in %s cannot be recognized. The first byte is %x",
            (const char*)bedfile, byte);
      byte = fgetc(fp);
      if(byte != 27)
         error("The PLINK format in %s cannot be recognized. The second byte is %x",
            (const char*)bedfile, byte);
      byte = fgetc(fp);
      if(byte != 1)
         error("Currently only SNP-major mode can be analyzed.");
      for(int m = 0;  m < (Bit64==64? longCount*64: shortCount*16); m++)
         G_SNP[m] = new unsigned char[SNP_spaces[c]];
      unsigned char *aSNP = new unsigned char[SNP_spaces[c]];
      int wholeblocks = ((SNP_spaces[c]-1)>>17);
      int leftoversize = SNP_spaces[c] & 0x1FFFF;
      int validCount = validMarkers[c].Length();
      for(int m = 0; m < validCount; m++){
         if(validMarkers[c][m]>-1){
            int pos = validMarkers[c][m];
            if(flips[c][m]){
               for(int i = 0; i < wholeblocks; i++)
                  if(fread(&aSNP[i<<17], 1, 0x20000, fp) < 0x20000)
                     error("Not enough genotypes at %dth marker\n", pos);
               if(fread(&aSNP[wholeblocks<<17], 1, leftoversize, fp) < leftoversize)
                  error("Not enough genotypes at %dth marker\n", pos);
               for(int i = 0; i < SNP_spaces[c]; i++)
                  G_SNP[pos][i] = flipByte[aSNP[i]];
            }else{   // no flip
               for(int i = 0; i < wholeblocks; i++)
                  if(fread(&G_SNP[pos][i<<17], 1, 0x20000, fp) < 0x20000)
                     error("Not enough genotypes at %dth marker\n", pos);
               if(fread(&G_SNP[pos][wholeblocks<<17], 1, leftoversize, fp) < leftoversize)
                  error("Not enough genotypes at %dth marker\n", pos);
            }  // end of if(flip)
         }else if(autosomeOnly || ((xvalidMarkers[c][m]==-1) && (yvalidMarkers[c][m]==-1) && (mtvalidMarkers[c][m]==-1)) ){
               fseek(fp, SNP_spaces[c], SEEK_CUR);
         }else{
            if(xvalidMarkers[c][m]>-1){ // X-chromosome
               newbyte = xvalidMarkers[c][m] / 16;
               newbit = xvalidMarkers[c][m] % 16;
               for(int i = 0; i < idCounts[c]; i++){
                  int k = i%4;
                  if(k==0){
                     byte = fgetc(fp);
                     if(flips[c][m]) byte = flipByte[byte];
                  }
                  bit1 = byte & base[k*2];
                  bit2 = byte & base[k*2+1];
                  if(bit1 && bit2){  // aa
                     XG[0][baseID+i][newbyte] |= shortbase[newbit];
                  }else if(!bit1 && !bit2){ // AA
                     XG[0][baseID+i][newbyte] |= shortbase[newbit];
                     XG[1][baseID+i][newbyte] |= shortbase[newbit];
                  }else if(!bit1 && bit2) // Aa
                     XG[1][baseID+i][newbyte] |= shortbase[newbit];
               }
            }else if(yvalidMarkers[c][m]>-1){ // Y-chromosome
               newbyte = yvalidMarkers[c][m] / 16;
               newbit = yvalidMarkers[c][m] % 16;
               for(int i = 0; i < idCounts[c]; i++){
                  int k = i%4;
                  if(k==0){
                     byte = fgetc(fp);
                     if(flips[c][m]) byte = flipByte[byte];
                  }
                  bit1 = byte & base[k*2];
                  bit2 = byte & base[k*2+1];
                  if(bit1 && bit2){  // aa
                     YG[0][baseID+i][newbyte] |= shortbase[newbit];
                  }else if(!bit1 && !bit2){ // AA
                     YG[0][baseID+i][newbyte] |= shortbase[newbit];
                     YG[1][baseID+i][newbyte] |= shortbase[newbit];
                  }else if(!bit1 && bit2) // Aa
                     YG[1][baseID+i][newbyte] |= shortbase[newbit];
               }
            }else if(mtvalidMarkers[c][m]>-1){ // mitochondrial
               newbyte = mtvalidMarkers[c][m] / 16;
               newbit = mtvalidMarkers[c][m] % 16;
               for(int i = 0; i < idCounts[c]; i++){
                  int k = i%4;
                  if(k==0){
                     byte = fgetc(fp);
                     if(flips[c][m]) byte = flipByte[byte];
                  }
                  bit1 = byte & base[k*2];
                  bit2 = byte & base[k*2+1];
                  if(bit1 && bit2){  // aa
                     MG[0][baseID+i][newbyte] |= shortbase[newbit];
                  }else if((!bit1) && (!bit2)){ // AA
                     MG[0][baseID+i][newbyte] |= shortbase[newbit];
                     MG[1][baseID+i][newbyte] |= shortbase[newbit];
                  }else if((!bit1) && bit2) // Aa
                     MG[1][baseID+i][newbyte] |= shortbase[newbit];
               }
            }// end of chromosomeType
         }  // end of if(autosome)
      }  // end of marker m
      fclose(fp);
      delete []aSNP;
      if(Bit64==64){
         unsigned long long int **localLG[2];
         for(int j = 0; j < 2; j++){
            localLG[j] = new unsigned long long int * [SNP_spaces[c]*4];
            for(int i = 0; i < SNP_spaces[c]*4; i++){
               localLG[j][i] = new unsigned long long int [longCount];
            }
         }
         ConvertMajorFromSNPtoIndividual64Bit(G_SNP, localLG, idCounts[c], 0, markerCount-1);
         for(int m = 0; m < longCount*64; m++)
            delete []G_SNP[m];
         int memsize = longCount<<3;
         for(int j = 0; j < 2; j++)
            for(int i = 0; i < idCounts[c]; i++)
               memcpy(&LG[j][baseID+i][0], &localLG[j][i][0], memsize);
         for(int j = 0; j < 2; j++){
            for(int i = 0; i < SNP_spaces[c]*4; i++)
               delete []localLG[j][i];
            delete []localLG[j];
         }
      }else{   // 32 bit
         unsigned short int **localGG[2];
         for(int j = 0; j < 2; j++){
            localGG[j] = new unsigned short int * [SNP_spaces[c]*4];
            for(int i = 0; i < SNP_spaces[c]*4; i++){
               localGG[j][i] = new unsigned short int [shortCount];
            }
         }
         ConvertMajorFromSNPtoIndividual(G_SNP, localGG, idCounts[c], 0, markerCount-1);
         for(int m = 0; m < (shortCount*16); m++)
            delete []G_SNP[m];
         int memsize = shortCount<<1;
         for(int j = 0; j < 2; j++)
            for(int i = 0; i < idCounts[c]; i++)
               memcpy(&GG[j][baseID+i][0], &localGG[j][i][0], memsize);
         for(int j = 0; j < 2; j++){
            for(int i = 0; i < SNP_spaces[c]*4; i++)
               delete []localGG[j][i];
            delete []localGG[j];
         }
      }
      baseID += idCounts[c];
   }  // end of cohort loop
   delete []G_SNP;
   printf("  PLINK binary genotypes loaded: %d samples\n", idCount);
   printf("  KING format genotype data successfully converted\n");
   for(int c = 0; c < cohortCount; c++){
      delete []validMarkers[c];
      delete []xvalidMarkers[c];
      delete []yvalidMarkers[c];
      delete []mtvalidMarkers[c];
      delete []valids[c];
      delete []flips[c];
   }
}

void KingEngine::ReadPlinkBinaryBigData(const char *filename, const char *famfile, const char *bimfile, const char *covfile, const char *phefile)
{
   printf("Loading genotype data in PLINK binary format...\n");
   plinkinput = filename;
   ReadPlinkPedFile(filename, (famfile==NULL || famfile[0]=='\0') ? NULL: famfile);
   ReadPlinkCovFile(filename, (covfile==NULL || covfile[0]=='\0') ? NULL: covfile);
   if(!genotypeOnly)
      ReadPlinkTraitFile(filename, (phefile==NULL || phefile[0]=='\0') ? NULL: phefile);
   IntArray *pvalidMarker = new IntArray;
   ReadPlinkMapFile(filename, *pvalidMarker, (bimfile==NULL || bimfile[0]=='\0') ? NULL: bimfile);
   for(int i = 0; i < ped.count; i++)
      if(ped[i].ngeno == 10000)
         ped[i].ngeno = markerCount+xmarkerCount+ymarkerCount+mtmarkerCount;
   if(markerCount){
      if(Bit64==64) longCount = (markerCount-1)/64+1;
      shortCount = (markerCount-1)/16+1;
   }// else longCount = shortCount = 0;
   if(markerCount == 0)// update later for sex-chr analysis
      error("No autosome SNPs are available. Please check your map file.");
   String bedfile = filename;
   bedfile.Add(".bed");
   FILE *fp = fopen(bedfile, "rb");
   if(fp==NULL) error("Cannot open %s to read", (const char*)bedfile);
   printf("Read in PLINK bed file %s...\n", (const char*)bedfile);
   unsigned char byte = fgetc(fp);
   if(byte != 108)
      error("The PLINK format in %s cannot be recognized. The first byte is %x",
         (const char*)bedfile, byte);
   byte = fgetc(fp);
   if(byte != 27)
      error("The PLINK format in %s cannot be recognized. The second byte is %x",
         (const char*)bedfile, byte);
   byte = fgetc(fp);
   if(byte != 1)
      error("Currently only SNP-major mode can be analyzed.");
   const int SNP_space = (idCount+3)>>2; // 1 for 1-4, 2 for 5-8 etc.
   if(Bit64==64){
      if(LG[0]==NULL)
         for(int j = 0; j < 2; j++){
            LG[j] = new unsigned long long int * [SNP_space*4];
            for(int i = 0; i < SNP_space*4; i++)
               LG[j][i] = new unsigned long long int [longCount];
         }
   }else{
      if(GG[0]==NULL)
         for(int j = 0; j < 2; j++){
            GG[j] = new unsigned short int * [SNP_space*4];
            for(int i = 0; i < SNP_space*4; i++)
               GG[j][i] = new unsigned short int [shortCount];
         }
   }
   int tenpercent = Bit64==64? ((((markerCount-1)>>10)+1)<<6): ((((markerCount-1)>>8)+1)<<4); // 6.25%
   if(tenpercent == 0) error("No autosome SNPs");
   if(lessmemFlag){
      int one = (Bit64==64? 64: 16);
      tenpercent = (tenpercent < one * defaultMaxCoreCount)? one: one*defaultMaxCoreCount;
   }
   unsigned char **localG_SNP = new unsigned char *[tenpercent];
   for(int m = 0;  m < tenpercent; m++)
      localG_SNP[m] = new unsigned char [SNP_space];
   int pos = 0;
   int posB = 0;
   const int wholeblocks = ((SNP_space-1) >> 17);
   const int leftoversize = SNP_space & 0x1FFFF;
   const int pvalidMarkerLength = (*pvalidMarker).Length();
   unsigned char **sexG_SNP[4];
   int sexPos[4];
   int sexCount[4];
   if(autosomeOnly){
      if(SNP_space > 0x20000){ // At least 128KB, or N > 512,000
         for(int m = 0; m < pvalidMarkerLength; m++)
            if((*pvalidMarker)[m]==1){
               for(int i = 0; i < wholeblocks; i++)
                  if(fread(&localG_SNP[posB][i<<17], 1, 0x20000, fp) < 0x20000)
                     error("Not enough genotypes at the %dth marker\n", m);
               if(fread(&localG_SNP[posB][wholeblocks<<17], 1, leftoversize, fp) < leftoversize)
                  error("Not enough genotypes at the %dth marker\n", m);
               pos ++;
               posB ++;
               if(posB == tenpercent){
                  if(!lessmemFlag) {
                     printf("%d%%\r", int((pos-posB)/tenpercent*6.25+0.5));
                     fflush(stdout);
                  }
                  if(Bit64==64)
                     ConvertMajorFromSNPtoIndividual64Bit(localG_SNP, LG, idCount, pos-tenpercent, pos-1);
                  else
                     ConvertMajorFromSNPtoIndividual(localG_SNP, GG, idCount, pos-tenpercent, pos-1);
                  posB = 0;
               }
            }else// non-autosome
               fseek(fp, SNP_space, SEEK_CUR);
      }else{   // N <= 512,000
         for(int m = 0; m < pvalidMarkerLength; m++)
            if((*pvalidMarker)[m]==1){
               if(fread(localG_SNP[posB], 1, SNP_space, fp)<SNP_space)
                  error("Not enough genotypes at the %dth marker\n", m);
               pos ++;
               posB ++;
               if(posB == tenpercent){
                  if(!lessmemFlag){
                     printf("%d%%\r", int((pos-posB)/tenpercent*6.25+0.5));
                     fflush(stdout);
                  }
                  if(Bit64==64)
                     ConvertMajorFromSNPtoIndividual64Bit(localG_SNP, LG, idCount, pos-tenpercent, pos-1);
                  else
                     ConvertMajorFromSNPtoIndividual(localG_SNP, GG, idCount, pos-tenpercent, pos-1);
                  posB = 0;
               }
            }else// non-autosome
               fseek(fp, SNP_space, SEEK_CUR);
      }
   }else{   // not autosomeOnly
      sexCount[0] = (((xmarkerCount-1)>>4)+1)<<4;
      sexCount[1] = (((ymarkerCount-1)>>4)+1)<<4;
      sexCount[2] = 0;
      sexCount[3] = (((mtmarkerCount-1)>>4)+1)<<4;
      for(int i = 0; i < 4; i++){
         if(sexCount[i]){
            sexG_SNP[i] = new unsigned char *[sexCount[i]];
            for(int j = 0; j < sexCount[i]; j++)
               sexG_SNP[i][j] = new unsigned char [SNP_space];
         }else
            sexG_SNP[i] = NULL;
         sexPos[i] = 0;
      }
      for(int m = 0; m < pvalidMarkerLength; m++)
         if((*pvalidMarker)[m]==1){
            for(int i = 0; i < wholeblocks; i++) // 128KB or N=512,000 at a time
               if(fread(&localG_SNP[posB][i<<17], 1, 0x20000, fp) < 0x20000)
                  error("Not enough genotypes at the %dth marker\n", m);
            if(fread(&localG_SNP[posB][wholeblocks<<17], 1, leftoversize, fp) < leftoversize)
               error("Not enough genotypes at the %dth marker\n", m);
            pos ++;
            posB ++;
            if(posB == tenpercent){
               if(!lessmemFlag){
                  printf("%d%%\r", int((pos-posB)/tenpercent*6.25+0.5));
                  fflush(stdout);
               }
               if(Bit64==64)
                  ConvertMajorFromSNPtoIndividual64Bit(localG_SNP, LG, idCount, pos-tenpercent, pos-1);
               else
                  ConvertMajorFromSNPtoIndividual(localG_SNP, GG, idCount, pos-tenpercent, pos-1);
               posB = 0;
            }
         }else{// non-autosome
            if((((*pvalidMarker)[m]!=SEXCHR) && ((*pvalidMarker)[m]!=SEXCHR+1) && ((*pvalidMarker)[m]!=SEXCHR+3)))
               fseek(fp, SNP_space, SEEK_CUR);
            else{
               int sex = (*pvalidMarker)[m] - SEXCHR;
               for(int i = 0; i < wholeblocks; i++) // 128KB or N=512,000 at a time
                  if(fread(&sexG_SNP[sex][sexPos[sex]][i<<17], 1, 0x20000, fp) < 0x20000)
                  error("Not enough genotypes at the %dth marker on sex chr %d\n", m, SEXCHR+sex);
               if(fread(&sexG_SNP[sex][sexPos[sex]][wholeblocks<<17], 1, leftoversize, fp) < leftoversize)
                  error("Not enough genotypes at the %dth marker on sex chr %d\n", m, SEXCHR+sex);
               sexPos[sex] ++;
            }
         }
   }  // end of if autosome
   fclose(fp);
   printf("  PLINK binary genotypes loaded.\n");
   delete pvalidMarker;
   if((idCount < 100 || !bigdataFlag) && (SaveFormat != "KING")){
      bigdataFlag = false;
      bigdataIdx.Dimension(markerCount);
      for(int m = 0; m < markerCount; m++)
         bigdataIdx[m] = markerCount - 1 - m;
   }
   if(posB){
      if(!lessmemFlag){
         printf("%d%%\r", int((pos-posB)/tenpercent*6.25+0.5));
         fflush(stdout);
      }
      if(Bit64==64)
         ConvertMajorFromSNPtoIndividual64Bit(localG_SNP, LG, idCount, pos-posB, pos-1);
      else
         ConvertMajorFromSNPtoIndividual(localG_SNP, GG, idCount, pos-posB, pos-1);
   }
   printf("  KING format genotype data successfully converted.\n");
   for(int m = 0; m < tenpercent; m++)
      delete []localG_SNP[m];
   delete []localG_SNP;
   if(!autosomeOnly){
      int samplesize = SNP_space<<2;
      if(xmarkerCount){
         xshortCount = (xmarkerCount-1)/16+1;
         for(int j = 0; j < 2; j++){
            XG[j] = new unsigned short int * [samplesize];
            for(int i = 0; i < samplesize; i++)
               XG[j][i] = new unsigned short int [xshortCount];
         }
         ConvertMajorFromSNPtoIndividual(sexG_SNP[0], XG, idCount, 0, xmarkerCount-1);
         for(int m = 0; m < sexCount[0]; m++)
            delete []sexG_SNP[0][m];
         delete []sexG_SNP[0];
      }else
         xshortCount = 0;
      if(ymarkerCount){
         yshortCount = (ymarkerCount-1)/16+1;
         for(int j = 0; j < 2; j++){
            YG[j] = new unsigned short int * [samplesize];
            for(int i = 0; i < samplesize; i++)
               YG[j][i] = new unsigned short int [yshortCount];
         }
         ConvertMajorFromSNPtoIndividual(sexG_SNP[1], YG, idCount, 0, ymarkerCount-1);
         for(int m = 0; m < sexCount[1]; m++)
            delete []sexG_SNP[1][m];
         delete []sexG_SNP[1];
      }else
         yshortCount = 0;
      if(mtmarkerCount){
         mtshortCount = (mtmarkerCount-1)/16+1;
         for(int j = 0; j < 2; j++){
            MG[j] = new unsigned short int * [samplesize];
            for(int i = 0; i < samplesize; i++)
               MG[j][i] = new unsigned short int [mtshortCount];
         }
         ConvertMajorFromSNPtoIndividual(sexG_SNP[3], MG, idCount, 0, mtmarkerCount-1);
         for(int m = 0; m < sexCount[3]; m++)
            delete []sexG_SNP[3][m];
         delete []sexG_SNP[3];
      }else
         mtshortCount = 0;
   }
}

void KingEngine::ConvertLGtoSLG(unsigned long long int **oldG[2], int oldCount, unsigned long long int **newG[2], int newCount)
{
   int oldlongCount = ((oldCount-1)>>6)+1; // 1-64 -> 1
   int newlongCount = ((newCount-1)>>6)+1; // 1-64 -> 1
   for(int j = 0; j < 2; j++){
      newG[j] = new unsigned long long int * [idCount];
      for(int i = 0; i < idCount; i++)
         newG[j][i] = new unsigned long long int [newlongCount];
   }
   Vector HetBySample(oldCount);
   HetBySample.Zero();
   const double rareHetCount = idCount*0.095;
   int aaCount;
   IntArray AACounts, AaCounts, missingCounts;
   ComputeAlleleFrequency64Bit(AACounts, AaCounts, missingCounts);
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
      private(aaCount)
#endif
   for(int m = 0; m < oldCount; m++){
      aaCount = idCount - AACounts[m] - AaCounts[m] - missingCounts[m];
      if(AaCounts[m] > rareHetCount){   // common variants
         HetBySample[m] = (4.0 / AaCounts[m]) * AACounts[m] * aaCount;
         if(AaCounts[m] < HetBySample[m]) HetBySample[m] = AaCounts[m];   // not way too many heterozygotes
         double temp = (AaCounts[m]*0.5+AACounts[m])/(AACounts[m]+AaCounts[m]+aaCount)*(AaCounts[m]*0.5+aaCount)*2;
         if(temp < HetBySample[m]) HetBySample[m] = temp;
         if(HetBySample[m] < rareHetCount)
            HetBySample[m] = rareHetCount;
      }else   // rare variants
         HetBySample[m] = AaCounts[m];
      if(AaCounts[m]+AACounts[m] == 0 || AaCounts[m]+aaCount == 0)  //monomorphic
         HetBySample[m] = -1;
   }  // Marker bit loop ends
   bigdataIdx.Index(HetBySample);
   IntArray allSNPs(oldCount);
   allSNPs.Set(-1);
   for(int m = 0; m < newCount; m++)
      allSNPs[bigdataIdx[oldCount-1-m]] = m;
   unsigned int *oldwords = new unsigned int [newCount];
   unsigned int *oldbits = new unsigned int [newCount];
   unsigned int *newwords = new unsigned int [newCount];
   unsigned int *newbits = new unsigned int [newCount];
   int index = 0;
   for(int m = 0; m < oldCount; m++)
      if(allSNPs[m] > -1){
         oldwords[index] = (m>>6);
         oldbits[index] = m&0x3F;
         newwords[index] = (allSNPs[m]>>6);
         newbits[index] = allSNPs[m]&0x3F;
         index++;
      }
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount)
{
#endif
   unsigned long long int *aGG = new unsigned long long int [newlongCount];
#ifdef _OPENMP
   #pragma omp for
#endif
   for(int i = 0; i < idCount; i++){
      for(int k = 0; k < 2; k++){
         for(int w = 0; w < newlongCount; w++)
            aGG[w] = 0;
         for(int n = 0; n < newCount; n++)
            aGG[newwords[n]] |= (((oldG[k][i][oldwords[n]]>>oldbits[n])&1)<<newbits[n]);
         memcpy(newG[k][i], aGG, (newlongCount<<3));
      }
   }
   delete []aGG;
#ifdef _OPENMP
}
#endif
   if(newCount&0x3F){  // last few SNPs need to be masked
      unsigned long long int mask = ((unsigned long long int)1<<(newCount&0x3F))-1;
      for(int k = 0; k < 2; k++)
         for(int i = 0; i < idCount; i++)
            newG[k][i][newlongCount-1] &= mask;
   }
   delete []oldwords;
   delete []oldbits;
   delete []newwords;
   delete []newbits;
}

void KingEngine::ConvertMajorFromSNPtoIndividual64Bit(unsigned char **localG_SNP, unsigned long long int **localLG[], int localidCount, int startPos, int stopPos)
{
   if(startPos&0x3F) error("Start Position should be multiple of 64 for SNP2Individual conversion");
   int startWord = startPos >> 6;
   int stopWord = stopPos >> 6;
   int SNP_space = (localidCount+3)/4; // 1 for 1-4, 2 for 5-8 etc.
   int posCount = stopPos - startPos + 1;
   if(localidCount%4)// the last byte has extra AA's
      for(int i = 0; i < 4-(localidCount%4); i++){
         int maskbit = base[(3-i)*2];
         for(int m = 0;  m < posCount; m++)
            localG_SNP[m][SNP_space-1] |= maskbit;
      }
   char revbase[256];
   for(int i = 0; i < 8; i++)
      revbase[base[i]] = i;
   char rightmost[256];
   for(int i = 0; i < 256; i++)
      rightmost[i] = revbase[i&(-i)];
   const int BLOCKSIZE=128; // Blocking size for individual groups
   const int BLOCKSIZE_BIT=(BLOCKSIZE<<3);
   unsigned long long int aword[BLOCKSIZE_BIT];   // BLOCKSIZE*8=1024
   const int CACHESIZE=16; // Blocking size for both SNP words
   const char CACHESIZE_SHIFT=4;
   unsigned long long int tempword[2][8192]; // BLOCKSIZE*4 * CACHESIZE = 8192
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) private(tempword, aword)
#endif
   for(int word = startWord; word <= stopWord; word += CACHESIZE){
      int wMax = (word > stopWord+1-CACHESIZE) ? stopWord+1: word+CACHESIZE;
      int memsize = (wMax-word)<<3;
      for(int group = 0; group < SNP_space; group += BLOCKSIZE){
         int gMax = (group > SNP_space-BLOCKSIZE) ? SNP_space: group+BLOCKSIZE;
         int iMax = (gMax-group)<<2;
         for(int w = word; w < wMax; w++){
            int ww = w - word;
            for(int g = 0; g < BLOCKSIZE_BIT; g++)  // BLOCKSIZE * 8
               aword[g] = 0;
            for(int b = 0; b < 64; b++){
               int m = ((w<<6)|b) - startPos;
               unsigned long long int longword = ((unsigned long long int)1<<b);
               for(int g = group; g < gMax; g++)
                  for(unsigned char byte = (~localG_SNP[m][g])&255; byte; byte &= (byte-1))
                     aword[((g-group)<<3)|rightmost[byte]] |= longword;
            }
            for(int i = 0; i < iMax; i++){
               unsigned long long int w1 = aword[i<<1];
               unsigned long long int w2 = aword[(i<<1)|1];
               int pos = (i<<CACHESIZE_SHIFT)|ww;
               tempword[0][pos] = (w1 & w2) | (~(w1 | w2));
               tempword[1][pos] = w1;
            }  // end of individual
         }  // end of w
         for(int k = 0; k < 2; k++)
            for(int i = group*4; i < gMax*4; i++)
               memcpy(&localLG[k][i][word], &tempword[k][((i-(group<<2))<<CACHESIZE_SHIFT)], memsize);
      }  // end of SNP word block
   }  // end of individual block
   if((stopPos+1)%64){  // last few SNPs need to be masked
      unsigned long long int mask = ((unsigned long long int)1<<((stopPos+1)%64))-1;
      for(int j = 0; j < 2; j++)
         for(int i = 0; i < localidCount; i++)
            localLG[j][i][stopWord] &= mask;
   }
}

void KingEngine::ConvertGGtoSG(unsigned short int **oldG[2], int oldCount, unsigned short int **newG[2], int newCount)
{
   int oldshortCount = (oldCount-1)/16+1; // 1-16->1
   int newshortCount = (newCount-1)/16+1; // 1-16->1
   for(int j = 0; j < 2; j++){
      newG[j] = new unsigned short int * [idCount];
      for(int i = 0; i < idCount; i++)
         newG[j][i] = new unsigned short int [newshortCount];
   }
   int N_Aa, N_AA;
   int AA, Aa, missing;
   double *hetinfo = new double[oldshortCount*16];
   Vector HetBySample(oldCount);
   HetBySample.Zero();
   const double rareHetCount = idCount*0.095;

   char revbase[65536];
   for(int i = 0; i < 16; i++)
      revbase[shortbase[i]] = i;
   char rightmost[65536];
   for(int i = 0; i < 65536; i++)
      rightmost[i] = revbase[i&(-i)];

   unsigned short int word;
   const int CACHESIZE=256;   // blocking size for SNP words
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount)\
      private(AA, Aa, missing, word)
{
#endif
   int AACount[4096], AaCount[4096], missingCount[4096], aaCount;
   double info[4096];
   int wMax;
#ifdef _OPENMP
   #pragma omp for
#endif
   for(int blockw = 0; blockw < oldshortCount; blockw += CACHESIZE){
      wMax = (blockw > oldshortCount-CACHESIZE) ? oldshortCount: blockw+CACHESIZE;
      for(int j = 0; j < 4096; j++)
         AACount[j] = AaCount[j] = missingCount[j] = 0;
      for(int i = 0; i < idCount; i++)
         for(int w = blockw; w < wMax; w++){
            int ww = (w-blockw)*16;
            AA = oldG[0][i][w] & oldG[1][i][w];
            Aa = (~oldG[0][i][w]) & oldG[1][i][w];
            missing = (~oldG[0][i][w]) & (~oldG[1][i][w]) & 0xFFFF;
            for(word = AA; word; word &= (word-1))
               AACount[ww+rightmost[word]] ++;
            for(word = Aa; word; word &= (word-1))
               AaCount[ww+rightmost[word]] ++;
            for(word = missing; word; word &= (word-1))
               missingCount[ww+rightmost[word]] ++;
         }
      for(int j = 0; j < 4096; j++){
         int aaCount = idCount - AACount[j] - AaCount[j] - missingCount[j];
         if(AaCount[j] > rareHetCount){   // common variants
            info[j] = (4.0 / AaCount[j]) * AACount[j] * aaCount;
            if(AaCount[j] < info[j]) info[j] = AaCount[j];   // not way too many heterozygotes
            double temp = (AaCount[j]*0.5+AACount[j])/(AACount[j]+AaCount[j]+aaCount)*(AaCount[j]*0.5+aaCount)*2;
            if(temp < info[j]) info[j] = temp;
            if(info[j] < rareHetCount)
               info[j] = rareHetCount;
         }else{   // rare variants
            info[j] = AaCount[j];
         }
         if(AaCount[j]+AACount[j] == 0 || AaCount[j]+aaCount == 0)  //monomorphic
            info[j] = -1;
      }  // Marker bit loop ends
      for(int j = 0; j < 4096; j++)
         if(blockw*16+j < oldCount)
            hetinfo[blockw*16+j] = info[j];
   } // Marker blockw loop ends
#ifdef _OPENMP
}
#endif
   for(int m = 0; m < oldCount; m++)
      HetBySample[m] = hetinfo[m];
   bigdataIdx.Index(HetBySample);

   IntArray allSNPs(oldCount);
   allSNPs.Set(-1);
   for(int m = 0; m < newCount; m++)
      allSNPs[bigdataIdx[oldCount-1-m]] = m;
   unsigned int *oldwords = new unsigned int [newCount];
   unsigned int *oldbits = new unsigned int [newCount];
   unsigned int *newwords = new unsigned int [newCount];
   unsigned int *newbits = new unsigned int [newCount];

   int index = 0;
   for(int m = 0; m < oldCount; m++)
      if(allSNPs[m] > -1){
         oldwords[index] = m/16;
         oldbits[index] = m%16;
         newwords[index] = allSNPs[m]/16;
         newbits[index] = allSNPs[m]%16;
         index++;
      }
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount)
{
#endif
   unsigned short int *aGG = new unsigned short int [newshortCount];
#ifdef _OPENMP
   #pragma omp for
#endif
   for(int i = 0; i < idCount; i++){
      for(int k = 0; k < 2; k++){
         for(int w = 0; w < newshortCount; w++)
            aGG[w] = 0;
         for(int n = 0; n < newCount; n++)
            aGG[newwords[n]] |= (((oldG[k][i][oldwords[n]]>>oldbits[n])&1)<<newbits[n]);
         for(int w = 0; w < newshortCount; w++)
            newG[k][i][w] = aGG[w];
      }
   }
   delete []aGG;
#ifdef _OPENMP
}
#endif
   if(newCount%16){  // last few SNPs need to be masked
      int mask = (1<<(newCount%16))-1;
      for(int k = 0; k < 2; k++)
         for(int i = 0; i < idCount; i++)
            newG[k][i][newshortCount-1] &= mask;
   }
   delete []oldwords;
   delete []oldbits;
   delete []newwords;
   delete []newbits;
   delete []hetinfo;
}

void KingEngine::ConvertMajorFromSNPtoIndividual(unsigned char **localG_SNP, unsigned short int **localGG[], int localidCount, int startPos, int stopPos)
{
   if(startPos&0xF) error("Start Position should be multiple of 16 for SNP2Individual conversion");
   int startWord = startPos >> 4;
   int stopWord = stopPos >> 4;
   int SNP_space = (localidCount+3)/4; // 1 for 1-4, 2 for 5-8 etc.
   int posCount = stopPos - startPos + 1;
   if(localidCount%4)// the last byte has extra AA's
      for(int i = 0; i < 4-(localidCount%4); i++){
         int maskbit = base[(3-i)*2];
         for(int m = 0;  m < posCount; m++)
            localG_SNP[m][SNP_space-1] |= maskbit;
      }
   char revbase[256];
   for(int i = 0; i < 8; i++)
      revbase[base[i]] = i;
   char rightmost[256];
   for(int i = 0; i < 256; i++)
      rightmost[i] = revbase[i&(-i)];
   const int BLOCKSIZE=128;
   const int BLOCKSIZE_BIT=(BLOCKSIZE<<3);
   unsigned short int aword[BLOCKSIZE_BIT];  // BLOCKSIZE*8=1024
   const int CACHESIZE=16; // Blocking size for SNP words
   const char CACHESIZE_SHIFT=4;
   unsigned short int tempword[2][8192]; // BLOCKSIZE*4 * CACHESIZE = 8192
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) private(tempword, aword)
#endif
   for(int word = startWord; word <= stopWord; word += CACHESIZE){
      int wMax = (word > stopWord+1-CACHESIZE) ? stopWord+1: word+CACHESIZE;
      int memsize = (wMax - word) << 1;
      for(int group = 0; group < SNP_space; group += BLOCKSIZE){
         int gMax = (group > SNP_space-BLOCKSIZE) ? SNP_space: group+BLOCKSIZE;
         int iMax = (gMax - group) << 2;
         for(int w = word; w < wMax; w++){
            int ww = w - word;
            for(int g = 0; g < BLOCKSIZE_BIT; g++)  // CACHESIZE * 8
               aword[g] = 0;
            for(int b = 0; b < 16; b++){
               int m = (w<<4|b)-startPos;
               unsigned short int shortword = shortbase[b];
               for(int g = group; g < gMax; g++)
                  for(unsigned char byte = (~localG_SNP[m][g])&255; byte; byte &= (byte-1))
                     aword[((g-group)<<3)+rightmost[byte]] |= shortword;
            }
            for(int i = 0; i < iMax; i++){
               unsigned short int w1 = aword[i<<1];
               unsigned short int w2 = aword[(i<<1)|1];
               int pos = (i<<CACHESIZE_SHIFT)|ww;
               tempword[0][pos] = (w1 & w2) | (~(w1 | w2));
               tempword[1][pos] = w1;
            }  // end of individual
         }
         for(int k = 0; k < 2; k++)
            for(int i = group*4; i < gMax*4; i++)
               memcpy(&localGG[k][i][word], &tempword[k][((i-(group<<2))<<CACHESIZE_SHIFT)], memsize);
      }  // end of SNP word block
   }  // end of individual block
   if((stopPos+1)%16){  // last few SNPs need to be masked
      unsigned short int mask = ((unsigned short int)1<<((stopPos+1)%16))-1;
      for(int j = 0; j < 2; j++)
         for(int i = 0; i < localidCount; i++)
            localGG[j][i][stopWord] &= mask;
   }
}

void KingEngine::ReadPlinkBinaryBigDataSNPMajor(const char *filename, const char *snplist, const char *famfile, const char *bimfile, const char *covfile, const char *phefile)
{
   printf("Loading genotype data in PLINK binary format...\n");
   StringArray filenames;
   filenames.AddTokens(filename, ',');
   if(filenames.Length()>1)
      error("Multiple files cannot be read currently.");
   filenames.Clear();

   String name=filename;
   int p = name.Find(".bed");
   if(p==-1)
      plinkinput = name;
   else
      plinkinput = name.SubStr(0, p);
   ReadPlinkPedFile(plinkinput, (famfile==NULL || famfile[0]=='\0') ? NULL: famfile);
   ReadPlinkCovFile(plinkinput, (covfile==NULL || covfile[0]=='\0') ? NULL: covfile);
   if(!genotypeOnly)
      ReadPlinkTraitFile(plinkinput, (phefile==NULL || phefile[0]=='\0') ? NULL: phefile);
   IntArray *pvalidMarker = new IntArray;
   ReadPlinkMapFile(plinkinput, *pvalidMarker, (bimfile==NULL || bimfile[0]=='\0') ? NULL: bimfile);
   for(int i = 0; i < ped.count; i++)
      if(ped[i].ngeno == 10000)
         ped[i].ngeno = markerCount+xmarkerCount+ymarkerCount+mtmarkerCount;
   if(markerCount){
      if(Bit64==64) longCount = (markerCount-1)/64+1;
      shortCount = (markerCount-1)/16+1;
   }else
      longCount = shortCount = 0;

   long int SNP_space = (idCount+3)/4; // 1 for 1-4, 2 for 5-8 etc.
   StringIntHash markerLookup;
   markerLookup.Clear();
   for(int i = 0; i < markerCount; i++)
      markerLookup.SetInteger(snpName[i], i);

   IntArray snpidx(0);
   String line;
   StringArray tokens;
   IFILE input;
   input=ifopen((const char*)snplist, "rt");
   if(input==NULL) error("Cannot open %s to read", (const char*)snplist);
   printf("Read in SNP list from file %s...\n", (const char*)snplist);
   tokens.Clear();
   line.ReadLine(input);
   tokens.AddTokens(line);
   if(tokens.Length() < 1) error("The first line is empty");
   int proxyCol=0;
   if(tokens[0]!="SNP"){
      ifclose(input);
      input = ifopen((const char*)snplist, "rt");
   }else{
      for(int i = 1; i < tokens.Length(); i++)
         if(tokens[i]=="PROXY"){
            proxyCol=i;
            break;
         }
   }
   while(!ifeof(input)){
      tokens.Clear();
      line.ReadLine(input);
      tokens.AddTokens(line);
      if(tokens.Length() < 1) continue;
      int k = markerLookup.Integer(tokens[0]);
      if(k==-1){  // SNP not available
         if(proxyCol && tokens.Length()>proxyCol)
            k = markerLookup.Integer(tokens[proxyCol]);
         if(k==-1){
            printf("SNP %s from file %s cannot be found in data %s.bim\n",
            (const char*)tokens[0], (const char*)snplist, (const char*)plinkinput);
            continue;
         }
      }
      snpidx.Push(k);
   }
   ifclose(input);

   QuickIndex index;
   index.Index(snpidx);
   sortedSNP.Dimension(0);
   for(int i = 0; i < snpidx.Length(); i++)
      sortedSNP.Push(snpidx[index[i]]);

   int validMarkerCount = sortedSNP.Length();
   G_SNP = new unsigned char *[validMarkerCount];
   for(int m = 0;  m < validMarkerCount; m++)
      G_SNP[m] = new unsigned char[SNP_space];

   String bedfile = plinkinput;
   bedfile.Add(".bed");
   FILE *fp = fopen(bedfile, "rb");
   if(fp != NULL){
      printf("Read in PLINK bed file %s...\n", (const char*)bedfile);
      unsigned char byte;
      byte = fgetc(fp);
      if(byte != 108)
         error("The PLINK format in %s cannot be recognized. The first byte is %x",
            (const char*)bedfile, byte);
      byte = fgetc(fp);
      if(byte != 27)
         error("The PLINK format in %s cannot be recognized. The second byte is %x",
            (const char*)bedfile, byte);
      byte = fgetc(fp);
      if(byte != 1)
         error("Currently only SNP-major mode can be analyzed.");
      unsigned char bit1, bit2;
      unsigned short int intbit1, intbit2;
      int xmarkerIndex = 0;
      int ymarkerIndex = 0;
      int mtmarkerIndex = 0;
      int newbyte;
      unsigned long long int newbit;
      int oldbyte;

      if(snplist==NULL){
         for(int m = 0; m < validMarkerCount; m++)
            if(fread(G_SNP[m], sizeof(unsigned char), SNP_space, fp)<SNP_space)
               error("Not enough genotypes at %dth SNP %s\n",
                  sortedSNP[m], (const char*)snpName[sortedSNP[m]]);
      }else
      for(int m = 0; m < validMarkerCount; m++){
         if(m==0)
            fseek(fp, SNP_space * sortedSNP[m], SEEK_CUR);
         else
            fseek(fp, SNP_space * (sortedSNP[m] - sortedSNP[m-1] - 1), SEEK_CUR);
         if(fread(G_SNP[m], sizeof(unsigned char), SNP_space, fp)<SNP_space)
            error("Not enough genotypes at %dth SNP %s\n",
               sortedSNP[m], (const char*)snpName[sortedSNP[m]]);
      }
   }
   fclose(fp);
   if(idCount%4)// the last byte has extra AA's
      for(int i = 0; i < 4-(idCount%4); i++){
         int maskbit = base[(3-i)*2];
         for(int m = 0;  m < validMarkerCount; m++)
            G_SNP[m][SNP_space-1] |= maskbit;
      }
   printf("  PLINK binary genotypes at %d SNPs are loaded.\n", sortedSNP.Length());
}


void KingEngine::ReadPlinkPedFile(const char *filename, const char *famfile)
{
   String line, tempS;
   StringArray tokens;
   FID.Dimension(0);
   PID.Dimension(0);
   FA.Dimension(0);
   MO.Dimension(0);
   SEX.Dimension(0);
   sampleName.Dimension(0);

   affectionNames.Dimension(1);
   affectionNames[0] = "DISEASEKING";
   IntArray disease(0);
   String pedfile;
   if(famfile)
      pedfile = famfile;
   else{
      pedfile = filename;
      pedfile.Add(".fam");
   }
   IFILE input = ifopen(pedfile, "rt");
   if(input == NULL){
      String pedfile2 = pedfile;
      pedfile2.Add(".gz");
      input = ifopen(pedfile2, "rt");
      if(input == NULL)
         error("Pedigree file %s cannot be opened", (const char*)pedfile);
      printf("Read in PLINK fam file %s...\n", (const char*)pedfile2);
   }else
      printf("Read in PLINK fam file %s...\n", (const char*)pedfile);
   while(!ifeof(input)){
      tokens.Clear();
      line.ReadLine(input);
      tokens.AddTokens(line);
      if(tokens.Length() < 6) continue;
      FID.Push(tokens[0]);
      PID.Push(tokens[1]);
      tempS = tokens[0];
      tempS.Add("->"); // '->'
      tempS.Add(tokens[1]);
      sampleName.Push(tempS);
      FA.Push(tokens[2]);
      MO.Push(tokens[3]);
      if(tokens[4].AsInteger() == -9)
         SEX.Push("0");
      else
         SEX.Push(tokens[4]);
      if(tokens[5].AsInteger() == -9)
         disease.Push(0);
      else
         disease.Push(tokens[5].AsInteger());
   }
   idCount = sampleName.Length();
   ifclose(input);

   StringIntHash sampleLookup;
   sampleLookup.Clear();
   for(int i = 0; i < idCount; i++)
      sampleLookup.SetInteger((const char*)sampleName[i], i);

   // check if pedigrees are complete
   IntArray inclusionFA(0);
   IntArray inclusionMO(0);
   StringArray newName = sampleName;
   int newid = 1;
   for(int i = 0; i < idCount; i++){
      if(FA[i]!="0"){
         tempS = FID[i];
         tempS.Add("->"); // '->'
         tempS.Add(FA[i]);
         if(sampleLookup.Integer(tempS)==-1){
            inclusionFA.Push(i);
            sampleLookup.SetInteger(tempS, newName.Length());
            newName.Push(tempS);
         }
      }
      if(MO[i]!="0"){
         tempS = FID[i];
         tempS.Add("->"); // '->'
         tempS.Add(MO[i]);
         if(sampleLookup.Integer(tempS)==-1){
            inclusionMO.Push(i);
            sampleLookup.SetInteger(tempS, newName.Length());
            newName.Push(tempS);
         }
      }
      if(FA[i]!="0" && MO[i]=="0"){
         MO[i] = "KING";
         MO[i] += newid;
         tempS = FID[i];
         tempS.Add("->"); // '->'
         tempS.Add(MO[i]);
         newName.Push(tempS);
         inclusionMO.Push(i);
         newid++;
      }else if(FA[i]=="0" && MO[i]!="0"){
         FA[i] = "KING";
         FA[i] += newid;
         tempS = FID[i];
         tempS.Add("->"); // '->'
         tempS.Add(FA[i]);
         newName.Push(tempS);
         inclusionFA.Push(i);
         newid++;
      }
   }
   if(inclusionFA.Length() || inclusionMO.Length()) {
      for(int i = 0; i < inclusionFA.Length(); i++){
         FID.Push(FID[inclusionFA[i]]);
         PID.Push(FA[inclusionFA[i]]);
         FA.Push("0");
         MO.Push("0");
         SEX.Push("1");
         disease.Push(0);
      }
      for(int i = 0; i < inclusionMO.Length(); i++){
         FID.Push(FID[inclusionMO[i]]);
         PID.Push(MO[inclusionMO[i]]);
         FA.Push("0");
         MO.Push("0");
         SEX.Push("2");
         disease.Push(0);
      }
   }
   affections.Dimension(1, FID.Length());
   for(int i = 0; i < FID.Length(); i++)
      affections[0][i] = disease[i];

   markerCount = 10000;
   MakePed();
   markerCount = 0;
   printf("  PLINK pedigrees loaded: %d samples\n", idCount);
}

void KingEngine::ReadPlinkCovFile(const char *filename, const char *covfile)
{
   StringArray tokens;
   String line;
   String mycovfile;
   if(covfile)
      mycovfile = covfile;
   else{
      mycovfile = filename;
      mycovfile.Add(".cov");
   }
   IFILE input = ifopen(mycovfile, "rt");
   if(input == NULL){
      String mycovfile2 = mycovfile;
      mycovfile2.Add(".gz");
      input = ifopen(mycovfile2, "rt");
   }
   if(input){
      printf("Read in PLINK covariate file %s...\n", (const char*)mycovfile);
      tokens.Clear();
      line.ReadLine(input);
      tokens.AddTokens(line);
      if(tokens.Length() < 2 || tokens[0]!="FID" || tokens[1]!="IID"){
         warning("Covariate file %s not in PLINK format and not recognized.\n",
            (const char*)mycovfile);
      }else{
         StringIntHash FIDLookup;
         FIDLookup.Clear();
         for(int f = 0; f < ped.familyCount; f++)
            FIDLookup.SetInteger((const char*)ped.families[f]->famid, f);
         int twinPos = -1;
         ped.covariateNames.Dimension(tokens.Length()-2);
         ped.covariateCount = ped.covariateNames.Length();
         for(int i = 0; i < ped.covariateCount; i++)
            ped.covariateNames[i] = tokens[i+2];
         for(int i = 0; i < ped.count; i++){
            ped[i].covariates.Dimension(ped.covariateCount);
            ped[i].covariates.Set(_NAN_);
         }
         twinPos = ped.covariateNames.SlowFind("MZTWIN");
         if(twinPos > -1) ped.haveTwins = 1;
         while(!ifeof(input)){
            tokens.Clear();
            line.ReadLine(input);
            tokens.AddTokens(line);
            if(tokens.Length() < ped.covariateCount+2) continue;
            int f = FIDLookup.Integer(tokens[0]);
            if(f > -1){
               for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
                  if(ped[i].pid != tokens[1]) continue;
                  for(int j = 2; j < tokens.Length(); j++)
                     if(tokens[j] == "-9" || tokens[j] == ".")
                        ped[i].covariates[j-2] = _NAN_;
                     else
                        ped[i].covariates[j-2] = tokens[j].AsDouble();
                  if(twinPos > -1){
                     ped[i].zygosity = tokens[twinPos+2].AsInteger();
                     if(ped[i].zygosity == -9) ped[i].zygosity = 0;
                  }
               }
            }
         }
         if(ped.covariateCount){
            printf("  %d covariates loaded: ",ped.covariateCount);
            for(int t = 0; t < ped.covariateCount; t++)
               printf(" %s", (const char*)ped.covariateNames[t]);
            printf("\n");
            if(twinPos > -1) printf(" MZ twin status loaded\n");
         }
      }
      ifclose(input);
   }
}

void KingEngine::ReadPlinkTraitFile(const char *filename, const char *phefile)
{
   bool diseaseMissing = true;
   for(int i = 0; i < ped.count; i++)
      if(ped[i].affections[0]){
         diseaseMissing = false;
         break;
      }
   if(diseaseMissing)
      ped.affectionCount = 0;

   StringArray tokens;
   String line;
   String phenofile;
   if(phefile)
      phenofile = phefile;
   else{
      phenofile = filename;
      phenofile.Add(".phe");
   }
   IFILE input = ifopen(phenofile, "rt");
   if(input == NULL){
      String phenofile2 = phenofile;
      phenofile2.Add(".gz");
      input = ifopen(phenofile2, "rt");
   }
   if(input){
      Matrix pheno;
      printf("Read in PLINK phenotype file %s...\n", (const char*)phenofile);
      tokens.Clear();
      line.ReadLine(input);
      tokens.AddTokens(line);
      if(tokens.Length() < 2 || tokens[0]!="FID" || tokens[1]!="IID"){
         warning("Phenotype file %s not in PLINK format and not recognized.\n",
            (const char*)phenofile);
      }else{
         int myTraitCount = tokens.Length()-2;
         StringArray myTraitNames(myTraitCount);
         for(int i = 0; i < myTraitCount; i++)
            myTraitNames[i] = tokens[i+2];
         pheno.Dimension(ped.count, myTraitCount);
         IntArray isQuantitative(myTraitCount);
         pheno.Set(_NAN_);
         IntArray myindex(0);
         StringIntHash FIDLookup;
         FIDLookup.Clear();
         for(int f = 0; f < ped.familyCount; f++)
            FIDLookup.SetInteger((const char*)ped.families[f]->famid, f);
         while(!ifeof(input)){
            tokens.Clear();
            line.ReadLine(input);
            tokens.AddTokens(line);
            if(tokens.Length() < myTraitCount+2) continue;
            int f = FIDLookup.Integer(tokens[0]);
            if(f > -1){
               for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
                  if(ped[i].pid == tokens[1]){
                     for(int j = 2; j < tokens.Length(); j++)
                        if(tokens[j] == "-9" || tokens[j]==".")
                           pheno[i][j-2] = _NAN_;
                        else
                           pheno[i][j-2] = tokens[j].AsDouble();
                     break;
                  }
               }
            }
         }
         for(int t = 0; t < myTraitCount; t++){
            isQuantitative[t] = 0;
            for(int i = 0; i < ped.count; i++)
               if(pheno[i][t] != _NAN_ && pheno[i][t] != 0 && pheno[i][t] != 1 && pheno[i][t] != 2){
                  isQuantitative[t] = 1;
                  myindex.Push(t);
                  break;
               }
         }
         ped.traitCount = isQuantitative.Sum();
         ped.traitNames.Dimension(0);
         for(int t = 0; t < myTraitCount; t++)
            if(isQuantitative[t])
               ped.traitNames.Push(myTraitNames[t]);
         for(int i = 0; i < ped.count; i++){
            ped[i].traits.Dimension(ped.traitCount);
            ped[i].traits.Set(_NAN_);
            for(int t = 0; t < ped.traitCount; t++)
               ped[i].traits[t] = pheno[i][myindex[t]];
         }
         if(ped.traitCount){
            printf("  %d quantitative traits loaded: ",ped.traitCount);
            if(ped.traitCount < 1000){
               for(int t = 0; t < ped.traitCount; t++)
                  printf(" %s", (const char*)ped.traitNames[t]);
               printf("\n");
            }else{
               printf("(the first 10 traits) ");
               for(int t = 0; t < 10; t++)
                  printf(" %s", (const char*)ped.traitNames[t]);
               printf(" ...\n");
            }
         }
         if(myTraitCount > ped.traitCount){ // Binary traits
            int disease = 0;
            for(int i = 0; i < ped.count; i++){
               if(ped[i].affections[0]) disease = 1;
               break;
            }
            ped.affectionCount = disease + myTraitCount - ped.traitCount;
            ped.affectionNames.Dimension(ped.affectionCount);
            IntArray valueCount(4);
            for(int t = 0; t < myTraitCount; t++)
               if(!isQuantitative[t]){
                  valueCount.Zero();
                  for(int i = 0; i < ped.count; i++)
                     if(pheno[i][t] == -9) valueCount[0] ++;
                     else if(pheno[i][t] == 0) valueCount[1] ++;
                     else if(pheno[i][t] == 1) valueCount[2] ++;
                     else if(pheno[i][t] == 2) valueCount[3] ++;
                  if(!valueCount[3] && valueCount[1] && valueCount[2]) // 0/1 convert to 1/2
                     for(int i = 0; i < ped.count; i++){
                        if(pheno[i][t] == 1) ped[i].affections[disease] = 2;
                        else if(pheno[i][t] == 0) ped[i].affections[disease] = 1;
                        else if(pheno[i][t] == _NAN_) ped[i].affections[disease] = 0;
                     }
                  else // already 1/2
                     for(int i = 0; i < ped.count; i++)
                        ped[i].affections[disease] = pheno[i][t]==_NAN_? 0: int(pheno[i][t]);
                  ped.affectionNames[disease] = myTraitNames[t];
                  disease++;
               }
            printf("  %d affection traits loaded: ", ped.affectionCount);
            for(int t = 0; t < ped.affectionCount; t++)
               printf(" %s", (const char*)ped.affectionNames[t]);
            printf("\n");
         }
      }
      ifclose(input);
   }
}

void KingEngine::ReadPlinkMapFile(const char *filename, IntArray & validMarker, const char *bimfile)
{
   String line;
   StringArray tokens;
   String mapfile;
   if(bimfile)
      mapfile = bimfile;
   else{
      mapfile = filename;
      mapfile.Add(".bim");
   }
   FILE *fp = fopen(mapfile, "rt");

   chromosomes.Dimension(0);
   bp.Dimension(0);  //positions.Dimension(0); 
   snpName.Dimension(0);
   StringArray labels[2], xlabels[2], ylabels[2], mtlabels[2];
   labels[0].Dimension(0); labels[1].Dimension(0);
   validMarker.Dimension(0);
   int chr;
   int xymarkerCount = 0;

   if(fp){
      printf("Read in PLINK bim file %s...\n", (const char*)mapfile);
      if(genotypeOnly){
         while(!feof(fp)){
            tokens.Clear();
            line.ReadLine(fp);
            tokens.AddTokens(line);
            if(tokens.Length() < 6) continue;
            if(tokens[0]=="X" || tokens[0]=="x")
               chr = SEXCHR;
            else if(tokens[0]=="Y" || tokens[0]=="y")
               chr = SEXCHR+1;
            else if(tokens[0]=="XY" || tokens[0]=="xy")
               chr = SEXCHR+2;
            else if(tokens[0]=="MT" || tokens[0]=="mt")
               chr = SEXCHR+3;
            else
               chr = tokens[0].AsInteger();
            if(chr > 0 && chr < SEXCHR)
               validMarker.Push(1);
            else if(chr==SEXCHR+2){
               validMarker.Push(1);
               xymarkerCount++;
            }else if(chr==SEXCHR)
               validMarker.Push(SEXCHR);
            else if(chr==SEXCHR+1)
               validMarker.Push(SEXCHR+1);
            else if(chr==SEXCHR+3)
               validMarker.Push(SEXCHR+3);
            else
               validMarker.Push(0);
         }
      }else{   // if NOT genotypeOnly
         while(!feof(fp)){
            tokens.Clear();
            line.ReadLine(fp);
            tokens.AddTokens(line);
            if(tokens.Length() < 6) continue;
            if(tokens[0]=="X" || tokens[0]=="x")
               chr = SEXCHR;
            else if(tokens[0]=="Y" || tokens[0]=="y")
               chr = SEXCHR+1;
            else if(tokens[0]=="XY" || tokens[0]=="xy")
               chr = SEXCHR+2;
            else if(tokens[0]=="MT" || tokens[0]=="mt")
               chr = SEXCHR+3;
            else
               chr = tokens[0].AsInteger();
            if((chr < 1) || (chr > SEXCHR+3) // chromosome not 0-26
               || (chrList.Length() && (chrList.Find(chr)==-1 // chr not specified
               || (start != _NAN_ && tokens[3].AsDouble()*0.000001 < start) // pos not specified
               || (stop != _NAN_ && tokens[3].AsDouble()*0.000001 > stop) ) ) ){ // skip
               validMarker.Push(0);
            }else if(chr == SEXCHR){ // chromosome X
               validMarker.Push(SEXCHR);
               if(!autosomeOnly){
                  xsnpName.Push(tokens[1]);
                  xbp.Push(tokens[3].AsInteger()); //xpositions.Push(tokens[3].AsDouble()*0.000001);
                  xlabels[0].Push(tokens[4]);
                  xlabels[1].Push(tokens[5]);
               }
            }else if(chr == SEXCHR+1){ // chromosome Y
               validMarker.Push(SEXCHR+1);
               if(!autosomeOnly){
                  ysnpName.Push(tokens[1]);
                  ybp.Push(tokens[3].AsInteger()); //ypositions.Push(tokens[3].AsDouble()*0.000001); 
                  ylabels[0].Push(tokens[4]);
                  ylabels[1].Push(tokens[5]);
               }
            }else if(chr == SEXCHR+3){ // mitochondrion
               validMarker.Push(SEXCHR+3);
               if(!autosomeOnly){
                  mtsnpName.Push(tokens[1]);
                  mtbp.Push(tokens[3].AsInteger());   //mtpositions.Push(tokens[3].AsDouble()*0.000001); 
                  mtlabels[0].Push(tokens[4]);
                  mtlabels[1].Push(tokens[5]);
               }
            }else{
               if(chr == SEXCHR+2) xymarkerCount++;
               validMarker.Push(1);
               chromosomes.Push(chr);
               bp.Push(tokens[3].AsInteger());  //positions.Push(tokens[3].AsDouble()*0.000001); 
               snpName.Push(tokens[1]);
               labels[0].Push(tokens[4]);
               labels[1].Push(tokens[5]);
            }
         }
      }
      fclose(fp);
   }else{
      String mapfile2 = mapfile;
      mapfile2.Add(".gz");
      IFILE input = ifopen(mapfile2, "rt");
      if(input == NULL)
         error("Map file %s cannot be opened", (const char*)mapfile);
      printf("Read in PLINK bim file %s...\n", (const char*)mapfile2);
      if(!autosomeOnly){
      xbp.Dimension(0); //xpositions.Dimension(0);
      xsnpName.Dimension(0);
      xlabels[0].Dimension(0); xlabels[1].Dimension(0);
      ybp.Dimension(0); //ypositions.Dimension(0); 
      ysnpName.Dimension(0);
      ylabels[0].Dimension(0); ylabels[1].Dimension(0);
      mtbp.Dimension(0);   //mtpositions.Dimension(0); 
      mtsnpName.Dimension(0);
      mtlabels[0].Dimension(0); mtlabels[1].Dimension(0);
      }

      while(!ifeof(input)){
      tokens.Clear();
      line.ReadLine(input);
      tokens.AddTokens(line);
      if(tokens.Length() < 6) continue;

      if(tokens[0]=="X" || tokens[0]=="x")
         chr = SEXCHR;
      else if(tokens[0]=="Y" || tokens[0]=="y")
         chr = SEXCHR+1;
      else if(tokens[0]=="XY" || tokens[0]=="xy")
         chr = SEXCHR+2;
      else if(tokens[0]=="MT" || tokens[0]=="mt")
         chr = SEXCHR+3;
      else
         chr = tokens[0].AsInteger();
      if((chr < 1) || (chr > SEXCHR+3) // chromosome not 1-26
         || (chrList.Length() && (chrList.Find(chr)==-1 // chr not specified
         || (start != _NAN_ && tokens[3].AsDouble()*0.000001 < start) // pos not specified
         || (stop != _NAN_ && tokens[3].AsDouble()*0.000001 > stop) ) ) ){ // skip
         validMarker.Push(0);
      }else if(chr == SEXCHR){ // chromosome X
         validMarker.Push(SEXCHR);
         if(!autosomeOnly){
            xsnpName.Push(tokens[1]);
            xbp.Push(tokens[3].AsInteger()); //xpositions.Push(tokens[3].AsDouble()*0.000001);
            xlabels[0].Push(tokens[4]);
            xlabels[1].Push(tokens[5]);
         }
      }else if(chr == SEXCHR+1){ // chromosome Y
         validMarker.Push(SEXCHR+1);
         if(!autosomeOnly){
            ysnpName.Push(tokens[1]);
            ybp.Push(tokens[3].AsInteger()); //ypositions.Push(tokens[3].AsDouble()*0.000001); 
            ylabels[0].Push(tokens[4]);
            ylabels[1].Push(tokens[5]);
         }
      }else if(chr == SEXCHR+3){ // mitochondrion
         validMarker.Push(SEXCHR+3);
         if(!autosomeOnly){
            mtsnpName.Push(tokens[1]);
            mtbp.Push(tokens[3].AsInteger());   //mtpositions.Push(tokens[3].AsDouble()*0.000001); 
            mtlabels[0].Push(tokens[4]);
            mtlabels[1].Push(tokens[5]);
         }
      }else{
         if(chr == SEXCHR+2) xymarkerCount++;
         validMarker.Push(1);
         if(!genotypeOnly){
            chromosomes.Push(chr);
            bp.Push(tokens[3].AsInteger());  //positions.Push(tokens[3].AsDouble()*0.000001); 
            snpName.Push(tokens[1]);
            labels[0].Push(tokens[4]);
            labels[1].Push(tokens[5]);
         }
      }
      }
      ifclose(input);
   } // end of ifeof(input)
   markerCount = xmarkerCount = ymarkerCount = mtmarkerCount = 0;
   for(int m = 0; m < validMarker.Length(); m++)
      if(validMarker[m]==1) markerCount++;
      else if(validMarker[m]==SEXCHR) xmarkerCount++;
      else if(validMarker[m]==SEXCHR+1) ymarkerCount++;
      else if(validMarker[m]==SEXCHR+3) mtmarkerCount++;
   int totalmarkerCount = markerCount + xmarkerCount + ymarkerCount + mtmarkerCount;

   printf("  Genotype data consist of %d autosome SNPs", markerCount);
   if(xymarkerCount) printf(" (including %d XY SNPs)", xymarkerCount);
   if(xmarkerCount) printf(", %d X-chromosome SNPs", xmarkerCount);
   if(ymarkerCount) printf(", %d Y-chromosome SNPs", ymarkerCount);
   if(mtmarkerCount) printf(", %d mitochondrial SNPs", mtmarkerCount);
   printf("\n");

   if(totalmarkerCount  < validMarker.Length())
      printf("  %d other SNPs are removed.\n",
         validMarker.Length() - markerCount - xmarkerCount - ymarkerCount - mtmarkerCount);
   printf("  PLINK maps loaded: %d SNPs\n", markerCount + xmarkerCount + ymarkerCount + mtmarkerCount);

   if(!genotypeOnly)
      for(int i = 0; i < 2; i++)
         alleleLabel[i] = labels[i];
   if(!autosomeOnly && !genotypeOnly){ // allocate memory for sex-chr
      if(xmarkerCount)
         for(int i = 0; i < 2; i++)
            xalleleLabel[i] = xlabels[i];
      if(ymarkerCount)
         for(int i = 0; i < 2; i++)
            yalleleLabel[i] = ylabels[i];
      if(mtmarkerCount)
         for(int i = 0; i < 2; i++)
            mtalleleLabel[i] = mtlabels[i];
   } /// end of !autosomeOnly
}


