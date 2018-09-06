// phase.cpp
// 6/20/2011 Wei-Min Chen

#include "analysis.h"
#include "MathStats.h"
#include "Kinship.h"
#include "QuickIndex.h"

void Engine::Haplotyping(bool SaveFlag)
{
   int AA, Aa, missing;
   int b;
   int HHH, oneinformative, homhom;
   int ids[3];

   if(geno.Length()==0) BuildShortBinary();

   int AACount1, AaCount1, missingCount1;
   double frequency;
   rmarkers.Dimension(0);
   String labelInHap[2];
   labelInHap[0].Clear(); labelInHap[1].Clear();
   IntArray IsPhased(idCount);
   IsPhased.Zero();

   ExtractUnrelated();
   printf("A subset of %d unrelated individuals are used to calculate allele frequencies\n",
      unrelatedList.Length());

   int hom1, hom2, rareAllele, homCount[16];
   double myMAF = 0.5;
   if(rareMAF!=_NAN_) myMAF = rareMAF;
   for(int m = 0; m < markerCount; m++){
      if(chrList.Length() && chrList.Find(chromosomes[m])==-1) continue;
      if(start != _NAN_ && positions[m] < start) continue;
      if(stop != _NAN_ && positions[m] > stop) continue;

      b = m/16;
      int j = m%16;
      AACount1 = AaCount1 = missingCount1 = 0;
      frequency = 0.0;

      for(int i = 0; i < idCount; i++){
         AA = GG[0][i][b] & GG[1][i][b];
         Aa = (~GG[0][i][b]) & GG[1][i][b];
         missing = (~GG[0][i][b]) & (~GG[1][i][b]);
         if(AA & shortbase[j])
            AACount1 ++;
         else if(Aa & shortbase[j])
            AaCount1 ++;
         else if(missing & shortbase[j])
            missingCount1 ++;
      }
      if(missingCount1 < idCount){
         frequency = (AACount1 + AaCount1*0.5) / (idCount - missingCount1);
         if(frequency < myMAF || frequency > 1-myMAF)
            rmarkers.Push(m);
         if(frequency > 1-myMAF){ // flip
            for(int i = 0; i < idCount; i++)
               GG[1][i][b] ^= (GG[0][i][b] & shortbase[j]);
            labelInHap[0].Add(alleleLabel[1][m]);
            labelInHap[1].Add(alleleLabel[0][m]);
         }else if(frequency < myMAF){
            labelInHap[0].Add(alleleLabel[0][m]);
            labelInHap[1].Add(alleleLabel[1][m]);
         }
      }
   }
  IntArray thap[2][3];
  for(int h = 0; h < 2; h++){
      hap[h] = new IntArray[idCount];
      for(int i = 0; i < idCount; i++){
         hap[h][i].Dimension(rmarkers.Length());
         hap[h][i].Set(-10000);
      }
  }

   MakeTrioForTDT();
   if(Ltrio.Length()==0){
      warning("TDT analysis requires parent-affected-offspring trios.\n");
      return;
   }
   if(rmarkers.Length() > 10000){
      warning("Too many SNPs in this region.\n");
      return;
   }else if(rmarkers.Length()==0){
      warning("No SNPs are included in this region.\n");
      return;
   }else
      printf("%d SNPs are used to construct haplotype data.\n", rmarkers.Length());

   for(int i = 0; i < Ltrio.Length()/3; i++){
      ids[0] = geno[Ltrio[i*3+1]];
      ids[1] = geno[Ltrio[i*3+2]];
      ids[2] = geno[Ltrio[i*3]];   // ids[2] is the child
      for(int k = 0; k < 3; k++)
         if(IsPhased[ids[k]])
            for(int h = 0; h < 2; h++)
               thap[h][k] = hap[h][ids[k]];
      for(int l = 0; l < rmarkers.Length(); l++){
         int m = rmarkers[l];
         b = m/16;
         int j = m % 16;

         oneinformative = (GG[0][ids[2]][b] | GG[1][ids[2]][b]) & // offspring not missing
            ( (~GG[0][ids[0]][b] & GG[1][ids[0]][b] & GG[0][ids[1]][b]) |
              (~GG[0][ids[1]][b] & GG[1][ids[1]][b] & GG[0][ids[0]][b]) );

         HHH = ~GG[0][ids[0]][b] & GG[1][ids[0]][b] &
               ~GG[0][ids[1]][b] & GG[1][ids[1]][b] &
               ~GG[0][ids[2]][b] & GG[1][ids[2]][b];
         homhom = (GG[0][ids[2]][b] | GG[1][ids[2]][b]) & // offspring not missing
            (GG[0][ids[0]][b] & GG[0][ids[1]][b]);
         if(homhom & shortbase[j]){
            if(GG[1][ids[0]][b] & GG[1][ids[1]][b] & GG[1][ids[2]][b] & shortbase[j])
               // AA * AA -> AA
               for(int k = 0; k < 3; k ++)
                  for(int h = 0; h < 2; h++)
                     hap[h][ids[k]][l] = 1;
            else if((~GG[1][ids[0]][b]) & (~GG[1][ids[1]][b]) & (~GG[1][ids[2]][b]) & shortbase[j])
               // aa * aa -> aa
               for(int k = 0; k < 3; k ++)
                  for(int h = 0; h < 2; h++)
                     hap[h][ids[k]][l] = 0;
            else if(GG[1][ids[0]][b] & (~GG[1][ids[1]][b]) & GG[1][ids[2]][b] & shortbase[j])
               // AA * aa -> Aa
               for(int h = 0; h < 2; h++){
                  hap[h][ids[0]][l] = 1;
                  hap[h][ids[1]][l] = 0;
                  hap[h][ids[2]][l] = 1-h;
               }
            else if((~GG[1][ids[0]][b]) & GG[1][ids[1]][b] & GG[1][ids[2]][b] & shortbase[j])
               // aa * AA -> aA
               for(int h = 0; h < 2; h++){
                  hap[h][ids[0]][l] = 0;
                  hap[h][ids[1]][l] = 1;
                  hap[h][ids[2]][l] = h;
               }
            }
         else if(HHH & shortbase[j]){
            // Aa * Aa -> Aa
            for(int k = 0; k < 3; k ++)
               for(int h = 0; h < 2; h++)
                  hap[h][ids[k]][l] = -9999;
         }else if(oneinformative & shortbase[j]){
            for(int k = 0; k < 2; k++){
               if((~GG[0][ids[k]][b]) & GG[1][ids[k]][b] & GG[0][ids[2]][b] & GG[1][ids[2]][b] & shortbase[j])
                  // Aa->AA: Aa * AA -> AA
                  for(int h = 0; h < 2; h++){
                        hap[h][ids[k]][l] = 1-h;
                        hap[h][ids[1-k]][l] = hap[h][ids[2]][l] = 1;
                  }
               else if((~GG[0][ids[k]][b]) & GG[1][ids[k]][b] & GG[0][ids[2]][b] & (~GG[1][ids[2]][b]) & shortbase[j])
                  // Aa->aa: aA * aa -> aa
                  for(int h = 0; h < 2; h++){
                        hap[h][ids[k]][l] = h;
                        hap[h][ids[1-k]][l] = hap[h][ids[2]][l] = 0;
                  }
               else if(GG[0][ids[k]][b] & (~GG[1][ids[k]][b]) & (~GG[0][ids[2]][b]) & GG[1][ids[2]][b] & shortbase[j])
                  // aa->Aa: aa * Aa -> aA
                  for(int h = 0; h < 2; h++){
                        hap[h][ids[k]][l] = 0;
                        hap[h][ids[1-k]][l] = 1-h;
                        hap[h][ids[2]][l] = (k==0? h: 1-h);
                  }
               else if(GG[0][ids[k]][b] & (~GG[1][ids[k]][b]) & (~GG[0][ids[2]][b]) & GG[1][ids[2]][b] & shortbase[j])
                  // AA->Aa: AA * aA -> Aa
                  for(int h = 0; h < 2; h++){
                        hap[h][ids[k]][l] = 1;
                        hap[h][ids[1-k]][l] = h;
                        hap[h][ids[2]][l] = (k==0?1-h:h);
                  }
            } // end of parent
         }else{ // missing or MI error

         }
      }// end of marker
      if(IsPhased[ids[0]] || IsPhased[ids[1]] || IsPhased[ids[2]]){
         // compare the difference
         // fill in missing haplotypes
         // remove haplotypes if NOT consistent

      }

      for(int k = 0; k < 3; k++){
         IsPhased[ids[k]] = 0;
         for(int l = 0; l < rmarkers.Length(); l++)
            if(hap[0][ids[k]][l]>=0){
               IsPhased[ids[k]] = 1;
               break;
            }
      }
   }  // end of trio

   if(SaveFlag){
      String hapdatfile = prefix;
      hapdatfile.Add("hap.dat");
      FILE *fp = fopen((const char*)hapdatfile, "wt");
      if(fp==NULL)
         error("File %s cannot open to write", (const char*)hapdatfile);
      for(int l = 0; l < rmarkers.Length(); l++)
         fprintf(fp, "M %s\n", (const char*)snpName[rmarkers[l]]);

      String happedfile = prefix;
      happedfile.Add("hap.ped");
      fp = fopen((const char*)happedfile, "wt");
      if(fp==NULL)
         error("File %s cannot open to write", (const char*)happedfile);
      int count = 0;
      int missingcount = 0;
      for(int f = 0; f < ped.familyCount; f++){
         for(int k = 0; k < id[f].Length(); k++){
            int i = geno[id[f][k]];
            if(i > -1 && IsPhased[i]){
               fprintf(fp, "%s %s %s %s %d",
                     (const char*)ped[id[f][k]].famid,
                     (const char*)ped[id[f][k]].pid,
                     (const char*)ped[id[f][k]].fatid,
                     (const char*)ped[id[f][k]].motid,
                     ped[id[f][k]].sex);
               for(int l = 0; l < rmarkers.Length(); l++)
                  if(hap[0][i][l]>=0)
                     for(int h = 0; h < 2; h++)
                        fprintf(fp, " %c", labelInHap[1-hap[h][i][l]][l]);
                  else{
                     fprintf(fp, " X X");
                     missingcount++;
                  }
               fprintf(fp,"\n");
               count++;
            }
         }
      }
      printf("%d individuals are phased and the haplotypes are saved in %s and %s.\n",
         count, (const char*)hapdatfile, (const char*)happedfile);
      printf("Missing rate in the phased haplotype data is %d/%d=%.2lf%%\n",
         missingcount, count*rmarkers.Length(), missingcount*100.0/count/rmarkers.Length());
   }
}

