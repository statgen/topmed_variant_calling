// rareTDT.cpp
// 5/27/2011 Wei-Min Chen

#include "analysis.h"
#include "MathStats.h"
#include "Kinship.h"
#include "QuickIndex.h"

void Engine::rareTDT()
{
   int AA, Aa, missing;
   int b, id0, id1, id2, id3;
   int informative, T, NT, HHH;
   int TCount, NTCount, informativeCount;

   if(geno.Length()==0) BuildShortBinary();

   IntArray rmarkers(0);
   int AACount1, AaCount1, missingCount1;
   double frequency;
   for(int m = 0; m < markerCount; m++){
      if(chrList.Length() && chrList.Find(chromosomes[m])==-1) continue;
      if(start != _NAN_ && bp[m] < start) continue;
      if(stop != _NAN_ && bp[m] > stop) continue;

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
         if(frequency < rareMAF || frequency > 1-rareMAF)
            rmarkers.Push(m);
         if(frequency > 1-rareMAF) // flip
            for(int i = 0; i < idCount; i++)
               GG[1][i][b] ^= (GG[0][i][b] & shortbase[j]);
      }
   }

   printf("%d markers are collapsed in TDT analysis\n", rmarkers.Length());

   MakeTrioForTDT();
   if(Ltrio.Length()==0){
      warning("TDT analysis requires parent-affected-offspring trios.\n");
      return;
   }


   String pedfile=prefix;
   pedfile.Add("TDT.info");
   FILE *fp = fopen(pedfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
   fprintf(fp, "FAMID\tIID\tFA\tMO");
   for(int m = 0; m < rmarkers.Length(); m++){
      fprintf(fp, "\t%s",  (const char*)snpName[rmarkers[m]]);
   }

   IntArray allT(rmarkers.Length());
   IntArray allNT(rmarkers.Length());
   int tTCount, tNTCount, tcounts, ntcounts;
   tcounts = ntcounts = informativeCount = 0;
   for(int i = 0; i < Ltrio.Length()/3; i++){
      TCount = NTCount = 0;
      allT.Zero();
      allNT.Zero();
      id1 = geno[Ltrio[i*3+1]];
      id2 = geno[Ltrio[i*3+2]];
      id3 = geno[Ltrio[i*3]];   // id3 is the child
      for(int l = 0; l < rmarkers.Length(); l++){
         tTCount = tNTCount = 0;
         int m = rmarkers[l];
         b = m/16;
         int j = m % 16;                           
         informative = (GG[0][id3][b] | GG[1][id3][b]) & // offspring not missing
            ( (~GG[0][id1][b] & GG[1][id1][b] & (GG[0][id2][b] | GG[1][id2][b])) |
              (~GG[0][id2][b] & GG[1][id2][b] & (GG[0][id1][b] | GG[1][id1][b])) );
         if(!informative) continue; // no informative trios at any SNPs
         HHH = ~GG[0][id1][b] & GG[1][id1][b] &
               ~GG[0][id2][b] & GG[1][id2][b] &
               ~GG[0][id3][b] & GG[1][id3][b];
         if(HHH)
               if(HHH & shortbase[j]){
                  TCount++;
                  NTCount++;
                  informativeCount++;
               }
         for(int k = 0; k < 2; k++){
            id0 = geno[Ltrio[i*3+1+k]];
            T = informative &   // Aa->AA or aa->Aa
              ( (~GG[0][id0][b] & GG[1][id0][b] & GG[0][id3][b] & GG[1][id3][b]) |
                (GG[0][id0][b] & ~GG[1][id0][b] & ~GG[0][id3][b] & GG[1][id3][b]) );
            NT = informative &   // AA->Aa or Aa->aa
              ( (~GG[0][id3][b] & GG[1][id3][b] & GG[0][id0][b] & GG[1][id0][b]) |
                (GG[0][id3][b] & ~GG[1][id3][b] & ~GG[0][id0][b] & GG[1][id0][b]) );
            if(T & shortbase[j])
               tTCount++;
            if(NT & shortbase[j])
               tNTCount++;
         }
         // if TCount and NTCount consistent
         if((tTCount || tNTCount) && TCount!=-1 && NTCount!=-1){
            if(TCount == 0 && NTCount == 0){
               TCount = tTCount; NTCount = tNTCount;
            }else{
               if(tTCount != TCount || tNTCount != NTCount){
                  TCount = NTCount = -1;
                  break;
               }
            }
            allT[l] = TCount;
            allNT[l] = NTCount;
         }
      }// end of marker
      if(TCount == -1){ // uninformative

      }else if(TCount || NTCount){ // informative
         fprintf(fp, "\n%s\t%s\t%s\t%s",
            (const char*)ped[Ltrio[i*3]].famid,
            (const char*)ped[Ltrio[i*3]].pid,
            (const char*)ped[Ltrio[i*3+1]].pid,
            (const char*)ped[Ltrio[i*3+2]].pid);
         for(int m = 0; m < rmarkers.Length(); m++)
            fprintf(fp, "\t%d", allT[m] - allNT[m]);
         informativeCount++;
         tcounts += TCount;
         ntcounts += NTCount;

      }
   }  // end of trio
   fprintf(fp, "\n");
   fclose(fp);

   double den = tcounts + ntcounts;
   double statistic, temp;
   if(den > 0.5){
      temp = tcounts - ntcounts;
      statistic = temp*temp/den;
   }else
      statistic = 0;
   double p = chidist(statistic, 1);
   printf("T=%d NT=%d N=%d chisq=%.3lf P=%G\n",
      tcounts, ntcounts, informativeCount, statistic, p);

   printf("Family information saved in file %s\n", (const char*)pedfile);
}

/*
   unsigned long long *G0 = new unsigned long long [idCount];
   unsigned long long *G1 = new unsigned long long [idCount];
   for(int i = 0; i < idCount; i++)
      G0[i] = G1[i] = 0;
   int MAXSIZE = sizeof(unsigned long long)*8;
*/

/*
         if(frequency < rareMAF){ //
            // assign genotypes to G0[i], G1[i]
            // frequencies.push(frequency)
            int change = rmarkers.Length()-1-j;
            if(change > 0)
               for(int i = 0; i < idCount; i++){
                  G0[i] |= (GG[0][i][b] & shortbase[j]) << change;
                  G1[i] |= (GG[1][i][b] & shortbase[j]) << change;
               }
            else if(change < 0)
               for(int i = 0; i < idCount; i++){
                  G0[i] |= (GG[0][i][b] & shortbase[j]) >> (-change);
                  G1[i] |= (GG[1][i][b] & shortbase[j]) >> (-change);
               }
            else
               for(int i = 0; i < idCount; i++){
                  G0[i] |= (GG[0][i][b] & shortbase[j]);
                  G1[i] |= (GG[1][i][b] & shortbase[j]);
               }
         }else if(frequency > 1-rareMAF){ // flip
            // assign flipped genotypes to G0[i], G1[i]
            if(rmarkers.Length()==MAXSIZE) error("Too many rare variants: > %d", MAXSIZE);
            rmarkers.Push(m);

            int change = rmarkers.Length()-1-j;
            if(change > 0)
               for(int i = 0; i < idCount; i++){
                  G0[i] |= (GG[0][i][b] & shortbase[j]) << change;
                  G1[i] |= ((GG[0][i][b] ^ GG[1][i][b]) & shortbase[j]) << change;
               }
            else if(change < 0)
               for(int i = 0; i < idCount; i++){
                  G0[i] |= (GG[0][i][b] & shortbase[j]) >> (-change);
                  G1[i] |= ((GG[0][i][b] ^ GG[1][i][b]) & shortbase[j]) >> (-change);
               }
            else
               for(int i = 0; i < idCount; i++){
                  G0[i] |= (GG[0][i][b] & shortbase[j]);
                  G1[i] |= ((GG[0][i][b] ^ GG[1][i][b]) & shortbase[j]);
               }
            // frequencies.push(1-frequency)
         }      */

