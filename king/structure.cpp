////////////////////////////////////////////////////////////////////// 
// structure.cpp
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
// August 7, 2018

#include "analysis.h"
#include <math.h>
#include "Kinship.h"
#include "KinshipX.h"
#include "MathStats.h"
#include "MathSVD.h"
#include "QuickIndex.h"
#include "MathCholesky.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

#ifdef WITH_LAPACK
extern "C" void dgesdd_(char*, int*, int*, double*, int*, double*, double *, int*, double*, int*, double*, int*, int*, int*);
extern "C" void dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);
// SUBROUTINE DGESVD( JOBU, JOBVT, M, N,	A, LDA,	S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
#endif

void Engine::mds()
{
   if(xflag){
      if(xsnpName.Length() < MINSNPCOUNT) {
         if(xsnpName.Length()==0)
            printf("No X-chromosome SNPs were included in the data.\n");
         else
            printf("There are < %d SNPs in the data.\n", MINSNPCOUNT);
         return;
      }
      printf("X-chromosome MDS starts at %s", currentTime());
   }else
      printf("MDS starts at %s", currentTime());
   if(xflag)
      printf("X-chromosome genotypes stored in %d words for each of %d individuals.\n",
         xshortCount, idCount);
   else
      printf("Genotypes stored in %d words for each of %d individuals.\n",
         Bit64==64? longCount: shortCount, idCount);
   int id1, id2, count;
   IntArray ID(0);
   if(xflag)
      for(int i = 0; i < idCount; i++){
         count = 0;
         for(int m = 0; m < xshortCount; m++){
            int k = XG[0][i][m] | XG[1][i][m];
            count += oneCount[k&255] + oneCount[(k>>8)&255];
         }
         if(count >= MINSNPCOUNT)
            ID.Push(i);
      }
   else{
      if(Bit64==64)
         for(int i = 0; i < idCount; i++){
            count = 0;
            for(int m = 0; m < longCount; m++)
               count += popcount(LG[0][i][m] | LG[1][i][m]);
            if(count >= MINSNPCOUNT) ID.Push(i);
         }
      else
      for(int i = 0; i < idCount; i++){
         count = 0;
         for(int m = 0; m < shortCount; m++){
            int k = GG[0][i][m] | GG[1][i][m];
            count += oneCount[k&255] + oneCount[(k>>8)&255];
         }
         if(count >= MINSNPCOUNT)
            ID.Push(i);
      }
   }

   int dimN = ID.Length();
   Matrix D;
   if(xflag)
      error("To be implemented");
   if(Bit64==64)
      ComputeDistanceMatrix64Bit(D); // computationally intensive here
   else
      ComputeDistanceMatrix(D);// computationally intensive here

   Vector tempV(dimN);
   tempV.Zero();
   for(int j = 0; j < dimN; j++)
      for(int i = 0; i < dimN; i++)
         tempV[j] += D[i][j];
   tempV.Multiply(1.0/dimN);
   // (I-11'/N) * D
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimN; j++)
         D[i][j] -= tempV[j];
   tempV.Zero();
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimN; j++)
         tempV[i] += D[i][j];
   tempV.Multiply(1.0/dimN);
   // D * (I-11'/N)
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimN; j++)
         D[i][j] -= tempV[i];
   D.Multiply(-0.5);
   printf("SVD starts at %s", currentTime());

#ifdef WITH_LAPACK
   printf("  LAPACK is used.\n");
   char JOBZ = 'A';
   int info;
   double *A = new double[dimN*dimN];
   int dimLA = dimN;
   for(int i = 0; i < dimN; i++)
     for(int j = 0; j < dimN; j++)
       A[i*dimN+j] = D[j][i];
   double *S = new double[dimN];
   double *U = new double[dimN*dimN];
   double *VT = new double[dimN*dimN];
   int *IWORK = new int[dimN*8];
   int LWORK = 8*dimN + 4*dimN*dimN;
   double *WORK = new double[LWORK];
   dgesdd_(&JOBZ, &dimLA, &dimLA, A, &dimLA, S, U, &dimLA, VT, &dimLA, WORK, &LWORK, IWORK, &info);
// SUBROUTINE DGESDD(JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO)

   delete []U;
   delete []IWORK;
   delete []WORK;
   delete []A;

   int dimPC = dimN>20?20:dimN;
   printf("Largest %d eigenvalues:", dimPC);
   for(int i = 0; i < dimPC; i++)
      printf(" %.2lf", S[i]);
   printf("\n");
   double totalVariance = 0;
   double variance20 = 0;
   for(int i = 0; i < dimPC; i++)
      variance20 += S[i];
   for(int i = 0; (S[i] > 1E-10) && (i < dimN-1); i++)
      totalVariance += S[i];
   printf("The first %d PCs are able to explain %.2lf / %.2lf = %.1lf%% of total variance.\n",
      dimPC, variance20, totalVariance, variance20/totalVariance*100);
   printf("The proportion of total variance explained (%%) by each PC is:\n  ");
   for(int i = 0; i < dimPC; i++)
      printf(" %.1lf", S[i]/totalVariance*100);
   printf("\n");
   delete S;
#else
   printf("  Please re-compile KING with LAPACK library.\n");
   SVD svd;
   svd.Decompose(D);
   printf("done\n");
   if(svd.n == 0) return;
   QuickIndex idx;
   idx.Index(svd.w);

   int dimPC = dimN>20?20:dimN;
   printf("Largest %d eigenvalues:", dimPC);
   for(int i = 0; i < dimPC; i++)
      printf(" %.2lf", svd.w[idx[dimN-1-i]]);
   printf("\n");
#endif

   String datfile = prefix;
   datfile.Add("pc.dat");
   FILE *fp = fopen(datfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)datfile);
   fprintf(fp, "S related\n");
   for(int i = 0; i < dimPC; i++)
      fprintf(fp, "C PC%d\n", i+1);
   if(detailFlag){
      fprintf(fp, "C Heterozygosity\n");
      fprintf(fp, "C Quality\n");
   }
   fclose(fp);

   String pedfile = prefix;
   pedfile.Add("pc.ped");
   fp = fopen(pedfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
   for(int i = 0; i < ID.Length(); i++){
      int id = ID[i];
      fprintf(fp, "%s %s %s %s %d",
         (const char*)ped[phenoid[id]].famid, (const char*)ped[phenoid[id]].pid,
         (const char*)ped[phenoid[id]].fatid, (const char*)ped[phenoid[id]].motid,
         ped[phenoid[id]].sex);
         fprintf(fp, " 1");
         for(int j = 0; j < dimPC; j++)
#ifdef WITH_LAPACK
            fprintf(fp, " %.4lf", VT[i*dimN+j]);
#else
            fprintf(fp, " %.4lf", svd.v[i][idx[dimN-1-j]]);
#endif
      fprintf(fp, "\n");
   }
   fclose(fp);
   printf("MDS ends at %s", currentTime());
   printf("%d principal components saved in files %s and %s\n",
      dimPC, (const char*)datfile, (const char*)pedfile);
#ifdef WITH_LAPACK
   delete []VT;
#endif
}


/*
   if(allflags & (1<<IBSFLAG))
      printf("IBS Distance is used in the MDS analysis.\n");
   else
      printf("Euclidean Distance is used in the MDS analysis.\n");
*/


/*   if(detailFlag) countGenotype();
   if(geno.Length()==0) {
      individualInfo = true;
      if(shortFlag) BuildShortBinary();
      else BuildBinary();
   }
*/

   /*
   if(xflag)
   for(int n = 0; n < ped.count; n++){
      id1 = geno[n];
      count = 0;
      for(int m = 0; m < xshortCount; m++)
         count += oneoneCount[XG[0][id1][m] | XG[1][id1][m]];
      if(count >= MINSNPCOUNT)
         ID.Push(n);
   }
   else
   for(int n = 0; n < ped.count; n++)
      if(ped[n].ngeno >= MINSNPCOUNT)
         ID.Push(n);
     */

   /*
   int IBS0Count, IBS2Count, het1Count, notMissingCount;
   Matrix D(dimN, dimN);
   D.Zero();

  IntArray ompindex(dimN-1);
#ifdef _OPENMP
   for(int i = 0; i < (dimN-1)/2; i++){
      ompindex[2*i] = i;
      ompindex[2*i+1] = dimN-2-i;
   }
   if(dimN%2==0)
      ompindex[dimN-2] = (dimN-1)/2;
#else
   for(int i = 0; i < dimN-1; i++)
      ompindex[i] = i;
#endif
#ifdef _OPENMP
   printf("%d CPU cores are used to compute the distant matrix.\n", defaultMaxCoreCount);
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
         private(IBS0Count, IBS2Count, het1Count, notMissingCount, id1, id2)
#endif
   for(int tempi = 0; tempi < dimN-1; tempi++){
      int i = ompindex[tempi];
      id1 = geno[ID[i]];
      for(int j = i+1; j < dimN; j++){
         id2 = geno[ID[j]];
         if(id1 < 0 || id2 < 0) continue;
         if(ibsFlag){ // IBS Distance
            IBS0Count = IBS2Count = notMissingCount = 0;
            if(xflag)
               for(int m = 0; m < xshortCount; m++){
                  IBS0Count += oneoneCount[XG[0][id1][m] & XG[0][id2][m] & (XG[1][id1][m] ^ XG[1][id2][m])];
                  IBS2Count += oneoneCount[~(XG[0][id1][m]^XG[0][id2][m]) & ~(XG[1][id1][m]^XG[1][id2][m]) & (XG[0][id1][m] | XG[1][id1][m])];
                  notMissingCount += oneoneCount[(XG[0][id1][m] | XG[1][id1][m]) & (XG[0][id2][m] | XG[1][id2][m])];
               }
            else
               for(int m = 0; m < shortCount; m++){
                  IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
                  IBS2Count += oneoneCount[~(GG[0][id1][m]^GG[0][id2][m]) & ~(GG[1][id1][m]^GG[1][id2][m]) & (GG[0][id1][m] | GG[1][id1][m])];
                  notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
               }
            if(notMissingCount)
               D[i][j] = D[j][i] = 1-(IBS2Count-IBS0Count)*1.0/notMissingCount;
         }else{ // Euclidean Distance
            IBS0Count = het1Count = notMissingCount = 0;
            if(xflag)
               for(int m = 0; m < xshortCount; m++){
                  IBS0Count += oneoneCount[XG[0][id1][m] & XG[0][id2][m] & (XG[1][id1][m] ^ XG[1][id2][m])];
                  notMissingCount += oneoneCount[(XG[0][id1][m] | XG[1][id1][m]) & (XG[0][id2][m] | XG[1][id2][m])];
                  het1Count += oneoneCount[(XG[0][id2][m] & (~XG[0][id1][m]) & XG[1][id1][m]) |
                        (XG[0][id1][m] & (~XG[0][id2][m]) & XG[1][id2][m])];
               }
            else
               for(int m = 0; m < shortCount; m++){
                  IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
                  notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
                  het1Count += oneoneCount[(GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m]) |
                     (GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]) ];
               }
            if(notMissingCount)
               D[i][j] = D[j][i] = (het1Count + 4.0*IBS0Count) / notMissingCount;
         }
      }
   }
   */

void Engine::SemifamilyKinship()
{
   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   Kinship kin;
   int id1, id2;
   double temp;
   int HetHetCount, IBS0Count, het1Count;
   if(pedKin == NULL)
      pedKin = new Matrix[ped.familyCount];
   if(Bit64 == 64)
   for(int f = 0; f < ped.familyCount; f++){
      pedKin[f].Dimension(id[f].Length(), id[f].Length());
      pedKin[f].Set(0.5);
      if(id[f].Length() < 2) continue;
      kin.Setup(*ped.families[f]);
      for(int i = 0; i < id[f].Length(); i++)
         for(int j = i+1; j < id[f].Length(); j++){
            id1 = geno[id[f][i]]; id2 = geno[id[f][j]];
            pedKin[f][i][j] = pedKin[f][j][i] = kin(ped[id[f][i]], ped[id[f][j]]);
            if(!semifamilyFlag) continue;
            if(pedKin[f][i][j] > 0.005) continue;
            HetHetCount = IBS0Count = het1Count = 0;
            for(int m = 0; m < longCount; m++){
               HetHetCount += popcount((~LG[0][id1][m]) & (LG[1][id1][m]) & (~LG[0][id2][m]) & LG[1][id2][m]);
               het1Count += (popcount((LG[0][id2][m] | LG[1][id2][m]) & (~LG[0][id1][m]) & LG[1][id1][m])
                  + popcount((LG[0][id1][m] | LG[1][id1][m]) & (~LG[0][id2][m]) & LG[1][id2][m]));
               IBS0Count += popcount(LG[0][id1][m] & LG[0][id2][m] & (LG[1][id1][m] ^ LG[1][id2][m]));
            }
            if(het1Count){
               temp = (HetHetCount - IBS0Count*2)*1.0/het1Count;
               if(temp > 0.0221)
                  pedKin[f][i][j] = pedKin[f][j][i] = temp;
            }
         }
      }
   else
   for(int f = 0; f < ped.familyCount; f++){
      pedKin[f].Dimension(id[f].Length(), id[f].Length());
      pedKin[f].Set(0.5);
      if(id[f].Length() < 2) continue;
      kin.Setup(*ped.families[f]);
      for(int i = 0; i < id[f].Length(); i++)
         for(int j = i+1; j < id[f].Length(); j++){
            id1 = geno[id[f][i]]; id2 = geno[id[f][j]];
            pedKin[f][i][j] = pedKin[f][j][i] = kin(ped[id[f][i]], ped[id[f][j]]);
            if(!semifamilyFlag) continue;
            if(pedKin[f][i][j] > 0.005) continue;
            HetHetCount = IBS0Count = het1Count = 0;
            if(xflag)
            for(int m = 0; m < xshortCount; m++){
               HetHetCount += oneoneCount[(~XG[0][id1][m]) & (XG[1][id1][m]) & (~XG[0][id2][m]) & XG[1][id2][m]];
               het1Count += (oneoneCount[(XG[0][id2][m] | XG[1][id2][m]) & (~XG[0][id1][m]) & XG[1][id1][m]]
                  + oneoneCount[(XG[0][id1][m] | XG[1][id1][m]) & (~XG[0][id2][m]) & XG[1][id2][m]]);
               IBS0Count += oneoneCount[XG[0][id1][m] & XG[0][id2][m] & (XG[1][id1][m] ^ XG[1][id2][m])];
            }
            else
            for(int m = 0; m < shortCount; m++){
               HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
               het1Count += (oneoneCount[(GG[0][id2][m] | GG[1][id2][m]) & (~GG[0][id1][m]) & GG[1][id1][m]]
                  + oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]]);
               IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
            }
            if(het1Count){
               temp = (HetHetCount - IBS0Count*2)*1.0/het1Count;
               if(temp > 0.0221)
                  pedKin[f][i][j] = pedKin[f][j][i] = temp;
            }
         }
      }
}

void Engine::mds_family()
{
//   if(projectFlag || ibsFlag || clusterFlag){

   if(projectFlag || (allflags&(1<<IBSFLAG)) || (allflags&(1<<ClusterFLAG)) ){
      if(projectFlag)
         printf("Only unrelated individuals are used in SVD and PCs of their relatives are derived.\n");
      mds_family_projection();
      return;
   }
   if(detailFlag) countGenotype();
   if(geno.Length()==0) {
      individualInfo = true;
      BuildShortBinary();
   }
   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   int IBS0Count, IBS1Count, notMissingCount;
   int id1, id2;
   double kinship;

   IntArray vid(idCount);
   vid.Set(-1);
   int vidCount = 0;
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = 0; i < id[f].Length(); i++)
         if(ped[id[f][i]].ngeno > 0)
            vid[geno[id[f][i]]] = vidCount++;

   int dimN = vidCount;
   Matrix D(dimN, dimN);
   D.Zero();
   printf("Euclidean Distance is used in the MDS analysis.\n");
   if(xflag){  // X-chromosome
      if(xsnpName.Length() < MINSNPCOUNT) return;
      printf("MDS for X-chromosome family data starts at %s", currentTime());
      printf("X-chromosome genotypes stored in %d words for each of %d individuals.\n",
         xshortCount, idCount);
      KinshipX kinx;
      for(int f = 0; f < ped.familyCount; f++){
         kinx.Setup(*ped.families[f]);
         for(int i = 0; i < id[f].Length(); i++){
            if(ped[id[f][i]].ngeno==0) continue;   // untyped
            for(int j = i+1; j < id[f].Length(); j++){
               if(ped[id[f][j]].ngeno==0) continue;   // untyped
               kinship = kinx(ped[id[f][i]], ped[id[f][j]]);
               if(kinship >= 0.4999) continue; // Distance should be 0
               id1 = geno[id[f][i]]; id2 = geno[id[f][j]];
               IBS0Count = IBS1Count = notMissingCount = 0;
               for(int m = 0; m < xshortCount; m++){
                  IBS0Count += oneoneCount[XG[0][id1][m] & XG[0][id2][m] & (XG[1][id1][m] ^ XG[1][id2][m])];
                  notMissingCount += oneoneCount[(XG[0][id1][m] | XG[1][id1][m]) & (XG[0][id2][m] | XG[1][id2][m])];
                  IBS1Count += oneoneCount[ (XG[0][id2][m] & (~XG[0][id1][m]) & XG[1][id1][m]) |
                  (XG[0][id1][m] & (~XG[0][id2][m]) & XG[1][id2][m]) ];
               }
               if(notMissingCount)
                  D[vid[id1]][vid[id2]] = D[vid[id2]][vid[id1]] =
               (IBS1Count + 4.0*IBS0Count) / notMissingCount / (1-2*kinship);
            }
         }
      }
      for(int f1 = 0; f1 < ped.familyCount; f1++)
      for(int i = 0; i < id[f1].Length(); i++){
         if(ped[id[f1][i]].ngeno==0) continue;   // untyped
         for(int f2 = f1+1; f2 < ped.familyCount; f2++)
         for(int j = 0; j < id[f2].Length(); j++){
            if(ped[id[f2][j]].ngeno==0) continue;   // untyped
            id1 = geno[id[f1][i]]; id2 = geno[id[f2][j]];
            IBS0Count = IBS1Count = notMissingCount = 0;
            for(int m = 0; m < xshortCount; m++){
               IBS0Count += oneoneCount[XG[0][id1][m] & XG[0][id2][m] & (XG[1][id1][m] ^ XG[1][id2][m])];
               notMissingCount += oneoneCount[(XG[0][id1][m] | XG[1][id1][m]) & (XG[0][id2][m] | XG[1][id2][m])];
               IBS1Count += oneoneCount[ (XG[0][id2][m] & (~XG[0][id1][m]) & XG[1][id1][m]) |
                  (XG[0][id1][m] & (~XG[0][id2][m]) & XG[1][id2][m]) ];
            }
            if(notMissingCount)
               D[vid[id1]][vid[id2]] = D[vid[id2]][vid[id1]] = (IBS1Count + 4.0*IBS0Count) / notMissingCount;
         }
      }
   }else{ // autosome
      printf("MDS for family data starts at %s", currentTime());
      printf("Genotypes stored in %d words for each of %d individuals.\n",
         shortCount, idCount);
      Kinship kin;
      for(int f = 0; f < ped.familyCount; f++){
         kin.Setup(*ped.families[f]);
         for(int i = 0; i < id[f].Length(); i++){
            if(ped[id[f][i]].ngeno==0) continue;   // untyped
            for(int j = i+1; j < id[f].Length(); j++){
               if(ped[id[f][j]].ngeno==0) continue;   // untyped
               kinship = kin(ped[id[f][i]], ped[id[f][j]]);
               if(kinship >= 0.4999) continue; // Distance should be 0
               id1 = geno[id[f][i]]; id2 = geno[id[f][j]];
               IBS0Count = IBS1Count = notMissingCount = 0;
               for(int m = 0; m < shortCount; m++){
                  IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
                  notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
                  IBS1Count += oneoneCount[ (GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m]) |
                     (GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]) ];
               }
               if(notMissingCount)
                  D[vid[id1]][vid[id2]] = D[vid[id2]][vid[id1]] =
                  (IBS1Count + 4.0*IBS0Count) / notMissingCount / (1-2*kinship);
            }
         }
      }
  IntArray ompindex(idCount-1);
#ifdef _OPENMP
   for(int i = 0; i < (idCount-1)/2; i++){
      ompindex[2*i] = i;
      ompindex[2*i+1] = idCount-2-i;
   }
   if(idCount%2==0)
      ompindex[idCount-2] = (idCount-1)/2;
#else
   for(int i = 0; i < idCount-1; i++)
      ompindex[i] = i;
#endif
#ifdef _OPENMP
   printf("%d CPU cores are used to compute the distant matrix.\n", defaultMaxCoreCount);
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
         private(IBS0Count, IBS1Count, notMissingCount, id1, id2)
#endif
   for(int i = 0; i < idCount-1; i++){
      id1 = ompindex[i];
      if(ped[phenoid[id1]].ngeno==0) continue;
      for(id2 = id1 + 1; id2 < idCount; id2++){
         if(ped[phenoid[id2]].ngeno==0) continue;
         if(ped[phenoid[id1]].famid == ped[phenoid[id2]].famid) continue;
            IBS0Count = IBS1Count = notMissingCount = 0;
            for(int m = 0; m < shortCount; m++){
               IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
               notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
               IBS1Count += oneoneCount[ (GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m]) |
                  (GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]) ];
            }
            if(notMissingCount)
               D[vid[id1]][vid[id2]] = D[vid[id2]][vid[id1]] = (IBS1Count + 4.0*IBS0Count) / notMissingCount;
         }
      }
   }

   Vector tempV(dimN);
   tempV.Zero();
   for(int j = 0; j < dimN; j++)
      for(int i = 0; i < dimN; i++)
         tempV[j] += D[i][j];
   tempV.Multiply(1.0/dimN);
   // (I-11'/N) * D
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimN; j++)
         D[i][j] -= tempV[j];
   tempV.Zero();
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimN; j++)
         tempV[i] += D[i][j];
   tempV.Multiply(1.0/dimN);
   // D * (I-11'/N)
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimN; j++)
         D[i][j] -= tempV[i];
   D.Multiply(-0.5);

   printf("SVD of a %d x %d matrix starts at %s", dimN, dimN, currentTime());

#ifdef WITH_LAPACK
   printf("  LAPACK is used.\n");
   char JOBZ = 'A';
   int info;
   double *A = new double[dimN*dimN];
   int dimLA = dimN;
   for(int i = 0; i < dimN; i++)
     for(int j = 0; j < dimN; j++)
       A[i*dimN+j] = D[j][i];
   double *S = new double[dimN];
   double *U = new double[dimN*dimN];
   double *VT = new double[dimN*dimN];
   int *IWORK = new int[dimN*8];
   int LWORK = 8*dimN + 4*dimN*dimN;
   double *WORK = new double[LWORK];
   dgesdd_(&JOBZ, &dimLA, &dimLA, A, &dimLA, S, U, &dimLA, VT, &dimLA, WORK, &LWORK, IWORK, &info);
// SUBROUTINE DGESDD(JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO)

   delete []U;
   delete []IWORK;
   delete []WORK;
   delete []A;

   int dimPC = dimN>20?20:dimN;
//   int dimPC = every==1? 20: every;
//   if(dimN < dimPC) dimPC = dimN;
   printf("Largest %d eigenvalues:", dimPC);
   for(int i = 0; i < dimPC; i++)
      printf(" %.2lf", S[i]);
   printf("\n");
   double totalVariance = 0;
   double variance20 = 0;
   for(int i = 0; i < dimPC; i++)
      variance20 += S[i];

   for(int i = 0; (S[i] > 1E-10) && (i < dimN-1); i++)
      totalVariance += S[i];
   printf("The first %d PCs are able to explain %.2lf / %.2lf = %.1lf%% of total variance.\n",
      dimPC, variance20, totalVariance, variance20/totalVariance*100);
   printf("The proportion of total variance explained (%%) by each PC is:\n  ");
   for(int i = 0; i < dimPC; i++)
      printf(" %.1lf", S[i]/totalVariance*100);
   printf("\n");
   delete []S;
#else
   printf("  Please re-compile KING with LAPACK library.\n");
   SVD svd;
   svd.Decompose(D);
   if(svd.n == 0) return;
   QuickIndex idx;
   idx.Index(svd.w);

   int dimPC = dimN>20?20:dimN;
   printf("Largest %d eigenvalues:", dimPC);
   for(int i = 0; i < dimPC; i++)
      printf(" %.2lf", svd.w[idx[dimN-1-i]]);
   printf("\n");
   double totalVariance = 0;
   double variance20 = 0;
   for(int i = 0; i < dimPC; i++)
      variance20 += svd.w[idx[dimN-1-i]];
   for(int i = 0; (svd.w[idx[dimN-1-i]] > 1E-10) && (i < dimN-1); i++)
      totalVariance += svd.w[idx[dimN-1-i]];
   printf("The first %d PCs are able to explain %.2lf / %.2lf = %.1lf%% of total variance.\n",
      dimPC, variance20, totalVariance, variance20/totalVariance*100);
   printf("The proportion of total variance explained (%%) by each PC is:\n  ");
   for(int i = 0; i < dimPC; i++)
      printf(" %.1lf", svd.w[idx[dimN-1-i]]/totalVariance*100);
   printf("\n");

#endif

   Matrix EV(dimN, dimPC);
   EV.Zero();
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = 0; i < id[f].Length(); i++){
         if(ped[id[f][i]].ngeno==0) continue;
         for(int j = 0; j < dimPC; j++){
#ifdef WITH_LAPACK
            EV[vid[geno[id[f][i]]]][j] = VT[vid[geno[id[f][i]]]*dimN+j];
#else
            EV[vid[geno[id[f][i]]]][j] = svd.v[vid[geno[id[f][i]]]][idx[dimN-1-j]];
#endif
         }
      }

   String datfile = prefix;
   if(xflag)
      datfile.Add("Xpc.dat");
   else
      datfile.Add("pc.dat");
   FILE *fp = fopen(datfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)datfile);
   fprintf(fp, "S related\n");
   for(int i = 0; i < dimPC; i++)
      fprintf(fp, "C PC%d\n", i+1);
   if(detailFlag){
      fprintf(fp, "C Heterozygosity\n");
      fprintf(fp, "C Quality\n");
   }
   fclose(fp);

   String pedfile = prefix;
   if(xflag)
      pedfile.Add("Xpc.ped");
   else
      pedfile.Add("pc.ped");
   fp = fopen(pedfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
   for(int f = 0; f < ped.familyCount; f++)
      for(int i = 0; i < id[f].Length(); i++){
         if(ped[id[f][i]].ngeno==0) continue;
         fprintf(fp, "%s %s %s %s %d",
         (const char*)ped[id[f][i]].famid, (const char*)ped[id[f][i]].pid,
         (const char*)ped[id[f][i]].fatid, (const char*)ped[id[f][i]].motid,
         ped[id[f][i]].sex);
         fprintf(fp, " 1"); // used in SVD
         int k = geno[id[f][i]];
         for(int j = 0; j < dimPC; j++)
            fprintf(fp, " %.4lf", EV[vid[k]][j]);
/*         if(detailFlag)
            fprintf(fp, " %.3lf %.3lf",
            pAaCount[k] * 1.0 / (pAACount[k] + pAaCount[k] + paaCount[k]),
            GetQuality(pAACount[k], pAaCount[k], paaCount[k]));
*/
         fprintf(fp, "\n");
   }
   fclose(fp);
   printf("MDS for family data ends at %s", currentTime());
   printf("%d principal components saved in files %s and %s\n",
      dimPC, (const char*)datfile, (const char*)pedfile);
#ifdef WITH_LAPACK
   delete []VT;
#endif
}

void Engine::mds_family_projection()
{
   if(!unrelatedExtraction){
      if(xflag){
         if(xsnpName.Length() < MINSNPCOUNT) return;
         printf("MDS for X-chromosome family data starts at %s", currentTime());
      }else
         printf("MDS for family data starts at %s", currentTime());
   }

   StringArray tokens;
   QuickIndex idx;

   if(allflags&(1<<ClusterFLAG))
      for(int f = 0; f < ped.familyCount; f++){
//         int k = ped[ped.families[f]->first].pid.FastFindChar(specialChar);
         int k = ped[ped.families[f]->first].pid.Find("->");
         if(ped.families[f]->famid.SubStr(0, 4) == "KING" && k > -1)
            for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
//               k = ped[i].pid.FastFindChar(specialChar); // a family could have diff family ID
               k = ped[i].pid.Find("->"); // a family could have diff family ID
               ped[i].famid = ped[i].pid.SubStr(0,k);
               ped[i].pid = ped[i].pid.SubStr(k+2);
               if(ped[i].fatid != "0")
                  ped[i].fatid = ped[i].fatid.SubStr(k+2);
               if(ped[i].motid != "0")
                  ped[i].motid = ped[i].motid.SubStr(k+2);
            }
      }
   if(detailFlag) countGenotype();
   if(geno.Length()==0) {
      individualInfo = true;
      if(shortFlag) BuildShortBinary();
      else BuildBinary();
   }

   if(!unrelatedExtraction){
      if(xflag)
         printf("X-chromosome genotypes stored in %d words for each of %d individuals.\n",
         xshortCount, idCount);
      else
         printf("Genotypes stored in %d words for each of %d individuals.\n",
         shortCount, idCount);
      if(allflags&(1<<IBSFLAG))
         printf("IBS Distance is used in the MDS analysis.\n");
      else
         printf("Euclidean Distance is used in the MDS analysis.\n");
   }
   IntArray first(ped.familyCount);
   IntArray last(ped.familyCount);
   first.Set(-1); last.Set(-1);
   IntArray ID(0), ID_UN(0), ID_unrelated;
   IntArray unrelatedCount;

   SemifamilyKinship();
   for(int f = 0; f < ped.familyCount; f++){
      first[f] = ID_UN.Length();
      if(id[f].Length() < 2) {
         if(id[f].Length() == 1){
            ID_UN.Push(id[f][0]);
            ID.Stack(id[f]);
            last[f] = first[f];
         }
         continue;
      }
      unrelatedCount.Dimension(id[f].Length());
      unrelatedCount.Zero();
      for(int i = 0; i < id[f].Length(); i++){
         if(ped[id[f][i]].ngeno == 0) { // untyped person not included in analysis
            unrelatedCount[i] = -1;
            continue;
         }
         for(int j = i+1; j < id[f].Length(); j++)
            if( (semifamilyFlag && pedKin[f][i][j] < 0.022) || (pedKin[f][i][j] < 0.001) ){
               unrelatedCount[i] ++;
               unrelatedCount[j] ++;
            }
      }
      idx.Index(unrelatedCount);
      ID_unrelated.Dimension(0);
      ID_unrelated.Push(idx[unrelatedCount.Length()-1]);
      first[f] = ID_UN.Length();
      for(int i = unrelatedCount.Length()-2; (i >= 0) && (unrelatedCount[idx[i]] > 0); i--){
         bool isUnrelated = true;
         for(int j = 0; j < ID_unrelated.Length(); j++)
            if( (semifamilyFlag && pedKin[f][idx[i]][ID_unrelated[j]] >= 0.022)
               || (!semifamilyFlag && pedKin[f][idx[i]][ID_unrelated[j]] >= 0.001) ) {
               isUnrelated = false;
               break;
            }
         if(isUnrelated) ID_unrelated.Push(idx[i]);
      }
      last[f] = ID_unrelated.Length()-1+ID_UN.Length();
      for(int i = 0; i < ID_unrelated.Length(); i++)
         ID_UN.Push(id[f][ID_unrelated[i]]);
      ID.Stack(id[f]);
   }

   if(unrelatedExtraction){ // extract a subset of unrelated individuals
      String file = prefix;
      file.Add("unrelated.txt");
      FILE *fp = fopen((const char*)file, "wt");
      if(fp==NULL) error("Cannot open %s to write", (const char*)file);
      for(int i = 0; i < ID_UN.Length(); i++)
         fprintf(fp, "%s\t%s\n", (const char*)ped[ID_UN[i]].famid, (const char*)ped[ID_UN[i]].pid);
      fclose(fp);
      printf("\nA list of %d unrelated individuals saved in file %s\n",
         ID_UN.Length(), (const char*)file);
      // generate a list of samples to be removed
      IntArray ID_remove=ID;
      for(int i = 0; i < ID_UN.Length(); i++)
         ID_remove.Delete(ID_remove.Find(ID_UN[i]));
      file = prefix;
      file.Add("unrelated_toberemoved.txt");
      fp = fopen((const char*)file, "wt");
      if(fp==NULL) error("Cannot open %s to write", (const char*)file);
      for(int i = 0; i < ID_remove.Length(); i++)
         fprintf(fp, "%s\t%s\n", (const char*)ped[ID_remove[i]].famid, (const char*)ped[ID_remove[i]].pid);
      fclose(fp);
      printf("An alternative list of %d to-be-removed individuals saved in file %s\n",
         ID_remove.Length(), (const char*)file);
      printf("\nExtracting a subset of unrelated individuals ends at %s", currentTime());
      return;
   }

   int dimN = ID_UN.Length();
   Matrix D(dimN, dimN);
   D.Zero();

   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   int IBS0Count, IBS2Count, het1Count, notMissingCount;
   int id1, id2;

  IntArray ompindex(dimN-1);
#ifdef _OPENMP
   for(int i = 0; i < (dimN-1)/2; i++){
      ompindex[2*i] = i;
      ompindex[2*i+1] = dimN-2-i;
   }
   if(dimN%2==0)
      ompindex[dimN-2] = (dimN-1)/2;
#else
   for(int i = 0; i < dimN-1; i++)
      ompindex[i] = i;
#endif
#ifdef _OPENMP
   printf("%d CPU cores are used to compute the distant matrix.\n", defaultMaxCoreCount);
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
         private(IBS0Count, IBS2Count, het1Count, notMissingCount, id1, id2)
#endif
   for(int tempi = 0; tempi < dimN-1; tempi++){
      int i = ompindex[tempi];
      id1 = geno[ID_UN[i]];
      for(int j = i+1; j < dimN; j++){
         id2 = geno[ID_UN[j]];
         if(id1 < 0 || id2 < 0) continue;
         if(allflags&(1<<IBSFLAG)){ // IBS Distance
            IBS0Count = IBS2Count = notMissingCount = 0;
            if(xflag)
               for(int m = 0; m < xshortCount; m++){
                  IBS0Count += oneoneCount[XG[0][id1][m] & XG[0][id2][m] & (XG[1][id1][m] ^ XG[1][id2][m])];
                  IBS2Count += oneoneCount[~(XG[0][id1][m]^XG[0][id2][m]) & ~(XG[1][id1][m]^XG[1][id2][m]) & (XG[0][id1][m] | XG[1][id1][m])];
                  notMissingCount += oneoneCount[(XG[0][id1][m] | XG[1][id1][m]) & (XG[0][id2][m] | XG[1][id2][m])];
               }
            else
               for(int m = 0; m < shortCount; m++){
                  IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
                  IBS2Count += oneoneCount[~(GG[0][id1][m]^GG[0][id2][m]) & ~(GG[1][id1][m]^GG[1][id2][m]) & (GG[0][id1][m] | GG[1][id1][m])];
                  notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
               }
            if(notMissingCount)
               D[i][j] = D[j][i] = 1-(IBS2Count-IBS0Count)*1.0/notMissingCount;
         }else{ // Euclidean Distance
            IBS0Count = het1Count = notMissingCount = 0;
            if(xflag)
               for(int m = 0; m < xshortCount; m++){
                  IBS0Count += oneoneCount[XG[0][id1][m] & XG[0][id2][m] & (XG[1][id1][m] ^ XG[1][id2][m])];
                  notMissingCount += oneoneCount[(XG[0][id1][m] | XG[1][id1][m]) & (XG[0][id2][m] | XG[1][id2][m])];
                  het1Count += oneoneCount[ (XG[0][id2][m] & (~XG[0][id1][m]) & XG[1][id1][m]) |
                     (XG[0][id1][m] & (~XG[0][id2][m]) & XG[1][id2][m]) ];
               }
            else
               for(int m = 0; m < shortCount; m++){
                  IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
                  notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
                  het1Count += oneoneCount[ (GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m]) |
                  (GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m]) ];
               }
            if(notMissingCount)
               D[i][j] = D[j][i] = (het1Count + 4.0*IBS0Count) / notMissingCount;
         }
      }
   }
   Vector tempV(dimN);
   tempV.Zero();
   for(int j = 0; j < dimN; j++)
      for(int i = 0; i < dimN; i++)
         tempV[j] += D[i][j];
   tempV.Multiply(1.0/dimN);
   // (I-11'/N) * D
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimN; j++)
         D[i][j] -= tempV[j];
   tempV.Zero();
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimN; j++)
         tempV[i] += D[i][j];
   tempV.Multiply(1.0/dimN);
   // D * (I-11'/N)
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimN; j++)
         D[i][j] -= tempV[i];
   D.Multiply(-0.5);

   printf("SVD of a %d x %d matrix starts at %s", dimN, dimN, currentTime());

#ifdef WITH_LAPACK
   printf("  LAPACK is used.\n");
   char JOBZ = 'A';
   int info;
   double *A = new double[dimN*dimN];
   int dimLA = dimN;
   for(int i = 0; i < dimN; i++)
     for(int j = 0; j < dimN; j++)
       A[i*dimN+j] = D[j][i];
   double *S = new double[dimN];
   double *U = new double[dimN*dimN];
   double *VT = new double[dimN*dimN];
   int *IWORK = new int[dimN*8];
   int LWORK = 8*dimN + 4*dimN*dimN;
   double *WORK = new double[LWORK];
   dgesdd_(&JOBZ, &dimLA, &dimLA, A, &dimLA, S, U, &dimLA, VT, &dimLA, WORK, &LWORK, IWORK, &info);
// SUBROUTINE DGESDD(JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO)

   delete []U;
   delete []IWORK;
   delete []WORK;
   delete []A;

   int dimPC = dimN>20?20:dimN;
   printf("Largest %d eigenvalues:", dimPC);
   for(int i = 0; i < dimPC; i++)
      printf(" %.2lf", S[i]);
   printf("\n");
   double totalVariance = 0;
   double variance20 = 0;
   for(int i = 0; i < dimPC; i++)
      variance20 += S[i];

   for(int i = 0; (S[i] > 1E-10) && (i < dimN-1); i++)
      totalVariance += S[i];
   printf("The first %d PCs are able to explain %.2lf / %.2lf = %.1lf%% of total variance.\n",
      dimPC, variance20, totalVariance, variance20/totalVariance*100);
   printf("The proportion of total variance explained (%%) by each PC is:\n  ");
   for(int i = 0; i < dimPC; i++)
      printf(" %.1lf", S[i]/totalVariance*100);
   printf("\n");
   delete []S;
#else
   printf("  Please re-compile KING with LAPACK library.\n");
   SVD svd;
   svd.Decompose(D);
   if(svd.n == 0) return;
   idx.Index(svd.w);

   int dimPC = dimN>20?20:dimN;
   printf("Largest %d eigenvalues:", dimPC);
   for(int i = 0; i < dimPC; i++)
      printf(" %.2lf", svd.w[idx[dimN-1-i]]);
   printf("\n");
   double totalVariance = 0;
   double variance20 = 0;
   for(int i = 0; i < dimPC; i++)
      variance20 += svd.w[idx[dimN-1-i]];
   for(int i = 0; (svd.w[idx[dimN-1-i]] > 1E-10) && (i < dimN-1); i++)
      totalVariance += svd.w[idx[dimN-1-i]];
   printf("The first %d PCs are able to explain %.2lf / %.2lf = %.1lf%% of total variance.\n",
      dimPC, variance20, totalVariance, variance20/totalVariance*100);
   printf("The proportion of total variance explained (%%) by each PC is:\n  ");
   for(int i = 0; i < dimPC; i++)
      printf(" %.1lf", svd.w[idx[dimN-1-i]]/totalVariance*100);
   printf("\n");

#endif

   Matrix EV(ID.Length(), dimPC);
   EV.Zero();
   ID.Dimension(0);

   IntArray ID_related;
   IntArray flagged;
   Vector kinship;

   flagged.Dimension(ped.count);
   flagged.Zero();
   double temp;
   for(int i = 0; i < ID_UN.Length(); i++)
      flagged[ID_UN[i]] = 1;
   IntArray tIdx;
   for(int f = 0; f < ped.familyCount; f++){
      if(last[f] < 0) continue;
      ID_related.Dimension(0);
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(ped[i].ngeno >= MINSNPCOUNT && flagged[i]==0)
            ID_related.Push(i);
      for(int i = first[f]; i <= last[f]; i++)
         for(int j = 0; j < dimPC; j++)
#ifdef WITH_LAPACK
            EV[geno[ID_UN[i]]][j] = VT[i*dimN+j];
#else
            EV[geno[ID_UN[i]]][j] = svd.v[i][idx[dimN-1-j]];
#endif
      tIdx.Dimension(ped.families[f]->count);
      tIdx.Set(-1);
      for(int i = 0; i < id[f].Length(); i++)
         tIdx[id[f][i]-ped.families[f]->first] = i;

      for(int i = 0; i < ID_related.Length(); i++){
         kinship.Dimension(last[f]-first[f]+1);
         kinship.Zero();
         for(int j = first[f]; j <= last[f]; j++)
            kinship[j-first[f]] = pedKin[f][tIdx[ID_related[i]-ped.families[f]->first]][tIdx[ID_UN[j]-ped.families[f]->first]]*2;
         for(int j = 0; j < dimPC; j++){
            temp = 0;
            for(int k = first[f]; k <= last[f]; k++){
#ifdef WITH_LAPACK
               EV[geno[ID_related[i]]][j] += kinship[k-first[f]] * VT[k*dimN+j];
#else
               EV[geno[ID_related[i]]][j] += kinship[k-first[f]] * svd.v[k][idx[dimN-1-j]];
#endif
               temp += kinship[k-first[f]];
            }
            if(temp > 1E-10)
               EV[geno[ID_related[i]]][j] /= temp;
         }
      }
   }

   String datfile = prefix;
   if(xflag)
      datfile.Add("Xpc.dat");
   else
      datfile.Add("pc.dat");
   FILE *fp = fopen(datfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)datfile);
   fprintf(fp, "S related\n");
   for(int i = 0; i < dimPC; i++)
      fprintf(fp, "C PC%d\n", i+1);
   if(detailFlag){
      fprintf(fp, "C Heterozygosity\n");
      fprintf(fp, "C Quality\n");
   }
   fclose(fp);

   String pedfile = prefix;
   if(xflag)
      pedfile.Add("Xpc.ped");
   else
      pedfile.Add("pc.ped");
   fp = fopen(pedfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
   IntArray inSVD(ped.count);
   inSVD.Set(-1);
   for(int i = 0; i < ID_UN.Length(); i++)
      inSVD[ID_UN[i]] = i;
   for(int i = 0; i < ped.count; i++){
      int k = geno[i];
/*      if(k == -1){
         for(int j = 0; j < dimPC+1; j++)
            fprintf(fp, " X");
         if(detailFlag)
            fprintf(fp, " X X");
      }else*/
      if(k == -1) continue;
      fprintf(fp, "%s %s %s %s %d",
            (const char*)ped[i].famid, (const char*)ped[i].pid,
            (const char*)ped[i].fatid, (const char*)ped[i].motid,
            ped[i].sex);
      if(inSVD[i]!=-1) fprintf(fp, " 1"); // unrelated in SVD
      else fprintf(fp, " 2"); // related in PCA
      for(int j = 0; j < dimPC; j++)
            fprintf(fp, " %.4lf", EV[k][j]);
/*      if(detailFlag)
            fprintf(fp, " %.3lf %.3lf",
             pAaCount[k] * 1.0 / (pAACount[k] + pAaCount[k] + paaCount[k]),
            GetQuality(pAACount[k], pAaCount[k], paaCount[k]));
*/
      fprintf(fp, "\n");
   }
   fclose(fp);
   printf("MDS for family data ends at %s", currentTime());
   printf("%d principal components saved in files %s and %s\n",
      dimPC, (const char*)datfile, (const char*)pedfile);
#ifdef WITH_LAPACK
   delete []VT;
#endif
}

void Engine::pca_projection(int every)
{
   printf("PCA starts at %s", currentTime());
   individualInfo = true;
   if(geno.Length()==0)
      BuildShortBinary();
   printf("Genotypes stored in %d words for each of %d individuals.\n",
         shortCount, idCount);

   Kinship kin;
   IntArray ID(0), ID_UN(0);
   for(int i = 0; i < ped.count; i++)
      if(ped[i].ngeno >= MINSNPCOUNT){
         if(ped[i].affections[0]!=2)
            ID_UN.Push(i);
         ID.Push(i);
      }

   int dimN = ID_UN.Length();
   if(dimN < 2) {
      printf("The number of unaffected individuals is < 2.\n");
      return;
   }
//   int dimM = (markerCount-1) / every +1;
   int dimM = markerCount;
   Matrix X(dimM, dimN);
   Matrix Z(dimM, ID.Length());
   X.Zero(); Z.Zero();
   int minorAllele;
   int validSNP=0;
   int byte, offset;
   double p, g, mean, se;
   int k, freqCount;
   bool flipFlag;

   for(int m = 0; m < markerCount; m ++){
      byte = m/16;
      offset = m%16;

      p = 0.0;
      freqCount = 0;
      for(int i = 0; i < dimN; i++){
         k = geno[ID_UN[i]];
         if(GG[1][k][byte] & shortbase[offset]){ // AA or Aa
            p ++;
            if(GG[0][k][byte] & shortbase[offset]) p ++; // AA
            freqCount += 2;
         }else if(GG[0][k][byte] & shortbase[offset]) // aa
            freqCount += 2;
      }
      if(freqCount == 0) continue;
      p /= freqCount;
      flipFlag = false;
      if(p > 0.5) {
         flipFlag = true;
         p = 1-p;
      }
      if(p < 0.001) continue;
      mean = p*2;
      se = sqrt(2*p*(1-p));
      if(se < 1E-100) continue;

      for(int i = 0; i < dimN; i++){
         k = geno[ID_UN[i]];
         if(k == -1) continue;
         g = -1.0;
         if(GG[0][k][byte] & shortbase[offset])  // homozygote
            g = (GG[1][k][byte] & shortbase[offset])? 2.0: 0.0;
         else if(GG[1][k][byte] & shortbase[offset])   // Aa
            g = 1.0;
         if(g > -0.9){
            if(flipFlag)
               X[validSNP][i] = (2 - g - mean) / se;
            else
               X[validSNP][i] = (g - mean) / se;
         }
      }
      for(int i = 0; i < ID.Length(); i++){
         k = geno[ID[i]];
         if(k == -1) continue;
         g = -1.0;
         if(GG[0][k][byte] & shortbase[offset])  // homozygote
            g = (GG[1][k][byte] & shortbase[offset])? 2.0: 0.0;
         else if(GG[1][k][byte] & shortbase[offset])   // Aa
            g = 1.0;
         if(g > -0.9) {
            if(flipFlag)
               Z[validSNP][i] = (2 - g - mean) / se;
            else
               Z[validSNP][i] = (g - mean) / se;
         }
      }
      validSNP++;
   }
   dimM = validSNP;
   X.Dimension(dimM, dimN);
   Z.Dimension(dimM, ID.Length());

   printf("Dimension of PCA: %d %d\n", dimM, dimN);
   printf("SVD...");

#ifdef WITH_LAPACK
   printf("  LAPACK is used.\n");
   char JOBU = 'S';
   char JOBVT = 'O';
   int LDU = dimM;
   int LDVT = 1;
   int dimS = dimM>dimN? dimN: dimM;
   int info;
   double *A = (double *)malloc(dimM*dimN*sizeof(double));
   for(int i = 0; i < dimM; i++)
     for(int j = 0; j < dimN; j++)
       A[j*dimM+i] = X[i][j];
   double *S = (double *)malloc(dimS*sizeof(double));
   double *U = (double *)malloc(LDU*dimS*sizeof(double));
   double *VT = (double *)malloc(LDVT*sizeof(double));
   int LWORK = dimM>dimN? dimM*2+8*dimN: dimN+8*dimM;
   double *WORK = (double *)malloc(LWORK*sizeof(double));
   dgesvd_(&JOBU, &JOBVT, &dimM, &dimN, A, &dimM, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &info);
// SUBROUTINE DGESVD( JOBU, JOBVT, M, N,	A, LDA,	S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
   if(info!=0) error("SVD failed.");
   free(WORK);
   free(A);
   free(VT);

//   int dimPC = dimN>20?20:dimN; // change to every later
   int dimPC = every==1? 20: every;
   if(dimN < dimPC) dimPC = dimN;
   if(dimM < dimPC) dimPC = dimM;
   printf("Largest %d eigenvalues:", dimPC);
   for(int i = 0; i < dimPC; i++)
      printf(" %.2lf", S[i]);
   printf("\n");

   Matrix EV(ID.Length(), dimPC);
   EV.Zero();
   double temp;
#ifdef _OPENMP
   printf("%d CPU cores are used to compute the projection.\n", defaultMaxCoreCount);
   #pragma omp parallel for num_threads(defaultMaxCoreCount) private(temp)
#endif
   for(int i = 0; i < ID.Length(); i++)
      for(int j = 0; j < dimPC; j++){
         temp = 0;
         for(int m = 0; m < dimM; m++)
            temp += Z.data[m]->data[i] * U[j*dimM+m];
         if(S[j] > 1E-100)
            temp /= S[j];
          EV.data[i]->data[j] = temp;
      }
   free(U);
   free(S);
#else
   printf("  Please re-compile KING with LAPACK library.\n");
   SVD svd;
   svd.Decompose(X);
   printf("done\n");
   if(svd.n == 0) return;
   QuickIndex idx;
   idx.Index(svd.w);

   int dimPC = dimN>20?20:dimN;
   if(dimM < dimPC) dimPC = dimM;
   printf("Largest %d eigenvalues:", dimPC);
   for(int i = 0; i < dimPC; i++)
      printf(" %.2lf", svd.w[idx[dimN-1-i]]);
   printf("\n");

   Matrix EV(ID.Length(), dimPC);
   EV.Zero();
   for(int i = 0; i < ID.Length(); i++)
      for(int j = 0; j < dimPC; j++){
         for(int m = 0; m < dimM; m++)
            EV[i][j] += Z[m][i] * svd.u[m][idx[dimN-1-j]];
         if(svd.w[idx[dimN-1-j]] > 1E-10)
            EV[i][j] /= svd.w[idx[dimN-1-j]];
      }
#endif

   String datfile = prefix;
   datfile.Add("pc.dat");
   FILE *fp = fopen(datfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)datfile);
   fprintf(fp, "S related\n");
   for(int i = 0; i < dimPC; i++)
      fprintf(fp, "C PC%d\n", i+1);
   fclose(fp);

   String pedfile = prefix;
   pedfile.Add("pc.ped");
   fp = fopen(pedfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
   IntArray typed(ped.count);
   typed.Set(-1);
   for(int i = 0; i < ID.Length(); i++)
      typed[ID[i]] = i;
   IntArray inSVD(ped.count);
   inSVD.Set(-1);
   for(int i = 0; i < ID_UN.Length(); i++)
      inSVD[ID_UN[i]] = i;
   for(int i = 0; i < ped.count; i++){
      if(typed[i] == -1) continue;
      /*
         for(int j = 0; j < dimPC+1; j++)
            fprintf(fp, " X");
      else{  */
      fprintf(fp, "%s %s %s %s %d",
         (const char*)ped[i].famid, (const char*)ped[i].pid,
         (const char*)ped[i].fatid, (const char*)ped[i].motid,
         ped[i].sex);

      if(inSVD[i]!=-1) fprintf(fp, " 1"); // unrelated in SVD
      else fprintf(fp, " 2"); // related in PCA
      for(int j = 0; j < dimPC; j++)
         fprintf(fp, " %.4lf", EV[typed[i]][j]);
      // }
      fprintf(fp, "\n");
   }
   fclose(fp);
   printf("PCA ends at %s", currentTime());
   printf("%d principal components saved in files %s and %s\n",
      dimPC, (const char*)datfile, (const char*)pedfile);
}

void Engine::pca_family(int every)
{
   printf("PCA starts at %s", currentTime());
   individualInfo = true;
   if(geno.Length()==0){
      if(shortFlag) BuildShortBinary();
      else BuildBinary();
   }
//   if(shortCount)
      printf("Genotypes stored in %d words for each of %d individuals.\n",
         shortCount, idCount);

   Kinship kin;
   IntArray ID(0), ID_UN(0), ID_unrelated;
   for(int f = 0; f < ped.familyCount; f++){
      ID_unrelated.Dimension(0);
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(ped[i].ngeno >= MINSNPCOUNT){
            ID.Push(i);
            if(ped[i].isFounder())
               ID_unrelated.Push(i);
         }
      kin.Setup(*ped.families[f]);
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
         if(ped[i].ngeno < MINSNPCOUNT) continue;
         bool isUnrelated = true;
         for(int j = 0; j < ID_unrelated.Length(); j++)
            if(kin(ped[i], ped[ID_unrelated[j]]) > 0) {
               isUnrelated = false;
               break;
            }
         if(isUnrelated) ID_unrelated.Push(i);
      }
      ID_UN.Stack(ID_unrelated);
   }
   int dimN = ID_UN.Length();
   if(dimN < 2) {
      printf("The number of unrelated individuals is < 2.\n");
      return;
   }
//   int dimM = (markerCount-1) / every +1;
   int dimM = markerCount;
#ifndef WITH_LAPACK
   if(dimM > 100000)
      printf("Note: it may take too long to perform PCA on %d markers.\n", dimM);
#endif
   Matrix X(dimM, dimN);
   Matrix Z(dimM, ID.Length());
   X.Zero(); Z.Zero();
   int minorAllele;
   int validSNP=0;
   int byte, offset;
   double p, g, mean, se;
   int k, freqCount;
   bool flipFlag;

   for(int m = 0; m < markerCount; m ++){
      byte = m/16;
      offset = m%16;

      p = 0.0;
      freqCount = 0;
      for(int i = 0; i < dimN; i++){
         k = geno[ID_UN[i]];
         if(GG[1][k][byte] & shortbase[offset]){ // AA or Aa
            p ++;
            if(GG[0][k][byte] & shortbase[offset]) p ++; // AA
            freqCount += 2;
         }else if(GG[0][k][byte] & shortbase[offset]) // aa
            freqCount += 2;
      }
      if(freqCount == 0) continue;
      p /= freqCount;
      flipFlag = false;
      if(p > 0.5) {
         flipFlag = true;
         p = 1-p;
      }
      if(p < 0.001) continue;
      mean = p*2;
      se = sqrt(2*p*(1-p));
      if(se < 1E-100) continue;

      for(int i = 0; i < dimN; i++){
         k = geno[ID_UN[i]];
         if(k == -1) continue;
         g = -1.0;
         if(GG[0][k][byte] & shortbase[offset])  // homozygote
            g = (GG[1][k][byte] & shortbase[offset])? 2.0: 0.0;
         else if(GG[1][k][byte] & shortbase[offset])   // Aa
            g = 1.0;
         if(g > -0.9){
            if(flipFlag)
               X[validSNP][i] = (2 - g - mean) / se;
            else
               X[validSNP][i] = (g - mean) / se;
         }
      }
      for(int i = 0; i < ID.Length(); i++){
         k = geno[ID[i]];
         if(k == -1) continue;
         g = -1.0;
         if(GG[0][k][byte] & shortbase[offset])  // homozygote
            g = (GG[1][k][byte] & shortbase[offset])? 2.0: 0.0;
         else if(GG[1][k][byte] & shortbase[offset])   // Aa
            g = 1.0;
         if(g > -0.9) {
            if(flipFlag)
               Z[validSNP][i] = (2 - g - mean) / se;
            else
               Z[validSNP][i] = (g - mean) / se;
         }
      }
      validSNP++;
   }
   dimM = validSNP;
   X.Dimension(dimM, dimN);
   Z.Dimension(dimM, ID.Length());

   printf("Dimension of PCA: %d %d\n", dimM, dimN);
   printf("SVD...");

#ifdef WITH_LAPACK
   printf("  LAPACK is used.\n");
   char JOBU = 'S';
   char JOBVT = 'O';
   int LDU = dimM;
   int LDVT = 1;
   int dimS = dimM>dimN? dimN: dimM;
   int info;
   double *A = (double *)malloc(dimM*dimN*sizeof(double));
   for(int i = 0; i < dimM; i++)
     for(int j = 0; j < dimN; j++)
       A[j*dimM+i] = X[i][j];
   double *S = (double *)malloc(dimS*sizeof(double));
   double *U = (double *)malloc(LDU*dimS*sizeof(double));
   double *VT = (double *)malloc(LDVT*sizeof(double));
   int LWORK = dimM>dimN? dimM*2+8*dimN: dimN+8*dimM;
   double *WORK = (double *)malloc(LWORK*sizeof(double));
   dgesvd_(&JOBU, &JOBVT, &dimM, &dimN, A, &dimM, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &info);
// SUBROUTINE DGESVD( JOBU, JOBVT, M, N,	A, LDA,	S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
   if(info!=0) error("SVD failed.");
   free(WORK);
   free(A);
   free(VT);

//   int dimPC = dimN>20?20:dimN;
   int dimPC = every==1? 20: every;
   if(dimN < dimPC) dimPC = dimN;
   if(dimM < dimPC) dimPC = dimM;
   printf("Largest %d eigenvalues:", dimPC);
   for(int i = 0; i < dimPC; i++)
      printf(" %.2lf", S[i]);
   printf("\n");

   Matrix EV(ID.Length(), dimPC);
   EV.Zero();
   for(int i = 0; i < ID.Length(); i++)
      for(int j = 0; j < dimPC; j++){
         for(int m = 0; m < dimM; m++)
            EV.data[i]->data[j] += Z.data[m]->data[i] * U[j*dimM+m];
         if(S[j] > 1E-100)
            EV[i][j] /= S[j];
      }
   free(U);
   free(S);
#else
   printf("  Please re-compile KING with LAPACK library.\n");
   SVD svd;
   svd.Decompose(X);
   printf("done\n");
   if(svd.n == 0) return;
   QuickIndex idx;
   idx.Index(svd.w);

   int dimPC = dimN>20?20:dimN;
   if(dimM < dimPC) dimPC = dimM;
   printf("Largest %d eigenvalues:", dimPC);
   for(int i = 0; i < dimPC; i++)
      printf(" %.2lf", svd.w[idx[dimN-1-i]]);
   printf("\n");

   Matrix EV(ID.Length(), dimPC);
   EV.Zero();
   for(int i = 0; i < ID.Length(); i++)
      for(int j = 0; j < dimPC; j++){
         for(int m = 0; m < dimM; m++)
            EV[i][j] += Z[m][i] * svd.u[m][idx[dimN-1-j]];
         if(svd.w[idx[dimN-1-j]] > 1E-10)
            EV[i][j] /= svd.w[idx[dimN-1-j]];
      }
#endif

   String datfile = prefix;
   datfile.Add("pc.dat");
   FILE *fp = fopen(datfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)datfile);
   fprintf(fp, "S related\n");
   for(int i = 0; i < dimPC; i++)
      fprintf(fp, "C PC%d\n", i+1);
   fclose(fp);

   String pedfile = prefix;
   pedfile.Add("pc.ped");
   fp = fopen(pedfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
   IntArray typed(ped.count);
   typed.Set(-1);
   for(int i = 0; i < ID.Length(); i++)
      typed[ID[i]] = i;
   IntArray inSVD(ped.count);
   inSVD.Set(-1);
   for(int i = 0; i < ID_UN.Length(); i++)
      inSVD[ID_UN[i]] = i;
   for(int i = 0; i < ped.count; i++){
      if(typed[i] == -1) continue;
      fprintf(fp, "%s %s %s %s %d",
         (const char*)ped[i].famid, (const char*)ped[i].pid,
         (const char*)ped[i].fatid, (const char*)ped[i].motid,
         ped[i].sex);
      /*
         for(int j = 0; j < dimPC+1; j++)
            fprintf(fp, " X");
      else{  */

         if(inSVD[i]!=-1) fprintf(fp, " 1"); // unrelated in SVD
         else fprintf(fp, " 2"); // related in PCA
         for(int j = 0; j < dimPC; j++)
            fprintf(fp, " %.4lf", EV[typed[i]][j]);
      // }
      fprintf(fp, "\n");
   }
   fclose(fp);
   printf("PCA ends at %s", currentTime());
   printf("%d principal components saved in files %s and %s\n",
      dimPC, (const char*)datfile, (const char*)pedfile);
}

void Engine::pca_family_check(int every)
{
   printf("PCA starts at %s", currentTime());

   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   individualInfo = true;
   if(geno.Length()==0) {
      if(shortFlag) BuildShortBinary();
      else BuildBinary();
   }
//   if(shortCount)
      printf("Genotypes stored in %d words for each of %d individuals.\n",
         shortCount, idCount);

   Kinship kin;
   QuickIndex idx;
   int id1, id2;
   Matrix kinship, phi;
   bool flipFlag;
   IntArray unrelatedCount;
   IntArray ID(0), ID_UN(0), ID_unrelated;
   int HetHet, Het1, Het2, IBS0, notMissing;
   int HetHetCount, IBS0Count, het1Count, het2Count, notMissingCount;
   for(int f = 0; f < ped.familyCount; f++){
      if(id[f].Length() < 2) {
         if(id[f].Length() == 1)
            ID_UN.Push(id[f][0]);
         ID.Stack(id[f]);
         continue;
      }
      kin.Setup(*ped.families[f]);
      kinship.Dimension(id[f].Length(), id[f].Length());
      phi.Dimension(id[f].Length(), id[f].Length());
      kinship.Set(1);
      phi.Set(1);
      for(int i = 0; i < id[f].Length(); i++)
         for(int j = i+1; j < id[f].Length(); j++){
            id1 = geno[id[f][i]]; id2 = geno[id[f][j]];
            if(id1 < 0 || id2 < 0) continue;
            HetHetCount = IBS0Count = het1Count = 0;
            for(int m = 0; m < shortCount; m++){
               HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
               het1Count += (oneoneCount[(GG[0][id2][m] | GG[1][id2][m]) & (~GG[0][id1][m]) & GG[1][id1][m]]
                  + oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]]);
               IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
            }
            if(het1Count)
               kinship[i][j] = kinship[j][i] = (HetHetCount - IBS0Count*2)*1.0/het1Count;
            phi[i][j] = phi[j][i] = kin(ped[id[f][i]], ped[id[f][j]]);
         }

      unrelatedCount.Dimension(id[f].Length());
      unrelatedCount.Zero();
      for(int i = 0; i < id[f].Length(); i++)
         for(int j = i+1; j < id[f].Length(); j++)
            if(phi[i][j] < 0.005 && kinship[i][j] < 0.0221){
               unrelatedCount[i] ++;
               unrelatedCount[j] ++;
            }
      idx.Index(unrelatedCount);
      ID_unrelated.Dimension(0);
      ID_unrelated.Push(idx[unrelatedCount.Length()-1]);
      for(int i = unrelatedCount.Length()-2; (i >= 0) && (unrelatedCount[idx[i]] > 0); i--){
         bool isUnrelated = true;
         for(int j = 0; j < ID_unrelated.Length(); j++)
            if(phi[idx[i]][ID_unrelated[j]] >= 0.005 || kinship[idx[i]][ID_unrelated[j]] >= 0.0221){
               isUnrelated = false;
               break;
            }
         if(isUnrelated) ID_unrelated.Push(idx[i]);
      }
      for(int i = 0; i < ID_unrelated.Length(); i++)
         ID_UN.Push(id[f][ID_unrelated[i]]);
      ID.Stack(id[f]);
   }
   int dimN = ID_UN.Length();
//   int dimM = (markerCount-1) / every + 1;
   int dimM = markerCount;
#ifndef WITH_LAPACK
   if(dimM > 100000)
      printf("Note: it may take too long to perform PCA on %d markers.\n", dimM);
#endif
   Matrix X(dimM, dimN);
   Matrix Z(dimM, ID.Length());
   X.Zero(); Z.Zero();
   int minorAllele;
   int validSNP=0;
   int byte, offset;
   double p, g, mean, se;
   int k, freqCount;

   for(int m = 0; m < markerCount; m ++){
      byte = m/16;
      offset = m%16;
      p = 0.0;
      freqCount = 0;
      for(int i = 0; i < dimN; i++){
         k = geno[ID_UN[i]];
         if(k == -1) continue;
         if(GG[1][k][byte] & shortbase[offset]){ // AA or Aa
            p ++;
            if(GG[0][k][byte] & shortbase[offset]) p ++; // AA
            freqCount += 2;
         }else if(GG[0][k][byte] & shortbase[offset]) // aa
            freqCount += 2;
      }
      if(freqCount==0) continue;
      p /= freqCount;
      flipFlag = false;
      if(p > 0.5) {
         flipFlag = true;
         p = 1-p;
      }
      if(p < 0.001) continue;
      mean = p*2;
      se = sqrt(2*p*(1-p));
      if(se < 1E-100) continue;

      for(int i = 0; i < dimN; i++){
         k = geno[ID_UN[i]];
         if(k == -1) continue;
         g = -1.0;
         if(GG[0][k][byte] & shortbase[offset])  // homozygote
            g = (GG[1][k][byte] & shortbase[offset])? 2.0: 0.0;
         else if(GG[1][k][byte] & shortbase[offset])   // Aa
            g = 1.0;
         if(g > -0.9) {
            if(flipFlag)
               X[validSNP][i] = (2 - g - mean) / se;
            else
               X[validSNP][i] = (g - mean) / se;
         }
      }
      for(int i = 0; i < ID.Length(); i++){
         k = geno[ID[i]];
         if(k == -1) continue;
         g = -1.0;
         if(GG[0][k][byte] & shortbase[offset])  // homozygote
            g = (GG[1][k][byte] & shortbase[offset])? 2.0: 0.0;
         else if(GG[1][k][byte] & shortbase[offset])   // Aa
            g = 1.0;
         if(g > -0.9) {
            if(flipFlag)
               Z[validSNP][i] = (2 - g - mean) / se;
            else
               Z[validSNP][i] = (g - mean) / se;
         }
      }
      validSNP++;
   }
   dimM = validSNP;
   X.Dimension(dimM, dimN);
   Z.Dimension(dimM, ID.Length());

   printf("Dimension of PCA: %d %d\n", dimM, dimN);
   printf("SVD...");

#ifdef WITH_LAPACK
   printf("  LAPACK is used.\n");
   char JOBU = 'S';
   char JOBVT = 'O';
   int LDU = dimM;
   int LDVT = 1;
   int dimS = dimM>dimN? dimN: dimM;
   int info;
   double *A = (double *)malloc(dimM*dimN*sizeof(double));
   for(int i = 0; i < dimM; i++)
     for(int j = 0; j < dimN; j++)
       A[j*dimM+i] = X[i][j];
   double *S = (double *)malloc(dimS*sizeof(double));
   double *U = (double *)malloc(LDU*dimS*sizeof(double));
   double *VT = (double *)malloc(LDVT*sizeof(double));
   int LWORK = dimM>dimN? dimM*2+8*dimN: dimN+8*dimM;
   double *WORK = (double *)malloc(LWORK*sizeof(double));
   dgesvd_(&JOBU, &JOBVT, &dimM, &dimN, A, &dimM, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &info);
// SUBROUTINE DGESVD( JOBU, JOBVT, M, N,	A, LDA,	S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
   if(info!=0) error("SVD failed.");
   free(WORK);
   free(A);
   free(VT);

   int dimPC = dimN>20?20:dimN;
   if(dimM < dimPC) dimPC = dimM;
   printf("Largest %d eigenvalues:", dimPC);
   for(int i = 0; i < dimPC; i++)
      printf(" %.2lf", S[i]);
   printf("\n");

   Matrix EV(ID.Length(), dimPC);
   EV.Zero();
   for(int i = 0; i < ID.Length(); i++)
      for(int j = 0; j < dimPC; j++){
         for(int m = 0; m < dimM; m++)
            EV.data[i]->data[j] += Z.data[m]->data[i] * U[j*dimM+m];
         if(S[j] > 1E-100)
            EV[i][j] /= S[j];
      }
   free(U);
   free(S);
#else
   printf("  Please re-compile KING with LAPACK library.\n");
   SVD svd;
   svd.Decompose(X);
   printf("done\n");
   if(svd.n == 0) return;
   idx.Index(svd.w);

//   int dimPC = dimN>20?20:dimN;
   int dimPC = every==1? 20: every;
   if(dimN < dimPC) dimPC = dimN;
   if(dimM < dimPC) dimPC = dimM;
   printf("Largest %d eigenvalues:", dimPC);
   for(int i = 0; i < dimPC; i++)
      printf(" %.2lf", svd.w[idx[dimN-1-i]]);
   printf("\n");
   Matrix EV(ID.Length(), dimPC);
   EV.Zero();
   for(int i = 0; i < ID.Length(); i++)
      for(int j = 0; j < dimPC; j++){
         for(int m = 0; m < dimM; m++)
            EV[i][j] += Z[m][i] * svd.u[m][idx[dimN-1-j]];
         if(svd.w[idx[dimN-1-j]] > 1E-10)
            EV[i][j] /= svd.w[idx[dimN-1-j]];
      }
#endif
   String datfile = prefix;
   datfile.Add("pc.dat");
   FILE *fp = fopen(datfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)datfile);
   fprintf(fp, "S related\n");
   for(int i = 0; i < dimPC; i++)
      fprintf(fp, "C PC%d\n", i+1);
   fclose(fp);

   String pedfile = prefix;
   pedfile.Add("pc.ped");
   fp = fopen(pedfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
   IntArray typed(ped.count);
   typed.Set(-1);
   for(int i = 0; i < ID.Length(); i++)
      typed[ID[i]] = i;
   IntArray inSVD(ped.count);
   inSVD.Set(-1);
   for(int i = 0; i < ID_UN.Length(); i++)
      inSVD[ID_UN[i]] = i;
   for(int i = 0; i < ped.count; i++){
      if(typed[i] == -1) continue;
      fprintf(fp, "%s %s %s %s %d",
         (const char*)ped[i].famid, (const char*)ped[i].pid,
         (const char*)ped[i].fatid, (const char*)ped[i].motid,
         ped[i].sex);

      /*
         for(int j = 0; j < dimPC+1; j++)
            fprintf(fp, " X");
      else{ */
         if(inSVD[i]!=-1) fprintf(fp, " 1"); // unrelated in SVD
         else fprintf(fp, " 2"); // related in PCA
         for(int j = 0; j < dimPC; j++)
            fprintf(fp, " %.4lf", EV[typed[i]][j]);
      fprintf(fp, "\n");
   }
   fclose(fp);
   printf("PCA ends at %s", currentTime());
   printf("%d principal components saved in files %s and %s\n",
      dimPC, (const char*)datfile, (const char*)pedfile);
}

void Engine::mds_kin_family()
{
   printf("MDS with kinship incorporated starts at %s", currentTime());
   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   if(detailFlag) countGenotype();
   if(geno.Length()==0){
      individualInfo = true;
      if(shortFlag) BuildShortBinary();
      else BuildBinary();
   }
   if(xflag)
      printf("X-chromosome genotypes stored in %d words for each of %d individuals.\n",
         xshortCount, idCount);
   else
      printf("Genotypes stored in %d words for each of %d individuals.\n",
         shortCount, idCount);
   if(allflags&(1<<IBSFLAG))
      printf("--ibs is ignored.\n");
   printf("Euclidean Distance is used in the MDS analysis.\n");

   Kinship kin;
   IntArray ID(0), ID_UN(0), ID_unrelated;
   int HetHetCount, IBS0Count, het1Count, het2Count, notMissingCount;
   for(int f = 0; f < ped.familyCount; f++){
      ID_unrelated.Dimension(0);
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(ped[i].ngeno >= MINSNPCOUNT){
            ID.Push(i);
            if(ped[i].isFounder())
               ID_unrelated.Push(i);
         }
      kin.Setup(*ped.families[f]);
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
         if(ped[i].ngeno < MINSNPCOUNT) continue;
         bool isUnrelated = true;
         for(int j = 0; j < ID_unrelated.Length(); j++)
            if(kin(ped[i], ped[ID_unrelated[j]]) > 0) {
               isUnrelated = false;
               break;
            }
         if(isUnrelated) ID_unrelated.Push(i);
      }
      ID_UN.Stack(ID_unrelated);
   }

   int dimN = ID_UN.Length();
   int id1, id2;
   double dist;
   Matrix D(dimN, dimN);
   D.Zero();
   for(int i = 0; i < dimN; i++){
      id1 = geno[ID_UN[i]];
      for(int j = i+1; j < dimN; j++){
         id2 = geno[ID_UN[j]];
         IBS0Count = het1Count = het2Count = notMissingCount = 0;
         if(xflag)
         for(int m = 0; m < xshortCount; m++){
            IBS0Count += oneoneCount[XG[0][id1][m] & XG[0][id2][m] & (XG[1][id1][m] ^ XG[1][id2][m])];
            notMissingCount += oneoneCount[(XG[0][id1][m] | XG[1][id1][m]) & (XG[0][id2][m] | XG[1][id2][m])];
            het1Count += oneoneCount[(XG[0][id2][m] & (~XG[0][id1][m]) & XG[1][id1][m])];
            het2Count += oneoneCount[(XG[0][id1][m] & (~XG[0][id2][m]) & XG[1][id2][m])];
         }
         else
         for(int m = 0; m < shortCount; m++){
            IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
            notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
            het1Count += oneoneCount[(GG[0][id2][m] & (~GG[0][id1][m]) & GG[1][id1][m])];
            het2Count += oneoneCount[(GG[0][id1][m] & (~GG[0][id2][m]) & GG[1][id2][m])];
         }
         if(notMissingCount){
            dist = het1Count + het2Count + 4.0*IBS0Count;
//            if( (dist < het1Count*distcoef) && (dist < het2Count*distcoef) )
//               dist = het1Count + het2Count;
            D[i][j] = D[j][i] =  dist / notMissingCount;
         }
      }
   }
   Vector tempV(dimN);
   tempV.Zero();
   for(int j = 0; j < dimN; j++)
      for(int i = 0; i < dimN; i++)
         tempV[j] += D[i][j];
   tempV.Multiply(1.0/dimN);
   // (I-11'/N) * D
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimN; j++)
         D[i][j] -= tempV[j];
   tempV.Zero();
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimN; j++)
         tempV[i] += D[i][j];
   tempV.Multiply(1.0/dimN);
   // D * (I-11'/N)
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimN; j++)
         D[i][j] -= tempV[i];
   D.Multiply(-0.5);

   printf("SVD start...");
#ifdef WITH_LAPACK
   printf("  LAPACK is used.\n");
   char JOBZ = 'A';
   int info;
   double *A = new double[dimN*dimN];
   int dimLA = dimN;
   for(int i = 0; i < dimN; i++)
     for(int j = 0; j < dimN; j++)
       A[i*dimN+j] = D[j][i];
   double *S = new double[dimN];
   double *U = new double[dimN*dimN];
   double *VT = new double[dimN*dimN];
   int *IWORK = new int[dimN*8];
   int LWORK = 8*dimN + 4*dimN*dimN;
   double *WORK = new double[LWORK];
   dgesdd_(&JOBZ, &dimLA, &dimLA, A, &dimLA, S, U, &dimLA, VT, &dimLA, WORK, &LWORK, IWORK, &info);
// SUBROUTINE DGESDD(JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO)

   delete []U;
   delete []IWORK;
   delete []WORK;
   delete []A;

   int dimPC = dimN>20?20:dimN;
   printf("Largest %d eigenvalues:", dimPC);
   for(int i = 0; i < dimPC; i++)
      printf(" %.2lf", S[i]);
   printf("\n");
   double totalVariance = 0;
   double variance20 = 0;
   for(int i = 0; i < dimPC; i++)
      variance20 += S[i];
   for(int i = 0; (S[i] > 1E-10) && (i < dimN-1); i++)
      totalVariance += S[i];
   printf("The first %d PCs are able to explain %.2lf / %.2lf = %.1lf%% of total variance.\n",
      dimPC, variance20, totalVariance, variance20/totalVariance*100);
   printf("The proportion of total variance explained (%%) by each PC is:\n  ");
   for(int i = 0; i < dimPC; i++)
      printf(" %.1lf", S[i]/totalVariance*100);
   printf("\n");
   delete []S;
#else
   printf("  Please re-compile KING with LAPACK library.\n");
   SVD svd;
   svd.Decompose(D);
   printf("done\n");
   if(svd.n == 0) return;
   QuickIndex idx;
   idx.Index(svd.w);

   int dimPC = dimN>20?20:dimN;
   printf("Largest %d eigenvalues:", dimPC);
   for(int i = 0; i < dimPC; i++)
      printf(" %.2lf", svd.w[idx[dimN-1-i]]);
   printf("\n");
#endif

   Matrix EV(ID.Length(), dimPC);
   EV.Zero();
   ID.Dimension(0);

   IntArray ID_related;
   IntArray flagged;
   int unrelatedCount = 0;
   int totalCount = 0;
   Vector kinship;
   double temp;
   for(int f = 0; f < ped.familyCount; f++){
      ID_unrelated.Dimension(0);
      ID_related.Dimension(0);
      flagged.Dimension(ped.families[f]->last - ped.families[f]->first + 1);
      flagged.Zero();
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(ped[i].ngeno >= MINSNPCOUNT){
            if(ped[i].isFounder()){
               ID_unrelated.Push(i);
               flagged[i-ped.families[f]->first] = 1;
            }
         }
      kin.Setup(*ped.families[f]);
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
         if(ped[i].ngeno < MINSNPCOUNT) continue;
         bool isUnrelated = true;
         for(int j = 0; j < ID_unrelated.Length(); j++)
            if(kin(ped[i], ped[ID_unrelated[j]]) > 0) {
               isUnrelated = false;
               break;
            }
         if(isUnrelated) {
            ID_unrelated.Push(i);
            flagged[i-ped.families[f]->first] = 1;
         }
      }
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(ped[i].ngeno >= MINSNPCOUNT && flagged[i-ped.families[f]->first]==0)
            ID_related.Push(i);
      for(int i = 0; i < ID_unrelated.Length(); i++)
         for(int j = 0; j < dimPC; j++)
#ifdef WITH_LAPACK
            EV[totalCount+i][j] = VT[(unrelatedCount+i)*dimN+j];
#else
            EV[totalCount+i][j] = svd.v[unrelatedCount+i][idx[dimN-1-j]];
#endif
      for(int i = 0; i < ID_related.Length(); i++){
         kinship.Dimension(ID_unrelated.Length());
         for(int j = 0; j < ID_unrelated.Length(); j++)
            kinship[j] = kin(ped[ID_related[i]], ped[ID_unrelated[j]])*2;
         for(int j = 0; j < dimPC; j++){
            temp = 0;
            for(int k = 0; k < ID_unrelated.Length(); k++){
#ifdef WITH_LAPACK
               EV[totalCount+ID_unrelated.Length()+i][j] += kinship[k] * VT[(unrelatedCount+k)*dimN+j];
#else
               EV[totalCount+ID_unrelated.Length()+i][j] += kinship[k] * svd.v[unrelatedCount+k][idx[dimN-1-j]];
#endif
               temp += kinship[k];
            }
            if(temp > 1E-10)
               EV[totalCount+ID_unrelated.Length()+i][j] /= temp;
         }
      }

      unrelatedCount += ID_unrelated.Length();
      totalCount += ID_unrelated.Length() + ID_related.Length();
      ID.Stack(ID_unrelated);
      ID.Stack(ID_related);
   }
   String datfile = prefix;
   datfile.Add("pc.dat");
   FILE *fp = fopen(datfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)datfile);
   fprintf(fp, "S related\n");
   for(int i = 0; i < dimPC; i++)
      fprintf(fp, "C PC%d\n", i+1);
   if(detailFlag){
      fprintf(fp, "C Heterozygosity\n");
      fprintf(fp, "C Quality\n");
   }
   fclose(fp);

   String pedfile = prefix;
   pedfile.Add("pc.ped");
   fp = fopen(pedfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
   IntArray typed(ped.count);
   typed.Set(-1);
   for(int i = 0; i < ID.Length(); i++)
      typed[ID[i]] = i;
   IntArray inSVD(ped.count);
   inSVD.Set(-1);
   for(int i = 0; i < ID_UN.Length(); i++)
      inSVD[ID_UN[i]] = i;
   for(int i = 0; i < ped.count; i++){
      int k = typed[i];
      if(k == -1) continue;
      fprintf(fp, "%s %s %s %s %d",
         (const char*)ped[i].famid, (const char*)ped[i].pid,
         (const char*)ped[i].fatid, (const char*)ped[i].motid,
         ped[i].sex);
                           /*
         for(int j = 0; j < dimPC+1; j++)
            fprintf(fp, " X");
         if(detailFlag)
            fprintf(fp, " X X");
      }else{                 */
         if(inSVD[i]!=-1) fprintf(fp, " 1"); // unrelated in SVD
         else fprintf(fp, " 2"); // related in PCA
         for(int j = 0; j < dimPC; j++)
            fprintf(fp, " %.4lf", EV[k][j]);
/*         if(detailFlag)
            fprintf(fp, " %.3lf %.3lf",
               pAaCount[k] * 1.0 / (pAACount[k] + pAaCount[k] + paaCount[k]),
               GetQuality(pAACount[k], pAaCount[k], paaCount[k]));
*/
      fprintf(fp, "\n");
   }
   fclose(fp);
   printf("MDS with kinship incorporated ends at %s", currentTime());
   printf("%d principal components saved in files %s and %s\n",
      dimPC, (const char*)datfile, (const char*)pedfile);
#ifdef WITH_LAPACK
   delete []VT;
#endif
}

void Engine::mds_moving()
{
   // currently by chromosome
   IntArray segIndex[256];
   IntArray segPartialIndex[256];
   IntArray segMask[256];
   int segCount = 22;

   // setup segments
   if(chromosomes.Length() < markerCount || positions.Length() < markerCount){
      printf("Chromosome and position information is incomplete. Moving-window MDS analysis is aborted.\n");
      return;
   }
   for(int i = 0; i < SEXCHR; i++){
      segIndex[i].Dimension(0);
      segPartialIndex[i].Dimension(0);
      segMask[i].Dimension(0);
   }
   int offset;
   short int mask;

   IntArray chr(segCount);
   int b;
   for(int m = 0; m < markerCount; m+=16){
      b = m / 16;
      for(offset = 1; offset < 16 && m+offset < markerCount; offset++)
         if(chromosomes[m+offset] != chromosomes[m]) break;
      if(offset == 16 || m+offset == markerCount){ // all 16 SNPs on the same chromosome
         if(chromosomes[m] > 0 && chromosomes[m] < SEXCHR)
            segIndex[chromosomes[m]-1].Push(b);
         else
            printf(".");
      }else{ // 16 SNPs on different chromosomes
         chr.Zero();
         for(int i = 0; i < 16 && m+i < markerCount; i++){
            if(chromosomes[m] > 0 && chromosomes[m] < SEXCHR)
               chr[chromosomes[m+i]-1] = 1;
            else
               printf(",");
         }
         for(int i = 0; i < segCount; i++)
            if(chr[i]){ // chromosome i+1 is involved
               segPartialIndex[i].Push(b);
               mask = 0;
               for(int j = 0; j < 16 && m+j < markerCount; j++)
                  if(chromosomes[m+j]-1 == i)
                     mask |= (1<<j);
               segMask[i].Push(mask);
            }
      }
   }

   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   if(detailFlag) countGenotype();
   if(geno.Length()==0) {
      individualInfo = true;
      BuildShortBinary();
   }
   printf("Genotypes stored in %d words for each of %d individuals.\n",
         shortCount, idCount);
   if(allflags&(1<<IBSFLAG))
      printf("IBS Distance is used in the MDS analysis.\n");
   else
      printf("Euclidean Distance is used in the MDS analysis.\n");

   IntArray ID(0);
   for(int n = 0; n < ped.count; n++)
      if(ped[n].ngeno >= MINSNPCOUNT )
         ID.Push(n);

   IntArray typed(ped.count);
   typed.Set(-1);
   for(int i = 0; i < ID.Length(); i++)
      typed[ID[i]] = i;

   int dimN = ID.Length();
   int dimPC = dimN>20?20:dimN;
   int id1, id2;
   int HetHetCount, IBS0Count, IBS2Count, het1Count, notMissingCount;
   Matrix D(dimN, dimN);
   Vector tempV(dimN);

#ifdef WITH_LAPACK
   printf("  LAPACK is used.\n");
   char JOBZ = 'A';
   double *S = new double[dimN];
   double *U = new double[dimN*dimN];
   double *VT = new double[dimN*dimN];
   int *IWORK = new int[dimN*8];
   int LWORK = 8*dimN + 4*dimN*dimN;
   double *WORK = new double[LWORK];
   double *A = new double[dimN*dimN];
   int dimLA = dimN;
   int info;
#else
   printf("  Please re-compile KING with LAPACK library.\n");
   SVD svd;
   QuickIndex idx;
#endif
   Matrix *PC = new Matrix[dimN];
   for(int i = 0; i < dimN; i++){
      PC[i].Dimension(dimPC, segCount);
      PC[i].Zero();
   }

   for(int seg = 0; seg < segCount; seg++){
      printf("%d ", seg+1);
      D.Zero();
      for(int i = 0; i < dimN; i++)
         for(int j = i+1; j < dimN; j++){
            id1 = geno[ID[i]]; id2 = geno[ID[j]];
            if(id1 < 0 || id2 < 0) continue;
            if(allflags&(1<<IBSFLAG)){ // IBS Distance
               IBS0Count = IBS2Count = notMissingCount = 0;
               for(int s = 0; s < segIndex[seg].Length(); s++){
                  int m = segIndex[seg][s];
                  IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
                  IBS2Count += oneoneCount[~(GG[0][id1][m]^GG[0][id2][m]) & ~(GG[1][id1][m]^GG[1][id2][m]) & (GG[0][id1][m] | GG[1][id1][m])];
                  notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
               }
               if(notMissingCount)
                  D[i][j] = D[j][i] = 1-(IBS2Count-IBS0Count)*1.0/notMissingCount;
            }else{ // Euclidean Distance
               HetHetCount = IBS0Count = het1Count = notMissingCount = 0;
               for(int s = 0; s < segIndex[seg].Length(); s++){
                  int m = segIndex[seg][s];
                  HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
                  IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
                  notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
                  het1Count += (oneoneCount[(GG[0][id2][m] | GG[1][id2][m]) & (~GG[0][id1][m]) & GG[1][id1][m]]
                     + oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]]);
               }
               if(notMissingCount)
                  D[i][j] = D[j][i] = (het1Count - 2.0*HetHetCount + 4.0*IBS0Count) / notMissingCount;
            }
         }

      tempV.Zero();
      for(int j = 0; j < dimN; j++)
         for(int i = 0; i < dimN; i++)
            tempV[j] += D[i][j];
      tempV.Multiply(1.0/dimN);
      // (I-11'/N) * D
      for(int i = 0; i < dimN; i++)
         for(int j = 0; j < dimN; j++)
            D[i][j] -= tempV[j];
      tempV.Zero();
      for(int i = 0; i < dimN; i++)
         for(int j = 0; j < dimN; j++)
            tempV[i] += D[i][j];
      tempV.Multiply(1.0/dimN);
      // D * (I-11'/N)
      for(int i = 0; i < dimN; i++)
         for(int j = 0; j < dimN; j++)
            D[i][j] -= tempV[i];
      D.Multiply(-0.5);
#ifdef WITH_LAPACK
      for(int i = 0; i < dimN; i++)
         for(int j = 0; j < dimN; j++)
            A[i*dimN+j] = D[j][i];
      dgesdd_(&JOBZ, &dimLA, &dimLA, A, &dimLA, S, U, &dimLA, VT, &dimLA, WORK, &LWORK, IWORK, &info);
      for(int k = 0; k < dimN; k++)
         for(int j = 0; j < dimPC; j++)
            PC[k][j][seg] = VT[k*dimN+j];
#else
      svd.Decompose(D);
      if(svd.n == 0) return;
      idx.Index(svd.w);
      for(int k = 0; k < dimN; k++)
         for(int j = 0; j < dimPC; j++)
            PC[k][j][seg] = svd.v[k][idx[dimN-1-j]];
#endif
      printf("MDS ends at %s", currentTime());
   }

#ifdef WITH_LAPACK
   delete []U;
   delete []IWORK;
   delete []WORK;
   delete []A;
   delete []S;
   delete []VT;
#endif
   for(int p = 0; p < dimPC; p++){
      String pedfile = prefix;
      pedfile += p+1;
      pedfile.Add("pc.ped");
      FILE *fp = fopen(pedfile, "wt");
      if(fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
      for(int i = 0; i < ped.count; i++){
         int k = typed[i];
         if(k == -1) // continue;
         fprintf(fp, "%s %s %s %s %d",
         (const char*)ped[i].famid, (const char*)ped[i].pid,
         (const char*)ped[i].fatid, (const char*)ped[i].motid,
         ped[i].sex);
         /*
         if(k == -1){
            for(int j = 0; j < segCount+1; j++)
               fprintf(fp, " X");
         }else{  */

            fprintf(fp, " 1");
            for(int j = 0; j < segCount; j++)
               fprintf(fp, " %.4lf", PC[k][p][j]);
         //}
         fprintf(fp, "\n");
      }
      fclose(fp);

      String datfile = prefix;
      datfile += p+1;
      datfile.Add("pc.dat");
      fp = fopen(datfile, "wt");
      if(fp == NULL) error("Cannot open %s to write.", (const char*)datfile);
      fprintf(fp, "S related\n");
      for(int i = 0; i < segCount; i++)
         fprintf(fp, "C PC%d_%d\n", p+1, i+1);
      fclose(fp);

      printf("Principal component %d results saved in files %s and %s\n",
         p+1, (const char*)datfile, (const char*)pedfile);
   }

   delete []PC;
   printf("MDS ends at %s", currentTime());
}


void Engine::mds_kin()
{
   printf("MDS with kinship incorporated starts at %s", currentTime());
   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   if(detailFlag) countGenotype();
   if(geno.Length()==0) {
      individualInfo = true;
      if(shortFlag) BuildShortBinary();
      else BuildBinary();
   }
   if(xflag)
      printf("X-chromosome genotypes stored in %d words for each of %d individuals.\n",
         xshortCount, idCount);
   else
      printf("Genotypes stored in %d words for each of %d individuals.\n",
         shortCount, idCount);
   if(allflags&(1<<IBSFLAG))
      printf("--ibs is ignored.\n");
   printf("Euclidean Distance is used in the MDS analysis.\n");

   IntArray ID(0);
   for(int n = 0; n < ped.count; n++)
      if(ped[n].ngeno >= MINSNPCOUNT )
         ID.Push(n);

   int dimN = ID.Length();
   double dist;
   const double distcoef = 1.911612;//up to 4th-degree relatives removed//1.823223;
   int id1, id2, HetHetCount, IBS0Count, het1Count, het2Count, notMissingCount;
   Matrix D(dimN, dimN);
   D.Zero();
   for(int i = 0; i < dimN; i++){
      id1 = geno[ID[i]];
      for(int j = i+1; j < dimN; j++){
         id2 = geno[ID[j]];
         HetHetCount = IBS0Count = het1Count = het2Count = notMissingCount = 0;
         if(xflag)
         for(int m = 0; m < xshortCount; m++){
            HetHetCount += oneoneCount[(~XG[0][id1][m]) & (XG[1][id1][m]) & (~XG[0][id2][m]) & XG[1][id2][m]];
            IBS0Count += oneoneCount[XG[0][id1][m] & XG[0][id2][m] & (XG[1][id1][m] ^ XG[1][id2][m])];
            notMissingCount += oneoneCount[(XG[0][id1][m] | XG[1][id1][m]) & (XG[0][id2][m] | XG[1][id2][m])];
            het1Count += oneoneCount[(XG[0][id2][m] | XG[1][id2][m]) & (~XG[0][id1][m]) & XG[1][id1][m]];
            het2Count += oneoneCount[(XG[0][id1][m] | XG[1][id1][m]) & (~XG[0][id2][m]) & XG[1][id2][m]];
         }
         else
         for(int m = 0; m < shortCount; m++){
            HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
            IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
            notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
            het1Count += oneoneCount[(GG[0][id2][m] | GG[1][id2][m]) & (~GG[0][id1][m]) & GG[1][id1][m]];
            het2Count += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
         }
         if(notMissingCount){
            dist = het1Count + het2Count - 2.0*HetHetCount + 4.0*IBS0Count;
            if( (dist < het1Count*distcoef) && (dist < het2Count*distcoef) )
               dist = het1Count + het2Count;
            D[i][j] = D[j][i] =  dist / notMissingCount;
         }
      }
   }
   Vector tempV(dimN);
   tempV.Zero();
   for(int j = 0; j < dimN; j++)
      for(int i = 0; i < dimN; i++)
         tempV[j] += D[i][j];
   tempV.Multiply(1.0/dimN);
   // (I-11'/N) * D
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimN; j++)
         D[i][j] -= tempV[j];
   tempV.Zero();
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimN; j++)
         tempV[i] += D[i][j];
   tempV.Multiply(1.0/dimN);
   // D * (I-11'/N)
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimN; j++)
         D[i][j] -= tempV[i];
   D.Multiply(-0.5);
   printf("SVD start...");

 #ifdef WITH_LAPACK
   printf("  LAPACK is used.\n");
   char JOBZ = 'A';
   int info;
   double *A = new double[dimN*dimN];
   int dimLA = dimN;
   for(int i = 0; i < dimN; i++)
     for(int j = 0; j < dimN; j++)
       A[i*dimN+j] = D[j][i];
   double *S = new double[dimN];
   double *U = new double[dimN*dimN];
   double *VT = new double[dimN*dimN];
   int *IWORK = new int[dimN*8];
   int LWORK = 8*dimN + 4*dimN*dimN;
   double *WORK = new double[LWORK];
   dgesdd_(&JOBZ, &dimLA, &dimLA, A, &dimLA, S, U, &dimLA, VT, &dimLA, WORK, &LWORK, IWORK, &info);
// SUBROUTINE DGESDD(JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO)
   delete []U;
   delete []IWORK;
   delete []WORK;
   delete []A;

   int dimPC = dimN>20?20:dimN;
   printf("Largest %d eigenvalues:", dimPC);
   for(int i = 0; i < dimPC; i++)
      printf(" %.2lf", S[i]);
   printf("\n");
   double totalVariance = 0;
   double variance20 = 0;
   for(int i = 0; i < dimPC; i++)
      variance20 += S[i];
   for(int i = 0; (S[i] > 1E-10) && (i < dimN-1); i++)
      totalVariance += S[i];
   printf("The first %d PCs are able to explain %.2lf / %.2lf = %.1lf%% of total variance.\n",
      dimPC, variance20, totalVariance, variance20/totalVariance*100);
   printf("The proportion of total variance explained (%%) by each PC is:\n  ");
   for(int i = 0; i < dimPC; i++)
      printf(" %.1lf", S[i]/totalVariance*100);
   printf("\n");
   delete []S;
#else
   printf("  Please re-compile KING with LAPACK library.\n");
   SVD svd;
   svd.Decompose(D);
   printf("done\n");
   if(svd.n == 0) return;
   QuickIndex idx;
   idx.Index(svd.w);

   int dimPC = dimN>20?20:dimN;
   printf("Largest %d eigenvalues:", dimPC);
   for(int i = 0; i < dimPC; i++)
      printf(" %.2lf", svd.w[idx[dimN-1-i]]);
   printf("\n");
#endif

   String datfile = prefix;
   datfile.Add("pc.dat");
   FILE *fp = fopen(datfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)datfile);
   fprintf(fp, "S related\n");
   for(int i = 0; i < dimPC; i++)
      fprintf(fp, "C PC%d\n", i+1);
   if(detailFlag){
      fprintf(fp, "C Heterozygosity\n");
      fprintf(fp, "C Quality\n");
   }
   fclose(fp);

   String pedfile = prefix;
   pedfile.Add("pc.ped");
   fp = fopen(pedfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
   IntArray typed(ped.count);
   typed.Set(-1);
   for(int i = 0; i < ID.Length(); i++)
      typed[ID[i]] = i;
   for(int i = 0; i < ped.count; i++){
      int k = typed[i];
      if(k == -1) continue;
      fprintf(fp, "%s %s %s %s %d",
         (const char*)ped[i].famid, (const char*)ped[i].pid,
         (const char*)ped[i].fatid, (const char*)ped[i].motid,
         ped[i].sex);
/*
         for(int j = 0; j < dimPC+1; j++)
            fprintf(fp, " X");
         if(detailFlag)
            fprintf(fp, " X X");
      }else{  */
         fprintf(fp, " 1");
         for(int j = 0; j < dimPC; j++)
#ifdef WITH_LAPACK
            fprintf(fp, " %.4lf", VT[k*dimN+j]);
#else
            fprintf(fp, " %.4lf", svd.v[k][idx[dimN-1-j]]);
#endif
/*         if(detailFlag)
            fprintf(fp, " %.3lf %.3lf",
               pAaCount[k] * 1.0 / (pAACount[k] + pAaCount[k] + paaCount[k]),
               GetQuality(pAACount[k], pAaCount[k], paaCount[k]));
*/      //}
      fprintf(fp, "\n");
   }
   fclose(fp);
   printf("PCA with kinship incorporated ends at %s", currentTime());
   printf("%d principal components saved in files %s and %s\n",
      dimPC, (const char*)datfile, (const char*)pedfile);
#ifdef WITH_LAPACK
   delete []VT;
#endif
}

/*
   if(ped.affectionCount>0){
      for(int i = 0; i < ped.count; i++)
         if(ped[i].affections[0]==2){
            projectFlag = true;
            break;
         }
   }
*/


void Engine::pca(int every)
{
   if(ped.affectionCount && projectFlag){
         printf("Only unaffected individuals are used in PCA and affected individuals are projected.\n");
         pca_projection(every);
         return;
      }

   printf("PCA starts at %s", currentTime());
   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++)
      oneoneCount[i] = oneCount[i&255] + oneCount[(i>>8)&255];
   individualInfo = true;
   if(geno.Length()==0) {
      if(shortFlag) BuildShortBinary();
      else BuildBinary();
   }
//   if(shortFlag)
      printf("Genotypes stored in %d words for each of %d individuals.\n",
         shortCount, idCount);

   IntArray ID(0);
   for(int n = 0; n < ped.count; n++)
      if(ped[n].ngeno >= MINSNPCOUNT )
         ID.Push(n);
   int dimN = ID.Length();
//   int dimM = (markerCount-1) / every + 1;
   int dimM = markerCount;
#ifndef WITH_LAPACK
   if(dimM > 100000)
      printf("Note: it may take too long to perform PCA on %d markers.\n", dimM);
#endif
   Matrix X(dimM, dimN);
   X.Zero();
   int minorAllele;
   int validSNP=0;
   int byte, offset;
   double p, g, mean, se;
   int k, freqCount;
   bool flipFlag;

   for(int m = 0; m < markerCount; m ++){
      byte = m/16;
      offset = m%16;

      p = 0.0;
      freqCount = 0;
      for(int i = 0; i < dimN; i++){
         k = geno[ID[i]];
         if(GG[1][k][byte] & shortbase[offset]){ // AA or Aa
            p ++;
            if(GG[0][k][byte] & shortbase[offset]) p ++; // AA
            freqCount += 2;
         }else if(GG[0][k][byte] & shortbase[offset]) // aa
            freqCount += 2;
      }
      if(freqCount==0) continue;
      p /= freqCount;
      flipFlag = false;
      if(p > 0.5) {
         flipFlag = true;
         p = 1-p;
      }
      if(p < 0.001) continue;
      mean = p*2;
      se = sqrt(2*p*(1-p));
      if(se < 1E-100) continue;

      for(int i = 0; i < dimN; i++){
         k = geno[ID[i]];
         if(k == -1) continue;
         g = -1.0;
         if(GG[0][k][byte] & shortbase[offset])  // homozygote
            g = (GG[1][k][byte] & shortbase[offset])? 2.0: 0.0;
         else if(GG[1][k][byte] & shortbase[offset])   // Aa
            g = 1.0;
         if(g > -0.9) {
            if(flipFlag)
               X[validSNP][i] = (2 - g - mean) / se;
            else
               X[validSNP][i] = (g - mean) / se;
         }
      }
      validSNP++;
   }
   dimM = validSNP;
   X.Dimension(dimM, dimN);
   printf("Dimension of PCA: %d %d\n", dimM, dimN);
   printf("SVD...");

#ifdef WITH_LAPACK
   printf("  LAPACK is used.\n");
   char JOBU = 'O';
   char JOBVT = 'A';
   int LDU = 1;
   int LDVT = dimN;
   int dimS = dimM>dimN? dimN: dimM;
   int info;
   double *A = (double *)malloc(dimM*dimN*sizeof(double));
   for(int i = 0; i < dimM; i++)
     for(int j = 0; j < dimN; j++)
       A[j*dimM+i] = X[i][j];
   double *S = (double *)malloc(dimS*sizeof(double));
   double *U = (double *)malloc(LDU*sizeof(double));
   double *VT = (double *)malloc(LDVT*dimN*sizeof(double));
   int LWORK = dimM>dimN? dimM*2+8*dimN: dimN+8*dimM;
   double *WORK = (double *)malloc(LWORK*sizeof(double));
   dgesvd_(&JOBU, &JOBVT, &dimM, &dimN, A, &dimM, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &info);
   // A = ? S VT
// SUBROUTINE DGESVD( JOBU, JOBVT, M, N,	A, LDA,	S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
   if(info!=0) error("SVD failed.");
   free(WORK);
   free(A);

//   int dimPC = dimN>20?20:dimN;
   int dimPC = every==1? 20: every;
   if(dimN < dimPC) dimPC = dimN;
   if(dimM < dimPC) dimPC = dimM;
   printf("Largest %d eigenvalues:", dimPC);
   for(int i = 0; i < dimPC; i++)
      printf(" %.2lf", S[i]);
   printf("\n");
   free(U);
   free(S);
#else
   printf("  Please re-compile KING with LAPACK library.\n");
   SVD svd;
   svd.Decompose(X);
   printf("done\n");
   if(svd.n == 0) return;
   QuickIndex idx;
   idx.Index(svd.w);

   int dimPC = dimN>20?20:dimN;
   printf("Largest %d eigenvalues:", dimPC);
   for(int i = 0; i < dimPC; i++)
      printf(" %.2lf", svd.w[idx[dimN-1-i]]);
   printf("\n");
#endif
   String datfile = prefix;
   datfile.Add("pc.dat");
   FILE *fp = fopen(datfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)datfile);
   fprintf(fp, "S related\n");
   for(int i = 0; i < dimPC; i++)
      fprintf(fp, "C PC%d\n", i+1);
   fclose(fp);

   String pedfile = prefix;
   pedfile.Add("pc.ped");
   fp = fopen(pedfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
   IntArray typed(ped.count);
   typed.Set(-1);
   for(int i = 0; i < ID.Length(); i++)
      typed[ID[i]] = i;
   for(int i = 0; i < ped.count; i++){
      if(typed[i] == -1) continue;
      fprintf(fp, "%s %s %s %s %d",
         (const char*)ped[i].famid, (const char*)ped[i].pid,
         (const char*)ped[i].fatid, (const char*)ped[i].motid,
         ped[i].sex);
      /*
         for(int j = 0; j < dimPC+1; j++)
            fprintf(fp, " X");
      else{ */
         fprintf(fp, " 1");
         for(int j = 0; j < dimPC; j++)
#ifdef WITH_LAPACK
            fprintf(fp, " %.4lf", VT[typed[i]*dimN+j]);
#else
            fprintf(fp, " %.4lf", svd.v[typed[i]][idx[dimN-1-j]]);
#endif
      //}
      fprintf(fp, "\n");
   }
   fclose(fp);
   printf("PCA ends at %s", currentTime());
   printf("%d principal components saved in files %s and %s\n",
      dimPC, (const char*)datfile, (const char*)pedfile);
#ifdef WITH_LAPACK
   free(VT);
#endif
}






/*
void Engine::pca_family_noadmix(int every)
{
   Kinship kin;
   IntArray FID(0), ID_UN(0), ID_unrelated;
   for(int f = 0; f < ped.familyCount; f++){
      ID_unrelated.Dimension(0);
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(ped[i].ngeno >= MINSNPCOUNT){
            if(ped[i].isFounder())
               ID_unrelated.Push(i);
         }
      kin.Setup(*ped.families[f]);
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
         if(ped[i].ngeno < MINSNPCOUNT) continue;
         bool isUnrelated = true;
         for(int j = 0; j < ID_unrelated.Length(); j++)
            if(kin(ped[i], ped[ID_unrelated[j]]) > 0) {
               isUnrelated = false;
               break;
            }
         if(isUnrelated) ID_unrelated.Push(i);
      }
      for(int i = 0; i < ID_unrelated.Length(); i++)
         ID_UN.Push(ID_unrelated[i]);
      if(ID_unrelated.Length()) FID.Push(f);
   }
   int dimN = FID.Length();
   IntArray *pid=NULL;
   pid = new IntArray[dimN];
   for(int f = 0; f < dimN; f++){
      pid[f].Dimension(0);
      kin.Setup(*ped.families[f]);
      for(int i = ped.families[FID[f]]->first; i <= ped.families[FID[f]]->last; i++){
         if(ped[i].ngeno < MINSNPCOUNT) continue;
         bool skip = false;
         for(int j = ped.families[FID[f]]->first; j < i; j++)
            if(kin(ped[i], ped[j]) > 0.4){
               skip = true;
               break;
            }
         if(skip) continue;
         pid[f].Push(i);
      }
   }

   double p;
   if(!freqLoaded)
      for(int m = 0; m < ped.markerCount; m++){
         int alleleCount = ped.GetMarkerInfo(m)->alleleLabels.Length();
         ped.GetMarkerInfo(m)->freq.Dimension(alleleCount);
         if(alleleCount != 3) continue;
         double freq = 0;
         int freqCount = 0;
         for(int i = 0; i < dimN; i++)
            if(ped[ID_UN[i]].markers[m].isKnown()){
               freq += ped[ID_UN[i]].markers[m].countAlleles(2);
               freqCount += 2;
            }
         if(freqCount) freq /= freqCount;
         ped.GetMarkerInfo(m)->freq[2] = freq;
         ped.GetMarkerInfo(m)->freq[1] = 1-freq;
      }

   Cholesky chol;

   int dimM = 0;
   for(int m = 0; m < ped.markerCount; m += every){
      int alleleCount = ped.GetMarkerInfo(m)->freq.Length();
      if(alleleCount != 3) continue;
      p = ped.GetMarkerInfo(m)->freq[2];
      if(p < 0.001 || p > 0.999) continue;
      dimM ++;
   }
   printf("Dimension of PCA: %d %d\n", dimM, dimN);
   Matrix X(dimM, dimN);
   X.Zero();
   int minorAllele;
   int row=0;
   Matrix *phi = NULL;
   phi = new Matrix[dimN];
   Vector OnePhiOne(dimN);
   OnePhiOne.Zero();
   for(int f = 0; f < dimN; f++){
      kin.Setup(*ped.families[FID[f]]);
      phi[f].Dimension(pid[f].Length(),pid[f].Length());
      if(pid[f].Length()==1){
         phi[f][0][0] = 1;
         OnePhiOne[f] = 1;
      }else{
         phi[f].Zero();
         for(int i = 0; i < pid[f].Length(); i++)
            for(int j = i; j < pid[f].Length(); j++)
               phi[f][i][j] = phi[f][j][i] = 2*kin(ped[pid[f][i]], ped[pid[f][j]]);

         chol.Decompose(phi[f]);
         chol.Invert();
         for(int i = 0; i < pid[f].Length(); i++)
            for(int j = 0; j < pid[f].Length(); j++){
               phi[f][i][j] = chol.inv[i][j];
               OnePhiOne[f] += chol.inv[i][j];
            }
      }
   }
   Vector G;
   for(int m = 0; m < ped.markerCount; m += every){
      int alleleCount = ped.GetMarkerInfo(m)->freq.Length();
      if(alleleCount != 3) continue;
      minorAllele = ped.GetMarkerInfo(m)->freq[2]>0.5 ? 1: 2;
      p = ped.GetMarkerInfo(m)->freq[minorAllele];
      if(p < 0.001) continue;
      for(int f = 0; f < dimN; f++){
         G.Dimension(pid[f].Length());
         for(int i = 0; i < pid[f].Length(); i++)
            if(ped[pid[f][i]].isGenotyped(m))
               G[i] = ped[pid[f][i]].markers[m].countAlleles(minorAllele);
            else
               G[i] = p*2;
         for(int i = 0; i < pid[f].Length(); i++)
            for(int j = 0; j < pid[f].Length(); j++)
               X[row][f] += phi[f][i][j] * G[j];
         X[row][f] = (X[row][f] - OnePhiOne[f]*2*p) / sqrt(2*p*(1-p)*OnePhiOne[f]);
      }
      row++;
   }
   printf("SVD...");
   SVD svd;
   svd.Decompose(X);
   printf("done\n");

   if(pid) delete []pid;
   if(phi) delete []phi;
   if(svd.n == 0) return;
   QuickIndex idx;
   idx.Index(svd.w);

   int dimPC = dimN>20?20:dimN;
   printf("Largest %d eigenvalues:", dimPC);
   for(int i = 0; i < dimPC; i++)
      printf(" %.2lf", svd.w[idx[dimN-1-i]]);
   printf("\n");

   String filename = prefix;
   filename.Add(".pc");
   FILE *fp = fopen(filename, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)filename);
   for(int n = 0; n < dimN; n++){
      fprintf(fp, "%s", (const char*)ped.families[FID[n]]->famid);
      for(int i = 0; i < dimPC; i++)
         fprintf(fp, " %.4lf", svd.v[n][idx[dimN-1-i]]);
      fprintf(fp, "\n");
   }
   fclose(fp);
   printf("%d principal components saved in file %s\n",
      dimPC, (const char*)filename);
      }
      */
/*
void Engine::OutputIndividualInfo()
{
   double freq;
   const double MAFe = 0.00000001+minMAF;
   int freqIndex;
   Matrix var(ped.count, binCount);
   Matrix mean(ped.count, binCount);
   var.Zero();
   mean.Zero();
   for(int i = 0; i < ped.count; i++){
      if(ped[i].ngeno < MINSNPCOUNT ) continue;
      for(int k = 0; k < binCount; k++){
         var[i][k] = 2 * freqP[k] * freqQ[k];
         mean[i][k] = 2 * freqP[k];
      }
   }

      String filename = prefix;
      filename.Add(".var");
      FILE *fp = fopen(filename, "wt");
      if(fp == NULL) error("Cannot open %s to write.", (const char*)filename);
      String filename2 = prefix;
      filename2.Add(".mean");
      FILE *fp2 = fopen(filename2, "wt");
      if(fp2 == NULL) error("Cannot open %s to write.", (const char*)filename2);
      Vector geno[1000];
      Vector lambda(binCount);
      Vector lambda2(binCount);
      for(int i = 0; i < ped.count; i++){
         if(ped[i].ngeno < MINSNPCOUNT ) continue;
         for(int k = 0; k < binCount; k++)
            geno[k].Dimension(0);
         for(int m = 0; m < ped.markerCount; m++){
            if(ped.GetMarkerInfo(m)->freq.Length() != 3) continue;
            if(!ped[i].isGenotyped(m)) continue;
            int g = ped[i].markers[m].countAlleles(2);
            freq = ped.GetMarkerInfo(m)->freq[2];
            if(freq > 0.5) {
               freq = 1 - freq;
               g = 2-g;
            }
            if(freq < MAFe) continue;
            freqIndex = int((freq - MAFe)*binCount/(0.5-minMAF));
            geno[freqIndex].Push(g);
         }
         fprintf(fp, "%s\t%s", (const char*)ped[i].famid, (const char*)ped[i].pid);
         fprintf(fp2, "%s\t%s", (const char*)ped[i].famid, (const char*)ped[i].pid);
         lambda.Zero(); lambda2.Zero();
         for(int k = 0; k < binCount; k++){
            var[i][k] = geno[k].Var();
            mean[i][k] = geno[k].Average();
            for(int j = 0; j < geno[k].Length(); j++)
               if(geno[k][j]==2) lambda[k]++;
            lambda[k] = lambda[k] / geno[k].Length() - 0.25 * mean[i][k] * mean[i][k];
            lambda2[k] = mean[i][k] / (2*freqP[k]);
         }
         fprintf(fp, "\t%.4lf", lambda.Average());
         for(int k = 0; k < binCount; k++)
            fprintf(fp, "\t%.4lf", lambda[k]);
         fprintf(fp, "\n");
         fprintf(fp2, "\t%.4lf", lambda2.Average());
         for(int k = 0; k < binCount; k++)
            fprintf(fp2, "\t%.4lf", lambda2[k]);
         fprintf(fp2, "\n");

      }
      fclose(fp);
      printf("Variance for each individual is saved in file %s\n",
         (const char*)filename);
      fclose(fp2);
      printf("Mean for each individual is saved in file %s\n",
         (const char*)filename2);
         }
*/
/*
void Engine::mds_family_check()
{
   individualInfo = true;
   BuildBinary();
   printf("Genotypes stored in %d bytes for each of %d individuals.\n",
      byteCount, idCount);

   Kinship kin;
   QuickIndex idx;
   int id1, id2;
   Matrix kinship, phi;

   IntArray unrelatedCount;
   IntArray ID(0), ID_UN(0), ID_unrelated;
   for(int f = 0; f < ped.familyCount; f++){
      if(id[f].Length() < 2) {
         if(id[f].Length() == 1)
            ID_UN.Push(id[f][0]);
         ID.Stack(id[f]);
         continue;
      }
      kin.Setup(*ped.families[f]);
      kinship.Dimension(id[f].Length(), id[f].Length());
      phi.Dimension(id[f].Length(), id[f].Length());
      kinship.Set(1);
      phi.Set(1);
      for(int i = 0; i < id[f].Length(); i++)
         for(int j = i+1; j < id[f].Length(); j++){
            id1 = geno[id[f][i]]; id2 = geno[id[f][j]];
            if(id1 < 0 || id2 < 0) continue;
            HetHetCount = IBS0Count = het1Count = 0;
            for(int m = 0; m < byteCount; m++){
               HetHet = (~G[0][id1][m]) & (G[1][id1][m]) & (~G[0][id2][m]) & G[1][id2][m];
               Het1 = (G[0][id2][m] | G[1][id2][m]) & (~G[0][id1][m]) & G[1][id1][m];
               Het2 = (G[0][id1][m] | G[1][id1][m]) & (~G[0][id2][m]) & G[1][id2][m];
               IBS0 = G[0][id1][m] & G[0][id2][m] & (G[1][id1][m] ^ G[1][id2][m]);
               notMissing = (G[0][id1][m] | G[1][id1][m]) & (G[0][id2][m] | G[1][id2][m]);
               HetHetCount += oneCount[HetHet];
               het1Count += oneCount[Het1] + oneCount[Het2];
               IBS0Count += oneCount[IBS0];
               notMissingCount += oneCount[notMissing];
            }
            if(notMissingCount)
            if(het1Count)
               kinship[i][j] = kinship[j][i] = (HetHetCount - IBS0Count*2)*1.0/het1Count;
            phi[i][j] = phi[j][i] = kin(ped[id[f][i]], ped[id[f][j]]);
         }

      unrelatedCount.Dimension(id[f].Length());
      unrelatedCount.Zero();
      for(int i = 0; i < id[f].Length(); i++)
         for(int j = i+1; j < id[f].Length(); j++)
            if(phi[i][j] < 0.005 && kinship[i][j] < 0.0221){
               unrelatedCount[i] ++;
               unrelatedCount[j] ++;
            }
      idx.Index(unrelatedCount);
      ID_unrelated.Dimension(0);
      ID_unrelated.Push(idx[unrelatedCount.Length()-1]);
      for(int i = unrelatedCount.Length()-2; (i >= 0) && (unrelatedCount[idx[i]] > 0); i--){
         bool isUnrelated = true;
         for(int j = 0; j < ID_unrelated.Length(); j++)
            if(phi[idx[i]][ID_unrelated[j]] >= 0.005 || kinship[idx[i]][ID_unrelated[j]] >= 0.0221){
               isUnrelated = false;
               break;
            }
         if(isUnrelated) ID_unrelated.Push(idx[i]);
      }
      for(int i = 0; i < ID_unrelated.Length(); i++)
         ID_UN.Push(id[f][ID_unrelated[i]]);
      ID.Stack(id[f]);
   }
   int dimN = ID_UN.Length();

   int dimM = 0;
   for(int m = 0; m < ped.markerCount; m += every){
      int alleleCount = ped.GetMarkerInfo(m)->freq.Length();
      if(alleleCount != 3) continue;
      if(ped.GetMarkerInfo(m)->freq[1] < 0.001 || ped.GetMarkerInfo(m)->freq[1] > 0.999) continue;
      dimM ++;
   }
   Matrix X(dimM, dimN);
   Matrix Z(dimM, ID.Length());
   X.Zero(); Z.Zero();
   int minorAllele;
   int validSNP=0;
   int byte, offset;
   double p, g, mean, se;
   int k, freqCount;

   for(int m = 0; m < ped.markerCount; m += every){
      int alleleCount = ped.GetMarkerInfo(m)->freq.Length();
      if(alleleCount != 3) continue;
      byte = m/8;
      offset = m%8;

      p = 0.0;
      freqCount = 0;
      for(int i = 0; i < dimN; i++){
         k = geno[ID_UN[i]];
         if(k == -1) continue;
         if(G[1][k][byte] & base[offset]){ // AA or Aa
            p ++;
            if(G[0][k][byte] & base[offset]) p ++; // AA
            freqCount += 2;
         }else if(G[0][k][byte] & base[offset]) // aa
            freqCount += 2;
      }
      if(freqCount==0) continue;
      p /= freqCount;
      if(p > 0.5) p = 1-p;
      if(p < 0.001) continue;
      mean = p*2;
      se = sqrt(2*p*(1-p));
      if(se < 1E-100) continue;

      for(int i = 0; i < dimN; i++){
         k = geno[ID_UN[i]];
         if(k == -1) continue;
         g = -1.0;
         if(G[0][k][byte] & base[offset])  // homozygote
            g = (G[1][k][byte] & base[offset])? 2.0: 0.0;
         else if(G[1][k][byte] & base[offset])   // Aa
            g = 1.0;
         if(g > -0.9) X[validSNP][i] = (g - mean) / se;
      }
      for(int i = 0; i < ID.Length(); i++){
         k = geno[ID[i]];
         if(k == -1) continue;
         g = -1.0;
         if(G[0][k][byte] & base[offset])  // homozygote
            g = (G[1][k][byte] & base[offset])? 2.0: 0.0;
         else if(G[1][k][byte] & base[offset])   // Aa
            g = 1.0;
         if(g > -0.9) Z[validSNP][i] = (g - mean) / se;
      }
      validSNP++;
   }

   printf("Dimension of PCA: %d %d\n", dimM, dimN);
   printf("SVD...");
   SVD svd;
   svd.Decompose(X);
   printf("done\n");
   if(svd.n == 0) return;
   idx.Index(svd.w);

   int dimPC = dimN>20?20:dimN;
   printf("Largest %d eigenvalues:", dimPC);
   for(int i = 0; i < dimPC; i++)
      printf(" %.2lf", svd.w[idx[dimN-1-i]]);
   printf("\n");

   Matrix EV(ID.Length(), dimPC);
   EV.Zero();
   for(int i = 0; i < ID.Length(); i++)
      for(int j = 0; j < dimPC; j++){
         for(int m = 0; m < dimM; m++)
            EV[i][j] += Z[m][i] * svd.u[m][idx[dimN-1-j]];
         if(svd.w[idx[dimN-1-j]] > 1E-10)
            EV[i][j] /= svd.w[idx[dimN-1-j]];
      }

   String datfile = prefix;
   datfile.Add("pc.dat");
   FILE *fp = fopen(datfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)datfile);
   fprintf(fp, "S related\n");
   for(int i = 0; i < dimPC; i++)
      fprintf(fp, "C PC%d\n", i+1);
   fclose(fp);

   String pedfile = prefix;
   pedfile.Add("pc.ped");
   fp = fopen(pedfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
   IntArray typed(ped.count);
   typed.Set(-1);
   for(int i = 0; i < ID.Length(); i++)
      typed[ID[i]] = i;
   IntArray inSVD(ped.count);
   inSVD.Set(-1);
   for(int i = 0; i < ID_UN.Length(); i++)
      inSVD[ID_UN[i]] = i;
   for(int i = 0; i < ped.count; i++){
      fprintf(fp, "%s %s %s %s %d",
         (const char*)ped[i].famid, (const char*)ped[i].pid,
         (const char*)ped[i].fatid, (const char*)ped[i].motid,
         ped[i].sex);
      if(typed[i] == -1)
         for(int j = 0; j < dimPC+1; j++)
            fprintf(fp, " X");
      else{
         if(inSVD[i]!=-1) fprintf(fp, " 1"); // unrelated in SVD
         else fprintf(fp, " 2"); // related in PCA
         for(int j = 0; j < dimPC; j++)
            fprintf(fp, " %.4lf", EV[typed[i]][j]);
      }
      fprintf(fp, "\n");
   }
   fclose(fp);
   printf("%d principal components saved in files %s and %s\n",
      dimPC, (const char*)datfile, (const char*)pedfile);
} */

/*
   for(int f = 0; f < ped.familyCount; f++){
      ID_unrelated.Dimension(0);
      ID_related.Dimension(0);
      flagged.Dimension(ped.families[f]->last - ped.families[f]->first + 1);
      flagged.Zero();
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(ped[i].ngeno >= MINSNPCOUNT){
            if(ped[i].isFounder()){
               ID_unrelated.Push(i);
               flagged[i-ped.families[f]->first] = 1;
            }
         }
      kin.Setup(*ped.families[f]);
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
         if(ped[i].ngeno < MINSNPCOUNT) continue;
         bool isUnrelated = true;
         for(int j = 0; j < ID_unrelated.Length(); j++)
            if(kin(ped[i], ped[ID_unrelated[j]]) > 0) {
               isUnrelated = false;
               break;
            }
         if(isUnrelated) {
            ID_unrelated.Push(i);
            flagged[i-ped.families[f]->first] = 1;
         }
      }
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(ped[i].ngeno >= MINSNPCOUNT && flagged[i-ped.families[f]->first]==0)
            ID_related.Push(i);

      for(int i = 0; i < ID_unrelated.Length(); i++)
         for(int j = 0; j < dimPC; j++)
            EV[totalCount+i][j] = svd.v[unrelatedCount+i][idx[dimN-1-j]];

      for(int i = 0; i < ID_related.Length(); i++){
         kinship.Dimension(ID_unrelated.Length());
         for(int j = 0; j < ID_unrelated.Length(); j++)
            kinship[j] = kin(ped[ID_related[i]], ped[ID_unrelated[j]])*2;
         for(int j = 0; j < dimPC; j++){
            temp = 0;
            for(int k = 0; k < ID_unrelated.Length(); k++){
               EV[totalCount+ID_unrelated.Length()+i][j] += kinship[k] * svd.v[unrelatedCount+k][idx[dimN-1-j]];
               temp += kinship[k];
            }
            if(temp > 1E-10)
               EV[totalCount+ID_unrelated.Length()+i][j] /= temp;
         }
      }

      unrelatedCount += ID_unrelated.Length();
      totalCount += ID_unrelated.Length() + ID_related.Length();
      ID.Stack(ID_unrelated);
      ID.Stack(ID_related);
   }
*/
/*
void Engine::mds_family()
{
   printf("MDS for family data starts at %s", currentTime());
   char oneoneCount[65536];
   for(int i = 0; i < 65536; i++){
      oneoneCount[i] = 0;
      for(int j = 0; j < 16; j++)
         if(i & shortbase[j]) oneoneCount[i]++;
   }
   if(detailFlag) countGenotype();
   if(geno.Length()==0) {
      individualInfo = true;
      if(shortFlag) BuildShortBinary();
      else BuildBinary();
   }
   if(shortFlag)
      printf("Genotypes stored in %d short words for each of %d individuals.\n",
         shortCount, idCount);
   else
      printf("Genotypes stored in %d bytes for each of %d individuals.\n",
         byteCount, idCount);
   if(ibsFlag)
      printf("IBS Distance is used in the MDS analysis.\n");
   else
      printf("Euclidean Distance is used in the MDS analysis.\n");
   Kinship kin;
   IntArray first(ped.familyCount);
   IntArray last(ped.familyCount);
   first.Set(-1); last.Set(-1);
   IntArray ID(0), ID_UN(0), ID_unrelated;
   int HetHetCount, IBS0Count, IBS2Count, het1Count, notMissingCount;
   for(int f = 0; f < ped.familyCount; f++){
      ID_unrelated.Dimension(0);
      bool firstFlag = true;
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(ped[i].ngeno >= MINSNPCOUNT){
            ID.Push(i);
            if(ped[i].isFounder()){
               ID_unrelated.Push(i);
               if(firstFlag){
                  first[f] = ID_unrelated.Length()-1+ID_UN.Length();
                  firstFlag = false;
               }
            }
         }
      kin.Setup(*ped.families[f]);
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
         if(ped[i].ngeno < MINSNPCOUNT) continue;
         bool isUnrelated = true;
         for(int j = 0; j < ID_unrelated.Length(); j++)
            if(kin(ped[i], ped[ID_unrelated[j]]) > 0) {
               isUnrelated = false;
               break;
            }
         if(isUnrelated) {
            ID_unrelated.Push(i);
            if(firstFlag){
               first[f] = ID_unrelated.Length()-1+ID_UN.Length();
               firstFlag = false;
            }
         }
//         last[f] = ID_unrelated.Length()-1+ID_UN.Length();
      }
      last[f] = ID_unrelated.Length()-1+ID_UN.Length();
      ID_UN.Stack(ID_unrelated);
   }
   int dimN = ID_UN.Length();
   int id1, id2;
   Matrix D(dimN, dimN);
   D.Zero();
   for(int i = 0; i < dimN; i++){
      id1 = geno[ID_UN[i]];
      for(int j = i+1; j < dimN; j++){
         id2 = geno[ID_UN[j]];
         if(id1 < 0 || id2 < 0) continue;
         if(ibsFlag){ // IBS Distance
            IBS0Count = IBS2Count = notMissingCount = 0;
            if(shortFlag)
               for(int m = 0; m < shortCount; m++){
                  IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
                  IBS2Count += oneoneCount[~(GG[0][id1][m]^GG[0][id2][m]) & ~(GG[1][id1][m]^GG[1][id2][m]) & (GG[0][id1][m] | GG[1][id1][m])];
                  notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
               }
            else
               for(int m = 0; m < byteCount; m++){
                  IBS0Count += oneCount[G[0][id1][m] & G[0][id2][m] & (G[1][id1][m] ^ G[1][id2][m])];
                  IBS2Count += oneCount[~(G[0][id1][m]^G[0][id2][m]) & ~(G[1][id1][m]^G[1][id2][m]) & (G[0][id1][m] | G[1][id1][m])];
                  notMissingCount += oneCount[(G[0][id1][m] | G[1][id1][m]) & (G[0][id2][m] | G[1][id2][m])];
               }
            if(notMissingCount)
               D[i][j] = D[j][i] = 1-(IBS2Count-IBS0Count)*1.0/notMissingCount;
         }else{ // Euclidean Distance
            HetHetCount = IBS0Count = het1Count = notMissingCount = 0;
            if(shortFlag)
               for(int m = 0; m < shortCount; m++){
                  HetHetCount += oneoneCount[(~GG[0][id1][m]) & (GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]];
                  IBS0Count += oneoneCount[GG[0][id1][m] & GG[0][id2][m] & (GG[1][id1][m] ^ GG[1][id2][m])];
                  notMissingCount += oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (GG[0][id2][m] | GG[1][id2][m])];
                  het1Count += (oneoneCount[(GG[0][id2][m] | GG[1][id2][m]) & (~GG[0][id1][m]) & GG[1][id1][m]]
                     + oneoneCount[(GG[0][id1][m] | GG[1][id1][m]) & (~GG[0][id2][m]) & GG[1][id2][m]]);
               }
            else
               for(int m = 0; m < byteCount; m++){
                  HetHetCount += oneCount[(~G[0][id1][m]) & (G[1][id1][m]) & (~G[0][id2][m]) & G[1][id2][m]];
                  IBS0Count += oneCount[G[0][id1][m] & G[0][id2][m] & (G[1][id1][m] ^ G[1][id2][m])];
                  notMissingCount += oneCount[(G[0][id1][m] | G[1][id1][m]) & (G[0][id2][m] | G[1][id2][m])];
                  het1Count += (oneCount[(G[0][id2][m] | G[1][id2][m]) & (~G[0][id1][m]) & G[1][id1][m]]
                     + oneCount[(G[0][id1][m] | G[1][id1][m]) & (~G[0][id2][m]) & G[1][id2][m]]);
            }
            if(notMissingCount)
               D[i][j] = D[j][i] = (het1Count - 2.0*HetHetCount + 4.0*IBS0Count) / notMissingCount;
         }
      }
   }
   Vector tempV(dimN);
   tempV.Zero();
   for(int j = 0; j < dimN; j++)
      for(int i = 0; i < dimN; i++)
         tempV[j] += D[i][j];
   tempV.Multiply(1.0/dimN);
   // (I-11'/N) * D
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimN; j++)
         D[i][j] -= tempV[j];
   tempV.Zero();
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimN; j++)
         tempV[i] += D[i][j];
   tempV.Multiply(1.0/dimN);
   // D * (I-11'/N)
   for(int i = 0; i < dimN; i++)
      for(int j = 0; j < dimN; j++)
         D[i][j] -= tempV[i];
   D.Multiply(-0.5);

   printf("SVD of a %d x %d matrix start...", dimN, dimN);
   SVD svd;
   svd.Decompose(D);
   printf("done\n");
   if(svd.n == 0) return;
   QuickIndex idx;
   idx.Index(svd.w);

   int dimPC = dimN>20?20:dimN;
   printf("Largest %d eigenvalues:", dimPC);
   for(int i = 0; i < dimPC; i++)
      printf(" %.2lf", svd.w[idx[dimN-1-i]]);
   printf("\n");

   Matrix EV(ID.Length(), dimPC);
   EV.Zero();
   ID.Dimension(0);

   IntArray ID_related;
   IntArray flagged;
   Vector kinship;
   double temp;

   flagged.Dimension(ped.count);
   flagged.Zero();
   for(int i = 0; i < ID_UN.Length(); i++)
      flagged[ID_UN[i]] = 1;
   for(int f = 0; f < ped.familyCount; f++){
      if(last[f] < 0) continue;
      kin.Setup(*ped.families[f]);
      ID_related.Dimension(0);
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         if(ped[i].ngeno >= MINSNPCOUNT && flagged[i]==0)
            ID_related.Push(i);
      for(int i = first[f]; i <= last[f]; i++)
         for(int j = 0; j < dimPC; j++)
            EV[geno[ID_UN[i]]][j] = svd.v[i][idx[dimN-1-j]];
      for(int i = 0; i < ID_related.Length(); i++){
         kinship.Dimension(last[f]-first[f]+1);
         for(int j = first[f]; j <= last[f]; j++)
            kinship[j-first[f]] = kin(ped[ID_related[i]], ped[ID_UN[j]])*2;
         for(int j = 0; j < dimPC; j++){
            temp = 0;
            for(int k = first[f]; k <= last[f]; k++){
               EV[geno[ID_related[i]]][j] += kinship[k-first[f]] * svd.v[k][idx[dimN-1-j]];
                 temp += kinship[k-first[f]];
            }
            if(temp > 1E-10)
               EV[geno[ID_related[i]]][j] /= temp;
         }
      }
   }

   String datfile = prefix;
   datfile.Add("pc.dat");
   FILE *fp = fopen(datfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)datfile);
   fprintf(fp, "S related\n");
   for(int i = 0; i < dimPC; i++)
      fprintf(fp, "C PC%d\n", i+1);
   if(detailFlag){
      fprintf(fp, "C Heterozygosity\n");
      fprintf(fp, "C Quality\n");
   }
   fclose(fp);

   String pedfile = prefix;
   pedfile.Add("pc.ped");
   fp = fopen(pedfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
   IntArray inSVD(ped.count);
   inSVD.Set(-1);
   for(int i = 0; i < ID_UN.Length(); i++)
      inSVD[ID_UN[i]] = i;
   for(int i = 0; i < ped.count; i++){
      fprintf(fp, "%s %s %s %s %d",
         (const char*)ped[i].famid, (const char*)ped[i].pid,
         (const char*)ped[i].fatid, (const char*)ped[i].motid,
         ped[i].sex);
      int k = geno[i];
      if(k == -1){
         for(int j = 0; j < dimPC+1; j++)
            fprintf(fp, " X");
         if(detailFlag)
            fprintf(fp, " X X");
      }else{
         if(inSVD[i]!=-1) fprintf(fp, " 1"); // unrelated in SVD
         else fprintf(fp, " 2"); // related in PCA
         for(int j = 0; j < dimPC; j++)
            fprintf(fp, " %.4lf", EV[k][j]);
         if(detailFlag)
            fprintf(fp, " %.3lf %.3lf",
             pAaCount[k] * 1.0 / (pAACount[k] + pAaCount[k] + paaCount[k]),
            GetQuality(pAACount[k], pAaCount[k], paaCount[k]));
      }
      fprintf(fp, "\n");
   }
   fclose(fp);
   printf("MDS for family data ends at %s", currentTime());
   printf("%d principal components saved in files %s and %s\n",
      dimPC, (const char*)datfile, (const char*)pedfile);
}
  */
  

