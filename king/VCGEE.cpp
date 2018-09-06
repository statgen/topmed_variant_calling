//////////////////////////////////////////////////////////////////////
// VCGEE.cpp
// Author: Wei-Min Chen
// Oct 10, 2004

#include "VCGEE.h"
#include "Pedigree.h"
#include "Kinship.h"
#include "MathStats.h"
#include "QuickIndex.h"
#include <math.h>

GEE::GEE(Pedigree & pedigree):ped(pedigree)
{
   pheno = NULL; traits = NULL; covariates = NULL; OD = NULL;
   LoopCount = 20;
   Epsilon = .0002;
   deltaScale = 1.0;
/*   Epsilon = .000002;
   deltaScale = .5;
   LoopCount = 50;
*/
   AtBorder = 0;
   borderIndex.Dimension(0);
   saturatedMean = 0;
   meanPerFamily.Dimension(ped.familyCount);
   meanPerFamily.Zero();
   moreFlag = 0;
   polyFlag = 0;
   variances.SetLabel("Variances");
   coef.SetLabel("Coefficients");
   SEcoef.SetLabel("SE");
   SEcoef_R.SetLabel("Robust SE");
   SEvariances.SetLabel("SE");
   SEvariances_R.SetLabel("Robust SE");
   covariance.SetLabel("Covariance");
   covariance_R.SetLabel("Cov_R");
   CovCoef.SetLabel("CoefCov");
}

GEE::~GEE()
{
   if(pheno) delete []pheno;
   if(traits) delete []traits;
   if(covariates) delete []covariates;
   if(OD) delete []OD;
}

void GEE::GetPhi(int f)
{
   Kinship kin;
   Phi.Dimension(size, size);
   kin.Setup(*ped.families[f]);

   for(int i = 0; i < size; i++){
      Phi[i][i] = 2 * kin(ped[pheno[f][i]], ped[pheno[f][i]]);
      for(int j = i+1; j < size; j++)
         Phi[j][i] = Phi[i][j] = 2 * kin(ped[pheno[f][i]], ped[pheno[f][j]]);
   }
}

void GEE::Refresh(int f)
{
   GetPhi(f);
   RefreshO(f);
}

int GEE::constraint(void)
{
   int ret = 0;
   for(int i = 0; i < parCount; i++)
      if(variances[i] < 0){
         variances[i] = 0.0000001;
         ret = 1;
      }
   return ret;
}

void GEE::solve()
{
   InitCoef();

   String logfile;

   if(polyFlag){
      logfile = prefix;
      logfile += "poly.log";
      polyfp = fopen((const char *)logfile, "wt");
      if(polyfp == NULL) error("Cannot open %s to write", (const char*)logfile);
      fprintf(polyfp, "Likelihood maximization...\n\n");
      fprintf(polyfp, "Iteration 0...\n");
      variances.Print(polyfp);
      coef.Print(polyfp);
   }
   double sum;
   Matrix OS;
   Cholesky chol;

   Matrix OC;     // OC = O * C
   Vector OR;     // OR = O * R
   Vector buf;
   SEvariances.Dimension(parCount);
   SEcoef.Dimension(coefCount);
   DVD.Dimension(parCount, parCount);
   Matrix BOB(coefCount, coefCount);
   Vector DVS(parCount);
   Vector BOY(coefCount);
   Matrix DVS_FAM(parCount, ped.familyCount);
   Matrix BOY_FAM(coefCount, ped.familyCount);
   if(OD) delete[]OD;
   OD = new Matrix[parCount];
   double delta0 = -100;
   double delta1 = -100;
   for(int loop = 0; loop < LoopCount; loop++){
      loglik = 0;
      DVD.Zero();
      BOB.Zero();
      DVS_FAM.Zero();
      BOY_FAM.Zero();
      for (int f = 0; f < ped.familyCount; f++){
         size = pheno[f].Length();
         if(size==0) continue;
         if(saturatedMean && size == 1) continue;
         Refresh(f);

         OR.Dimension(size);
         OC.Dimension(coefCount, size);
         for(int i = 0; i < parCount; i++)
            OD[i].Dimension(size, size);
         OS.Dimension(size, size);

         buf = traits[f];
         for(int i = 0; i < size; i++){
            if(saturatedMean) buf[i] -= meanPerFamily[f];
            for(int j = 0; j < coefCount; j++)
                buf[i] -= coef[j] * covariates[f][j][i];
         }

         if(saturatedMean){
            Vector means(coefCount+1);
            means.Zero();
//            CholeskyOmega.Invert();
            double count = 0;
            for(int i = 0; i < size; i++)
               for(int j = 0; j < size; j++){
//                  count += CholeskyOmega.inv[i][j];
//                  means[coefCount] += traits[f][i]*CholeskyOmega.inv[i][j];
//                  for(int k = 0; k < coefCount; k++)
//                     means[k] += covariates[f][k][i]*CholeskyOmega.inv[i][j];
                }
            for(int k = 0; k < coefCount+1; k++)
               means[k] /= count;
            Vector temp;
            for(int k = 0; k < coefCount; k++){
               temp = covariates[f][k];
               for(int i = 0; i < size; i++)
                  temp[i] -= means[k];
//               CholeskyOmega.BackSubst0(temp);
//               OC[k] = CholeskyOmega.x;
            }
            temp = buf;
            for(int i = 0; i < size; i++)
               temp[i] -= means[coefCount];
//            CholeskyOmega.BackSubst0(temp);
//            OR = CholeskyOmega.x;
         }else{}
         RefreshOD(f);

         W2.Product(covariates[f], OmegaInv);
         for(int u = 0; u < size; u++)
            for(int i = 0; i < coefCount; i++){
               BOY_FAM.data[i]->data[f] += W2.data[i]->data[u] * buf[u];
               for(int j = 0; j < coefCount; j++)
                  BOB.data[i]->data[j] += W2.data[i]->data[u] * covariates[f].data[j]->data[u];
            }
         Vector tempV(size);
         tempV.Zero();
         for(int u = 0; u < size; u++)
            for(int v = 0; v < size; v++)
               tempV[u] += OmegaInv.data[u]->data[v] * buf[v];
         for(int u = 0; u < size; u++)
            for(int v = 0; v < size; v++)
               OS.data[u]->data[v] = tempV[u] * buf[v];
         for(int u = 0; u < size; u++)
            OS[u][u] --;

         // DVDij = Tr(ODiODj)/2
         for(int u = 0; u < size; u++)
            for(int v = 0; v < size; v++)
               for(int i = 0; i < parCount; i++){
                  for(int j = i; j < parCount; j++)
                     DVD.data[i]->data[j] += OD[i].data[u]->data[v] * OD[j].data[v]->data[u];
                  DVS_FAM.data[i]->data[f] += OD[i].data[u]->data[v] * OS.data[v]->data[u] * 0.5 ;
               }
      }  // end of f

      for(int i = 0; i < coefCount; i++) BOY[i] = BOY_FAM[i].Sum();
      for(int i = 0; i < parCount; i++) DVS[i] = DVS_FAM[i].Sum();
      CovCoef.CholeskyInvert(BOB);
      delta2.Dimension(coefCount);
      delta2.Zero();
      for(int i = 0; i < coefCount; i++)
         for(int j = 0; j < coefCount; j++)
            delta2[i] += CovCoef[i][j] * BOY[j];
      coef.Add(delta2);

      if(saturatedMean){
         Vector temp;
         for (int f = 0; f < ped.familyCount; f++){
            size = pheno[f].Length();
            if(size==0) continue;
         if(saturatedMean && size ==1) continue;
            GetPhi(f);
            RefreshO(f);
            temp.Dimension(size);
//            CholeskyOmega.Decompose(Omega);
//            CholeskyOmega.Invert();

            for(int i = 0; i < size; i++){
               temp[i] = traits[f][i];
               for(int j = 0; j < coefCount; j++)
                  temp[i] -= coef[j] * covariates[f][j][i];
            }
            double count=0;
            meanPerFamily[f] = 0;
            for(int i = 0; i < size; i++)
               for(int j = 0; j < size; j++){
//                  count += CholeskyOmega.inv[i][j];
//                  meanPerFamily[f] += CholeskyOmega.inv[i][j] * temp[j];
               }
            meanPerFamily[f] /= count;
         }
      }
   
/*      if(AtBorder){
         newDVD.Dimension(borderIndex.Length(), borderIndex.Length());
         for(int i = 0; i < borderIndex.Length(); i++)
            for(int j = i; j < borderIndex.Length(); j++)
               newDVD[i][j] = newDVD[j][i] = DVD[borderIndex[i]][borderIndex[j]] * 0.5;
         W2.CholeskyInvert(newDVD);
         covariance.Dimension(parCount, parCount);
         covariance.Zero();
         for(int i = 0; i < borderIndex.Length(); i++)
            for(int j = 0; j < borderIndex.Length(); j++)
               covariance[borderIndex[i]][borderIndex[j]] = W2[i][j];
      }else{   */
         for(int i = 0; i < parCount; i++)
            for(int j = i; j < parCount; j++)
               DVD[j][i] = DVD[i][j] = DVD[i][j] * 0.5;
         covariance.CholeskyInvert(DVD);
//      }
      delta.Dimension(parCount);
      delta.Zero();
      for(int i = 0; i < parCount; i++)
         for(int j = 0; j < parCount; j++)
            delta[i] += covariance[i][j] * DVS[j];
      if(deltaScale != 1.0)
         for(int i = 0; i < parCount; i++)
            delta[i] *= deltaScale;
      variances.Add(delta);
      if(constraint()) AtBorder = 1;
      else AtBorder = 0;

      if(polyFlag) {
         fprintf(polyfp, "Iteration %d...\n", loop+1);
         variances.Print(polyfp);
         coef.Print(polyfp);
      }
//      if(moreFlag){ variances.Print(); coef.Print(); }
//      variances.Print(); coef.Print();

      if(StopRule()) break;
      // dead loop: delta == delta0
//      if(AtBorder && delta0>0 && fabs(delta.SumSquares()-delta0) < Epsilon*Epsilon) break;
      if(AtBorder && delta0>0 && fabs(delta.SumSquares()-delta0) < Epsilon) break;
      delta0 = delta1;
      delta1 = AtBorder? delta.SumSquares(): -100;
   }  // end of loop

   Matrix BOB2(coefCount, coefCount);
   BOB2.Zero();
   for(int i = 0; i < coefCount; i++)
      for(int j = 0; j < coefCount; j++)
         for(int f = 0; f < ped.familyCount; f++)
              BOB2[i][j] += BOY_FAM[i][f] * BOY_FAM[j][f];
   CovCoef_R.Dimension(coefCount, coefCount);
   CovCoef_R.Zero();
   for(int i = 0; i < coefCount; i++)
      for(int j = i; j < coefCount; j++){
         for(int u = 0; u < coefCount; u++)
            for(int v = 0; v < coefCount; v++)
               CovCoef_R[i][j] += CovCoef[i][u] * BOB2[u][v] * CovCoef[v][j];
         CovCoef_R[j][i] = CovCoef_R[i][j];
      }

   SEcoef_R.Dimension(coefCount);
   for(int i = 0; i < coefCount; i++){
      if(CovCoef[i][i]>0) SEcoef[i] = sqrt(CovCoef[i][i]);
      else continue;
      if(CovCoef_R[i][i]>0) SEcoef_R[i] = sqrt(CovCoef_R[i][i]);
      else continue;
   }

   Matrix DVD2(parCount, parCount);
   DVD2.Zero();
   for(int i = 0; i < parCount; i++)
      for(int j = 0; j < parCount; j++)
         for(int f = 0; f < ped.familyCount; f++)
              DVD2[i][j] += DVS_FAM[i][f] * DVS_FAM[j][f];
   covariance_R.Dimension(parCount, parCount);
   covariance_R.Zero();
   for(int i = 0; i < parCount; i++)
      for(int j = i; j < parCount; j++){
         for(int u = 0; u < parCount; u++)
            for(int v = 0; v < parCount; v++)
               covariance_R[i][j] += covariance[i][u] * DVD2[u][v] * covariance[v][j];
         covariance_R[j][i] = covariance_R[i][j];
      }
   SEvariances_R.Dimension(parCount);
   for(int i = 0; i < parCount; i++){
      if(covariance[i][i]>0) SEvariances[i] = sqrt(covariance[i][i]);
      else continue;
      if(covariance_R[i][i]>0) SEvariances_R[i] = sqrt(covariance_R[i][i]);
      else continue;
   }


  // MLE stuff
   loglik = 0;
   sum=0;
   for (int f = 0; f < ped.familyCount; f++){
      size = pheno[f].Length();
      if(size==0) continue;
      sum += size;
      GetPhi(f);
      RefreshO(f);
      buf = traits[f];
      for(int i = 0; i < size; i++)
         for(int j = 0; j < coefCount; j++)
            buf[i] -= coef[j] * covariates[f][j][i];
      chol.Decompose(Omega);
      chol.BackSubst0(buf);
      for(int i = 0; i < size; i++) {
         loglik -= chol.x[i] * chol.x[i] / 2;
         loglik -= log(fabs(chol.L[i][i]));
      }
   }

//   printf("loglik (no constant): %lf\n\n", loglik);

   if(polyFlag){
      fprintf(polyfp, "\nSUMMARY\nloglik (no constant): %lf\n", loglik);
      fprintf(polyfp, "loglik (with constant): %lf\n", loglik - log(2*M_PI)*sum*.5);
      fprintf(polyfp, "\n");
      variances.Print(polyfp);
      SEvariances.Print(polyfp);
      SEvariances_R.Print(polyfp);
      fprintf(polyfp, "\n");
      coef.Print(polyfp);
      SEcoef.Print(polyfp);
      SEcoef_R.Print(polyfp);
      covariance.Print(polyfp);
      CovCoef.Print(polyfp);
   }
   summary();
   if(polyFlag) {
      fclose(polyfp);
      printf("loglikelihoods and SEs are in %s\n", (const char*)logfile);
   }
}

int GEE::StopRule(void)
{
   int stop = 0;
   if(delta.SumSquares() < Epsilon*Epsilon*variances.Sum()*variances.Sum()) stop = 1;
   return stop;
}








/*
         for(int i = 0; i < parCount; i++)
            for(int j = i; j < parCount; j++)
               for(int u = 0; u < size; u++){
                  DVD.data[i]->data[j] += OD[i].data[u]->data[u] * OD[j].data[u]->data[u];
                  for(int v = u+1; v < size; v++)
                     DVD.data[i]->data[j] += 2 * OD[i].data[u]->data[v] * OD[j].data[v]->data[u];
               }
         for(int i = 0; i < parCount; i++)
               for(int u = 0; u < size; u++)
                  for(int v = 0; v < size; v++)
                     DVS_FAM.data[i]->data[f] += OD[i].data[u]->data[v] * OS.data[v]->data[u] * .5;
*/

/*
   // MLE stuff
   loglik = 0;
   for (int f = 0; f < ped.familyCount; f++){
      size = pheno[f].Length();
      if(size==0) continue;
      GetPhi(f);
      RefreshO(f);
      buf = traits[f];
      for(int i = 0; i < size; i++)
         for(int j = 0; j < coefCount; j++)
            buf[i] -= coef[j] * covariates[f][j][i];
      CholeskyOmega.Decompose(Omega);
      CholeskyOmega.BackSubst0(buf);
      for(int i = 0; i < size; i++) {
         loglik -= CholeskyOmega.x[i] * CholeskyOmega.x[i] / 2;
         loglik -= log(fabs(CholeskyOmega.L[i][i]));
      }
   } */







/*
Matrix GEE::BlockInverse(Matrix & M, int blockCount, int extra)
{
   int size = M.rows;
   Matrix W;
   if(extra==0){  // balanced
      if(blockCount==1) W.CholeskyInvert(M);
      else if(blockCount==2) W.Block2Invert(M);
      else //return BlockInverse(M, blockCount-1, size / blockCount);
         W.Block2Invert(M, size/blockCount*(blockCount-2));
   }else{
         W.Block2Invert(M, (size-extra)/blockCount*(blockCount-2)+extra);
   }
   return W;
} */
/*      Matrix W(size, size);
      Matrix B[5];
      int count = size - extra;
      B[0].Dimension(count, count);
      B[1].Dimension(count, extra);
      B[2].Dimension(extra, count);
      B[3].Dimension(extra, extra);

      for(int i = 0; i < count; i++)
         for(int j = 0; j < count; j++)
            B[0].data[i]->data[j] = M.data[i]->data[j];
      for(int i = 0; i < count; i++)
         for(int j = 0; j < extra; j++)
            B[1].data[i]->data[j] = B[2].data[j]->data[i] = M.data[i]->data[count+j];
      for(int i = 0; i < extra; i++)
         for(int j = 0; j < extra; j++)
            B[3].data[i]->data[j] = M.data[count+i]->data[count+j];
      B[0] = BlockInverse(B[0], blockCount, 0); // A^-1
      Matrix E, F;
      F.Product(B[0], B[1]);  // F = A^-1 B
      B[4].Product(B[2], F);
      B[1] = B[3];
      B[1].AddMultiple(-1, B[4]);  // E = D - B' A^-1 B
//      E = OmegaInverse(E); // 2, 2
      E.CholeskyInvert(B[1]);
      B[1].Product(F, E);         // 1, 2
      B[4].Transpose(F);
      B[2].Product(B[1], B[4]);
      B[0].Add(B[2]);               // 1, 1
      for(int i = 0; i < count; i++)
         for(int j = 0; j < count; j++)
            W.data[i]->data[j] = B[0].data[i]->data[j];
      for(int i = 0; i < count; i++)
         for(int j = 0; j < extra; j++)
            W.data[i]->data[count+j] = W.data[count+j]->data[i] = -B[1].data[i]->data[j];
      for(int i = 0; i < extra; i++)
         for(int j = 0; j < extra; j++)
            W.data[count+i]->data[count+j] = E.data[i]->data[j];
      return W;
*/
          /*
void GEE::Block2Inverse(Matrix & M)
{
   int size = M.rows;
   int count = size / 2;
   W2.Dimension(size, size);
   for(int b = 0; b < 5; b++)
      B[b].Dimension(count, count);
   for(int i = 0; i < count; i++)
      for(int j = 0; j < count; j++){
         B[0].data[i]->data[j] = M.data[i]->data[j];
         B[1].data[i]->data[j] = M.data[i]->data[count+j];
         B[2].data[i]->data[j] = M.data[count+i]->data[count+j];
      }
   B[3].Product(B[0], B[2]);
   B[4].Product(B[1], B[1]);
   B[3].AddMultiple(-1, B[4]);
   B[4] = OmegaInverse(B[3]);
   B[3].Product(B[2], B[4]);    // 1, 1
   B[2].Product(B[0], B[4]);    // 2, 2
   B[0].Product(B[1], B[4]);    // 1, 2
   for(int i = 0; i < count; i++)
      for(int j = 0; j < count; j++){
         W2.data[i]->data[j] = B[3].data[i]->data[j];
         W2.data[count+i]->data[count+j] = B[2].data[i]->data[j];
         W2.data[i]->data[count+j] = W2.data[count+j]->data[i] = -B[0].data[i]->data[j];
      }
}           */
/*
Matrix GEE::OmegaInverse(Matrix & M)
{
   int N = M.rows;
   if(M.cols != N) return NULL;
   if(N > 3){
      chol.Decompose(M);
      chol.Invert();
      return chol.inv;
   }else{
      Matrix W(N, N);
      double det;
      if(N==1) W[0][0] = 1/M[0][0];
      else if(N==2){
         det = M[0][0]*M[1][1] - M[0][1]*M[0][1];
         W[0][0] = M[1][1] / det;
         W[1][1] = M[0][0] / det;
         W[0][1] = W[1][0] = -M[0][1] / det;
      }else if(N==3){
         det = M[0][0] * M[1][1] * M[2][2]
         + 2 * M[0][1] * M[0][2] * M[1][2]
         - M[1][1] * M[0][2] * M[0][2]
         - M[0][0] * M[1][2] * M[1][2]
         - M[2][2] * M[0][1] * M[0][1];
         W[0][0] = (M[1][1]*M[2][2] - M[1][2]*M[1][2]) / det;
         W[1][1] = (M[0][0]*M[2][2] - M[0][2]*M[0][2]) / det;
         W[2][2] = (M[0][0]*M[1][1] - M[0][1]*M[0][1]) / det;
         W[1][0] = W[0][1] = (M[0][2]*M[1][2] - M[0][1]*M[2][2]) / det;
         W[2][0] = W[0][2] = (M[0][1]*M[1][2] - M[0][2]*M[1][1]) / det;
         W[2][1] = W[1][2] = (M[0][1]*M[0][2] - M[0][0]*M[1][2]) / det;
      }
      return W;
   }
} */
/*
void GEE::RefreshOD(int f)
{
   OmegaInv.CholeskyInvert(Omega);
   for(int i = 0; i < parCount; i++){
      Matrix MD(size, size);
      for(int u = 0; u < size; u++)
         for(int v = 0; v < size; v++)
            MD[u][v] = D[Index(v,u)][i];
      OD[i].Product(OmegaInv, MD);
   }
}   */
   /*
   Vector TempVector(size);
   for(int i = 0; i < parCount; i++)
      for(int u = 0; u < size; u++){      // column u
         for(int v = 0; v < size; v++)
            TempVector[v] = D[Index(v,u)][i];
            CholeskyOmega.BackSubst(TempVector); // most intensive computation
            for(int v = 0; v < size; v++)    // row v
               OD[i][v][u] = CholeskyOmega.x[v];
      }*/


