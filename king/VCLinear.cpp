//////////////////////////////////////////////////////////////////////
// VCLINEAR.cpp
// Author: Wei-Min Chen
// April 2, 2008

#include "VCLinear.h"
#include "Kinship.h"
#include "KinshipX.h"
#include "MathStats.h"
#include <math.h>

// ***************** GEEVC_LINEAR *********************
GEEVC_LINEAR::GEEVC_LINEAR(Pedigree & pedigree):GEE(pedigree)
{
   varianceComponents = NULL;
   personValid.Dimension(pedigree.count);
   personValid.Set(1);
   mCovariate.Dimension(0);
}

GEEVC_LINEAR::~GEEVC_LINEAR()
{
   if(varianceComponents) delete []varianceComponents;
}

double GEEVC_LINEAR::residual(int p)
{
   if(!ped[p].isPhenotyped(trait)) return _NAN_;
   for(int k = 0; k < mCovariate.Length(); k++)
      if(!ped[p].isControlled(mCovariate[k])) return _NAN_;
   double r = ped[p].traits[trait] - coef[0];
   for(int i = 1; i < coef.Length(); i++)
      r -= coef[i] * ped[p].covariates[mCovariate[i-1]];
   return r;
}
 
void GEEVC_LINEAR::init()
{}

void GEEVC_LINEAR::summary()
{
   totalVariance = variances.Sum();
   H2 = variances[1] / totalVariance;
//   seH2 = sqrt((1 - 2 * H2) * covariance[1][1] + H2*H2*(covariance[0][0]+covariance[1][1])) / totalVariance;
   borderIndex.Dimension(0);
   if(polyFlag){
      seH2 = sqrt(variances[1] * variances[1] * covariance[0][0]
         + variances[0] * variances[0] * covariance[1][1]
         - 2 * variances[0] * variances[1] * covariance[0][1]) / totalVariance / totalVariance;
      fprintf(polyfp, "\nH2: %.4f (%.4f)\n", H2, seH2);
   }
}

void GEEVC_LINEAR::print()
{}

void GEEVC_LINEAR::RefreshOD(int f)
{
   OmegaInv.CholeskyInvert(Omega);
   for(int i = 0; i < parCount; i++)
      OD[i].Product(OmegaInv, varianceComponents[i]);
}

void GEEVC_LINEAR::InitCoef()
{
   if(pheno==NULL) pheno = new IntArray[ped.familyCount];
   if(traits==NULL) traits = new Vector[ped.familyCount];
   if(covariates==NULL) covariates = new Matrix[ped.familyCount];
   isNuclear.Dimension(ped.familyCount);
   isNuclear.Zero();
   nuclearP1.Dimension(ped.familyCount);
   nuclearP1.Set(-1);
   nuclearP2.Dimension(ped.familyCount);
   nuclearP2.Set(-1);
   for(int f = 0; f < ped.familyCount; f++) {
      pheno[f].Dimension(0);
      for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
         int missing = 0;
         if(!personValid[i]) continue;
         for(int k = 0; k < mCovariate.Length(); k++){
            if(!ped[i].isControlled(mCovariate[k])) {
               missing=1;
               break;
            }
         }
         if(!missing && ped[i].isPhenotyped(trait))
            pheno[f].Push(i);
      }
      if(ped.families[f]->isNuclear()) {
         isNuclear[f] = 1;
         for(int i = 0; i < pheno[f].Length(); i++)
            if(ped[pheno[f][i]].father == NULL && ped[pheno[f][i]].mother == NULL){
               if(nuclearP1[f] == -1) nuclearP1[f] = i;
               else nuclearP2[f] = i;
            }
      }
//      if(isNuclear[f]) printf("%d %d %d\t", nuclearP1[f], nuclearP2[f], sibshipSize[f]);
   }

   if(saturatedMean) coef.Dimension(mCovariate.Length());
   else coef.Dimension(mCovariate.Length()+1);
   coefCount = coef.Length();
   ValidFamilies = ValidPersons = 0;
   for (int f = 0; f < ped.familyCount; f++){
      int size = pheno[f].Length();
      traits[f].Dimension(size);
      covariates[f].Dimension(coefCount, size);
      if(size == 0) continue;
      ValidFamilies++;
      ValidPersons += size;
      for(int u = 0; u < size; u++){
         traits[f][u] = ped[pheno[f][u]].traits[trait];
         if(saturatedMean)
            for(int i = 0; i < coefCount; i++)
               covariates[f][i][u] = ped[pheno[f][u]].covariates[mCovariate[i]];
         else{
            covariates[f][0][u] = 1;
            for(int i = 1; i < coefCount; i++)
               covariates[f][i][u] = ped[pheno[f][u]].covariates[mCovariate[i-1]];
         }
      }
   }
   parCount = variances.Length();
}

// *********************** POLY ***************************

void POLY::InitCoef()
{
   GEEVC_LINEAR::InitCoef();
//   variances.SetLabel("Variances");
//   coef.SetLabel("Coefficients");
   variances.Dimension(2);
   variances.Zero();
   parCount = variances.Length();
   if(varianceComponents==NULL) varianceComponents = new Matrix[parCount];

   double IVY, PVYY, IVI, PVP;
   IVY = PVYY = IVI = PVP = 0;
   Vector XVY(coefCount);
   XVY.Zero();
   Matrix XVX(coefCount, coefCount);
   XVX.Zero();

   if(saturatedMean){
      Vector means;
      for(int f = 0; f < ped.familyCount; f++){
         int size = pheno[f].Length();
         if(size == 0) continue;
         if(saturatedMean && size==1) continue;
         means.Dimension(coefCount+1);
         means.Zero();
         for(int i = 0; i < size; i++){
            means[coefCount] += traits[f][i];
            for(int j = 0; j < coefCount; j++)
               means[j] += covariates[f][j][i];
         }
         for(int j = 0; j < coefCount+1; j++)
            means[j] /= size;
         for(int u = 0; u < size; u++)
            for(int i = 0; i < coefCount; i++) {
               XVY[i] += (covariates[f][i][u]-means[i]) * (traits[f][u]-means[coefCount]);
               for(int j = 0; j < coefCount; j++)
                  XVX[i][j] += (covariates[f][i][u]-means[i]) * (covariates[f][j][u]-means[j]);
            }
      }
   }else
      for (int f = 0; f < ped.familyCount; f++){
         int size = pheno[f].Length();
         if(size == 0) continue;
         for(int u = 0; u < size; u++)
            for(int i = 0; i < coefCount; i++) {
               XVY[i] += covariates[f][i][u] * traits[f][u];
               for(int j = 0; j < coefCount; j++)
                  XVX[i][j] += covariates[f][i][u] * covariates[f][j][u];
            }
      }
   if(coefCount == 1) coef[0] = XVY[0] / XVX[0][0];
   else if(coefCount > 1){
      Matrix tempM;
      tempM.CholeskyInvert(XVX);
      for(int i = 0; i < coefCount; i++){
         coef[i] = 0;
         for(int u = 0; u < coefCount; u++)
            coef[i] += tempM[i][u] * XVY[u];
      }
   }
   if(saturatedMean)
      for (int f = 0; f < ped.familyCount; f++){
         int size = pheno[f].Length();
         if(size==0) continue;
         if(saturatedMean && size==1) continue;
         meanPerFamily[f] = 0;
         for(int i = 0; i < size; i++){
            meanPerFamily[f] += traits[f][i];
            for(int j = 0; j < coefCount; j++)
               meanPerFamily[f] -= coef[j] * covariates[f][j][i];
         }
         meanPerFamily[f] /= size;
      }
   Kinship kin;
   for(int f = 0; f < ped.familyCount; f++){
      int size = pheno[f].Length();
      if(size==0) continue;
      if(saturatedMean && size==1) continue;
      Matrix Phi;
      Phi.Dimension(size, size);
      kin.Setup(*ped.families[f]);
      for(int i = 0; i < size; i++)
         for(int j = i; j < size; j++)
            Phi[j][i] = Phi[i][j] = 2 * kin(ped[pheno[f][i]], ped[pheno[f][j]]);

      Vector buf;
      buf.Dimension(0);
      for(int i = 0; i < size; i++){
         double temp = traits[f][i];
         if(saturatedMean) temp -= meanPerFamily[f];
         for(int k = 0; k < coefCount; k++)
            temp -= coef[k] * covariates[f][k][i];
         buf.Push(temp);
      }
       for(int i = 0; i < size; i++){
         IVY += buf[i]*buf[i];
         IVI ++;
      }
      for(int u = 0; u < size; u++)
         for(int v = u+1; v < size; v++) {
            PVYY += Phi[u][v] * buf[u] * buf[v];
            PVP += Phi[u][v] * Phi[u][v];
         }
   }
   variances[0] = (IVY - IVI * PVYY / PVP) / IVI;
   if(variances[0]<0) variances[0] =  0.000001;
   variances[1] = PVYY / PVP;
   if(variances[1]<0) variances[1] =  0.000001;
   if(moreFlag) {variances.Print(); coef.Print();}
}

void POLY::RefreshO(int f)
{
   Omega.Dimension(size, size);
   for(int i = 0; i < size; i++){
      Omega[i][i] = variances[0] + variances[1];
      for(int j = i+1; j < size; j++)
         Omega[i][j] = Omega[j][i] = Phi[i][j]*variances[1];
   }
}

void POLY::RefreshOD(int f)
{
  if(isNuclear[f]){
      OmegaInv.Dimension(size, size);
      if(size == 1) {
         OmegaInv[0][0] = 1/Omega[0][0];
         OD[1].Dimension(size, size);
         OD[1][0][0] = OmegaInv[0][0];
      }else if(nuclearP2[f]==-1){
         double r = Omega.data[0]->data[size-1]/Omega.data[0]->data[0];
         double s = (1+(size-1)*r) * (1-r) * Omega.data[0]->data[0];
         double m1 = (1+(size-2)*r) / s;
         double m2 = -r / s;
         s = m1 + m2*(size-1)*.5;
         r = (m1 + m2*size)/2;
         OmegaInv.Set(m2);
         OD[1].Dimension(size, size);
         OD[1].Set(r);
         for(int i = 0; i < size; i++){
            OmegaInv.data[i]->data[i] = m1;
            OD[1].data[i]->data[i] = s;
         }
      }else{
         double r = Omega.data[0]->data[size-1]/Omega.data[0]->data[0];
         if(Omega.data[0]->data[size-1] == 0) r = Omega.data[0]->data[1] / Omega.data[0]->data[0];
         double s = (1+(size-3)*r-2*r*r*(size-2)) * Omega.data[0]->data[0];
         for(int i = 0; i < size; i++){
            if(i == nuclearP1[f] || i == nuclearP2[f]){
               OmegaInv.data[i]->data[i] = (1-r)*(1+(size-2)*r) / s;
               for(int j = i+1; j < size; j++)
                  if(j == nuclearP1[f] || j == nuclearP2[f])
                     OmegaInv.data[i]->data[j] = OmegaInv.data[j]->data[i] = r*r*(size-2)/ s;
                  else OmegaInv.data[i]->data[j] = OmegaInv.data[j]->data[i] = -r/s;
            }else{
               OmegaInv.data[i]->data[i] = ((size-3)*r*(1-2*r)/(1-r) + 1) / s;
               for(int j = i+1; j < size; j++)
                  if(j == nuclearP1[f] || j == nuclearP2[f])
                     OmegaInv.data[i]->data[j] = OmegaInv.data[j]->data[i] = -r/s;
                  else OmegaInv.data[i]->data[j] = OmegaInv.data[j]->data[i] = -r*(1-2*r)/(1-r)/s;
            }
        }
      }
      OD[1].sProductPhi(OmegaInv, Phi);   // intensive
   }else{
      OmegaInv.CholeskyInvert(Omega); // intensive
      OD[1].sProductPhi(OmegaInv, Phi);    // intensive
   }
   OD[0] = OmegaInv;
   /*
   OmegaInv.CholeskyInvert(Omega);
   OD[1].Product(OmegaInv, Phi);
   OD[0] = OmegaInv;
     */
}

int POLY::StopRule(void)
{
   int stop = 0;
   if(fabs(delta.Max()) < variances.Sum() * Epsilon) stop = 1;
   return stop;
}

// ************************** GEEVC_LINKAGE ******************************

void GEEVC_LINKAGE::InitCoef()
{
   GEEVC_LINEAR::InitCoef();
   if(varianceComponents==NULL) varianceComponents = new Matrix[parCount];
}

void GEEVC_LINKAGE::RefreshO(int f)
{
   int k;
   Omega.Dimension(size, size);

   for(int i = 0; i < parCount; i++)
      varianceComponents[i].Dimension(size, size);
   varianceComponents[0].Identity();
   varianceComponents[1] = Phi;

   varianceComponents[2].Identity();
   for(int i = 0; i < size; i++)
      for(int j = 0; j < i; j++)
         varianceComponents[2][i][j] =
            varianceComponents[2][j][i] =
               ibd[f][ped[pheno[f][i]].traverse * ped.families[f]->count + ped[pheno[f][j]].traverse];
   Omega.Zero();
   for(int i = 0; i < size; i++)
      for(int j = 0; j < size; j++)
         for(int k = 0; k < parCount; k++)
            Omega[i][j] += variances[k] * varianceComponents[k][i][j];
}

void GEEVC_LINKAGE::summary()
{
   GEEVC_LINEAR::summary();
   h2 = variances[2] / totalVariance;
}

// *************************** GEEVC_ASSOC *********************

void GEEVC_ASSOC::InitCoef()
{
   personValid.Set(1);
   for(int i = 0; i < ped.count; i++)
      if(IBS[i] == _NAN_){
         personValid[i] = 0;
//         printf("%d ", i);
      }
      
   GEEVC_LINKAGE::InitCoef();

//   coef.Dimension(ped.covariateCount+2);
   coef.Dimension(mCovariate.Length()+2);
   coef[coef.Length()-1] = 0.001;
   coefCount = coef.Length();
   for (int f = 0; f < ped.familyCount; f++){
      int size = pheno[f].Length();
      if(size == 0) continue;
      covariates[f].Dimension(coefCount, size);
      for(int u = 0; u < size; u++){
         covariates[f][0][u] = 1;
         for(int i = 1; i < coefCount-1; i++)
            covariates[f][i][u] = ped[pheno[f][u]].covariates[mCovariate[i-1]];
         for(int i = coefCount-1; i < coefCount; i++){
            covariates[f][i][u] = IBS[pheno[f][u]];
         }
      }
   }
}

void GEEVC_ASSOC::summary()
{
   GEEVC_LINKAGE::summary();
}


 /*
    int k;
   for(int i = 0; i < parCount; i++)
      varianceComponents[i].Dimension(size, size);
   varianceComponents[0].Identity();
   varianceComponents[1] = Phi;
   Omega.Zero();
   for(int i = 0; i < size; i++)
      for(int j = 0; j < size; j++)
         for(int k = 0; k < parCount; k++)
            Omega[i][j] += variances[k] * varianceComponents[k][i][j];
   */
//   int count1, count2; count1=count2=0;
/*      for(int i = 0; i < size; i++){
         variances[0] += buf[i]*buf[i];
         count1 ++;
      }
      for(int i = 0; i < size; i++)
         for(int j = i+1; j < size; j++){
            if(Phi[i][j] > 0)
               variances[1] += buf[i]*buf[j] / Phi[i][j];
            count2 ++;
         }       */
/*
   variances[1] /= count2;
   if(variances[1] < 0) variances[1] = .000001;
   variances[0] = variances[0] / count1 - variances[1];
   if(variances[0] < 0) variances[0] = .000001;  */

