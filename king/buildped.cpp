//////////////////////////////////////////////////////////////////////
// buildped.cpp
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

#include <math.h>
#include "analysis.h"
#include "Kinship.h"
#include "Intervals.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

bool Engine::SplitPedigree()
{
   Kinship kin;
   IntArray queue, splitLength;
   bool splitFlag = false;
   bool pedigreeFlag = false;
   String updateidsfile = prefix;
   FILE *fp;
   for(int f = 0; f < ped.familyCount; f++){
      Family *pf = ped.families[f];
      if(pf->count < 2) continue;
      if(!pedigreeFlag){
         updateidsfile.Add("splitped.txt");
         fp = fopen((const char*)updateidsfile, "wt");
         pedigreeFlag = true;
      }
      queue.Dimension(0);
      splitLength.Dimension(0);
      int p = 0;
      kin.Setup(*pf);
      for(int i = pf->first; i <= pf->last; i++){
         int exist = queue.Find(i);
         if(exist > -1) continue;
         for(queue.Push(i); p < queue.Length(); p++)
            for(int j = pf->first; j <= pf->last; j++)
               if(j!=queue[p] && queue.Find(j)==-1 && kin(ped[queue[p]], ped[j]) > 0.0000001)
                  queue.Push(j);
         splitLength.Push(queue.Length());
      }
      int splitCount = splitLength.Length();
      if(splitCount > 1){   // split family f
         for(int s = 0; s < splitCount; s++)
            for(int i = (s? splitLength[s-1]: 0); i < splitLength[s]; i++){
               int k = queue[i];
               fprintf(fp, "%s %s %s_S%d %s %s %s %d %d %d\n",
                  (const char*)pf->famid, (const char*)ped[k].pid,
                  (const char*)pf->famid, s+1, (const char*)ped[k].pid,
                  (const char*)ped[k].fatid, (const char*)ped[k].motid,
                  ped[k].sex, ped[k].affections[0], ped[k].ngeno>0? 0: 1);
            }
         splitFlag = true;
      }else for(int i = pf->first; i <= pf->last; i++)
            fprintf(fp, "%s %s %s %s %s %s %d %d %d\n",
               (const char*)pf->famid, (const char*)ped[i].pid,
               (const char*)pf->famid, (const char*)ped[i].pid,
               (const char*)ped[i].fatid, (const char*)ped[i].motid,
               ped[i].sex, ped.affectionCount? ped[i].affections[0]: 0, ped[i].ngeno>0? 0: 1);
   }
   if(pedigreeFlag) fclose(fp);
   bool valid = true;
//   if(splitFlag) printf("  Pedigrees are split into connected ones as shown in %s\n", (const char*)updateidsfile);
//   else if(pedigreeFlag) printf("  Pedigrees are reformated for plotting as shown in %s\n", (const char*)updateidsfile);
//   else {printf("  No pedigrees found in the dataset.\n"); valid = false;}
   if(!pedigreeFlag){printf("  No pedigrees found in the dataset.\n"); valid = false;}
   return valid;
}

int Engine::BuildOneFamily(int f, IntArray chrSeg, double totalLength, String & message)
{
   int idfCount = id[f].Length();
   if(idfCount < 2) return 0;
   L0.Dimension(0);
   Lfs.Dimension(0);
   Lpo.Dimension(0);
   L2.Dimension(0);
   Family *pf = ped.families[f];
   IntArray valid(pf->count); valid.Set(1);
   IntArray *poConnection = new IntArray[pf->count];
   IntArray *d2Connection = new IntArray[pf->count];
   IntArray *d123Connection = new IntArray[pf->count];
   for(int i = 0; i < pf->count; i++){
      poConnection[i].Dimension(0);
      d2Connection[i].Dimension(0);
      d123Connection[i].Dimension(0);
   }
   IntArray pairs(0);
   Matrix relationship(pf->count, pf->count); relationship.Zero();
   IntArray *pairIndex;
   Kinship kin; kin.Setup(*ped.families[f]);
   for(int i = 0; i < idfCount; i++)
      for(int j = i+1; j < idfCount; j++){
         if(kin(ped[id[f][i]], ped[id[f][j]]) > 0)
            relationship[id[f][i]-pf->first][id[f][j]-pf->first] =
            relationship[id[f][j]-pf->first][id[f][i]-pf->first] = kin(ped[id[f][i]], ped[id[f][j]]);
         pairs.Push(geno[id[f][i]]); pairs.Push(geno[id[f][j]]);
      }
   int pairCount = pairs.Length()/2;
   if(pairCount==0) return 0;
   IntArray HetHetCounts, IBS0Counts, het1Counts, het2Counts, HomHomCounts, IBSCounts;
   //****************************Compute kinship coefficient*****************************
   if(Bit64==64)
      KinshipInSubset64Bit(pairs, HetHetCounts, IBS0Counts, het1Counts, het2Counts, HomHomCounts, IBSCounts);
   else
      KinshipInSubset(pairs, HetHetCounts, IBS0Counts, het1Counts, het2Counts, HomHomCounts, IBSCounts);
   //****************************Compute kinship coefficient*****************************
   Vector ibdprops, maxLengths, ibd2props, maxLengths2;
   IntArray *ibd1segs = new IntArray [pairCount];
   bool IBDvalidFlag = chrSeg.Length()>0;
   if(IBDvalidFlag && Bit64==64){
      IBDSegInSubset64Bit(pairs, ibdprops, maxLengths, ibd2props, maxLengths2, ibd1segs);
      pairIndex = new IntArray[pf->count];
      for(int i = 0; i < pf->count; i++){
         pairIndex[i].Dimension(pf->count);
         pairIndex[i].Set(-1);
      }
   }
   int traverse[2];
   for(int p = 0; p < pairCount; p++){
      if(!het1Counts[p] && !het2Counts[p] && !HetHetCounts[p]) continue;
      for(int s = 0; s < 2; s++)
         traverse[s] = phenoid[pairs[p*2+s]] - pf->first;
      double CHet = HetHetCounts[p] * 1.0 / (HetHetCounts[p]+het1Counts[p]+het2Counts[p]);
      if(IBDvalidFlag && Bit64==64){
         double ibdprop = ibdprops[p];
         double ibd2prop = ibd2props[p];
         double pi = ibd2prop + ibdprop * 0.5;
         pairIndex[traverse[0]][traverse[1]] = pairIndex[traverse[1]][traverse[0]] = p;
         if(CHet<0.8){
            if(pi > 0.3535534){  // 1st-degree
               if(ibdprop + ibd2prop > 0.96 || (ibdprop + ibd2prop > 0.9 && ibd2prop <= 0.08)){
                  if(relationship[traverse[0]][traverse[1]] < 0.177){
                     Lpo.Push(traverse[0]); Lpo.Push(traverse[1]);
                  }
               }else if(((ibd2prop > 0.15 && ibd2prop < 0.5) || (pi > 0.4 && pi < 0.8)) && (relationship[traverse[0]][traverse[1]] < 0.177)){
                  Lfs.Push(traverse[0]); Lfs.Push(traverse[1]);
               }
            }else if(pi > 0.2 && ibd2prop < 0.08){  // very strict 2nd-degree
               if(relationship[traverse[0]][traverse[1]] < 0.08){
                  L2.Push(traverse[0]); L2.Push(traverse[1]);
               }
               d2Connection[traverse[0]].Push(traverse[1]);
               d2Connection[traverse[1]].Push(traverse[0]);
            }
            if(pi > 0.04419417){// closer than 4th-degree
               d123Connection[traverse[0]].Push(traverse[1]);
               d123Connection[traverse[1]].Push(traverse[0]);
            }
         }else{ // Duplicate
            L0.Push(traverse[0]); L0.Push(traverse[1]);
         }
         if(relationship[traverse[0]][traverse[1]] < 0.0442)
            relationship[traverse[0]][traverse[1]] = relationship[traverse[1]][traverse[0]] = pi * 0.5;
         else if(relationship[traverse[0]][traverse[1]] > 0.177 && pi * 0.5 < 0.177){
            printf("Warning: (%s %s) does not look like 1st-degree relatives.\n",
               (const char*)ped[traverse[0]+pf->first].pid,
               (const char*)ped[traverse[1]+pf->first].pid);
            printf("please fix within-family errors first before pedigree recontruction.\n");
            relationship[traverse[0]][traverse[1]] = relationship[traverse[1]][traverse[0]] = pi * 0.5;
         }
      }else{   // IBDseg unavailable
         double kinship = (HetHetCounts[p] - IBS0Counts[p]*2.0) / (HetHetCounts[p]*2+het1Counts[p]+het2Counts[p]);
         if(relationship[traverse[0]][traverse[1]] > 0.0221) continue;
         relationship[traverse[0]][traverse[1]] = relationship[traverse[1]][traverse[0]] = kinship;
         if(kinship <= 0.0221) continue;
         if(markerCount < 6000 && kinship <= 0.088388) continue;
         int deg = int(-log(kinship)/0.6931472-0.5);
         if(deg == 0) {L0.Push(traverse[0]); L0.Push(traverse[1]);}
         else if(deg == 1 && kinship > 0.18){
            int IBS0Count = IBS0Counts[p];
            int notMissingCount = het1Counts[p] + het2Counts[p] + HomHomCounts[p] + HetHetCounts[p];
            if(IBS0Count > errorrateCutoff * notMissingCount) {
               Lfs.Push(traverse[0]); Lfs.Push(traverse[1]);
            }else{   // PO
               Lpo.Push(traverse[0]); Lpo.Push(traverse[1]);
            }
         }else if(deg == 2){
            L2.Push(traverse[0]); L2.Push(traverse[1]);
            d2Connection[traverse[0]].Push(traverse[1]);
            d2Connection[traverse[1]].Push(traverse[0]);
            d123Connection[traverse[0]].Push(traverse[1]);
            d123Connection[traverse[1]].Push(traverse[0]);
         }
      }
   }  // end of pairs
   int L0Count = L0.Length()/2;
   int fsCount = Lfs.Length()/2;
   int poCount = Lpo.Length()/2;
   if(!L0Count && !poCount && !fsCount) return 0;
   char buffer[256];
   message.Clear();
   sprintf(buffer, "Family %s:\n", (const char*)pf->famid);
   message += buffer;
   int smaller, larger, ids[2];

   if(L0Count){   // If dup exists, remove one of them
      IntArray childCount(pf->count); childCount.Zero();
      for(int i = pf->first; i <= pf->last; i++){
         if(ped[i].father) childCount[ped[i].father->serial-pf->first] ++;
         if(ped[i].mother) childCount[ped[i].mother->serial-pf->first] ++;
      }
      for(int j = 0; j < L0Count; j++){
         for(int s = 0; s < 2; s++) ids[s] = L0[j*2+s];
         if(!valid[ids[0]] || !valid[ids[1]]) continue;
         if(childCount[ids[0]] < childCount[ids[1]]){
            smaller = ids[0]; larger = ids[1];
         }else if(childCount[ids[1]] < childCount[ids[0]]){
            smaller = ids[1]; larger = ids[0];
         }else if(ped[ids[0]+pf->first].sibCount < ped[ids[1]+pf->first].sibCount){
            smaller = ids[0]; larger = ids[1];
         }else if(ped[ids[1]+pf->first].sibCount < ped[ids[0]+pf->first].sibCount){
            smaller = ids[1]; larger = ids[0];
         }else{smaller = ids[0]; larger = ids[1];}
         valid[smaller] = 0;  // Dup with smaller family size to be removed
         Person *plarger = &ped[larger+pf->first];
         Person *psmaller = &ped[smaller+pf->first];
         sprintf(buffer, "  Duplicate %s (of %s) is removed.\n",
            (const char*)psmaller->pid, (const char*)plarger->pid);
         message += buffer;
         if(childCount[smaller])// All smaller's children now have larger as parent
            for(int i = pf->first; i <= pf->last; i++){
               if(valid[i-pf->first] && psmaller->sex == 1 && ped[i].father && ped[i].father->serial==smaller+pf->first){
                  ped[i].father = &ped[larger+pf->first];
                  ped[i].fatid = plarger->pid;
               }else if(valid[i-pf->first] && psmaller->sex == 2 && ped[i].mother && ped[i].mother->serial==smaller+pf->first){
                  ped[i].mother = &ped[larger+pf->first];
                  ped[i].motid = plarger->pid;
               }
            }
         if(!plarger->father && !plarger->mother && psmaller->father && psmaller->mother){
            plarger->father = psmaller->father;
            plarger->fatid = psmaller->fatid;
            plarger->mother = psmaller->mother;
            plarger->motid = psmaller->motid;
         }
      }  // End of j for L0 loop
   }  // End of if L0 not empty

   int ss[2], pa[2], parentid[2][2];
   Person *fss[2];
   bool parentvalid[2][2];
   int newparent = 0;
   String tempS;
   // Reconstruct FS
   IntArray parents[2];
   for(int s = 0; s < 2; s++) parents[s].Dimension(0);
   IntArray sibship[10000];
   int sibshipCount = 0;   // Hopefully sibshipCount is never > 10000
   for(int i = pf->first; i <= pf->last; i++)
      if(ped[i].sibCount > 1 && ped[i].sibs[0]->serial == i){
         sibship[sibshipCount].Dimension(0);
         for(int s = 0; s < ped[i].sibCount; s++)
            if(valid[ped[i].sibs[s]->serial-pf->first])
               sibship[sibshipCount].Push(ped[i].sibs[s]->serial-pf->first);
         ped[i].sibCount = sibship[sibshipCount].Length();
         if(ped[i].sibCount < 2) continue;
         parents[0].Push(ped[i].father->serial-pf->first);
         parents[1].Push(ped[i].mother->serial-pf->first);
         sibshipCount++;
      }
   for(int j = 0; j < fsCount; j++){
      for(int s = 0; s < 2; s++) ids[s] = Lfs[j*2+s];
      if(valid[ids[0]] == 0 || valid[ids[1]] == 0) continue;
      for(int s = 0; s < 2; s++) {
         ss[s] = -1;
         for(int k = 0; k < sibshipCount; k++)
            if(sibship[k].Find(ids[s]) > -1) ss[s] = k;
         pa[s] = -1;
         fss[s] = &ped[ids[s] + pf->first];
         parentvalid[s][0] = (fss[s]->father && valid[fss[s]->father->serial-pf->first]);
         parentvalid[s][1] = (fss[s]->mother && valid[fss[s]->mother->serial-pf->first]);
         parentid[s][0] = fss[s]->father? fss[s]->father->serial: -1;
         parentid[s][1] = fss[s]->mother? fss[s]->mother->serial: -1;
      }
      for(int p = 0; p < 2; p++)// get ID for parents
         if(parentvalid[0][p] && parentvalid[1][p]){
            if((p==0 && fss[0]->father->ngeno > fss[1]->father->ngeno) ||
               (p==1 && fss[0]->mother->ngeno > fss[1]->mother->ngeno))
               pa[p] = parentid[0][p];
            else
               pa[p] = parentid[1][p];
         }else if(parentvalid[0][p] || parentvalid[1][p])
            pa[p] = parentid[parentvalid[0][p]? 0: 1][p];
      if(ss[0] > -1 && ss[1] > -1){ // RULE FS2
         if(ss[0] == ss[1]) continue;  // pass if belong to the same sibship
         else if(ss[0] < ss[1]) {smaller=ss[0]; larger=ss[1];}
         else{smaller=ss[1]; larger=ss[0];}  // combine two sibships
         sprintf(buffer, "  RULE FS2: Sibship (%s", (const char*)ped[sibship[ss[0]][0]+pf->first].pid);
         message += buffer;
         for(int s = 1; s < sibship[ss[0]].Length(); s++){
            sprintf(buffer, " %s", (const char*)ped[sibship[ss[0]][s]+pf->first].pid);
            message += buffer;
         }
         sprintf(buffer, ") and sibship (%s", (const char*)ped[sibship[ss[1]][0]+pf->first].pid);
         message += buffer;
         for(int s = 1; s < sibship[ss[1]].Length(); s++){
            sprintf(buffer, " %s", (const char*)ped[sibship[ss[1]][s]+pf->first].pid);
            message += buffer;
         }
         sprintf(buffer, ") are combined\n");
         message += buffer;
         sibship[smaller].Stack(sibship[larger]);
         for(int s = 0; s < 2; s++)
            if(pa[s] > -1) parents[smaller] = pa[s] - pf->first;
         if(larger < sibshipCount-1){
            sibship[larger] = sibship[sibshipCount-1];
            for(int s = 0; s < 2; s++)
               parents[s][larger] = parents[s][sibshipCount-1];
         }
         sibshipCount--;
         parents[0].Delete(sibshipCount);
         parents[1].Delete(sibshipCount);
      }else if(ss[0] == -1 && ss[1] == -1){  // RULE FS0
         sibship[sibshipCount].Dimension(0);
         for(int s = 0; s < 2; s++)
            sibship[sibshipCount].Push(ids[s]);
         sibshipCount++;
         sprintf(buffer, "  RULE FS0: Sibship (%s", (const char*)ped[sibship[sibshipCount-1][0]+pf->first].pid);
         message += buffer;
         for(int s = 1; s < sibship[sibshipCount-1].Length(); s++){
            sprintf(buffer, " %s", (const char*)ped[sibship[sibshipCount-1][s]+pf->first].pid);
            message += buffer;
         }
         sprintf(buffer, ")'s parents are (");
         message += buffer;
         for(int s = 0; s < 2; s++)
            if(pa[s] > -1){
               parents[s].Push(pa[s] - pf->first);
               sprintf(buffer, s==0?"%s":" %s)\n", (const char*)ped[pa[s]].pid);
               message += buffer;
            }else{
               parents[s].Push(missingBase + 1000000 + newparent);
               inclusionList[0].Push(ped.families[f]->famid);
               tempS = missingBase + newparent;
               sprintf(buffer, s==0?"%s":" %s)\n", (const char*)tempS);
               message += buffer;
               newparent++;
               inclusionList[1].Push(tempS);
               tempS = s+1;
               inclusionList[2].Push(tempS);
               inclusionList[3].Push("");
               inclusionList[4].Push("");
            }
      }else{   // RULE FS1
         int ls = ss[0] > -1? 0: 1; // 1-ls join ls's sibship
         sibship[ss[ls]].Push(ids[1-ls]);
         for(int s = 0; s < 2; s++)
            if(pa[s] > -1) parents[ls] = pa[s] - pf->first;
         sprintf(buffer, "  RULE FS1: %s joins in sibship (%s",
            (const char*)fss[1-ls]->pid, (const char*)ped[sibship[ss[ls]][0]+pf->first].pid);
         message += buffer;
         for(int s = 1; s < sibship[ss[ls]].Length()-1; s++){
            sprintf(buffer, " %s", (const char*)ped[sibship[ss[ls]][s]+pf->first].pid);
            message += buffer;
         }
         sprintf(buffer, ")\n");
         message += buffer;
      }  // End of if rules
   }  // End of j for FS loop
   for(int j = 0; j < sibshipCount; j++){// KING -> newparent
      if(parents[0][j] < missingBase+1000000 && parents[1][j] < missingBase+1000000 && parents[0][j] > -1 && parents[1][j] > -1)
         for(int s = 0; s < 2; s++)
            if(ped[parents[s][j]+pf->first].pid.SubStr(0, 4)=="KING" && ped[parents[s][j]+pf->first].ngeno==0){
               ped[parents[s][j]+pf->first].pid = missingBase + newparent;
               newparent ++;
            }
      for(int k = 0; k < sibship[j].Length(); k++){
         Person *pj = &ped[sibship[j][k]+pf->first];
         if(parents[0][j] < missingBase+1000000 && parents[1][j] < missingBase+1000000 && parents[0][j] > -1 && parents[1][j] > -1){
            pj->father = &ped[parents[0][j]+pf->first];
            pj->mother = &ped[parents[1][j]+pf->first];
            pj->fatid = ped[parents[0][j]+pf->first].pid;
            pj->motid = ped[parents[1][j]+pf->first].pid;
         }else{
            pj->fatid = parents[0][j]-1000000;
            pj->motid = parents[1][j]-1000000;
         }
      }  // End of k for sibs loop
   }  // End of j for sibship loop

   // Reconstruct parent-offspring
   Person *pos[2];
   int parent, offspring;
   IntArray ibd1seg, tempArray, tempArray2, relatives;
   StringArray missing_id(0);
   for(int j = 0; j < poCount; j++){
      for(int s = 0; s < 2; s++) {
         ids[s] = Lpo[j*2+s];
         pos[s] = &ped[ids[s]+pf->first];
         ss[s] = -1;
         for(int k = 0; k < sibshipCount; k++)
            if(sibship[k].Find(ids[s]) > -1) ss[s] = k;
      }
      if(valid[ids[0]] == 0 || valid[ids[1]] == 0) continue;
      if(pos[0]->fatid == pos[1]->pid || pos[0]->motid == pos[1]->pid ||
         pos[1]->fatid == pos[0]->pid || pos[1]->motid == pos[0]->pid) continue;
      sprintf(buffer, "  Reconstruct parent-offspring pair (%s, %s)...\n",
         (const char*)pos[0]->pid, (const char*)pos[1]->pid);
      message += buffer;
      parent = offspring = -1;
      int sibshiplist = -1;
      String rule("");
      if(ss[0] == -1 && ss[1] == -1){ // singletons
         for(int s = 0; s < 2; s++){
            parentvalid[s][0] = (pos[s]->father && valid[pos[s]->father->serial-pf->first] && pos[s]->father->ngeno >= MINSNPCOUNT);
            parentvalid[s][1] = (pos[s]->mother && valid[pos[s]->mother->serial-pf->first] && pos[s]->mother->ngeno >= MINSNPCOUNT);
         }
         if(parentvalid[0][0] || parentvalid[0][1] || parentvalid[1][0] || parentvalid[1][1]){
            for(int s = 0; s < 2; s++) // s: one or the other of the PO pair
               for(int t = 0; t < 2; t++){   // t: one's parent
                  if(!parentvalid[s][t]) continue;
                  if(pos[1-s]->sex == t+1){  // The other has the same sex as the parent
                     parent = ids[s];  // then the other cannot be the parent
                     offspring = ids[1-s];
                     rule = "PO.PA.NA";
                  }else{   // the other has differnt or unknown sex
                     int pa_traverse = (t?pos[s]->mother:pos[s]->father)->serial - pf->first;
                     if(relationship[ids[1-s]][pa_traverse] < 0.0625){   // The other is unrelated to the parent
                        parent = ids[1-s];   // then the other is a parent
                        offspring = ids[s];
                        rule="PO.PA.UN";
                     }else if(relationship[ids[1-s]][pa_traverse] > 0.0884){ // The other is D2 to the parent
                        parent = ids[s];  // then the other is an offspring
                        offspring = ids[1-s];
                        rule="PO.PA.D2";
                     }
                  }
               }  // End of loop s,t for pair & its parent
         }else{   // none of the PO has known parents
            for(int s = 0; s < 2; s++) // RULE PO.D2
               for(int k = 0; k < d123Connection[ids[s]].Length(); k++)
                  if(IBDvalidFlag && Bit64==64){
                     int relative = d123Connection[ids[s]][k];
                     if(relative == ids[1-s]) continue;
                     tempArray = ibd1segs[pairIndex[ids[s]][relative]]; // RO
                     double tempLength = SegmentLength(tempArray, bp);
                     double joinLength = JoinLength(ibd1segs[pairIndex[ids[1-s]][relative]], tempArray, bp);   // RP & RO
                     if(joinLength < tempLength * 0.125){ // Would be 0.5 or 1 if related
                        parent = ids[1-s]; // One's relative unrelated to the other
                        offspring = ids[s];  // then the other is the parent
                        rule="PO.REL.UN";
                        sprintf(buffer, "    (P=%s, O=%s, R=%s): PR|OR = %.1lfMb / %.1lfMb\n",
                           (const char*)ped[parent+pf->first].pid,
                           (const char*)ped[offspring+pf->first].pid,
                           (const char*)ped[relative+pf->first].pid,
                           joinLength*0.000001, tempLength*0.000001);
                        message += buffer;
                        if(relationship[ids[s]][relative] > 0.177){  // relative is the other parent
                           if(ped[parent+pf->first].sex == 1){// the other parent is mother
                              ped[offspring+pf->first].motid = ped[relative+pf->first].pid;
                              ped[offspring+pf->first].mother = &ped[relative+pf->first];
                              sprintf(buffer, "    %s is now mother of %s\n",
                                 (const char*)ped[relative+pf->first].pid, (const char*)ped[offspring+pf->first].pid);
                           }else{   // the other parent is father
                              ped[offspring+pf->first].fatid = ped[relative+pf->first].pid;
                              ped[offspring+pf->first].father = &ped[relative+pf->first];
                              sprintf(buffer, "    %s is now father of %s\n",
                                 (const char*)ped[relative+pf->first].pid, (const char*)ped[offspring+pf->first].pid);
                           }
                           message += buffer;
                        }
                        break;
                     }
                  }else if(relationship[ids[1-s]][d2Connection[ids[s]][k]] < 0.03125){// integrated inference unavailable
                     parent = ids[1-s]; // One's relative unrelated to the other
                     offspring = ids[s];  // then the other is the parent
                     rule="PO.D2.UN";
                  }
            if(parent==-1 && IBDvalidFlag && Bit64==64){ // Apply rule PO.RELS.UN
               for(int s = 0; s < 2; s++){
                  relatives.Dimension(0);
                  for(int k = 0; k < d123Connection[ids[s]].Length(); k++){ // RULE PO.D2
                     int relative = d123Connection[ids[s]][k];
                     if(relative == ids[1-s]) continue;
                     tempArray = ibd1segs[pairIndex[ids[1-s]][relative]];  // OR
                     double tempLength = SegmentLength(tempArray, bp);   // OR length in Mb
                     if(tempLength > 10000000){  // OR Length > 10Mb
                        double joinLength = JoinLength(ibd1segs[pairIndex[ids[s]][relative]], tempArray, bp);
                        if(joinLength > tempLength * 0.8 &&  // ids[s] is the parent
                        joinLength > SegmentLength(ibd1segs[pairIndex[ids[s]][relative]], bp) * 0.125)  // not Unrelated
                           relatives.Push(relative);
                     }
                  }
                  int relativeCount = relatives.Length();
                  if(relativeCount < 2) continue;
                  for(int k1 = 0; k1 < relativeCount; k1++){
                     int relative = relatives[k1];
                     tempArray = ibd1segs[pairIndex[ids[s]][relative]]; // PR1
                     for(int k2 = k1+1; k2 < relativeCount; k2++){
                        int relative2 = relatives[k2];
                        double joinRatio = RoRP(tempArray,
                           ibd1segs[pairIndex[ids[s]][relative2]],
                           ibd1segs[pairIndex[relative][relative2]], bp);
                        if(joinRatio < 0.125 && joinRatio >= 0){   // R1 R2 unrelated
                           parent = ids[s];
                           offspring = ids[1-s];
                           rule="PO.RELS.UN";
                           sprintf(buffer, "    (P=%s, O=%s, R1=%s, R2=%s): R1R2 relatedness = %.2G\n",
                              (const char*)ped[parent+pf->first].pid,
                              (const char*)ped[offspring+pf->first].pid,
                              (const char*)ped[relative+pf->first].pid,
                              (const char*)ped[relative2+pf->first].pid,
                              joinRatio);
                           message += buffer;
                        }
                     }  // End of k2 for the k2-th relative
                  }  // End of k1 for the k1-th relative
               }  // End of s for either of the pair order
            }  // End of if condition for PO.RELS.UN
            if(!(IBDvalidFlag && Bit64==64) && offspring>-1){
               int len = offspring>-1? poConnection[offspring].Length():0;
               int parentindex = -1;
               for(int i = 0; i < len; i++)
                  if(relationship[parent][poConnection[offspring][i]] < 0.03125){
                     parentindex = i;
                     break;
                  }
               if(parentindex > -1){
                  int otherparent = poConnection[offspring][parentindex];
                  if(ped[parent+pf->first].sex == 1){// the other parent is mother
                     ped[offspring+pf->first].motid = ped[otherparent+pf->first].pid;
                     ped[offspring+pf->first].mother = &ped[otherparent+pf->first];
                     sprintf(buffer, "    %s is now mother of %s\n",
                        (const char*)ped[otherparent+pf->first].pid, (const char*)ped[offspring+pf->first].pid);
                  }else{   // the other parent is father
                     ped[offspring+pf->first].fatid = ped[otherparent+pf->first].pid;
                     ped[offspring+pf->first].father = &ped[otherparent+pf->first];
                     sprintf(buffer, "    %s is now father of %s\n",
                        (const char*)ped[otherparent+pf->first].pid, (const char*)ped[offspring+pf->first].pid);
                  }
                  message += buffer;
               }
            }
            if(parent==-1){
               if(cAge == -1 || pos[0]->covariates[cAge] == _NAN_ || pos[1]->covariates[cAge] == _NAN_){
                  poConnection[ids[0]].Push(ids[1]);
                  poConnection[ids[1]].Push(ids[0]);
               }else{
                  smaller = pos[0]->covariates[cAge] < pos[1]->covariates[cAge]? 0: 1;
                  larger = 1-smaller;
                  if(pos[larger]->covariates[cAge] > pos[smaller]->covariates[cAge]+10){
                     parent = ids[larger];
                     offspring = ids[smaller];
                     sprintf(buffer, "  Age information is used\n");
                     message += buffer;
                  }else{// Whoever older is the parent if age known
                     printf("    Warning: parent-offspring age difference < 10.\n");
                     continue;
                  }
               }
            }
         }
      }else if(ss[0]>-1 || ss[1]>-1){  // one or both have sibs
         rule = "PO.S";
         for(int s = 0; s < 2; s++){
            if(ss[s] == -1) continue;  // now s has sibs
            int r = sibship[ss[s]].Find(ids[s]);
            if(relationship[ids[1-s]][sibship[ss[s]][r==0?1:0]] > 0.2){   // The other is a parent
               parent = ids[1-s];
               sibshiplist = ss[s];
            }else{   // The other is second degree
               parent = ids[s];
               offspring = ids[1-s];
            }
            sprintf(buffer, "  %s's sibship is used to determine the parent/offspring\n",
               (const char*)pos[s]->pid);
            message += buffer;
            break;
         }
      }
      if(parent == -1) {   // parent/offspring order still unknown
         if(cAge == -1 || pos[0]->covariates[cAge] == _NAN_ || pos[1]->covariates[cAge] == _NAN_) continue;
         smaller = pos[0]->covariates[cAge] < pos[1]->covariates[cAge]? 0: 1;
         larger = 1 - smaller;
         if(pos[larger]->covariates[cAge] > pos[smaller]->covariates[cAge]+10){
            parent = ids[larger]; // Whoever older is the parent if age known
            offspring = ids[smaller];
            sprintf(buffer, "  Age information is used\n");
            message += buffer;
         }else{
            printf("    Warning: parent-offspring age difference < 10.\n");
            continue;
         }
      }
      int parentsex = ped[parent+pf->first].sex;
      if(parentsex != 1 && parentsex != 2) continue;
      String string_pa = parentsex==1? "father": "mother";
      if(offspring > -1){
         tempS = parentsex==1? ped[offspring+pf->first].fatid: ped[offspring+pf->first].motid;
         if(tempS == "0")
            sprintf(buffer, "  RULE %s: %s is now %s of %s\n",
               (const char*)rule, (const char*)ped[parent+pf->first].pid,
               (const char*)string_pa, (const char*)ped[offspring+pf->first].pid);
         else
            sprintf(buffer, "  RULE %s: %s is now %s of %s, replacing %s (%d genotypes)\n",
               (const char*)rule, (const char*)ped[parent+pf->first].pid,
               (const char*)string_pa, (const char*)ped[offspring+pf->first].pid, 
               (const char*)tempS,
               (parentsex == 1)? (ped[offspring+pf->first].father ?
               ped[ped[offspring+pf->first].father->serial].ngeno: 0):
               (ped[offspring+pf->first].mother? ped[ped[offspring+pf->first].mother->serial].ngeno: 0));
         message += buffer;
         if(parentsex==1){
            ped[offspring+pf->first].fatid = ped[parent+pf->first].pid;
            ped[offspring+pf->first].father = &ped[parent+pf->first];
         }else{
            ped[offspring+pf->first].motid = ped[parent+pf->first].pid;
            ped[offspring+pf->first].mother = &ped[parent+pf->first];
         }
         if((parentsex==1 && ped[offspring+pf->first].motid == "0") ||
            (parentsex==2 && ped[offspring+pf->first].fatid == "0")){
               inclusionList[0].Push(ped.families[f]->famid);
               tempS = missingBase + newparent;
               if(parentsex==1){
                  sprintf(buffer, "    %s is created as %s's mother.\n",
                     (const char*)tempS, (const char*)ped[offspring+pf->first].pid);
                  ped[offspring+pf->first].motid = tempS;
               }else{
                  sprintf(buffer, "    %s is created as %s's father.\n",
                     (const char*)tempS, (const char*)ped[offspring+pf->first].pid);
                  ped[offspring+pf->first].fatid = tempS;
               }
               message += buffer;
               newparent++;
               inclusionList[1].Push(tempS);
               tempS = 3-parentsex;
               inclusionList[2].Push(tempS);
               inclusionList[3].Push("");
               inclusionList[4].Push("");
         }
      }else{   // all sibs are offspring
         tempS = (parentsex==1? ped[sibship[sibshiplist][0]+pf->first].fatid:
            ped[sibship[sibshiplist][0]+pf->first].motid);
         if(tempS == "0")
            sprintf(buffer, "  RULE %s: %s is now %s of %s's sibship\n",
            (const char*)rule, (const char*)ped[parent+pf->first].pid,
            (const char*)string_pa, (const char*)ped[sibship[sibshiplist][0]+pf->first].pid);
         else
            sprintf(buffer, "  RULE %s: %s is now %s of %s's sibship, replacing %s (%d genotypes)\n",
            (const char*)rule, (const char*)ped[parent+pf->first].pid, (const char*)string_pa,
            (const char*)ped[sibship[sibshiplist][0]+pf->first].pid, (const char*)tempS,
            (parentsex == 1)? (ped[sibship[sibshiplist][0]+pf->first].father ?
               ped[ped[sibship[sibshiplist][0]+pf->first].father->serial].ngeno: 0):
               (ped[sibship[sibshiplist][0]+pf->first].mother? ped[ped[sibship[sibshiplist][0]+pf->first].mother->serial].ngeno: 0));
         message += buffer;
         for(int k = 0; k < sibship[sibshiplist].Length(); k++)
            if(parentsex==1){
               ped[sibship[sibshiplist][k]+pf->first].fatid = ped[parent+pf->first].pid;
               ped[sibship[sibshiplist][k]+pf->first].father = &ped[parent+pf->first];
            }else{
               ped[sibship[sibshiplist][k]+pf->first].motid = ped[parent+pf->first].pid;
               ped[sibship[sibshiplist][k]+pf->first].mother = &ped[parent+pf->first];
            }
      }
   }  // End of j loop for PO pairs
   for(int i = 0; i < pf->count; i++){
      int couple1 = -1; int couple2 = -1;
      for(int j = 0; j < poConnection[i].Length(); j++)
         for(int k = j+1; k < poConnection[i].Length(); k++)
            if(relationship[poConnection[i][j]][poConnection[i][k]] < 0.0625){// two PO are unrelated
               if(ped[poConnection[i][j]+pf->first].sex==1 || ped[poConnection[i][k]+pf->first].sex==2){
                  couple1 = poConnection[i][j];
                  couple2 = poConnection[i][k];
               }else if(ped[poConnection[i][j]+pf->first].sex==2 || ped[poConnection[i][k]+pf->first].sex==1){
                  couple1 = poConnection[i][k];
                  couple2 = poConnection[i][j];
               }
            }
      if(couple1 > -1){ // two parents identified
         ped[i+pf->first].fatid = ped[couple1+pf->first].pid;
         ped[i+pf->first].motid = ped[couple2+pf->first].pid;
         sprintf(buffer, "  RULE PO.2P: %s and %s are %s's parents\n",
            (const char*)ped[i+pf->first].fatid,
            (const char*)ped[i+pf->first].motid,
            (const char*)ped[i+pf->first].pid);
         message += buffer;
         for(int k = 0; k < poConnection[i].Length(); k++){
            if(poConnection[i][k] == couple1 || poConnection[i][k] == couple2) continue;
            if(ped[i+pf->first].sex == 1)
               ped[poConnection[i][k] + pf->first].fatid = i;
            else
               ped[poConnection[i][k] + pf->first].motid = i;
         }
      }else if(d2Connection[i].Length()){
         for(int k = 0; k < poConnection[i].Length(); k++)
            if(relationship[d2Connection[i][0]][poConnection[i][k]] > 0.03125 &&
               relationship[d2Connection[i][0]][poConnection[i][k]] < 0.0884){
               parent = i; // One's relative is D3 with the other
               offspring = poConnection[i][k];  // then the other is offspring
               int sex = ped[i+pf->first].sex;
               if(sex != 1 && sex != 2) continue;
               if(sex==1){
                  ped[offspring + pf->first].fatid = ped[parent+pf->first].pid;
                  sprintf(buffer, "  RULE PO.D2.D3: %s is now father of %s\n",
                     (const char*)ped[parent+pf->first].pid,
                     (const char*)ped[offspring+pf->first].pid);
               }else{
                  ped[offspring + pf->first].motid = ped[parent+pf->first].pid;
                  sprintf(buffer, "  RULE PO.D2.D3: %s is now mother of %s\n",
                     (const char*)ped[parent+pf->first].pid,
                     (const char*)ped[offspring+pf->first].pid);
               }
               message += buffer;
               if((sex==1 && ped[poConnection[i][k] + pf->first].motid=="0") ||
                  (sex==2 && ped[poConnection[i][k] + pf->first].fatid=="0") ){
                     inclusionList[0].Push(ped.families[f]->famid);
                     tempS = missingBase + newparent;
                     if(sex==1){
                        sprintf(buffer, "    %s is created as %s's mother.\n",
                           (const char*)tempS, (const char*)ped[offspring+pf->first].pid);
                        ped[offspring+pf->first].motid = tempS;
                     }else{
                        sprintf(buffer, "    %s is created as %s's father.\n",
                           (const char*)tempS, (const char*)ped[offspring+pf->first].pid);
                        ped[offspring+pf->first].fatid = tempS;
                     }
                     message += buffer;
                     newparent++;
                     inclusionList[1].Push(tempS);
                     tempS = 3-sex;
                     inclusionList[2].Push(tempS);
                     inclusionList[3].Push("");
                     inclusionList[4].Push("");
                  }
            }
      }  // End of if there are two parents or D2
   }  // End of loop over each person

   if(IBDvalidFlag && Bit64==64){   // Check 2nd-degree
      sprintf(buffer, "\n");
      message += buffer;
      IntArray relativeType[3];
      String stringType[3] = {"HS", "AV", "GG"};
      IntArray unrelatedTo[2];
      int L2Count = L2.Length()/2;
      for(int j = 0; j < L2Count; j++){
         for(int s = 0; s < 2; s++) ids[s] = L2[j*2+s];
         if(valid[ids[0]] == 0 || valid[ids[1]] == 0) continue;
         tempArray = ibd1segs[pairIndex[ids[0]][ids[1]]]; //S1S2
         for(int s = 0; s < 2; s++) unrelatedTo[s].Dimension(0);
         for(int s = 0; s < 2; s++)
            for(int k = 0; k < d123Connection[ids[s]].Length(); k++){
               int relative = d123Connection[ids[s]][k];
               if(relative == ids[1-s]) continue;
               double joinRatio = RoRP(tempArray,
                  ibd1segs[pairIndex[ids[s]][relative]],
                  ibd1segs[pairIndex[ids[1-s]][relative]], bp);
               if(joinRatio < 0.125 && joinRatio >= 0) unrelatedTo[1-s].Push(relative);
            }
         if(unrelatedTo[0].Length() && unrelatedTo[1].Length()){  // D2.HS.UN2
            for(int s = 0; s < 2; s++){
               sprintf(buffer, "    HS %s unrelated to", (const char*)ped[ids[1-s]+pf->first].pid);
               message += buffer;
               for(int k = 0; k < unrelatedTo[1-s].Length(); k++){
                  sprintf(buffer, " %s", (const char*)ped[unrelatedTo[1-s][k]+pf->first].pid);
                  message += buffer;
               }
               sprintf(buffer, "\n");
               message += buffer;
            }
            sprintf(buffer, "  INFERENCE D2.HS.UN2: %s and %s are HS\n",
               (const char*)ped[ids[0]+pf->first].pid,
               (const char*)ped[ids[1]+pf->first].pid);
            message += buffer;
            bool validFlag = true;
            int parentSex[2]; parentSex[0]=parentSex[1]=0;
            for(int s = 0; s < 2; s++){
               if(ped[ids[s]+pf->first].father) {
                  parentSex[s] |= 1;
                  if(ped[ids[1-s]+pf->first].father == ped[ids[s]+pf->first].father)
                     validFlag = false;
               }
               if(ped[ids[s]+pf->first].mother) {
                  parentSex[s] |= 2;
                  if(ped[ids[1-s]+pf->first].mother == ped[ids[s]+pf->first].mother)
                     validFlag = false;
               }
            }
            if(validFlag && (parentSex[0] || parentSex[1]) && parentSex[0]!=3 && parentSex[1]!=3){   // not reconstructed yet
               int s = parentSex[0]? 0: 1;
               if(parentSex[s] == 2){  // father is missing
                  ped[ids[1-s]+pf->first].fatid = ped[ids[s]+pf->first].fatid;
                  sprintf(buffer, "  RULE D2.HS.UN2: %s's father is now %s\n",
                     (const char*)ped[ids[1-s]+pf->first].pid,
                     (const char*)ped[ids[s]+pf->first].fatid);
               }else{   // mother is missing
                  ped[ids[1-s]+pf->first].motid = ped[ids[s]+pf->first].motid;
                  sprintf(buffer, "  RULE D2.HS.UN2: %s's mother is now %s\n",
                     (const char*)ped[ids[1-s]+pf->first].pid,
                     (const char*)ped[ids[s]+pf->first].motid);
               }
               message += buffer;
            }
            for(int k = 0; k < d2Connection[ids[0]].Length(); k++){  // Check AV
               int sharedD2 = d2Connection[ids[0]][k];
               int m = d2Connection[ids[1]].Find(sharedD2);
               if(m==-1) continue;
               double joinRatio = RoRP(ibd1segs[pairIndex[ids[0]][sharedD2]],
                  ibd1segs[pairIndex[ids[1]][sharedD2]],
                  ibd1segs[pairIndex[ids[0]][ids[1]]], bp);
               if(joinRatio < 0.85 && joinRatio > 0.6){
                  sprintf(buffer, "    Join3/Join2 = %.2lf\n", joinRatio);
                  message += buffer;
                  sprintf(buffer, "  INFERENCE D2.AV.HS: %s is %s of %s and %s\n",
                     (const char*)ped[sharedD2+pf->first].pid,
                     ped[sharedD2+pf->first].sex == 1? "uncle": "Aunt",
                     (const char*)ped[ids[0]+pf->first].pid,
                     (const char*)ped[ids[1]+pf->first].pid);
                  message += buffer;
                  for(int s = 0; s < 2; s++){
                     int pa_sex = 0;
                     if(ped[ids[s]+pf->first].father && !ped[ids[s]+pf->first].mother &&
                        relationship[ped[ids[s]+pf->first].father->serial-pf->first][sharedD2] < 0.177)
                        pa_sex = 1;
                     if(ped[ids[s]+pf->first].mother && !ped[ids[s]+pf->first].father &&
                        relationship[ped[ids[s]+pf->first].mother->serial-pf->first][sharedD2] < 0.177)
                        pa_sex = 2;
                     if(pa_sex){
                        int m = inclusionList[1].Find((pa_sex==1)?
                           ped[ids[s]+pf->first].motid: ped[ids[s]+pf->first].fatid);
                        if(m==-1) continue;
                        if(ped[sharedD2+pf->first].fatid == "0" && ped[sharedD2+pf->first].motid == "0"){
                           inclusionList[0].Push(ped[sharedD2+pf->first].famid);
                           tempS = missingBase + newparent;
                           sprintf(buffer, "    %s is created as %s's father.\n",
                              (const char*)tempS, (const char*)ped[sharedD2+pf->first].pid);
                           message += buffer;
                           ped[sharedD2+pf->first].fatid = tempS;
                           newparent++;
                           inclusionList[3].Push("");
                           inclusionList[4].Push("");
                           inclusionList[3][m] = tempS;
                           inclusionList[2].Push("1");
                           inclusionList[1].Push(tempS);

                           inclusionList[0].Push(ped[sharedD2+pf->first].famid);
                           tempS = missingBase + newparent;
                           sprintf(buffer, "    %s is created as %s's mother.\n",
                              (const char*)tempS, (const char*)ped[sharedD2+pf->first].pid);
                           message += buffer;
                           ped[sharedD2+pf->first].motid = tempS;
                           newparent++;
                           inclusionList[3].Push("");
                           inclusionList[4].Push("");
                           inclusionList[4][m] = tempS;
                           inclusionList[2].Push("2");
                           inclusionList[1].Push(tempS);
                        }else{   // uncle/aunt's parents known
                           inclusionList[3][m] = ped[sharedD2+pf->first].fatid;
                           inclusionList[4][m] = ped[sharedD2+pf->first].motid;
                        }
                        sprintf(buffer, "  RULE D2.AV.HS: %s of %s and %s (%s) now has parents (%s, %s)\n",
                           pa_sex==1? "mother": "father",
                           (const char*)ped[ids[s]+pf->first].pid,
                           (const char*)ped[ids[1-s]+pf->first].pid,
                           (const char*)ped[ids[s]+pf->first].motid,
                           (const char*)ped[sharedD2+pf->first].fatid,
                           (const char*)ped[sharedD2+pf->first].motid);
                        message += buffer;
                        break;
                     }
                  }
               }else{
                  sprintf(buffer, "    %s is HS or %s with %s and %s\n",
                     (const char*)ped[sharedD2+pf->first].pid,
                     ped[sharedD2+pf->first].sex == 1? "grandfather": "grandmother",
                     (const char*)ped[ids[0]+pf->first].pid,
                     (const char*)ped[ids[1]+pf->first].pid);
                  message += buffer;
               }
            }
            continue;
         }
         double longestLength = 0.0, longestRatio;
         int longestR1, longestR2, longestType;
         for(int s = 0; s < 2; s++){
            for(int k = 1; k < 3; k++) relativeType[k].Dimension(0);
            if(s==0) relativeType[0].Dimension(0);
            for(int k = 0; k < d123Connection[ids[s]].Length(); k++){
               int relative = d123Connection[ids[s]][k];
               if(relative == ids[1-s] || relationship[ids[s]][relative] < relationship[ids[1-s]][relative]) continue;
               tempArray = ibd1segs[pairIndex[ids[s]][relative]];  // PR
               tempArray2 = ibd1segs[pairIndex[ids[1-s]][relative]];  // OR
               SegmentIntersect(tempArray, tempArray2, ibd1seg);
               double intersectLength = SegmentLength(ibd1seg, bp);
               double temp2Length = SegmentLength(tempArray2, bp);   // OR length in Mb
               double tempLength = SegmentLength(tempArray, bp);   // OR length in Mb
               if(temp2Length > 100000000){  // OR Length > 100Mb
                  double ratio = intersectLength * 1.0 / temp2Length;
                  int type1 = int(0.5-log(ratio)/log(2));
                  if(type1 < 0) type1 = 0;
                  else if(type1 > 3) continue;
                  else if(type1 > 2) type1 = 2;
                  double ratio2 = intersectLength * 1.0 / tempLength;
                  int type2 = int(0.5-log(ratio2)/log(2));
                  if(type2 < 0) type2 = 0;
                  else if(type2 > 3) continue;
                  else if(type2 > 2) type2 = 2;
                  int type = type2 - type1;
                  if(type==0 && type1==1 && type2==1 && ratio2 < ratio*1.414214 || // HS
                     type==1 && type1==1 && type2==2 || // AV
                     type==2 && type1==0 && type2==2 && ratio > 0.8 && ratio2 > ratio * 2.828427)// GG
                  relativeType[type].Push(relative);
               }
            }
            for(int k = (s==0)?1:0; k < 3; k++){
               int relativeCount = relativeType[k].Length();
               if(relativeCount < 2) continue;
               for(int k1 = 0; k1 < relativeCount; k1++){
                  int relative = relativeType[k][k1];
                  tempArray = ibd1segs[pairIndex[ids[s]][relative]];
                  for(int k2 = k1+1; k2 < relativeCount; k2++){
                     int relative2 = relativeType[k][k2];
                     double joinRatio = RoRP(tempArray,
                        ibd1segs[pairIndex[ids[s]][relative2]],
                        ibd1segs[pairIndex[relative][relative2]], bp);
                     if(joinRatio < 0.125 && joinRatio >= 0){   // R1 R2 unrelated
                        double length = SegmentLength(tempArray, bp)+SegmentLength(ibd1segs[pairIndex[ids[s]][relative2]], bp);
                        if(length > longestLength){
                           longestType = k;
                           parent = ids[s];
                           offspring = ids[1-s];
                           longestR1 = relative;
                           longestR2 = relative2;
                           longestLength = length;
                           longestRatio = joinRatio;
                        }
                     }
                  }  // End of k2 for k2-th relative
               }  // End of k1 for k1-th relative
            }  // End of k for 3 types
         }  // End of s for either of the pair order
         if(longestLength > 10){
            sprintf(buffer, "    %s: (P1=%s, P2=%s, R1=%s, R2=%s): R1R2 relatedness = %.2G, P1R length=%.1lfMb\n",
               (const char*)stringType[longestType],
               (const char*)ped[parent+pf->first].pid,
               (const char*)ped[offspring+pf->first].pid,
               (const char*)ped[longestR1+pf->first].pid,
               (const char*)ped[longestR2+pf->first].pid,
               longestRatio, longestLength*0.000001);
            message += buffer;
            if(longestType==1){//AV
               if(ped[offspring+pf->first].father && !ped[offspring+pf->first].mother){
                  sprintf(buffer, "  RULE D2.AV.RELS.UN: %s's mother %s's parents are (%s, %s)\n",
                     (const char*)ped[offspring+pf->first].pid,
                     (const char*)ped[offspring+pf->first].motid,
                     (const char*)ped[parent+pf->first].fatid,
                     (const char*)ped[parent+pf->first].motid);
                  message += buffer;
                  String informativeMissing = ped[offspring+pf->first].motid;
                  int m = inclusionList[1].Find(informativeMissing);
                  if(m>-1){
                     inclusionList[3][m] = ped[parent+pf->first].fatid;
                     inclusionList[4][m] = ped[parent+pf->first].motid;
                  }
               }else if(ped[offspring+pf->first].mother && !ped[offspring+pf->first].father){
                  sprintf(buffer, "  RULE D2.AV.RELS.UN: %s's father %s's parents are (%s, %s)\n",
                     (const char*)ped[offspring+pf->first].pid,
                     (const char*)ped[offspring+pf->first].fatid,
                     (const char*)ped[parent+pf->first].fatid,
                     (const char*)ped[parent+pf->first].motid);
                  message += buffer;
                  String informativeMissing = ped[offspring+pf->first].fatid;
                  int m = inclusionList[1].Find(informativeMissing);
                  if(m>-1){
                     inclusionList[3][m] = ped[parent+pf->first].fatid;
                     inclusionList[4][m] = ped[parent+pf->first].motid;
                  }
               }
            }else if(longestType==2){   // GG

            }else if(longestType==0){   // HS

            }
         }
      }  // End of j for L2
   }  // End of if IBDseg
   for(int i = pf->first; i <= pf->last; i++)
      if(!valid[i-pf->first]){
         tempS = ped[i].famid;
         tempS += "->";
         tempS += ped[i].pid;
         exclusionList.Push(tempS);
      }
   delete []poConnection;
   delete []d2Connection;
   delete []d123Connection;
   if(IBDvalidFlag && Bit64==64) delete []pairIndex;
   missingBase += newparent;
   return 1;
}

void Engine::internalKING(int degree)
{
   if(shortCount==0) error("No genotype data");
   printf("Autosome genotypes stored in %d", Bit64==64? longCount:shortCount);
   printf(" words for each of %d individuals.\n", idCount);
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
   }else{
      if(SG[0]==NULL){
         if(degree==1)
            ConvertGGtoSG(GG, markerCount, SG, (stop2 < shortCount)? (stop2<<4): markerCount);
         else if(degree==2)
            ConvertGGtoSG(GG, markerCount, SG, (stop1 < shortCount)? (stop1<<4): markerCount);
         else
            error("Degree of relatedness not defined.");
      }
   }
   bool IBDvalidFlag = false;
   if(Bit64==64)
      IBDvalidFlag = PreSegment();
   if(!IBDvalidFlag){
      printf("%s\n", (const char*)segmessage);
      printf("  Inference will be based on kinship estimation only.\n");
   }
#ifdef _OPENMP
   printf("%d CPU cores are used to compute the pairwise kinship coefficients...\n",
      defaultMaxCoreCount);
#endif
   relativedegree = degree;
   IntArray *allpairs0 = new IntArray [defaultMaxCoreCount];
   IntArray *allpairs = new IntArray [defaultMaxCoreCount];
   // most computationally intensive here: SCREENING RELATIVES
   if(Bit64==64)
      ScreenCloseRelativesInSubset64Bit(allpairs0);
   else
      ScreenCloseRelativesInSubset(allpairs0);
   int id1, id2;
   long long int midrelativeCount = 0;
   for(int t = 0; t < defaultMaxCoreCount; t++)
      allpairs[t].Dimension(0);
   for(int t = 0; t < defaultMaxCoreCount; t++){
      int pairCount = allpairs0[t].Length()/2;
      for(int i = 0; i < pairCount; i++){
         id1 = allpairs0[t][i*2];
         id2 = allpairs0[t][i*2+1];
         if(ped[phenoid[id1]].famid != ped[phenoid[id2]].famid){
            allpairs[t].Push(id1);
            allpairs[t].Push(id2);
         }
      }
      midrelativeCount += allpairs[t].Length()/2;
   }
   delete []allpairs0;
   if(midrelativeCount==0) return;

   double lowerbound = pow(2.0, -(relativedegree+1.5));
   IntArray HetHetCounts, IBS0Counts, het1Counts, het2Counts, HomHomCounts, IBSCounts;
   double smaller, kinship;
   Vector ibdprops, maxLengths, ibd2props, maxLengths2;
   double ibdprop, ibd2prop;
   L0.Dimension(0);
   Lfs.Dimension(0);
   Lpo.Dimension(0);
   L2.Dimension(0);
   IntArray L1(0);
   Vector IBS0L1(0);
   for(int t = 0; t < defaultMaxCoreCount; t++){
      int pairCount = allpairs[t].Length()/2;
      if(pairCount == 0) continue;

      //****************************Compute kinship coefficient*****************************
      if(Bit64==64)
         KinshipInSubset64Bit(allpairs[t], HetHetCounts, IBS0Counts, het1Counts, het2Counts, HomHomCounts, IBSCounts);
      else
         KinshipInSubset(allpairs[t], HetHetCounts, IBS0Counts, het1Counts, het2Counts, HomHomCounts, IBSCounts);
      //****************************Compute kinship coefficient*****************************

      if(Bit64==32 || !IBDvalidFlag){
         for(int p = 0; p < pairCount; p++){
            id1 = allpairs[t][p*2];
            id2 = allpairs[t][p*2+1];
            smaller = HetHetCounts[p] + (het1Counts[p] < het2Counts[p]? het1Counts[p]: het2Counts[p]);
            kinship = 0.5 - ((het1Counts[p]+het2Counts[p])*0.25+IBS0Counts[p])/smaller;
            if(kinship <= 0) continue;
            int deg = int(-log(kinship)/0.6931472-0.5);;
            if(deg == 0){
               L0.Push(phenoid[id1]);
               L0.Push(phenoid[id2]);
            }else if(deg == 1){
               int IBS0Count = IBS0Counts[p];
               int notMissingCount = het1Counts[p] + het2Counts[p] + HomHomCounts[p] + HetHetCounts[p];
               double ibs0 = IBS0Count*1.0/notMissingCount;
               IBS0L1.Push(ibs0);
               L1.Push(phenoid[id1]);
               L1.Push(phenoid[id2]);
            }else if(deg == 2){
               L2.Push(phenoid[id1]); L2.Push(phenoid[id2]);
            }
         }
      }else{   // 64 bit && IBDvalidFlag
         IBDSegInSubset64Bit(allpairs[t], ibdprops, maxLengths, ibd2props, maxLengths2);
      for(int p = 0; p < pairCount; p++){
         if(!het1Counts[p] && !het2Counts[p] && !HetHetCounts[p]) continue;
         id1 = allpairs[t][p*2];
         id2 = allpairs[t][p*2+1];
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
   }  // end of t loop for thread

   if(Bit64==32 || !IBDvalidFlag){
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
         printf("Cutoff value between FS and PO is set at %.4f\n", errorrateCutoff);
      }
      for(int i = 0; i < IBS0L1.Length(); i++){
         if(IBS0L1[i] > errorrateCutoff) {
            Lfs.Push(L1[i*2]);
            Lfs.Push(L1[i*2+1]);
         }else{
            Lpo.Push(L1[i*2]);
            Lpo.Push(L1[i*2+1]);
         }
      }
   }
   delete []allpairs;
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
   if(idCount >= 10)  // fast computation
      internalKING(degree);
   else if(Bit64 == 64){
      printf("This function is currently disabled for tiny dataset with sample size < 10.\n");
      return 0;
   }else
      runKING();
   printf("Clustering up to %d%s-degree relatives in families...\n",
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
         if(Bit64 == 64 && totalLength > 10000000){
            IntArray ids, pairs, pairIndex, HetHetCounts, IBS0Counts, het1Counts, het2Counts, HomHomCounts, IBSCounts;
            Vector ibdprops, maxLengths, ibd2props, maxLengths2;
            String type, outfile;
            outfile.Copy(prefix);
            outfile.Add("cluster.kin");
            fp = fopen((const char*)outfile, "wt");
            fprintf(fp, "FID\tID1\tID2\tSex1\tSex2\tN_SNP\tHetHet\tIBS0\tHetConc\tHomIBS0\tKinship\tIBD1Seg\tIBD2Seg\tPropIBD\tInfType\n");
            for(int c = 0; c < clusterCount; c++){
               ids.Dimension(0);
               for(int k = 0; k < cluster[c].Length(); k++){
                  int f = cluster[c][k];
                  for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
                     if(geno[i]!=-1) ids.Push(geno[i]);
               }
               pairs.Dimension(0);
               int idsCount = ids.Length();
               for(int i = 0; i < idsCount; i++)
                  for(int j = i+1; j < idsCount; j++){
                     pairs.Push(ids[i]);
                     pairs.Push(ids[j]);
                  }
               int pairCount = pairs.Length()/2;
               if(pairCount==0) continue;
               KinshipInSubset64Bit(pairs, HetHetCounts, IBS0Counts, het1Counts, het2Counts, HomHomCounts, IBSCounts);
               IBDSegInSubset64Bit(pairs, ibdprops, maxLengths, ibd2props, maxLengths2);
               for(int p = 0; p < pairCount; p++){
                  if(!het1Counts[p] && !het2Counts[p] && !HetHetCounts[p]) continue;
                  int id1 = pairs[p*2]; int id2 = pairs[p*2+1];
                  double kinship = (HetHetCounts[p] - IBS0Counts[p]*2.0) / (HetHetCounts[p]*2+het1Counts[p]+het2Counts[p]);
                  int notMissingCount = HetHetCounts[p]+het1Counts[p]+het2Counts[p]+HomHomCounts[p];
                  double ibs0 = IBS0Counts[p] * 1.0/notMissingCount;
                  double CHet = HetHetCounts[p] * 1.0 / (HetHetCounts[p]+het1Counts[p]+het2Counts[p]);
                  double ibdprop = ibdprops[p];
                  double ibd2prop = ibd2props[p];
                  double pi = ibd2prop + ibdprop * 0.5;
                  if(CHet<0.8){  // not MZ/Dup
                     if(pi > 0.3535534){  // 1st-degree
                        if(ibdprop + ibd2prop > 0.96 || (ibdprop + ibd2prop > 0.9 && ibd2prop <= 0.08))
                           type = "PO";
                        else if(ibd2prop > 0.08)
                           type = "FS";
                        else
                           type = "2nd";
                     }else if(pi > 0.1767767){  // 2nd-degree
                        if(pi > 0.32 && ibd2prop > 0.15)
                           type = "FS";
                        else
                           type = "2nd";
                     }else if(pi > 0.08838835)
                        type = "3rd";
                     else if(pi > 0.04419417)
                        type = "4th";
                     else
                        type = "UN";
                  }else // Duplicate
                     type="Dup/MZ";
                  int phenoid1 = phenoid[id1];
                  int phenoid2 = phenoid[id2];
                  fprintf(fp, "KING%d\t%s\t%s\t%d\t%d\t%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%s\n",
                     c+1, (const char*)ped[phenoid1].pid, (const char*)ped[phenoid2].pid,
                     ped[phenoid1].sex, ped[phenoid2].sex, notMissingCount,
                     HetHetCounts[p]*1.0/notMissingCount, ibs0, CHet, // CHet
                     ibs0/(ibs0+(IBSCounts[p]-HetHetCounts[p])*1.0/notMissingCount), kinship,
                     ibdprops[p], ibd2props[p], ibd2props[p] + ibdprops[p]*0.5, (const char*)type);       
               }  // End of pairs
            }  // End of c for clusters
            fclose(fp);
            printf("Pair-wise relatedness in newly clustered families saved in %s.\n", (const char*)outfile);
         }  // End of if IBDseg                                                                               
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

bool Engine::rebuild(int id_added)
{
   bool built = false;
   printf("Pedigree reconstruction starts at %s", currentTime());
   printf("Reconstructing pedigree...\n");
   cAge = ped.covariateNames.SlowFind("AGE");
   if(cAge > -1)
      printf("Covariate %s is used for pedigree reconstruction.\n",
         (const char*)ped.covariateNames[cAge]);
   else
      printf("Age information not provided.\n");
   // Hopefully missingBase is not in the range of IID; default 100000
   for(missingBase=id_added; ;missingBase += 100){
      bool overlapFlag = false;
      for(int i = 0; i < ped.count; i++)
         if(int(ped[i].pid) >= missingBase && int(ped[i].pid) < missingBase+100){
            overlapFlag = true;
            break;
         }
      if(!overlapFlag) break;
   }
   if(Bit64==64)PreSegment(/*chrSeg, totalLength, segmessage*/);
   String logfile = prefix; logfile.Add("build.log");
   FILE *fp2 = fopen((const char*)logfile, "wt");
   String updatefile = prefix; updatefile.Add("updateparents.txt");
   FILE *fp = fopen((const char*)updatefile, "wt");
   String message;
   for(int f = 0; f < ped.familyCount; f++) // ready to be parallalized later
      if(BuildOneFamily(f, chrSeg, totalLength, message)){ // pedigree reconstruction in this family
         built = true;
         if(rplotFlag)
            for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
//               if(ped[i].pid.SubStr(0,4) == "KING" && ped[i].ngeno == 0) continue;
               fprintf(fp, "%s\t%s\t%s\t%s\t%d\t%d\t%d\n",
                  (const char*)ped[i].famid, (const char*)ped[i].pid,
                  (const char*)ped[i].fatid, (const char*)ped[i].motid,
                  ped[i].sex, ped.affectionCount? ped[i].affections[0]: 0,
                  ped[i].ngeno? 0: 1);
            }
         else
            for(int i = ped.families[f]->first; i <= ped.families[f]->last; i++){
               if(ped[i].pid.SubStr(0,4) == "KING" && ped[i].ngeno == 0) continue;
               fprintf(fp, "%s\t%s\t%s\t%s\n",
                  (const char*)ped[i].famid, (const char*)ped[i].pid,
                  (ped[i].fatid.SubStr(0, 4) == "KING" && ped[i].father && ped[i].father->ngeno==0)? "0": (const char*)ped[i].fatid,
                  (ped[i].motid.SubStr(0, 4) == "KING" && ped[i].mother && ped[i].mother->ngeno==0)? "0": (const char*)ped[i].motid);
            }
         printf("%s", (const char*)message);
         fprintf(fp2, "%s", (const char*)message);
      }
   bool WriteFlag = false;
   int inclusionCount = inclusionList[1].Length();
   if(rplotFlag)
      for(int i = 0; i < inclusionCount; i++){
         fprintf(fp, "%s\t%s\t%s\t%s\t%d\t0\t1\n",
            (const char*)inclusionList[0][i],
            (const char*)inclusionList[1][i],
            inclusionList[3][i]==""? "0": (const char*)inclusionList[3][i],
            inclusionList[4][i]==""? "0": (const char*)inclusionList[4][i],
            (int)inclusionList[2][i]);
         if(inclusionList[3][i]!="" || inclusionList[4][i]!="") WriteFlag = true;
      }
   else
      for(int i = 0; i < inclusionCount; i++)
         if(inclusionList[3][i]!="" || inclusionList[4][i]!=""){
            fprintf(fp, "%s\t%s\t%s\t%s\n",
               (const char*)inclusionList[0][i],
               (const char*)inclusionList[1][i],
               (const char*)inclusionList[3][i],
               (const char*)inclusionList[4][i]);
            WriteFlag = true;
         }
   fclose(fp);
   fclose(fp2);
   printf("\n");
   String updateidsfile = prefix;
   printf("Details of pedigree reconstruction are also available in log file %s\n", (const char*)logfile);
   updateidsfile.Add("updateids.txt");
   printf("Update-ID information is saved in file %s\n", (const char*)updateidsfile);
   if(built)
      printf("Update-parent information is saved in file %s\n", (const char*)updatefile);
   else
      printf("No pedigrees can be reconstructed.\n");
   if(WriteFlag) printf("PLINK format genotypes are created for allowing untyped individuals in pedigrees.\n");
   String temp = prefix;
   if(SaveFormat == "MERLIN"){
      printf("Start writing reconstructed pedigrees in MERLIN format...\n");
      WriteMerlin();
   }
   if(SaveFormat == "PLINK" || WriteFlag){
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
   return built;
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
   for(missingBase=100000; ;missingBase += 100000){
      bool overlapFlag = false;
      for(int i = 0; i < ped.count; i++)
         if(int(ped[i].pid) >= missingBase && int(ped[i].pid) < missingBase+100000){
            overlapFlag = true;
            break;
         }
      if(!overlapFlag) break;
   }
   if(Bit64==64)
      PreSegment(/*chrSeg, totalLength, segmessage*/);
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
   printf("Update-parent information is saved in file %s\n", (const char*)updatefile);

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

