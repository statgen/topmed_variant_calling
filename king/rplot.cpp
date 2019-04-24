//////////////////////////////////////////////////////////////////////
// rplot.cpp
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

#include "rplot.h"
#include "MathMatrix.h"
#include "StringArray.h"

int CheckRout(const char *scriptfile);
void plotIBD1vsIBD2(const char *prefix, FILE *fp);

void plotMIerror(const char *prefix)
{
   String scriptfile=prefix;
   scriptfile.Add("_MIerrorplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING MI error plot, by Wei-Min Chen and Zhennan Zhu\n", (const char*)scriptfile);
   fprintf(fp, "library(kinship2)\n");
   fprintf(fp, "library(igraph)\n");
   fprintf(fp, "ped <- read.table(file=\"%ssplitped.txt\", stringsAsFactors=FALSE)[,c(1,2,5,6,7,8,9)]\n", (const char*)prefix);
   fprintf(fp, "ped$V8[ped$V8==-9 | ped$V8==0 | ped$V8==1] <- 0\n");
   fprintf(fp, "ped$V8[ped$V8==2] <- 1\n");
   fprintf(fp, "colnames(ped) <- c(\"FID\", \"ID\", \"FA\", \"MO\", \"Sex\", \"Affected\", \"Status\")\n");
   fprintf(fp, "data <- read.table(\"%s.kin\",header = TRUE, stringsAsFactors = FALSE)\n", (const char*)prefix);
   fprintf(fp, "mi.err <- data[data[, \"Error\"]==1, \"FID\" ]\n");
   fprintf(fp, "Inf.color <- c(\"purple\", \"red\", \"green\", \"blue\", \"yellow\")\n");
   fprintf(fp, "Inf.type <- c(\"Dup/MZ\", \"PO\", \"FS\", \"2nd\", \"3rd\")\n");
   fprintf(fp, "data.fam <- merge(data, ped, by.x = c(\"FID\", \"ID1\"), by.y = c(\"FID\", \"ID\"))\n");
   fprintf(fp, "data.all <- merge(data.fam, ped, by.x = c(\"FID\", \"ID2\"), by.y = c(\"FID\", \"ID\"))[, c(\"FID\", \"ID1\", \"ID2\", \"Sex.x\", \"Sex.y\", \"InfType\")]\n");
   fprintf(fp, "data.all <- data.all[data.all[, \"InfType\"] %%in%% Inf.type, ]\n");
   fprintf(fp, "for (i in 1:5) data.all[data.all$InfType == Inf.type[i], 6] <- i\n");
   fprintf(fp, "postscript(\"%s_MIerrorplot.ps\", paper=\"letter\", horizontal=T, fonts=c(\"serif\", \"Palatino\"))\n", (const char*)prefix);
   fprintf(fp, "par(mfrow=c(1, 2))\n");
   fprintf(fp, "for (famid in unique(mi.err)){\n");
   fprintf(fp, "  fam <- ped[ped[, \"FID\"]==famid,]\n");
   fprintf(fp, "if (all(fam[, c(\"FA\", \"MO\")]==0)){\n");
   fprintf(fp, "    g.empty <- make_empty_graph(n = nrow(fam), directed = FALSE)\n");
   fprintf(fp, "    plot(g.empty, vertex.size=27, vertex.color=NA, vertex.label.cex=.5, vertex.label.dist=1.6,\n");
   fprintf(fp, "         vertex.label.degree= pi/2, vertex.label.color=\"black\", vertex.label= fam[,\"ID\"],\n");
   fprintf(fp, "         edge.color=NA, layout=layout_with_fr(g.empty, grid=\"nogrid\"), asp=0,\n");
   fprintf(fp, "         vertex.shape=c(\"none\", \"square\", \"circle\")[1+fam[,\"Sex\"]])}else{\n");
   fprintf(fp, "  pedplot <- pedigree(id = fam$ID, dadid = fam$FA, momid = fam$MO, sex = as.numeric(fam$Sex),\n");
   fprintf(fp, "    affected = as.numeric(fam$Affected), status = as.numeric(fam$Status), famid = fam$FID, missid = 0)\n");
   fprintf(fp, "  plot(pedplot[toString(famid)], cex = 0.5, symbolsize = 2.8)}\n");
   fprintf(fp, "  fam.sub <- data.all[data.all$FID==famid,][, 2:6]\n");
   fprintf(fp, "  id <- unique(mapply(c, fam.sub[,c(1,3)], fam.sub[,c(2, 4)]))\n");
   fprintf(fp, "  g <- graph_from_data_frame(d=fam.sub, vertices=id[, 1], directed=FALSE)\n");
   fprintf(fp, "  plot(g, edge.width=1.5, vertex.size=27, vertex.color=NA, vertex.label.cex=0.5,\n");
   fprintf(fp, "       edge.color=Inf.color[as.numeric(fam.sub$InfType)], layout=layout_with_fr(g, grid=\"nogrid\"), asp=0,\n");
   fprintf(fp, "       vertex.shape=c(\"none\", \"square\", \"circle\")[1+as.numeric(id[, 2])], margin=0.3)\n");
   fprintf(fp, "  legend(\"bottomright\", Inf.type, lty=1, col=Inf.color, text.col=Inf.color, cex=0.7, bty=\"n\")\n");
   fprintf(fp, "  mtext(paste(\"Documented versus Inferred Family\", famid, \"in %s\"), side = 3, line = -2, outer = TRUE)}\n", (const char*)prefix);
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   switch(error){
      case 3:  printf("R plot failed for missing R libraries kinship2/igraph .\n");
               printf("  Please rerun R code %s (or KING) after kinship2/igraph is installed.\n", (const char*)scriptfile);
               break;
      case 2:  printf("R code %s failed.\n", (const char*)scriptfile);
               break;
      case 1:  printf("R code %s failed.\n", (const char*)scriptfile);
               printf("  Please check %sout for details.\n", (const char*)scriptfile);
               break;
      case 0:  sprintf(command, "ps2pdf %s_MIerrorplot.ps", (const char*)prefix);
               system(command);
               printf("MI error plots are generated in %s_MIerrorplot.pdf\n", (const char*)prefix);
               break;
      default: printf("Unexpected error. Please contact KING authors.\n");
   }
}

void plotUniqueFamily(const char *prefix, int degree, const char *analysis)
{
   if(analysis != "related" && analysis != "ibdseg") return;
   String inputfile = prefix;
   if(analysis == "related") inputfile.Add(".kin0");
   else inputfile.Add(".seg");
   String scriptfile=prefix;
   scriptfile.Add("_uniqfamplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --%s, by Wei-Min Chen and Zhennan Zhu\n", (const char*)scriptfile, (const char*)analysis);
   fprintf(fp, "library(igraph)\n");
   fprintf(fp, "data <- read.table(file=\"%s\", header=TRUE, stringsAsFactors=FALSE)[,c(\"ID1\", \"ID2\", \"InfType\")]\n", (const char*)inputfile);
   fprintf(fp, "Inf.color <- c(\"purple\", \"red\", \"green\", \"blue\")\n");
   fprintf(fp, "Inf.type <- c(\"Dup/MZ\", \"PO\", \"FS\", \"2nd\")\n");
   fprintf(fp, "relatives <- data[data$InfType %%in%% Inf.type, ]\n");
   fprintf(fp, "for (i in 1:length(Inf.type)) relatives[relatives$InfType == Inf.type[i], \"InfType\"] <- i\n");
   fprintf(fp, "colnames(relatives)[3] <- \"color\"\n");
   fprintf(fp, "BuildPed <- function(v){\n");
   fprintf(fp, " buildType <- \"\"\n");
   fprintf(fp, " if(v[1]==3 && v[2]==3 && v[4]==2 && v[6]==1) buildType <- \"2_GG/HalfSibs\"\n");
   fprintf(fp, " else if(v[1]==3 && v[2]==3 && v[3]==1 && v[4]==2) buildType <- \"2_MZtwins+1_Parent/Child\"\n");
   fprintf(fp, " else if(v[1]==4 && v[2]==5 && v[3]==1 && v[4]==4) buildType <- \"2_Parents+2_MZtwins\"\n");
   fprintf(fp, " if(v[5]>0){\n");
   fprintf(fp, "  s <- floor(sqrt(v[5]*2))+1\n");
   fprintf(fp, "  if(s*(s-1)/2==v[5]){\n");
   fprintf(fp, "   if(v[2]==v[5] && v[4]==0) buildType <- paste(s, \"_FullSiblings\", sep=\"\")\n");
   fprintf(fp, "   else if(v[2]==v[5]+s && v[4]==s) buildType <- paste(\"1_Parent+\", s, \"_FullSiblings\", sep=\"\")\n");
   fprintf(fp, "   else if(v[2]==v[5]+s*2 && v[4]==s*2) buildType <- paste(\"2_Parents+\", s, \"_FullSiblings\", sep=\"\")\n");
   fprintf(fp, "   else buildType <- paste(s, \"_FullSiblings+\", v[1]-s, \"_Relatives\", sep=\"\")\n");
   fprintf(fp, "  }else if(v[3]==0 && v[4]==0 && v[6]==0) buildType <- \"Undetermined_FS_Only\"\n");
   fprintf(fp, "  else buildType <- \"Undetermined_FS\"\n");
   fprintf(fp, " }else if(v[4]==0 && v[6]==0){\n");
   fprintf(fp, "  buildType <- ifelse(v[3]==1, \"2_MZtwins\", \"Undetermined_MZ_Only\")\n");
   fprintf(fp, " }else if(v[3]==0 && v[4]==0){\n");
   fprintf(fp, "  buildType <- ifelse(v[6]==1, \"2_Second-Degree_Relatives\", \"Undetermined_2nd_Only\")\n");
   fprintf(fp, " }else if(v[3]==0 && v[6]==0){\n");
   fprintf(fp, "  if(v[4]==1) buildType <- \"1_Parent+1_Child\"\n");
   fprintf(fp, "  else if(v[4]==2) buildType <- \"2_Parents+1_Child(Trio)\"\n");
   fprintf(fp, "  else buildType <- \"Undetermined_PO_Only\"}\n");
   fprintf(fp, "if(buildType==\"\"){\n");
   fprintf(fp, " if(v[3]>0) buildType <- \"Undetermined_MZ\"\n");
   fprintf(fp, " else if(v[4]>0) buildType <- \"Undetermined_PO\"}\n");
   fprintf(fp, " return(buildType)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "g <- graph_from_data_frame(d = relatives, vertices = unique(c(relatives$ID1, relatives$ID2)), directed = FALSE)\n");
   fprintf(fp, "imc <- cluster_infomap(g)\n");
   fprintf(fp, "imc.membership <- membership(imc)\n");
   fprintf(fp, "community.info <- cbind(imc$membership,imc$names)\n");
   fprintf(fp, "colnames(community.info) <- c(\"membership\", \"ID\")\n");
   fprintf(fp, "community.rel <- merge(relatives, community.info, by.x = c(\"ID1\"), by.y = c(\"ID\"))\n");
   fprintf(fp, "community.rel$membership <- as.numeric(as.character(community.rel$membership))\n");
   fprintf(fp, "color_counts <- function(i) {\n");
   fprintf(fp, "  sub.colors <- community.rel[community.rel[, \"membership\"]==i, \"color\"]\n");
   fprintf(fp, "  return(sapply(1:4, function(y) sum(sub.colors==y)))}\n");
   fprintf(fp, "colors <- t(sapply(1:length(imc), color_counts))\n");
   fprintf(fp, "fam.info <- cbind(community.index=1:length(imc), sizes.node = sizes(imc), edge.colors=apply(colors, 1, sum), colors)\n");
   fprintf(fp, "fam.info <- as.data.frame(fam.info)\n");
   fprintf(fp, "all.counts <- aggregate(cbind(fam.info[0],counts=1),fam.info[,-1], length)\n");
   fprintf(fp, "fam.info.counts <- merge(fam.info, all.counts, by= c(\"sizes.node\",\"edge.colors\",\"V4\",\"V5\",\"V6\",\"V7\"))\n");
   fprintf(fp, "uniq.fam.info <- fam.info.counts[!duplicated(fam.info.counts[, c(1:6,8)]), ][,c(7, 1:6, 8)]\n");
   fprintf(fp, "uniq.fam.info <- uniq.fam.info[order(uniq.fam.info$counts,decreasing=TRUE),]\n");
   fprintf(fp, "all.counts <- uniq.fam.info[,\"counts\"]\n");
   fprintf(fp, "uniq.cluster <- induced_subgraph(g, (1:length(V(g)))[imc.membership %%in%% uniq.fam.info[, \"community.index\"]])\n");
   fprintf(fp, "all.names <- sapply(V(uniq.cluster)$name, function(x) membership(imc)[names(imc.membership)%%in%%x])\n");
   fprintf(fp, "all.builds <- apply((uniq.fam.info[, -1]),1,BuildPed)\n");
   fprintf(fp, "lo <- layout_(uniq.cluster, with_fr(), normalize())\n");
   fprintf(fp, "LocationForACluster<-function(x){\n");
   fprintf(fp, "  lo.local <- lo[all.names==x,]\n");
   fprintf(fp, "  return(c(min(as.numeric(lo.local[,1])), max(as.numeric(lo.local[,1])), max(as.numeric(lo.local[,2]))))}\n");
   fprintf(fp, "locations <- sapply(uniq.fam.info[, 1], LocationForACluster)\n");
   fprintf(fp, "postscript(\"%s_uniqfamplot.ps\", paper=\"letter\", horizontal=T)\n", prefix);
   fprintf(fp, "plot(uniq.cluster, vertex.color=NA, vertex.size=1, vertex.label=NA, layout=lo, asp=0,\n");
   fprintf(fp, "  edge.color=Inf.color[as.numeric(E(uniq.cluster)$color)], main=\"All Unique Family Configurations in %s\")\n", prefix);
   fprintf(fp, "text((locations[1,]+locations[2,])/2, locations[3,]+0.04, all.counts)\n");
   fprintf(fp, "legend(\"bottomright\", Inf.type, lty = 1, col = Inf.color, text.col = Inf.color, cex = 0.7, bty = \"n\")\n");
   if(degree>1){
      fprintf(fp, "par(mfrow=c(2,3))\n");
      fprintf(fp, "is.built <- all.builds!=\"\"\n");
      fprintf(fp, "for(i in uniq.fam.info[, 1][is.built]){\n");
      fprintf(fp, "  index <- (1:length(uniq.fam.info[, 1]))[uniq.fam.info[, 1]==i]\n");
      fprintf(fp, "  g1 <- induced_subgraph(g, (1:length(V(g)))[imc.membership == i])\n");
      fprintf(fp, "  if(substr(all.builds[index],1,12)!=\"Undetermined\"&&all.counts[index]>1||all.counts[index]>10)\n");
      fprintf(fp, "  plot(g1, vertex.color=NA, vertex.size=3, vertex.label=NA, layout=layout_with_fr, asp=0,\n");
      fprintf(fp, "  edge.color=Inf.color[as.numeric(edge.attributes(g1)$color)], main=paste(all.counts[index], all.builds[index], \"Families\"))}\n");
   }
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   if(!CheckRout(scriptfile)){
      sprintf(command, "ps2pdf %s_uniqfamplot.ps", prefix);
      system(command);
      printf("Unique family plot is generated in %s_uniqfamplot.pdf\n\n", prefix);
   }
}

void plotDuplicate(const char *prefix)
{
   String scriptfile=prefix;
   scriptfile.Add("_duplicateplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --duplicate, by Wei-Min Chen and Zhennan Zhu\n", (const char*)scriptfile);
   fprintf(fp, "library(igraph)\n");
   fprintf(fp, "data <- read.table(file=\"%s.con\", header=TRUE, stringsAsFactors=FALSE)[,c(2,4)]\n", prefix);
   fprintf(fp, "if(dim(data)[1]==0) q()\n");
   fprintf(fp, "postscript(\"%s_duplicateplot.ps\", paper=\"letter\", horizontal=T, fonts=c(\"serif\", \"Palatino\"))\n", prefix);
   fprintf(fp, "if(substr(data[1,1],1,4)==\"QRY_\"||substr(data[1,1],1,4)==\"REF_\"||substr(data[1,2],1,4)==\"QRY_\"||substr(data[1,2],1,4)==\"REF_\"){\n");
   fprintf(fp, " ordered <- rbind(data[substr(data[,1],1,3)==\"QRY\" & substr(data[,2],1,3)==\"REF\", c(1,2)],\n");
   fprintf(fp, "  data[substr(data[,1],1,3)==\"REF\" & substr(data[,2],1,3)==\"QRY\", c(2,1)])\n");
   fprintf(fp, " ordered <- cbind(substr(ordered[,1],5,1000),substr(ordered[,2],5,1000))\n");
   fprintf(fp, " mismatched <- ordered[ordered[,1]!=ordered[,2],]\n");
   fprintf(fp, " if(dim(mismatched)[1]==0) q()\n");
   fprintf(fp, " g <- graph_from_data_frame(d=mismatched, vertices=unique(c(mismatched[,1],mismatched[,2])), directed=TRUE)\n");
   fprintf(fp, " plot(g, vertex.shape=\"none\", vertex.label.cex=0.5, edge.arrow.size=0.3,\n");
   fprintf(fp, "  layout=layout_with_fr, asp=0, main=paste(dim(mismatched)[1], \"Pairs Are Matched (QUERY<-REF) in %s\"))\n", prefix);
   fprintf(fp, "}else{\n");
   fprintf(fp, "g <- graph_from_data_frame(d=data, vertices=unique(c(data$ID1,data$ID2)), directed=FALSE)\n");
   fprintf(fp, "if(dim(data)[1]<500) {plot(g, vertex.shape=\"none\", vertex.label.cex=0.5,\n");
   fprintf(fp, "  layout=layout_with_fr, asp=0, main=paste(dim(data)[1], \"Duplicate Pairs in %s\"))\n", prefix);
   fprintf(fp, "}else{g2 <- induced_subgraph(g, c((1:length(V(g)))[degree(g)>1],(1:length(V(g)))[degree(g)==1])[1:200])\n");
   fprintf(fp, "plot(g2, vertex.shape=\"none\", vertex.label.cex=0.5, layout=layout_with_fr, asp=0, main=paste(\"200 Duplicates in %s\"))}}\n", prefix);
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   switch(error){
      case 3:  printf("--duplicate is done but R plot failed for missing R library igraph.\n");
               printf("  Please intall igraph and run R code %s (or KING) again.\n\n", (const char*)scriptfile);
               break;
      case 2:  printf("--duplicate is done but R code %s failed.\n\n", (const char*)scriptfile);
               break;
      case 1:  printf("--duplicate is done but R code %s failed.\n", (const char*)scriptfile);
               printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
               break;
      case 0:  sprintf(command, "ps2pdf %s_duplicateplot.ps", prefix);
               system(command);
               printf("Duplicate plot is generated in %s_duplicateplot.pdf\n\n", prefix);
               break;
      default: printf("Unexpected error. Please contact KING authors.\n");
   }
}

void plotCluster(const char *prefix)
{
   String scriptfile=prefix;
   scriptfile.Add("_clusterplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --cluster, by Wei-Min Chen and Zhennan Zhu\n", (const char*)scriptfile);
   fprintf(fp, "library(igraph)\n");
   fprintf(fp, "data <- read.table(file=\"%scluster.kin\", header=TRUE, stringsAsFactors=FALSE)\n", (const char*)prefix);
   fprintf(fp, "postscript(\"%s_clusterplot.ps\", paper=\"letter\", horizontal=T, fonts=c(\"serif\", \"Palatino\"))\n", (const char*)prefix);
   fprintf(fp, "Inf.color <- c(\"purple\", \"red\", \"green\", \"blue\", \"yellow\")\n");
   fprintf(fp, "Inf.type <- c(\"Dup/MZ\", \"PO\", \"FS\", \"2nd\", \"3rd\")\n");
   fprintf(fp, "for(famid in unique(data$FID)){\n");
   fprintf(fp, "  fam <- data[data$FID==famid & data$InfType %%in%% Inf.type,][,c(2,3,4,5,15)]\n");
   fprintf(fp, "  id <- unique(mapply(c, fam[,c(1,3)], fam[,c(2,4)]))\n");
   fprintf(fp, "  for(i in 1:5) fam[fam$InfType==Inf.type[i],5] <- i\n");
   fprintf(fp, "  g <- graph_from_data_frame(d=fam, vertices=id[,1], directed=FALSE)\n");
   fprintf(fp, "  plot(g, edge.width=1.5, vertex.size=4, vertex.color=NA, vertex.label.cex=0.5,\n");
   fprintf(fp, "    edge.color=Inf.color[as.numeric(fam$InfType)], layout=layout_with_fr(g,grid=\"nogrid\"), asp=0,\n");
   fprintf(fp, "    vertex.shape=c(\"none\",\"square\",\"circle\")[1+as.numeric(id[,2])])\n");
   fprintf(fp, "  legend(\"bottomright\", Inf.type, lty=1, col=Inf.color, text.col=Inf.color, pt.cex=2, cex=0.8, bty=\"n\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   switch(error){
      case 3:  printf("--cluster is done but R plot failed for missing R library igraph.\n");
               printf("  Please intall igraph and run R code %s (or KING) again.\n\n", (const char*)scriptfile);
               break;
      case 2:  printf("--cluster is done but R code %s failed.\n\n", (const char*)scriptfile);
               break;
      case 1:  printf("--cluster is done but R code %s failed.\n", (const char*)scriptfile);
               printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
               break;
      case 0:  sprintf(command, "ps2pdf %s_clusterplot.ps", (const char*)prefix);
               system(command);
               printf("Plots of newly clustered families are generated in %s_clusterplot.pdf\n\n", (const char*)prefix);
               break;
      default: printf("Unexpected error. Please contact KING authors.\n");
   }
}

void plotSplitped(const char *prefix)
{
   String scriptfile=prefix;
   scriptfile.Add("_pedplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING pedigree plot, by Wei-Min Chen and Zhennan Zhu\n", (const char*)scriptfile);
   fprintf(fp, "library(kinship2)\n");
   fprintf(fp, "ped <- read.table(file=\"%ssplitped.txt\", stringsAsFactors=FALSE)[,3:9]\n", (const char*)prefix);
   fprintf(fp, "postscript(\"%s_pedplot.ps\", paper=\"letter\", horizontal=T)\n",(const char*)prefix);
//   fprintf(fp, "ped$V8[ped$V8==-9 | ped$V8==0] <- NA\n");
//   fprintf(fp, "ped$V8[ped$V8==1] <- 0\n");
   fprintf(fp, "ped$V8[ped$V8==-9 | ped$V8==0 | ped$V8==1] <- 0\n");
   fprintf(fp, "ped$V8[ped$V8==2] <- 1\n");
   fprintf(fp, "pedAll <- pedigree(id = ped$V4, dadid = ped$V5, momid = ped$V6, sex = as.numeric(ped$V7), affected = as.numeric(ped$V8), status = as.numeric(ped$V9), famid = ped$V3, missid = 0)\n");
   fprintf(fp, "for(f in unique(ped$V3))\n");
   fprintf(fp, "  if(any(ped$V5[ped$V3 == f] != 0 | ped$V6[ped$V3 == f] != 0))\n");
   fprintf(fp, "    plot(pedAll[toString(f)], cex=0.5, symbolsize = 2.8)\n");
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   switch(error){
      case 3:  printf("R plot failed for missing R library kinship2.\n");
               printf("  Please rerun R code %s (or KING) after kinship2 is installed.\n\n", (const char*)scriptfile);
               break;
      case 2:  printf("R code %s failed.\n\n", (const char*)scriptfile);
               break;
      case 1:  printf("R code %s failed.\n", (const char*)scriptfile);
               printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
               break;
      case 0:  sprintf(command, "ps2pdf %s_pedplot.ps", (const char*)prefix);
               system(command);
               printf("Pedigree plots are generated in %s_pedplot.pdf\n\n", (const char*)prefix);
               break;
      default: printf("Unexpected error. Please contact KING authors.\n");
   }
}

void plotBuild(const char *prefix)
{
   String scriptfile=prefix;
   scriptfile.Add("_buildplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --build, by Wei-Min Chen and Zhennan Zhu\n", (const char*)scriptfile);
   fprintf(fp, "library(kinship2)\n");
   fprintf(fp, "ped <- read.table(file=\"%supdateparents.txt\", stringsAsFactors=FALSE)\n", (const char*)prefix);
   fprintf(fp, "ped$V6[ped$V6==-9 | ped$V6==0 | ped$V6==1] <- 0\n");
   fprintf(fp, "ped$V6[ped$V6==2] <- 1\n");
   fprintf(fp, "pedAll <- pedigree(id = ped$V2, dadid = ped$V3, momid = ped$V4, sex = as.numeric(ped$V5), affected = as.numeric(ped$V6), status = as.numeric(ped$V7),  famid = ped$V1, missid = 0)\n");
   fprintf(fp, "postscript(\"%s_buildplot.ps\", paper=\"letter\", horizontal=T)\n",(const char*)prefix);
   fprintf(fp, "for(f in unique(ped$V1))\n");
   fprintf(fp, "  if(any(ped$V3[ped$V1 == f] != 0 | ped$V4[ped$V1 == f] != 0))\n");
   fprintf(fp, "    plot(pedAll[toString(f)], cex=0.5)\n");
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   switch(error){
      case 3:  printf("--build is done but R plot failed for missing R library kinship2.\n");
               printf("  Please install kinship2 and run R code %s (or KING) again.\n\n", (const char*)scriptfile);
               break;
      case 2:  printf("--build is done but R code %s failed.\n\n", (const char*)scriptfile);
               break;
      case 1:  printf("--build is done but R code %s failed.\n", (const char*)scriptfile);
               printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
               break;
      case 0:  sprintf(command, "ps2pdf %s_buildplot.ps", (const char*)prefix);
               system(command);
               printf("Plots of newly reconstruction pedigrees are generated in %s_buildplot.pdf\n\n", (const char*)prefix);
               break;
      default: printf("Unexpected error. Please contact KING authors.\n");
   }
}

void plotHEreg(const char *prefix, int SEXCHR)
{
   String scriptfile=prefix;
   scriptfile.Add("_herplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --hereg, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_herplot.ps\", paper=\"letter\", horizontal=T)\n",(const char*)prefix);
   fprintf(fp, "data<-c()\n");
   fprintf(fp, "if(file.exists(\"%s.her\")) data <- read.table(file=\"%s.her\", header=T)\n",
      (const char*)prefix,(const char*)prefix);
   fprintf(fp, "if(length(data$LOD)>0){\n");
   fprintf(fp, "alltraits <- as.character(data$Trait)\n");
   fprintf(fp, "uniqtraits <- unique(alltraits)\n");
   fprintf(fp, "for(trait in uniqtraits){\n");
   fprintf(fp, "localdata <- data[alltraits==trait,]\n");
   fprintf(fp, "Pos <- localdata$Pos\n");
   fprintf(fp, "baseStop <- c()\n");
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "  Pos[localdata$Chr==i] <- localdata$Pos[localdata$Chr==i] + base\n");
   fprintf(fp, "  base <- base + max(0, localdata$Pos[localdata$Chr==i])\n");
   fprintf(fp, "  baseStop <- c(baseStop, base)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "baseStart <- c(0, baseStop[-%d])\n", SEXCHR-1);
   fprintf(fp, "baseMiddle <- (baseStart + baseStop)/2\n");
   fprintf(fp, "LOD=localdata$LOD\n");
   fprintf(fp, "plot(Pos, LOD, type=\"l\", xlab=\"Chromosome\", ylab = \"LOD Score\", xaxt = 'n',\n");
   fprintf(fp, "  cex = 1.5, cex.axis=1.5, cex.main=1.5, cex.lab = 1.5, cex.sub = 1.5,\n");
   fprintf(fp, "  ylim=c(0, min(max(LOD),10)),\n");
   fprintf(fp, "  lwd = 1.5, main = paste(\"Haseman-Elston Regression for Trait\" , trait))\n", (const char*)prefix);
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "oldlocaldata <- localdata\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "localdata <- oldlocaldata[oldlocaldata$Chr==i,]\n");
   fprintf(fp, "repeat{\n");
   fprintf(fp, "maxLOD <- max(0,localdata$LOD)\n");
   fprintf(fp, "if(maxLOD >= 2.2){\n");
   fprintf(fp, "localpos.max <- median(localdata$Pos[localdata$LOD==maxLOD])\n");
   fprintf(fp, "text(base+localpos.max, min(maxLOD,10)-0.1, paste(i, \":\", localpos.max, sep=\"\"), col=\"red\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "if(maxLOD >= 2.2) localdata$LOD[localdata$Pos > localpos.max - 10 & localdata$Pos < localpos.max + 10] <- 0;\n");
   fprintf(fp, "if(maxLOD < 2.2) break\n");
   fprintf(fp, "}\n");
   fprintf(fp, "  base <- base + max(0, localdata$Pos)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "axis(1, labels=FALSE, line = 0, at = c(0,baseStop))\n");
   fprintf(fp, "axis(1, labels=c(1:%d), line = -0.5, lty=\"blank\", at = baseMiddle)\n", SEXCHR-1);
   fprintf(fp, "abline(h = 3.6, lty = 2, col=\"red\")\n");
   fprintf(fp, "abline(h = 2.2, lty = 3, col=\"green\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "dev.off()\n");
   fprintf(fp, "}\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("--hereg is done but R code %s failed.\n\n", (const char*)scriptfile);
   else if(error){
      printf("--hereg is done but R code %s failed.\n", (const char*)scriptfile);
      printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
   }else{
      sprintf(command, "ps2pdf %s_herplot.ps", (const char*)prefix);
      system(command);
      printf("Haseman-Elston regression plots are generated in %s_herplot.pdf\n", (const char*)prefix);
   }
}

void plotPopROH(const char *prefix, int SEXCHR)
{
   String scriptfile=prefix;
   scriptfile.Add("_poprohplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --poproh, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_poprohplot.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix);
   fprintf(fp, "data<-c()\n");
   fprintf(fp, "if(file.exists(\"%s.rohdiff\")) data <- read.table(file=\"%s.rohdiff\", header=T)\n",
      (const char*)prefix,(const char*)prefix);
   fprintf(fp, "if(length(data$LOD)>0){\n");
   fprintf(fp, "N <- max(data$Pop_Pos)\n");
   fprintf(fp, "label.pop <- 1:N\n");
   fprintf(fp, "for(p1 in 1:(N-1)){\n");
   fprintf(fp, "  for(p2 in (p1+1):N){\n");
   fprintf(fp, "localdata <- data[data$Pop_Neg == p1 & data$Pop_Pos == p2,]\n");
   fprintf(fp, "Pos <- localdata$Pos\n");
   fprintf(fp, "baseStop <- c()\n");
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "  Pos[localdata$Chr==i] <- localdata$Pos[localdata$Chr==i] + base\n");
   fprintf(fp, "  base <- base + max(0, localdata$Pos[localdata$Chr==i])\n");
   fprintf(fp, "  baseStop <- c(baseStop, base)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "baseStart <- c(0, baseStop[-%d])\n", SEXCHR-1);
   fprintf(fp, "baseMiddle <- (baseStart + baseStop)/2\n");
   fprintf(fp, "LOD=localdata$LOD\n");
   fprintf(fp, "plot(Pos, LOD, type=\"l\", xlab=\"Chromosome\", ylab = \"LOD Score\", xaxt = 'n',\n");
   fprintf(fp, "  cex = 1.5, cex.axis=1.5, cex.main=1.5, cex.lab = 1.5, cex.sub = 1.5, lwd = 1.5,\n");
   fprintf(fp, "  main = paste(\"ROH difference in \", label.pop[p1], \" (-) and \", label.pop[p2], \" (+)\", sep=\"\"),\n");
   fprintf(fp, "  ylim=c(max(min(LOD),-10), min(max(LOD),10)))\n");
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "maxLOD <- max(0,localdata$LOD[localdata$Chr==i])\n");
   fprintf(fp, "if(maxLOD >= 3.6){\n");
   fprintf(fp, "localpos <- median(localdata$Pos[localdata$Chr==i][localdata$LOD[localdata$Chr==i]==maxLOD])\n");
   fprintf(fp, "text(base+localpos, min(maxLOD,10)-0.1, paste(i, \":\", localpos, sep=\"\"), col=\"red\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "minLOD <- min(0,localdata$LOD[localdata$Chr==i])\n");
   fprintf(fp, "if(minLOD <= -3.6){\n");
   fprintf(fp, "localpos <- median(localdata$Pos[localdata$Chr==i][localdata$LOD[localdata$Chr==i]==minLOD])\n");
   fprintf(fp, "text(base+localpos, max(minLOD,-10)+0.1, paste(i, \":\", localpos, sep=\"\"), col=\"red\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "  base <- base + max(0, localdata$Pos[localdata$Chr==i])\n");
   fprintf(fp, "}\n");
   fprintf(fp, "axis(1, labels=FALSE, line = 0, at = c(0,baseStop))\n");
   fprintf(fp, "axis(1, labels=c(1:%d), line = -0.5, lty=\"blank\", at = baseMiddle)\n", SEXCHR-1);
   fprintf(fp, "abline(h = 3.6, lty = 2, col=\"red\")\n");
   fprintf(fp, "abline(h = -3.6, lty = 2, col=\"red\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "}\n");
   fprintf(fp, "dev.off()\n");
   fprintf(fp, "}\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("--poproh is done but R code %s failed.\n\n", (const char*)scriptfile);
   else if(error){
      printf("--poproh is done but R code %s failed.\n", (const char*)scriptfile);
      printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
   }else{
      sprintf(command, "ps2pdf %s_poprohplot.ps", (const char*)prefix);
      system(command);
      printf("ROH difference between populations is generated in %s_poprohplot.pdf\n", (const char*)prefix);
   }
}

void plotROHforQT(const char *prefix, int SEXCHR)
{
   String scriptfile=prefix;
   scriptfile.Add("_mthomoplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --mthomo, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_mthomoplot.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix);
   fprintf(fp, "data<-c()\n");
   fprintf(fp, "if(file.exists(\"%s.mthomo\")) data <- read.table(file=\"%s.mthomo\", header=T)\n",
      (const char*)prefix,(const char*)prefix);
   fprintf(fp, "if(length(data$LOD)>0){\n");
   fprintf(fp, "alltraits <- as.character(data$Trait)\n");
   fprintf(fp, "uniqtraits <- unique(alltraits)\n");
   fprintf(fp, "for(trait in uniqtraits){\n");
   fprintf(fp, "localdata <- data[alltraits==trait,]\n");
   fprintf(fp, "Pos <- localdata$Pos\n");
   fprintf(fp, "baseStop <- c()\n");
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "  Pos[localdata$Chr==i] <- localdata$Pos[localdata$Chr==i] + base\n");
   fprintf(fp, "  base <- base + max(0, localdata$Pos[localdata$Chr==i])\n");
   fprintf(fp, "  baseStop <- c(baseStop, base)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "baseStart <- c(0, baseStop[-%d])\n", SEXCHR-1);
   fprintf(fp, "baseMiddle <- (baseStart + baseStop)/2\n");
   fprintf(fp, "LOD=localdata$LOD\n");
   fprintf(fp, "plot(Pos, LOD, type=\"l\", xlab=\"Chromosome\", ylab = \"LOD Score\", xaxt = 'n',\n");
   fprintf(fp, "  cex = 1.5, cex.axis=1.5, cex.main=1.5, cex.lab = 1.5, cex.sub = 1.5,\n");
   fprintf(fp, "  ylim=c(max(min(LOD),-10), min(max(LOD),10)),\n");
   fprintf(fp, "  lwd = 1.5, main = paste(\"Homozygosity Mapping of\" , trait))\n", (const char*)prefix);
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "oldlocaldata <- localdata\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "localdata <- oldlocaldata[oldlocaldata$Chr==i,]\n");
   fprintf(fp, "repeat{\n");
   fprintf(fp, "maxLOD <- max(0,localdata$LOD)\n");
   fprintf(fp, "if(maxLOD >= 2.2){\n");
   fprintf(fp, "localpos.max <- median(localdata$Pos[localdata$LOD==maxLOD])\n");
   fprintf(fp, "text(base+localpos.max, min(maxLOD,10)-0.1, paste(i, \":\", localpos.max, sep=\"\"), col=\"red\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "minLOD <- min(0,localdata$LOD)\n");
   fprintf(fp, "if(minLOD <= -2.2){\n");
   fprintf(fp, "localpos.min <- median(localdata$Pos[localdata$LOD==minLOD])\n");
   fprintf(fp, "text(base+localpos.min, max(minLOD,-10)+0.1, paste(i, \":\", localpos.min, sep=\"\"), col=\"red\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "if(maxLOD >= 2.2) localdata$LOD[localdata$Pos > localpos.max - 10 & localdata$Pos < localpos.max + 10] <- 0;\n");
   fprintf(fp, "if(minLOD <= -2.2) localdata$LOD[localdata$Pos > localpos.min - 10 & localdata$Pos < localpos.min + 10] <- 0;\n");
   fprintf(fp, "if(maxLOD < 2.2 && minLOD > -2.2) break\n");
   fprintf(fp, "}\n");
   fprintf(fp, "  base <- base + max(0, localdata$Pos)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "axis(1, labels=FALSE, line = 0, at = c(0,baseStop))\n");
   fprintf(fp, "axis(1, labels=c(1:%d), line = -0.5, lty=\"blank\", at = baseMiddle)\n", SEXCHR-1);
   fprintf(fp, "abline(h = 3.6, lty = 2, col=\"red\")\n");
   fprintf(fp, "abline(h = 2.2, lty = 3, col=\"green\")\n");
   fprintf(fp, "abline(h = -3.6, lty = 2, col=\"red\")\n");
   fprintf(fp, "abline(h = -2.2, lty = 3, col=\"green\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "dev.off()\n");
   fprintf(fp, "}\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("--mthomo is done but R code %s failed.\n\n", (const char*)scriptfile);
   else if(error){
      printf("--mthomo is done but R code %s failed.\n", (const char*)scriptfile);
      printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
   }else{
      sprintf(command, "ps2pdf %s_mthomoplot.ps", (const char*)prefix);
      system(command);
      printf("Homozygosity mapping plots are generated in %s_mthomoplot.pdf\n", (const char*)prefix);
   }
}

void plotROHmapping(const char *prefix, const char *stratName, int SEXCHR)
{
   String postfix = stratName[0]=='\0'? "homomap": "homomapMH";
   String scriptfile=prefix;
   scriptfile.Add("_");
   scriptfile.Add(postfix);
   scriptfile.Add("plot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --homomap, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_%splot.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix, (const char*)postfix);
   fprintf(fp, "data<-c()\n");
   fprintf(fp, "if(file.exists(\"%s.%s\")) data <- read.table(file=\"%s.%s\", header=T)\n",
      (const char*)prefix, (const char*)postfix,
      (const char*)prefix, (const char*)postfix);
   fprintf(fp, "if(length(data$LOD)>0){\n");
   fprintf(fp, "Pos <- data$Pos\n");
   fprintf(fp, "baseStop <- c()\n");
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "  Pos[data$Chr==i] <- data$Pos[data$Chr==i] + base\n");
   fprintf(fp, "  base <- base + max(0, data$Pos[data$Chr==i])\n");
   fprintf(fp, "  baseStop <- c(baseStop, base)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "baseStart <- c(0, baseStop[-%d])\n", SEXCHR-1);
   fprintf(fp, "baseMiddle <- (baseStart + baseStop)/2\n");
   fprintf(fp, "plot(Pos, data$LOD, type=\"l\", xlab=\"Chromosome\", ylab = \"LOD Score\", xaxt = 'n',\n");
   fprintf(fp, "  cex = 1.5, cex.axis=1.5, cex.main=1.5, cex.lab = 1.5, cex.sub = 1.5,\n");
   fprintf(fp, "  ylim=c(max(min(data$LOD),-5), min(max(data$LOD),10)),\n");
   fprintf(fp, "  lwd = 1.5, main = \"Homozygosity Mapping in %s\")\n", (const char*)prefix);
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "localdata <- data[data$Chr==i,]\n");
   fprintf(fp, "repeat{\n");
   fprintf(fp, "maxLOD <- max(0,localdata$LOD)\n");
   fprintf(fp, "if(maxLOD >= 3.6){\n");
   fprintf(fp, "localpos.max <- median(localdata$Pos[localdata$LOD==maxLOD])\n");
   fprintf(fp, "text(base+localpos.max, min(maxLOD,10)-0.1, paste(i, \":\", localpos.max, sep=\"\"), col=\"red\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "minLOD <- min(0,localdata$LOD)\n");
   fprintf(fp, "if(minLOD <= -3.6){\n");
   fprintf(fp, "localpos.min <- median(localdata$Pos[localdata$LOD==minLOD])\n");
   fprintf(fp, "text(base+localpos.min, max(minLOD,-5)+0.1, paste(i, \":\", localpos.min, sep=\"\"), col=\"red\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "if(maxLOD >= 3.6) localdata$LOD[localdata$Pos > localpos.max - 10 & localdata$Pos < localpos.max + 10] <- 0;\n");
   fprintf(fp, "if(minLOD <= -3.6) localdata$LOD[localdata$Pos > localpos.min - 10 & localdata$Pos < localpos.min + 10] <- 0;\n");
   fprintf(fp, "if(maxLOD < 3.6 && minLOD > -3.6) break\n");
   fprintf(fp, "}\n");
   fprintf(fp, "  base <- base + max(0, localdata$Pos)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "axis(1, labels=FALSE, line = 0, at = c(0,baseStop))\n");
   fprintf(fp, "axis(1, labels=c(1:%d), line = -0.5, lty=\"blank\", at = baseMiddle)\n", SEXCHR-1);
   fprintf(fp, "abline(h = 3.6, lty = 2, col=\"red\")\n");
   fprintf(fp, "abline(h = 2.2, lty = 3, col=\"green\")\n");
   fprintf(fp, "abline(h = -3.6, lty = 2, col=\"red\")\n");
   fprintf(fp, "dev.off()\n");
   fprintf(fp, "}\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("--homomap is done but R code %s failed.\n\n", (const char*)scriptfile);
   else if(error){
      printf("--homomap is done but R code %s failed.\n", (const char*)scriptfile);
      printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
   }else{
      sprintf(command, "ps2pdf %s_%splot.ps",
         (const char*)prefix, (const char*)postfix);
      system(command);
      printf("Homozygosity mapping plot is generated in %s_%splot.pdf\n",
         (const char*)prefix, (const char*)postfix);
   }
}

void plotPopDist(const char *prefix)
{
   String scriptfile=prefix;
   scriptfile.Add("_popdistplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --popdist, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_popdistplot.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix);
   fprintf(fp, "data<-c()\n");
   fprintf(fp, "if(file.exists(\"%s.dst\")) data <- read.table(file=\"%s.dst\", header=T)\n",
      (const char*)prefix,(const char*)prefix);
   fprintf(fp, "if(length(data$DistIBD)>0){\n");
   fprintf(fp, "N <- sqrt(2*length(data$DistIBD)+0.25)-0.5\n");
   fprintf(fp, "M <- matrix(0,N,N)\n");
   fprintf(fp, "M[upper.tri(M, diag=TRUE)] <- data$DistIBD\n");
   fprintf(fp, "M[lower.tri(M, diag=TRUE)] <- data$DistIBD\n");
   fprintf(fp, "plot(hclust(as.dist(M)),main=\"Clustering of %s Populations By Average IBD\",xlab=\"\",ylab=\"1 - IBD Proportion\")\n",
      (const char*)prefix);
   fprintf(fp, "mds <- cmdscale(M)\n");
   fprintf(fp, "plot(mds[,1],mds[,2],type=\"n\",xlab=\"\",ylab=\"\", axes=FALSE,main=\"MDS of %s Populations Using Average IBD\")\n",
      (const char*)prefix);
   fprintf(fp, "text(mds[,1],mds[,2],1:N)\n");
   fprintf(fp, "dev.off()\n");
   fprintf(fp, "}\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("--popdist is done but R code %s failed.\n\n", (const char*)scriptfile);
   else if(error){
      printf("--popdist is done but R code %s failed.\n", (const char*)scriptfile);
      printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
   }else{
      sprintf(command, "ps2pdf %s_popdistplot.ps", (const char*)prefix);
      system(command);
      printf("Population distance plot is generated in %s_popdistplot.pdf\n", (const char*)prefix);
   }
}

void plotIBDmapping(const char *prefix, int SEXCHR)
{
   String scriptfile=prefix;
   scriptfile.Add("_ibdmapplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --ibdmap, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_ibdmapplot.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix);
   fprintf(fp, "data<-c()\n");
   fprintf(fp, "if(file.exists(\"%s.ibdmap\")) data <- read.table(file=\"%s.ibdmap\", header=T)\n",
      (const char*)prefix,(const char*)prefix);
   fprintf(fp, "if(length(data$P)>0){\n");
   fprintf(fp, "Pos <- data$PosMb\n");
   fprintf(fp, "baseStop <- c()\n");
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "  Pos[data$Chr==i] <- data$PosMb[data$Chr==i] + base\n");
   fprintf(fp, "  base <- base + max(0, data$PosMb[data$Chr==i])\n");
   fprintf(fp, "  baseStop <- c(baseStop, base)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "baseStart <- c(0, baseStop[-%d])\n", SEXCHR-1);
   fprintf(fp, "baseMiddle <- (baseStart + baseStop)/2\n");
   fprintf(fp, "logP <- rep(10, length(data$P))\n");
   fprintf(fp, "logP[data$P>0] <- -log(data$P[data$P>0])/log(10)\n");
   fprintf(fp, "plot(Pos, logP, type=\"l\", xlab=\"Chromosome\", ylab = expression(paste(\"-\", log[10],\"P\",sep=\"\")), xaxt = 'n',\n");
   fprintf(fp, "  cex = 1.2, cex.axis=1.2, cex.main=1.5, cex.lab = 1.2, cex.sub = 1.2, ylim=c(0,max(c(logP[logP<10]),4)),\n");
   fprintf(fp, "  lwd = 1.2, main = \"Permutation P Values of IBD Mapping in %s\")\n", (const char*)prefix);
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "minP <- min(1,data$P[data$Chr==i])\n");
   fprintf(fp, "if(minP < 0.0001){\n");
   fprintf(fp, "localpos <- median(data$Pos[data$Chr==i][data$P[data$Chr==i]==minP])\n");
   fprintf(fp, "text(base+localpos, -log(minP, 10)-0.1, paste(i, \":\", localpos, sep=\"\"), col=\"red\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "  base <- base + max(0, data$Pos[data$Chr==i])\n");
   fprintf(fp, "}\n");
   fprintf(fp, "axis(1, labels=FALSE, line = 0, at = c(0,baseStop))\n");
   fprintf(fp, "axis(1, labels=c(1:%d), line = -0.5, lty=\"blank\", at = baseMiddle)\n", SEXCHR-1);
   fprintf(fp, "abline(h = 4, lty = 2, col=\"red\")\n");
   fprintf(fp, "dev.off()\n");
   fprintf(fp, "}\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("--ibdmap is done but R code %s failed.\n\n", (const char*)scriptfile);
   else if(error){
      printf("--ibdmap is done but R code %s failed.\n", (const char*)scriptfile);
      printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
   }else{
      sprintf(command, "ps2pdf %s_ibdmapplot.ps", (const char*)prefix);
      system(command);
      printf("IBD mapping plot is generated in %s_ibdmapplot.pdf\n", (const char*)prefix);
   }
}

void plotAncestry(const char *prefix)
{
   String scriptfile=prefix;
   scriptfile.Add("_ancestryplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --ancestry, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_ancestryplot.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix);
   fprintf(fp, "data<-c()\n");
   fprintf(fp, "if(file.exists(\"%s.anc\")) data <- read.table(file=\"%s.anc\", header=T)\n",
      (const char*)prefix,(const char*)prefix);
   fprintf(fp, "if(length(data$Admix)>0){\n");
   fprintf(fp, "  plot(data$Admix, data$Ancestry, xlab=\"Admixed Proportion\", ylab = \"Ancestry\",\n");
   fprintf(fp, "  cex = 1.5, cex.axis=1.5, cex.main=1.5, cex.lab = 1.5, cex.sub = 1.5,\n");
   fprintf(fp, "  lwd = 1.5, main = \"Ancestry in %s\")\n", (const char*)prefix);
   fprintf(fp, "  abline(a=0, b=0.5, col=\"red\")\n");
   fprintf(fp, "  abline(a=1, b=-0.5, col=\"red\")\n");
   fprintf(fp, "  y <- seq(0,1,0.001)\n");
   fprintf(fp, "  x <- y*(1-y)*2\n");
   fprintf(fp, "  lines(x, y, col=\"red\")\n");
   fprintf(fp, "  plot(data$Anc_P1, data$Anc_P2, xlab=\"Ancestry of Parent 1\", ylab = \"Ancestry of Parent 2\",\n");
   fprintf(fp, "  cex = 1.5, cex.axis=1.5, cex.main=1.5, cex.lab = 1.5, cex.sub = 1.5,\n");
   fprintf(fp, "  lwd = 1.5, main = \"Parental Ancestry in %s\")\n", (const char*)prefix);
   fprintf(fp, "}\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("--ancestry is done but R code %s failed.\n\n", (const char*)scriptfile);
   else if(error){
      printf("--ancestry is done but R code %s failed.\n", (const char*)scriptfile);
      printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
   }else{
      sprintf(command, "ps2pdf %s_ancestryplot.ps", (const char*)prefix);
      system(command);
      printf("Ancestry plot is generated in %s_ancestryplot.pdf\n", (const char*)prefix);
   }
}

void plotNPL(const char *prefix, int SEXCHR)
{
   String scriptfile=prefix;
   scriptfile.Add("_nplplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --npl, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_nplplot.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix);
   fprintf(fp, "data<-c()\n");
   fprintf(fp, "if(file.exists(\"%s.npl\")) data <- read.table(file=\"%s.npl\", header=T)\n",
      (const char*)prefix,(const char*)prefix);
   fprintf(fp, "if(length(data$LODwDSP)>0){\n");
   fprintf(fp, "Pos <- data$Pos\n");
   fprintf(fp, "baseStop <- c()\n");
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "  Pos[data$Chr==i] <- data$Pos[data$Chr==i] + base\n");
   fprintf(fp, "  base <- base + max(0,data$Pos[data$Chr==i])\n");
   fprintf(fp, "  baseStop <- c(baseStop, base)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "baseStart <- c(0, baseStop[-%d])\n", SEXCHR-1);
   fprintf(fp, "baseMiddle <- (baseStart + baseStop)/2\n");
   fprintf(fp, "plot(Pos, data$LODwDSP, type=\"l\", xlab=\"Chromosome\", ylab = \"LOD Score\", xaxt = 'n',\n");
   fprintf(fp, "  cex = 1.5, cex.axis=1.5, cex.main=1.5, cex.lab = 1.5, cex.sub = 1.5,\n");
   fprintf(fp, "  lwd = 1.5, main = \"ASP/DSP NPL Scan in %s\", ylim=c(0,min(max(data$LODwDSP),10)))\n", (const char*)prefix);
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "maxLOD <- max(0,data$LODwDSP[data$Chr==i])\n");
   fprintf(fp, "if(maxLOD >= 3.6){\n");
   fprintf(fp, "localpos <- median(data$Pos[data$Chr==i][data$LODwDSP[data$Chr==i]==maxLOD])\n");
   fprintf(fp, "text(base+localpos, min(maxLOD,10)-0.1, paste(i, \":\", localpos, sep=\"\"), col=\"red\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "  base <- base + max(0, data$Pos[data$Chr==i])\n");
   fprintf(fp, "}\n");
   fprintf(fp, "axis(1, labels=FALSE, line = 0, at = c(0,baseStop))\n");
   fprintf(fp, "axis(1, labels=c(1:%d), line = -0.5, lty=\"blank\", at = baseMiddle)\n", SEXCHR-1);
   fprintf(fp, "abline(h = 3.6, lty = 2, col=\"red\")\n");
   fprintf(fp, "abline(h = 2.2, lty = 3, col=\"green\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "Pos <- data$Pos\n");
   fprintf(fp, "baseStop <- c()\n");
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "  Pos[data$Chr==i] <- data$Pos[data$Chr==i] + base\n");
   fprintf(fp, "  base <- base + max(0,data$Pos[data$Chr==i])\n");
   fprintf(fp, "  baseStop <- c(baseStop, base)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "baseStart <- c(0, baseStop[-%d])\n", SEXCHR-1);
   fprintf(fp, "baseMiddle <- (baseStart + baseStop)/2\n");
   fprintf(fp, "plot(Pos, data$LOD_ASP, type=\"l\", xlab=\"Chromosome\", ylab = \"LOD Score\", xaxt = 'n',\n");
   fprintf(fp, "  cex = 1.5, cex.axis=1.5, cex.main=1.5, cex.lab = 1.5, cex.sub = 1.5,\n");
   fprintf(fp, "  lwd = 1.5, main = \"ASP NPL Scan in %s\", ylim=c(0,min(max(data$LOD_ASP),10)))\n", (const char*)prefix);
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "maxLOD <- max(0,data$LOD_ASP[data$Chr==i])\n");
   fprintf(fp, "if(maxLOD >= 3.6){\n");
   fprintf(fp, "localpos <- median(data$Pos[data$Chr==i][data$LOD_ASP[data$Chr==i]==maxLOD])\n");
   fprintf(fp, "text(base+localpos, min(maxLOD,10)-0.1, paste(i, \":\", localpos, sep=\"\"), col=\"red\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "  base <- base + max(0, data$Pos[data$Chr==i])\n");
   fprintf(fp, "}\n");
   fprintf(fp, "axis(1, labels=FALSE, line = 0, at = c(0,baseStop))\n");
   fprintf(fp, "axis(1, labels=c(1:%d), line = -0.5, lty=\"blank\", at = baseMiddle)\n", SEXCHR-1);
   fprintf(fp, "abline(h = 3.6, lty = 2, col=\"red\")\n");
   fprintf(fp, "abline(h = 2.2, lty = 3, col=\"green\")\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("--npl is done but R code %s failed.\n\n", (const char*)scriptfile);
   else if(error){
      printf("--npl is done but R code %s failed.\n", (const char*)scriptfile);
      printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
   }else{
      sprintf(command, "ps2pdf %s_nplplot.ps", (const char*)prefix);
      system(command);
      printf("ASP NPL plot is generated in %s_nplplot.pdf\n", (const char*)prefix);
   }
}

void plotAUCmapping(const char *prefix, int SEXCHR)
{
   String scriptfile=prefix;
   scriptfile.Add("_aucmapplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --aucmap, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_aucmapplot.ps\", paper=\"letter\", horizontal=T)\n",(const char*)prefix);
   fprintf(fp, "data<-c()\n");
   fprintf(fp, "if(file.exists(\"%s.aucmap\")) data <- read.table(file=\"%s.aucmap\", header=T)\n",
      (const char*)prefix,(const char*)prefix);
   fprintf(fp, "if(length(data$AUC)>0){\n");
   fprintf(fp, "Pos <- data$Pos\n");
   fprintf(fp, "baseStop <- c()\n");
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "  Pos[data$Chr==i] <- data$Pos[data$Chr==i] + base\n");
   fprintf(fp, "  base <- base + max(data$Pos[data$Chr==i])\n");
   fprintf(fp, "  baseStop <- c(baseStop, base)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "baseStart <- c(0, baseStop[-%d])\n", SEXCHR-1);
   fprintf(fp, "baseMiddle <- (baseStart + baseStop)/2\n");
   fprintf(fp, "plot(Pos[data$Success>0.9], data$AUC[data$Success>0.9], type=\"l\", xlab=\"Chromosome\", ylab = \"AUC\", xaxt = 'n',\n");
   fprintf(fp, "  cex = 1.5, cex.axis=1.5, cex.main=1.5, cex.lab = 1.5, cex.sub = 1.5, ylim=c(0.5,max(data$AUC[data$Success>0.9])),\n");
   fprintf(fp, "  lwd = 1.5, main = \"Risk Prediction Using Position-Specific IBD Relatives in %s\")\n",
      (const char*)prefix);
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "maxAUC <- max(0,data$AUC[data$Chr==i])\n");
   fprintf(fp, "if(maxAUC > 0.6){\n");
   fprintf(fp, "localpos <- median(data$Pos[data$Chr==i][data$AUC[data$Chr==i]==maxAUC])\n");
   fprintf(fp, "text(base+localpos, maxAUC-0.01, paste(i, \":\", localpos, sep=\"\"), col=\"red\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "  base <- base + max(0, data$Pos[data$Chr==i])\n");
   fprintf(fp, "}\n");
   fprintf(fp, "axis(1, labels=FALSE, line = 0, at = c(0,baseStop))\n");
   fprintf(fp, "axis(1, labels=c(1:%d), line = -0.5, lty=\"blank\", at = baseMiddle)\n", SEXCHR-1);
   fprintf(fp, "abline(h = 0.6, lty = 2, col=\"red\")\n");
   fprintf(fp, "dev.off()\n");
   fprintf(fp, "}\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("--aucmap is done but R code %s failed.\n\n", (const char*)scriptfile);
   else if(error){
      printf("--aucmap is done but R code %s failed.\n", (const char*)scriptfile);
      printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
   }else{
      sprintf(command, "ps2pdf %s_aucmapplot.ps", (const char*)prefix);
      system(command);
      printf("AUC mapping plot is generated in %s_aucmapplot.pdf\n", (const char*)prefix);
   }
}

void plotRelationship(const char *prefix)
{
   String scriptfile=prefix;
   scriptfile.Add("_relplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --related, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_relplot.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix);
   fprintf(fp, "data<-c()\n");
   fprintf(fp, "if(file.exists(\"%s.kin\")) data <- read.table(file=\"%s.kin\", header=T)\n",
      (const char*)prefix,(const char*)prefix);
   fprintf(fp, "allcolors <- c(\"purple\", \"red\", \"green\", \"blue\", \"magenta\", \"gold\", \"black\")\n");
   fprintf(fp, "allrelatives <- c(\"MZ Twin\", \"Parent-Offspring\", \"Full Siblings\", \"2nd-Degree\", \"3rd-Degree\", \"More Distant\", \"Unrelated\")\n");
   fprintf(fp, "if(length(data$IBD1Seg)>0 & length(data$IBD2Seg)>0){\n");
   // Page 1 plot: IBD1 vs IBD2
   fprintf(fp, "allpair <- data$PropIBD>0 | data$Kinship>0.04419\n");
   fprintf(fp, "d0 <- data$Phi==0.5\n");
   fprintf(fp, "d1.PO <- data$Phi==0.25 & data$Z0==0\n");
   fprintf(fp, "d1.FS <- data$Phi==0.25 & data$Z0>0\n");
   fprintf(fp, "d2 <- data$Phi>0.08839 & data$Phi<=0.17678\n");
   fprintf(fp, "d3 <- data$Phi>0.04419 & data$Phi<=0.08839\n");
   fprintf(fp, "dO <- data$Phi>0 & data$Phi<=0.04419\n");
   fprintf(fp, "dU <- data$Phi==0 & allpair\n");
   fprintf(fp, "plot(data$IBD1Seg[dU], data$IBD2Seg[dU], type=\"p\", col = \"black\", cex.lab=1.2,\n");
   fprintf(fp, "xlim=c(min(data$IBD1Seg[allpair]), max(data$IBD1Seg[allpair])),\n");
   fprintf(fp, "ylim=c(min(data$IBD2Seg[allpair]), max(data$IBD2Seg[allpair])),\n");
   fprintf(fp, "main = \"IBD Segments In %s Families\",\n", (const char*)prefix);
   fprintf(fp, "xlab=expression(paste(\"Length Proportion of IBD1 Segments (\", pi[1], \")\", sep=\"\")),\n");
   fprintf(fp, "ylab=expression(paste(\"Length Proportion of IBD2 Segments (\", pi[2], \")\", sep=\"\")))\n");
   fprintf(fp, "points(data$IBD1Seg[d0], data$IBD2Seg[d0], col=\"purple\")\n");
   fprintf(fp, "points(data$IBD1Seg[d1.PO], data$IBD2Seg[d1.PO], col=\"red\")\n");
   fprintf(fp, "points(data$IBD1Seg[d1.FS], data$IBD2Seg[d1.FS], col=\"green\")\n");
   fprintf(fp, "points(data$IBD1Seg[d2], data$IBD2Seg[d2], col=\"blue\")\n");
   fprintf(fp, "points(data$IBD1Seg[d3], data$IBD2Seg[d3], col=\"magenta\")\n");
   fprintf(fp, "points(data$IBD1Seg[dO], data$IBD2Seg[dO], col=\"gold\")\n");
   fprintf(fp, "points(data$IBD1Seg[dU & data$PropIBD>0.08838835], data$IBD2Seg[dU & data$PropIBD>0.08838835], col=\"black\")\n");
   fprintf(fp, "points(data$IBD1Seg[d1.FS & data$IBD2Seg<0.08], data$IBD2Seg[d1.FS & data$IBD2Seg<0.08], col=\"green\")\n");
   fprintf(fp, "points(data$IBD1Seg[d1.PO & data$IBD1Seg+data$IBD2Seg<0.9], data$IBD2Seg[d1.PO & data$IBD1Seg+data$IBD2Seg<0.9], col=\"red\")\n");
   fprintf(fp, "abline(h = 0.08, col = \"green\", lty = 3, lwd = 2)\n");   // FS && MZ
   fprintf(fp, "abline(a = 0.96, b = -1, col = \"red\", lty = 3, lwd = 2)\n");   // PO
   fprintf(fp, "abline(a = 0.3535534, b = -0.5, col = \"green\", lty = 3, lwd = 2)\n");// FS
   fprintf(fp, "abline(a = 0.1767767, b = -0.5, col = \"blue\", lty = 3, lwd = 2)\n");// 2nd/3rd
   fprintf(fp, "abline(a = 0.08838835, b = -0.5, col = \"magenta\", lty = 3, lwd = 2)\n");// 3rd/4th
   fprintf(fp, "abline(a = 0.04419, b = -0.5, col = \"gold\", lty = 3, lwd = 2)\n");// 4th/UN
   fprintf(fp, "legend(\"topright\", allrelatives, col=allcolors, text.col=allcolors, pch=19, cex=1.2)\n");

   // Page 2 plot: Kinship vs. PropIBD
   fprintf(fp, "allpair <- data$PropIBD>0 | data$Kinship>0.04419\n");
   fprintf(fp, "d0 <- data$Phi==0.5\n");
   fprintf(fp, "d1.PO <- data$Phi==0.25 & data$Z0==0\n");
   fprintf(fp, "d1.FS <- data$Phi==0.25 & data$Z0>0\n");
   fprintf(fp, "d2 <- data$Phi>0.08839 & data$Phi<=0.17678\n");
   fprintf(fp, "d3 <- data$Phi>0.04419 & data$Phi<=0.08839\n");
   fprintf(fp, "dO <- data$Phi>0 & data$Phi<=0.04419\n");
   fprintf(fp, "dU <- data$Phi==0 & allpair\n");
   fprintf(fp, "plot(data$PropIBD[dU], data$Kinship[dU], type=\"p\", col = \"black\", cex.lab=1.2,\n");
   fprintf(fp, "xlim=c(min(data$PropIBD[allpair]), max(data$PropIBD[allpair])),\n");
   fprintf(fp, "ylim=c(min(data$Kinship[allpair]), max(data$Kinship[allpair])),\n");
   fprintf(fp, "main = paste(\"Kinship vs Proportion IBD (Corr=\", round(cor(data$Kinship[data$PropIBD>0], data$PropIBD[data$PropIBD>0]),digit=3),\") in %s Families\",sep=\"\"),\n",
      (const char*)prefix);
   fprintf(fp, "xlab=expression(paste(\"Proportion of Genomes IBD (\", pi,\"=\",pi[2]+pi[1]/2,\")\",sep=\"\")),\n");
   fprintf(fp, "ylab = expression(paste(\"Estimated Kinship Coefficient (\", phi, \")\", sep=\"\")))\n");
   fprintf(fp, "points(data$PropIBD[d0], data$Kinship[d0], col=\"purple\")\n");
   fprintf(fp, "points(data$PropIBD[d1.FS], data$Kinship[d1.FS], col=\"green\")\n");
   fprintf(fp, "points(data$PropIBD[d1.PO], data$Kinship[d1.PO], col=\"red\")\n");
   fprintf(fp, "points(data$PropIBD[d2], data$Kinship[d2], col=\"blue\")\n");
   fprintf(fp, "points(data$PropIBD[d3], data$Kinship[d3], col=\"magenta\")\n");
   fprintf(fp, "points(data$PropIBD[dO], data$Kinship[dO], col=\"gold\")\n");
   fprintf(fp, "points(data$PropIBD[dU & data$Kinship>0.088], data$Kinship[dU & data$Kinship>0.088], col=\"black\")\n");
   fprintf(fp, "abline(h = 0.35355, col = \"purple\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.17678, col = \"green\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.08838, col = \"blue\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.04419, col = \"magenta\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.02210, col = \"gold\", lty = 3)\n");
   fprintf(fp, "abline(a = 0, b = 0.5, lty = 1)\n");
   fprintf(fp, "abline(a = 0, b = 0.7071068, lty = 3)\n");
   fprintf(fp, "abline(a = 0, b = 0.3535534, lty = 3)\n");
   fprintf(fp, "abline(v = 0.70711, col = \"purple\", lty = 3)\n");
   fprintf(fp, "abline(v = 0.35355, col = \"green\", lty = 3)\n");
   fprintf(fp, "abline(v = 0.17678, col = \"blue\", lty = 3)\n");
   fprintf(fp, "abline(v = 0.08838, col = \"magenta\", lty = 3)\n");
   fprintf(fp, "abline(v = 0.04419, col = \"gold\", lty = 3)\n");
   fprintf(fp, "text(x=0.35355, y=0.35355, \"1st\", adj=c(0,1), col=\"green\")\n");
   fprintf(fp, "text(x=0.17678, y=0.17678, \"2nd\", adj=c(0,1), col=\"blue\")\n");
   fprintf(fp, "text(x=0.08839, y=0.08839, \"3rd\", adj=c(0,1), col=\"magenta\")\n");
   fprintf(fp, "text(x=0.04419, y=0.04419, \"4th\", adj=c(0,1), col=\"gold\")\n");
   fprintf(fp, "legend(\"bottomright\", allrelatives, col=allcolors, text.col=allcolors, pch=19, cex=1.2)\n");

   fprintf(fp, "}else if(length(data$Kinship)>0){\n");
   // In absence of IBDSeg, plot 1: Kinship vs HomIBS0
   fprintf(fp, "d0 <- data$Phi==0.5\n");
   fprintf(fp, "d1.PO <- data$Phi==0.25 & data$Z0==0\n");
   fprintf(fp, "d1.FS <- data$Phi==0.25 & data$Z0>0\n");
   fprintf(fp, "d2 <- data$Phi==0.125\n");
   fprintf(fp, "dU <- data$Phi==0\n");
   fprintf(fp, "dO <- !d0 & !d1.PO & !d1.FS & !d2 & !dU\n");
   fprintf(fp, "plot(data$HomIBS0[dU], data$Kinship[dU], type=\"p\", col = \"black\", cex.lab=1.3,\n");
   fprintf(fp, "xlim=c(min(data$HomIBS0), max(data$HomIBS0)),\n");
   fprintf(fp, "ylim=c(min(data$Kinship), max(data$Kinship)),\n");
   fprintf(fp, "main = \"Kinship vs IBS0 in %s Families\",\n", (const char*)prefix);
   fprintf(fp, "xlab=\"Proportion of Zero IBS In Minor Homozygote Pairs\", ylab = \"Estimated Kinship Coefficient\")\n");
   fprintf(fp, "points(data$HomIBS0[d0], data$Kinship[d0], col=\"purple\")\n");
   fprintf(fp, "points(data$HomIBS0[d1.PO], data$Kinship[d1.PO], col=\"red\")\n");
   fprintf(fp, "points(data$HomIBS0[d1.FS], data$Kinship[d1.FS], col=\"green\")\n");
   fprintf(fp, "points(data$HomIBS0[d2], data$Kinship[d2], col=\"blue\")\n");
   fprintf(fp, "points(data$HomIBS0[dO], data$Kinship[dO], col=\"gold\")\n");
   fprintf(fp, "points(data$HomIBS0[dU & data$Kinship>0.088], data$Kinship[dU & data$Kinship>0.088], col=\"black\")\n");
   fprintf(fp, "abline(h = 0.3536, col = \"red\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.1768, col = \"green\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.0884, col = \"blue\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.0442, col = \"black\", lty = 3)\n");
   fprintf(fp, "legend(\"topright\", allrelatives, col=allcolors, text.col=allcolors, pch=19, cex=1.2)\n");
   // in absence of IBDSeg, plot 2: Kinship vs. HetConc
   fprintf(fp, "plot(data$HetConc[dU], data$Kinship[dU], type=\"p\",\n");
   fprintf(fp, "col = \"black\", cex.lab=1.3,\n");
   fprintf(fp, "xlim = c(min(data$HetConc), max(data$HetConc[data$HetConc<0.8])),\n");
   fprintf(fp, "ylim = c(min(data$Kinship), max(data$Kinship[data$HetConc<0.8])),\n");
   fprintf(fp, "main = \"Kinship vs Heterozygote Concordance In %s Families\",\n", (const char*)prefix);
   fprintf(fp, "xlab=\"Heterozygote Concordance Rate\", ylab = \"Estimated Kinship Coefficient\")\n");
   fprintf(fp, "points(data$HetConc[d1.PO], data$Kinship[d1.PO], col=\"red\")\n");
   fprintf(fp, "points(data$HetConc[d1.FS], data$Kinship[d1.FS], col=\"green\")\n");
   fprintf(fp, "points(data$HetConc[d2], data$Kinship[d2], col=\"blue\")\n");
   fprintf(fp, "points(data$HetConc[dO], data$Kinship[dO], col=\"gold\")\n");
   fprintf(fp, "points(data$HetConc[dU & data$Kinship>0.088], data$Kinship[dU & data$Kinship>0.088], col=\"black\")\n");
   fprintf(fp, "abline(h = 0.3536, col = \"red\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.1768, col = \"green\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.0884, col = \"blue\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.0442, col = \"black\", lty = 3)\n");
   fprintf(fp, "legend(\"bottomright\", allrelatives, col=allcolors, text.col=allcolors, pch=19, cex=1.2)\n");
   fprintf(fp, "if(sum(d1.FS)>20){\n");
   fprintf(fp, "y.cut <- 2^-2.5\n");
   fprintf(fp, "x.FS <- quantile(data$HetConc[d1.FS], probs=c(0.1, 0.9))\n");
   fprintf(fp, "y.FS <- quantile(data$Kinship[d1.FS], probs=c(0.1, 0.9))\n");
   fprintf(fp, "slope.FS <- (y.FS[2]-y.FS[1]) / (x.FS[2]-x.FS[1])\n");
   fprintf(fp, "a.FS <- y.FS[1]-x.FS[1]*slope.FS\n");
   fprintf(fp, "abline(a=a.FS, b=slope.FS, col=\"green\")\n");
   fprintf(fp, "x.cut.FS <- (y.cut - a.FS)/slope.FS\n");
   fprintf(fp, "}\n");
   fprintf(fp, "if(sum(d2)>20){\n");
   fprintf(fp, "x.d2 <- quantile(data$HetConc[d2], probs=c(0.9, 0.1))\n");
   fprintf(fp, "y.d2 <- quantile(data$Kinship[d2], probs=c(0.9, 0.1))\n");
   fprintf(fp, "slope.d2 <- (y.d2[2]-y.d2[1]) / (x.d2[2]-x.d2[1])\n");
   fprintf(fp, "a.d2 <- y.d2[1]-x.d2[1]*slope.d2\n");
   fprintf(fp, "abline(a=a.d2, b=slope.d2, col=\"blue\")\n");
   fprintf(fp, "x.cut.d2 <- (y.cut - a.d2)/slope.d2\n");
   fprintf(fp, "x.cut <- sqrt(x.cut.d2 * x.cut.FS)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "if(sum(d1.FS)>20 & sum(d2)>20){\n");
   fprintf(fp, "print(c(x.cut.d2, x.cut.FS, x.cut))\n");
   fprintf(fp, "abline(v=x.cut, col=\"purple\", lty=3)\n");
   fprintf(fp, "text(x=x.cut, y=0.3, labels=round(x.cut,digit=3),col=\"purple\")\n");
   fprintf(fp, "points(data$HetConc[d1.FS & data$HetConc < x.cut], data$Kinship[d1.FS & data$HetConc < x.cut], col=\"green\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "}\n");
   fprintf(fp, "\nif(!file.exists(\"%s.kin0\")) quit()\n", (const char*)prefix);
   fprintf(fp, "data <- read.table(file = \"%s.kin0\", header=T)\n", (const char*)prefix);
   plotIBD1vsIBD2(prefix, fp);

   // Page 2 plot: Kinship vs. PropIBD
   fprintf(fp, "if(length(data$IBD1Seg)>0 & length(data$IBD2Seg)>0){\n");
   fprintf(fp, "d0 <- data$IBD2Seg>0.7\n");
   fprintf(fp, "d1.PO <- (!d0) & data$IBD1Seg+data$IBD2Seg>0.96 | (data$IBD1Seg+data$IBD2Seg>0.9 & data$IBD2Seg<0.08)\n");
   fprintf(fp, "d1.FS <- (!d0) & (!d1.PO) & data$PropIBD>0.35355 & data$IBD2Seg>=0.08\n");
   fprintf(fp, "d2 <- data$PropIBD>0.17678 & data$IBD1Seg+data$IBD2Seg<=0.9 & (!d1.FS)\n");
   fprintf(fp, "d3 <- data$PropIBD>0.08839 & data$PropIBD<=0.17678\n");
   fprintf(fp, "d4 <- data$PropIBD>0.04419 & data$PropIBD<=0.08839\n");
   fprintf(fp, "dU <- data$PropIBD<=0.04419\n");
   fprintf(fp, "plot(data$PropIBD[d1.FS], data$Kinship[d1.FS], type=\"p\", col=\"green\", cex.lab=1.2,\n");
   fprintf(fp, "xlim=c(min(data$PropIBD), max(data$PropIBD)),\n");
   fprintf(fp, "ylim=c(min(data$Kinship), max(data$Kinship)),\n");
   fprintf(fp, "main = paste(\"Kinship vs Proportion IBD (Corr=\", round(cor(data$Kinship[data$PropIBD>0], data$PropIBD[data$PropIBD>0]),digit=3),\") in Inferred %s Relatives\",sep=\"\"),\n",
      (const char*)prefix);
   fprintf(fp, "xlab=expression(paste(\"Proportion of Genomes IBD (\", pi,\"=\",pi[2]+pi[1]/2,\")\",sep=\"\")), \n");
   fprintf(fp, "ylab = expression(paste(\"Estimated Kinship Coefficient (\", phi, \")\", sep=\"\")))\n");
   fprintf(fp, "points(data$PropIBD[d1.PO], data$Kinship[d1.PO], col=\"red\")\n");
   fprintf(fp, "points(data$PropIBD[d0], data$Kinship[d0], col=\"purple\")\n");
   fprintf(fp, "points(data$PropIBD[d2], data$Kinship[d2], col=\"blue\")\n");
   fprintf(fp, "points(data$PropIBD[d3], data$Kinship[d3], col=\"magenta\")\n");
   fprintf(fp, "points(data$PropIBD[d4], data$Kinship[d4], col=\"gold\")\n");
   fprintf(fp, "points(data$PropIBD[dU], data$Kinship[dU], col=\"black\")\n");
   fprintf(fp, "abline(h = 0.35355, col = \"purple\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.17678, col = \"green\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.08839, col = \"blue\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.04419, col = \"magenta\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.02210, col = \"gold\", lty = 3)\n");
   fprintf(fp, "abline(a = 0, b = 0.5, lty = 1)\n");
   fprintf(fp, "abline(a = 0, b = 0.7071068, lty = 3)\n");
   fprintf(fp, "abline(a = 0, b = 0.3535534, lty = 3)\n");
   fprintf(fp, "abline(v = 0.70711, col = \"purple\", lty = 3)\n");
   fprintf(fp, "abline(v = 0.35355, col = \"green\", lty = 3)\n");
   fprintf(fp, "abline(v = 0.17678, col = \"blue\", lty = 3)\n");
   fprintf(fp, "abline(v = 0.08839, col = \"magenta\", lty = 3)\n");
   fprintf(fp, "abline(v = 0.04419, col = \"gold\", lty = 3)\n");
   fprintf(fp, "text(x=0.35355, y=0.35355, \"1st\", adj=c(0,1), col=\"green\")\n");
   fprintf(fp, "text(x=0.17678, y=0.17678, \"2nd\", adj=c(0,1), col=\"blue\")\n");
   fprintf(fp, "text(x=0.08839, y=0.08839, \"3rd\", adj=c(0,1), col=\"magenta\")\n");
   fprintf(fp, "text(x=0.04419, y=0.04419, \"4th\", adj=c(0,1), col=\"gold\")\n");
   fprintf(fp, "legend(\"bottomright\", c(\"Inferred MZ\", \"Inferred PO\", \"Inferred FS\", \"Inferred 2nd\", \"Inferred 3rd\", \"Inferred 4th\", \"Inferred UN\"),\n");
   fprintf(fp, "col=allcolors, text.col = allcolors, pch = 19, cex = 1.2)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("--related is done but R code %s failed.\n\n", (const char*)scriptfile);
   else if(error){
      printf("--related is done but R code %s failed.\n", (const char*)scriptfile);
      printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
   }else{
      sprintf(command, "ps2pdf %s_relplot.ps", (const char*)prefix);
      system(command);
      printf("Relationship plot is generated in %s_relplot.pdf\n", (const char*)prefix);
   }
}

void plotIBDSeg(const char *prefix)
{
   String scriptfile=prefix;
   scriptfile.Add("_ibd1vsibd2.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --ibdseg, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_ibd1vsibd2.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix);
   fprintf(fp, "data<-c()\n");
   fprintf(fp, "if(file.exists(\"%s.seg\")) data <- read.table(file=\"%s.seg\", header=T)\n",
      (const char*)prefix,(const char*)prefix);
   plotIBD1vsIBD2(prefix, fp);
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("--ibdseg is done but R code %s failed.\n\n", (const char*)scriptfile);
   else if(error){
      printf("--ibdseg is done but R code %s failed.\n", (const char*)scriptfile);
      printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
   }else{
      sprintf(command, "ps2pdf %s_ibd1vsibd2.ps", (const char*)prefix);
      system(command);
      printf("IBD1 vs IBD2 plot is generated in %s_ibd1vsibd2.pdf\n", (const char*)prefix);
   }
}

void plotGenderError(const char *prefix, IntArray & plotx, Vector & ploty, IntArray & plotz, double xHeterozygosity, int gendererrorCount)
{
   String genderplotdata=prefix;
   genderplotdata.Add("_gender_autodata.txt");
   FILE *fp = fopen(genderplotdata, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)genderplotdata);
   fprintf(fp, "YCount\txHeterozygosity\tSEX\n");
   for(int i = 0; i < plotx.Length(); i++)
      fprintf(fp, "%d\t%.5lf\t%d\n", plotx[i], ploty[i], plotz[i]);
   fclose(fp);
   String genderplot=prefix;
   genderplot.Add("_gender_autoplot.ps");
   String scriptfile=prefix;
   scriptfile.Add("_gender_autoplot.R");
   fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --autoQC, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "data <- read.table(file=\"%s\", header=T)\n", (const char*)genderplotdata);
   fprintf(fp, "postscript(\"%s\", paper=\"letter\", horizontal=T)\n", (const char*)genderplot);
   fprintf(fp, "isFemale<-data$SEX==2\n");
   fprintf(fp, "isMale<-data$SEX==1\n");
   fprintf(fp, "isUnknown<-data$SEX==0\n");
   fprintf(fp, "cutoff<-max(data$YCount)/2\n");
   fprintf(fp, "plot(data$YCount[isFemale],data$xHeterozygosity[isFemale], type=\"p\",\n");
   fprintf(fp, "  col=\"red\", xlim=c(0,max(data$YCount)), ylim=c(0,max(data$xHeterozygosity)),\n");
   if(gendererrorCount){
      fprintf(fp, "  main=\"Gender Checking in %s Samples (%d Samples Mislabeled)\", \n",
         (const char*)prefix, gendererrorCount);
      fprintf(fp, "  xlab=\"# Y-Chr SNPs\", ylab=\"X-Chr Heterozygosity\")\n");
   }else{
      fprintf(fp, "  main=\"Gender Checking in %s Samples\", \n", (const char*)prefix);
      fprintf(fp, "  xlab=\"# Y-Chr SNPs\", ylab=\"X-Chr Heterozygosity\")\n");
   }
   fprintf(fp, "points(data$YCount[isMale], data$xHeterozygosity[isMale], col=\"blue\")\n");
   fprintf(fp, "points(data$YCount[isUnknown], data$xHeterozygosity[isUnknown], col=\"black\")\n");
   fprintf(fp, "points(data$YCount[isFemale&data$YCount>cutoff], data$xHeterozygosity[isFemale&data$YCount>cutoff], col=\"red\")\n");
   fprintf(fp, "points(data$YCount[isMale&data$YCount<cutoff], data$xHeterozygosity[isMale&data$YCount<cutoff], col=\"blue\")\n");
   fprintf(fp, "abline(v=cutoff, col=\"purple\", lty=2)\n");
   fprintf(fp, "abline(v=cutoff*2/3, col=\"purple\")\n");
   fprintf(fp, "abline(v=cutoff*4/3, col=\"purple\")\n");
   fprintf(fp, "abline(h=%.4lf, col=\"purple\")\n", xHeterozygosity);
   fprintf(fp, "legend(\"topright\", c(\"Female\", \"Male\", \"Unknown\"), col=c(\"red\", \"blue\", \"black\"),\n");
   fprintf(fp, "  text.col = c(\"red\", \"blue\", \"black\"), pch = 19, cex = 1.2)\n");
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("--autoQC is done but R code %s failed.\n\n", (const char*)scriptfile);
   else if(error){
      printf("--autoQC is done but R code %s failed.\n", (const char*)scriptfile);
      printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
   }else{
      sprintf(command, "ps2pdf %s", (const char*)genderplot);
      system(command);
      printf("Gender plot are generated in %s_gender_autoplot.pdf\n", (const char*)prefix);
   }
}

void plotPopStructure(const char *prefix, int projectFlag)
{
   String scriptfile=prefix;
   scriptfile.Add("_pcplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --mds or --pca, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_pcplot.ps\", paper=\"letter\", horizontal=T)\n",(const char*)prefix);
   fprintf(fp, "data <- read.table(file=\"%spc.txt\", header=T)\n", (const char*)prefix);
   fprintf(fp, "plot(data$PC1, data$PC2, type=\"p\", xlab=\"PC1\", ylab=\"PC2\", main = \"Population Structure in %s\")\n",(const char*)prefix);
   if(projectFlag){
      fprintf(fp, "points(data[data[,6]==2,7], data[data[,6]==2,8], col = \"red\")\n");
      fprintf(fp, "legend(\"topright\", c(\"Reference Sample to Generate PCs\", \"Study Sample Projected to Reference PC Space\"),\n");
      fprintf(fp, "col=c(\"black\",\"red\"),text.col = c(\"black\", \"red\"), pch = 19, cex = 1.2)\n");
   }else
      fprintf(fp, "plot(data$PC3, data$PC4, type=\"p\", xlab=\"PC3\", ylab=\"PC4\", main = \"Population Structure in %s\")\n",(const char*)prefix);
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("--mds or --pca is done but R code %s failed.\n\n", (const char*)scriptfile);
   else if(error){
      printf("--ibdmap or --pca is done but R code %s failed.\n", (const char*)scriptfile);
      printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
   }else{
      sprintf(command, "ps2pdf %s_pcplot.ps", (const char*)prefix);
      system(command);
      printf("Population structure plot is generated in %s_pcplot.pdf\n", (const char*)prefix);
   }
}


int CheckRout(const char *scriptfile)
{
   String outfile=scriptfile;
   outfile.Add("out");
   int errorFlag = 1;
   char buff[256];   // define the buffer and allocate the length
   FILE *fp = fopen((const char*)outfile, "rb");
   if(fp != NULL){
      fseek(fp, -17, SEEK_END);// set pointer to the end of file minus the length you need. Presumably there can be more than one new line caracter
      fread(buff, 16, 1, fp); // read the contents of the file starting from where fseek() positioned us
      buff[16] = '\0';        // close the string
      String lastline=buff;
      if(lastline.Find("Execution halted")==-1) errorFlag = 0;
      else{
         fseek(fp, -100, SEEK_END);  // set pointer to the end of file minus the length you need. Presumably there can be more than one new line caracter
            fread(buff, 80, 1, fp); // read the contents of the file starting from where fseek() positioned us
         buff[80] = '\0';           // close the string
         lastline=buff;
         if(lastline.Find("Error in library")>-1) errorFlag = 3;
      }
      fclose(fp);                   // close the file
   }else errorFlag = 2;             // unable to open .Rout file
   return errorFlag;
}

void plotIBD1vsIBD2(const char *prefix, FILE *fp)
{
   fprintf(fp, "if(length(data$IBD1Seg)>0 & length(data$IBD2Seg)>0){\n");
   fprintf(fp, "d0 <- data$IBD2Seg>0.7\n");
   fprintf(fp, "d1.PO <- (!d0) & data$IBD1Seg+data$IBD2Seg>0.96 | (data$IBD1Seg+data$IBD2Seg>0.9 & data$IBD2Seg<0.08)\n");
   fprintf(fp, "d1.FS <- (!d0) & (!d1.PO) & data$PropIBD>0.35355 & data$IBD2Seg>=0.08\n");
   fprintf(fp, "d2 <- data$PropIBD>0.17678 & data$IBD1Seg+data$IBD2Seg<=0.9 & (!d1.FS)\n");
   fprintf(fp, "d3 <- data$PropIBD>0.08839 & data$PropIBD<=0.17678\n");
   fprintf(fp, "d4 <- data$PropIBD>0.04419 & data$PropIBD<=0.08839\n");
   fprintf(fp, "dU <- data$PropIBD>0 & data$PropIBD<=0.04419\n");
   fprintf(fp, "plot(data$IBD1Seg[dU], data$IBD2Seg[dU], type=\"p\", col = \"black\", cex.lab=1.2,\n");
   fprintf(fp, "xlim=c(min(data$IBD1Seg), max(data$IBD1Seg)),\n");
   fprintf(fp, "ylim=c(min(data$IBD2Seg), max(data$IBD2Seg)),\n");
   fprintf(fp, "main = \"IBD Segments In Inferred %s Relatives\",\n", (const char*)prefix);
   fprintf(fp, "xlab=expression(paste(\"Length Proportion of IBD1 Segments (\", pi[1], \")\", sep=\"\")),\n");
   fprintf(fp, "ylab=expression(paste(\"Length Proportion of IBD2 Segments (\", pi[2], \")\", sep=\"\")))\n");
   fprintf(fp, "points(data$IBD1Seg[d0], data$IBD2Seg[d0], col=\"purple\")\n");
   fprintf(fp, "points(data$IBD1Seg[d1.PO], data$IBD2Seg[d1.PO], col=\"red\")\n");
   fprintf(fp, "points(data$IBD1Seg[d1.FS], data$IBD2Seg[d1.FS], col=\"green\")\n");
   fprintf(fp, "points(data$IBD1Seg[d2], data$IBD2Seg[d2], col=\"blue\")\n");
   fprintf(fp, "points(data$IBD1Seg[d3], data$IBD2Seg[d3], col=\"magenta\")\n");
   fprintf(fp, "points(data$IBD1Seg[d4], data$IBD2Seg[d4], col=\"gold\")\n");
   fprintf(fp, "abline(h = 0.08, col = \"green\", lty = 3, lwd = 2)\n");   // FS && MZ
   fprintf(fp, "abline(a = 0.96, b = -1, col = \"red\", lty = 3, lwd = 2)\n");   // PO
   fprintf(fp, "abline(a = 0.3535534, b = -0.5, col = \"green\", lty = 3, lwd = 2)\n");// FS
   fprintf(fp, "abline(a = 0.1767767, b = -0.5, col = \"blue\", lty = 3, lwd = 2)\n");// 2nd/3rd
   fprintf(fp, "abline(a = 0.08838835, b = -0.5, col = \"magenta\", lty = 3, lwd = 2)\n");// 3rd/4th
   fprintf(fp, "abline(a = 0.04419, b = -0.5, col = \"gold\", lty = 3, lwd = 2)\n");// 4th/UN
   fprintf(fp, "allcolors <- c(\"purple\", \"red\", \"green\", \"blue\", \"magenta\", \"gold\", \"black\")\n");
   fprintf(fp, "legend(\"topright\", c(\"Inferred MZ\", \"Inferred PO\", \"Inferred FS\", \"Inferred 2nd\", \"Inferred 3rd\", \"Inferred 4th\", \"Inferred UN\"),\n");
   fprintf(fp, "  col=allcolors, text.col = allcolors, pch = 19, cex = 1.2)\n");
   fprintf(fp, "}\n");
}

/*
void plotHetConcvsIBD2(const char *prefix)
{
   String scriptfile=prefix;
   scriptfile.Add("_ibd2plot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "#%s for KING, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_ibd2plot.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix);
   fprintf(fp, "data <- read.table(file=\"%s.seg2\", header=T)\n", (const char*)prefix);
   fprintf(fp, "valid <- data$Pr_IBD2>0.005\n");
   fprintf(fp, "if(sum(valid)>0){\n");
   fprintf(fp, "plot(data$Pr_IBD2[valid], data$HetConc[valid], type=\"p\",\n");
   fprintf(fp, "col = \"black\", cex.lab=1.3,\n");
   fprintf(fp, "main = \"Relationships In %s Families\",\n", (const char*)prefix);
   fprintf(fp, "xlab=\"Proportion of IBD2 Segments\", ylab = \"Heterozygote Concordance Rate\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("R code %s cannot be run.\n", (const char*)scriptfile);
   else if(error)
      printf("Errors found in R code %s. Please check %sout for details.\n\n", (const char*)scriptfile, (const char*)scriptfile);
   else{
      sprintf(command, "ps2pdf %s_ibd2plot.ps", (const char*)prefix);
      system(command);
      printf("  Relationship plot is generated in %s_ibd2plot.pdf\n", (const char*)prefix);
   }
}

void plotIBD2(const char *prefix)
{
   String scriptfile=prefix;
   scriptfile.Add("_ibd2plot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_ibd2plot.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix);
   fprintf(fp, "data <- read.table(file=\"%s.ibs\", header=T)\n", (const char*)prefix);
   fprintf(fp, "if(dim(data)[1]>0){\n");
   fprintf(fp, "d0 <- data$Phi==0.5\n");
   fprintf(fp, "d1.PO <- data$Phi==0.25 & data$Z0==0\n");
   fprintf(fp, "d1.FS <- data$Phi==0.25 & data$Z0>0\n");
   fprintf(fp, "d2 <- data$Phi==0.125\n");
   fprintf(fp, "dU <- data$Phi==0\n");
   fprintf(fp, "dO <- !d0 & !d1.PO & !d1.FS & !d2 & !dU\n");
   fprintf(fp, "plot(data$Pr_IBD2[dU], data$Kinship[dU], type=\"p\",\n");
   fprintf(fp, "col = \"black\", cex.lab=1.3,\n");
   fprintf(fp, "xlim = c(min(data$Pr_IBD2), max(data$Pr_IBD2)),\n");
   fprintf(fp, "ylim = c(min(data$Kinship), max(data$Kinship)),\n");
   fprintf(fp, "main = \"Relationships In %s Families\",\n", (const char*)prefix);
   fprintf(fp, "xlab=\"Proportion of IBD2 Segments\", ylab = \"Estimated Kinship Coefficient\")\n");
   fprintf(fp, "points(data$Pr_IBD2[d0], data$Kinship[d0], col=\"purple\")\n");
   fprintf(fp, "points(data$Pr_IBD2[d1.PO], data$Kinship[d1.PO], col=\"red\")\n");
   fprintf(fp, "points(data$Pr_IBD2[d1.FS], data$Kinship[d1.FS], col=\"green\")\n");
   fprintf(fp, "points(data$Pr_IBD2[d2], data$Kinship[d2], col=\"blue\")\n");
   fprintf(fp, "points(data$Pr_IBD2[dO], data$Kinship[dO], col=\"gold\")\n");
   fprintf(fp, "points(data$Pr_IBD2[dU & data$Kinship>0.088], data$Kinship[dU & data$Kinship>0.088], col=\"black\")\n");
   fprintf(fp, "points(data$Pr_IBD2[d1.FS & data$Pr_IBD2<0.05], data$Kinship[d1.FS & data$Pr_IBD2<0.05], col=\"green\")\n");
   fprintf(fp, "abline(h = 0.3536, col = \"red\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.1768, col = \"green\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.0884, col = \"blue\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.0442, col = \"black\", lty = 3)\n");
   fprintf(fp, "legend(\"bottomright\", c(\"MZ Twin\", \"Parent-Offspring\", \"Full Siblings\", \"2nd-Degree\", \"More Distant\", \"Unrelated\"),\n");
   fprintf(fp, "col=c(\"purple\", \"red\", \"green\", \"blue\", \"gold\", \"black\"),\n");
   fprintf(fp, "text.col = c(\"purple\", \"red\", \"green\", \"blue\", \"gold\", \"black\"), pch = 19, cex = 1.2)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("R code %s cannot be run.\n", (const char*)scriptfile);
   else if(error)
      printf("Errors found in R code %s. Please check %sout for details.\n\n", (const char*)scriptfile, (const char*)scriptfile);
   else{
      sprintf(command, "ps2pdf %s_ibd2plot.ps", (const char*)prefix);
      system(command);
      printf("  Relationship plot is generated in %s_ibd2plot.pdf\n", (const char*)prefix);
   }
}
*/

/*
   fprintf(fp, "data <- read.table(file=\"%s\", nrows=2000, header=TRUE, stringsAsFactors=FALSE)[,c(\"ID1\", \"ID2\", \"InfType\")]\n", (const char*)inputfile);
   fprintf(fp, "postscript(\"%s_relativeplot.ps\", paper=\"letter\", horizontal=T)\n", prefix);
   if(degree == 1){
      fprintf(fp, "Inf.color <- c(\"purple\", \"red\", \"green\")\n");
      fprintf(fp, "Inf.type <- c(\"DUP/MZ\", \"PO\", \"FS\")\n");
   }else if(degree == 2){
      fprintf(fp, "Inf.color <- c(\"purple\", \"red\", \"green\", \"blue\")\n");
      fprintf(fp, "Inf.type <- c(\"DUP/MZ\", \"PO\", \"FS\", \"2nd\")\n");
   }else{
      fprintf(fp, "Inf.color <- c(\"purple\", \"red\", \"green\", \"blue\", \"yellow\")\n");
      fprintf(fp, "Inf.type <- c(\"DUP/MZ\", \"PO\", \"FS\", \"2nd\", \"3rd\")\n");
   }
   fprintf(fp, "relatives <- data[data$InfType %%in%% Inf.type,]\n");
   fprintf(fp, "count.rp <- min(dim(relatives)[1],1000)\n");
   fprintf(fp, "relatives <- relatives[1:count.rp,]\n");
   fprintf(fp, "for(i in 1:length(Inf.type)) relatives[relatives$InfType==Inf.type[i],\"InfType\"] <- i\n");
   fprintf(fp, "g <- graph_from_data_frame(d=relatives, vertices=unique(c(relatives$ID1,relatives$ID2)), directed=FALSE)\n");
   fprintf(fp, "plot(g, vertex.color=NA, vertex.size=1, vertex.label=NA, layout=layout_with_fr, asp=0,\n");
   fprintf(fp, "edge.color=Inf.color[as.numeric(relatives$InfType)], main=paste(count.rp, \"Relative Pairs in %s\"))\n", prefix);
   fprintf(fp, "legend(\"bottomright\", Inf.type, lty=1, col=Inf.color, text.col=Inf.color, pt.cex=2, cex=0.8, bty=\"n\")\n");
*/
/*
   fprintf(fp, " }else if(v[1]==4){\n");
   fprintf(fp, "  if(v[2]==5 && v[4]==4 && v[5]==1) buildType <- \"2_Parents+2_FullSiblings\"\n");
   fprintf(fp, "  else if(v[2]==6 && v[4]==3 && v[5]==3) buildType <- \"1_Parent+3_FullSiblings\"\n");
   fprintf(fp, "  else if(v[2]==6 && v[5]==6) buildType <- \"4_FullSiblings\"\n");
   fprintf(fp, " }else if(v[1]==5){\n");
   fprintf(fp, "  if(v[2]==9 && v[4]==6 && v[5]==3) buildType <- \"2_Parents+3_FullSiblings\"\n");
   fprintf(fp, "  else if(v[2]==10 && v[4]==4 && v[5]==6) buildType <- \"1_Parents+4_FullSiblings\"\n");
   fprintf(fp, "  else if(v[2]==10 && v[5]==10) buildType <- \"5_FullSiblings\"\n");
   fprintf(fp, " }else if(v[1]==6){\n");
   fprintf(fp, "  if(v[2]==14 && v[4]==8 && v[5]==6) buildType <- \"2_Parents+4_FullSiblings\"\n");
   fprintf(fp, "  else if(v[2]==15 && v[4]==5 && v[5]==10) buildType <- \"1_Parent+5_FullSiblings\"\n");
   fprintf(fp, "  else if(v[2]==15 && v[5]==15) buildType <- \"6_FullSiblings\"\n");
   fprintf(fp, " }else if(v[1]==7){\n");
   fprintf(fp, "  if(v[2]==20 && v[4]==10 && v[5]==10) buildType <- \"2_Parents+5_FullSiblings\"\n");
   fprintf(fp, "  else if(v[2]==21 && v[4]==6 && v[5]==15) buildType <- \"1_Parent+6_FullSiblings\"\n");
   fprintf(fp, "  else if(v[2]==21 && v[5]==21) buildType <- \"7_FullSiblings\"\n");
   fprintf(fp, " }\n");
*/

