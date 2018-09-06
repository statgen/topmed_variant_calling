//////////////////////////////////////////////////////////////////////
// rplot.cpp
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

void Engine::plotPopROH()
{
   String scriptfile=prefix;
   scriptfile.Add("_poprohplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
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
//   fprintf(fp, "abline(h = 2.2, lty = 3, col=\"green\")\n");
   fprintf(fp, "abline(h = -3.6, lty = 2, col=\"red\")\n");
//   fprintf(fp, "abline(h = -2.2, lty = 3, col=\"green\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "}\n");
   fprintf(fp, "dev.off()\n");
   fprintf(fp, "}\n");
   fclose(fp);

   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   sprintf(command, "ps2pdf %s_poprohplot.ps", (const char*)prefix);
   system(command);
   printf("  ROH difference between populations is generated in %s_poprohplot.pdf\n", (const char*)prefix);
}

void Engine::plotROHforQT()
{
   String scriptfile=prefix;
   scriptfile.Add("_mthomoplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
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
   sprintf(command, "ps2pdf %s_mthomoplot.ps", (const char*)prefix);
   system(command);
   printf("  Homozygosity mapping plots are generated in %s_mthomoplot.pdf\n", (const char*)prefix);
}

void Engine::plotROHmapping(const char *stratName)
{
   String postfix = stratName[0]=='\0'? "homomap": "homomapMH";
   String scriptfile=prefix;
   scriptfile.Add("_");
   scriptfile.Add(postfix);
   scriptfile.Add("plot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
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
   sprintf(command, "ps2pdf %s_%splot.ps",
      (const char*)prefix, (const char*)postfix);
   system(command);
   printf("  Homozygosity mapping plot is generated in %s_%splot.pdf\n",
      (const char*)prefix, (const char*)postfix);
}

void Engine::plotPopDist()
{
   String scriptfile=prefix;
   scriptfile.Add("_popdistplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
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
   sprintf(command, "ps2pdf %s_popdistplot.ps", (const char*)prefix);
   system(command);
   printf("  Population distance plot is generated in %s_popdistplot.pdf\n", (const char*)prefix);
}

void Engine::plotIBDmapping()
{
   String scriptfile=prefix;
   scriptfile.Add("_ibdmapplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
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
   sprintf(command, "ps2pdf %s_ibdmapplot.ps", (const char*)prefix);
   system(command);
   printf("  IBD mapping plot is generated in %s_ibdmapplot.pdf\n", (const char*)prefix);
}

void Engine::plotAncestry()
{
   String scriptfile=prefix;
   scriptfile.Add("_ancestryplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
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
   sprintf(command, "ps2pdf %s_ancestryplot.ps", (const char*)prefix);
   system(command);
   printf("  Ancestry plot is generated in %s_ancestryplot.pdf\n", (const char*)prefix);
}

void Engine::plotNPL()
{
   String scriptfile=prefix;
   scriptfile.Add("_nplplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_nplplot.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix);
   fprintf(fp, "data<-c()\n");
   fprintf(fp, "if(file.exists(\"%s.npl\")) data <- read.table(file=\"%s.npl\", header=T)\n",
      (const char*)prefix,(const char*)prefix);
   fprintf(fp, "if(length(data$LOD_DSP)>0){\n");
//   fprintf(fp, "relativePos <- (data$Start + data$Stop) / 2\n");
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
   fprintf(fp, "plot(Pos, data$LOD_DSP, type=\"l\", xlab=\"Chromosome\", ylab = \"LOD Score\", xaxt = 'n',\n");
   fprintf(fp, "  cex = 1.5, cex.axis=1.5, cex.main=1.5, cex.lab = 1.5, cex.sub = 1.5,\n");
   fprintf(fp, "  lwd = 1.5, main = \"ASP NPL Analysis in %s\", ylim=c(0,min(max(data$LOD_DSP),10)))\n", (const char*)prefix);

   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "maxLOD <- max(0,data$LOD_DSP[data$Chr==i])\n");
   fprintf(fp, "if(maxLOD >= 3.6){\n");
   fprintf(fp, "localpos <- median(data$Pos[data$Chr==i][data$LOD_DSP[data$Chr==i]==maxLOD])\n");
   fprintf(fp, "text(base+localpos, min(maxLOD,10)-0.1, paste(i, \":\", localpos, sep=\"\"), col=\"red\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "  base <- base + max(0, data$Pos[data$Chr==i])\n");
   fprintf(fp, "}\n");

   fprintf(fp, "axis(1, labels=FALSE, line = 0, at = c(0,baseStop))\n");
   fprintf(fp, "axis(1, labels=c(1:%d), line = -0.5, lty=\"blank\", at = baseMiddle)\n", SEXCHR-1);
   fprintf(fp, "abline(h = 3.6, lty = 2, col=\"red\")\n");
   fprintf(fp, "abline(h = 2.2, lty = 3, col=\"green\")\n");
   fprintf(fp, "}\n");
   fclose(fp);

   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   sprintf(command, "ps2pdf %s_nplplot.ps", (const char*)prefix);
   system(command);
   printf("  ASP NPL plot is generated in %s_nplplot.pdf\n", (const char*)prefix);
}

void Engine::plotAUCmapping()
{
   String scriptfile=prefix;
   scriptfile.Add("_aucmapplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_aucmapplot.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix);

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
   sprintf(command, "ps2pdf %s_aucmapplot.ps", (const char*)prefix);
   system(command);
   printf("  AUC mapping plot is generated in %s_aucmapplot.pdf\n", (const char*)prefix);
}

void Engine::plotRelationship()
{
   String scriptfile=prefix;
   scriptfile.Add("_relplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_relplot.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix);
//   fprintf(fp, "png(\"%s_relplot%03d.png\", paper=\"letter\", horizontal=T)\n", (const char*)prefix);

   fprintf(fp, "data<-c()\n");
   fprintf(fp, "if(file.exists(\"%s.kin\")) data <- read.table(file=\"%s.kin\", header=T)\n",
      (const char*)prefix,(const char*)prefix);
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
   fprintf(fp, "u <- par(\"usr\")\n");
   fprintf(fp, "legend(u[2], u[4], xjust=1, yjust=1,\n");
   fprintf(fp, "legend=c(\"MZ Twin\", \"Parent-Offspring\", \"Full Siblings\", \"2nd-Degree\", \"3rd-Degree\", \"More Distant\", \"Unrelated\"),\n");
   fprintf(fp, "col=c(\"purple\", \"red\", \"green\", \"blue\", \"magenta\", \"gold\", \"black\"),\n");
   fprintf(fp, "text.col = c(\"purple\", \"red\", \"green\", \"blue\", \"magenta\", \"gold\", \"black\"), pch = 19, cex = 1.2)\n");

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
   fprintf(fp, "u <- par(\"usr\")\n");
   fprintf(fp, "legend(u[2], u[3], xjust=1, yjust=0,\n");
   fprintf(fp, "legend=c(\"MZ Twin\", \"Parent-Offspring\", \"Full Siblings\", \"2nd-Degree\", \"3rd-Degree\", \"More Distant\", \"Unrelated\"),\n");
   fprintf(fp, "col=c(\"purple\", \"red\", \"green\", \"blue\", \"magenta\", \"gold\", \"black\"),\n");
   fprintf(fp, "text.col = c(\"purple\", \"red\", \"green\", \"blue\", \"magenta\", \"gold\", \"black\"), pch = 19, cex = 1.2)\n");

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
   fprintf(fp, "u <- par(\"usr\")\n");
   fprintf(fp, "legend(u[2], u[4], xjust=1, yjust=1,\n");
   fprintf(fp, "legend=c(\"MZ Twin\", \"Parent-Offspring\", \"Full Siblings\", \"2nd-Degree\", \"More Distant\", \"Unrelated\"),\n");
   fprintf(fp, "col=c(\"purple\", \"red\", \"green\", \"blue\", \"gold\", \"black\"),\n");
   fprintf(fp, "text.col = c(\"purple\", \"red\", \"green\", \"blue\", \"gold\", \"black\"), pch = 19, cex = 1.2)\n");
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
   fprintf(fp, "u <- par(\"usr\")\n");
   fprintf(fp, "legend(u[2], u[3], xjust=1, yjust=0,\n");
   fprintf(fp, "legend=c(\"Parent-Offspring\", \"Full Siblings\", \"2nd-Degree\", \"More Distant\", \"Unrelated\"),\n");
   fprintf(fp, "col=c(\"red\", \"green\", \"blue\", \"gold\", \"black\"),\n");
   fprintf(fp, "text.col = c(\"red\", \"green\", \"blue\", \"gold\", \"black\"), pch = 19, cex = 1.2)\n");
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
   plotIBD1vsIBD2(fp);

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
   fprintf(fp, "u <- par(\"usr\")\n");
   fprintf(fp, "legend(u[2], u[3], xjust=1, yjust=0,\n");
   fprintf(fp, "legend=c(\"Inferred MZ\", \"Inferred PO\", \"Inferred FS\", \"Inferred 2nd\", \"Inferred 3rd\", \"Inferred 4th\", \"Inferred UN\"),\n");
   fprintf(fp, "col=c(\"purple\", \"red\", \"green\", \"blue\", \"magenta\", \"gold\", \"black\"),\n");
   fprintf(fp, "text.col = c(\"purple\", \"red\", \"green\", \"blue\", \"magenta\", \"gold\", \"black\"), pch = 19, cex = 1.2)\n");
   fprintf(fp, "}\n");

   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   sprintf(command, "ps2pdf %s_relplot.ps", (const char*)prefix);
   system(command);
   printf("  Relationship plot is generated in %s_relplot.pdf\n", (const char*)prefix);
}

void Engine::plotIBDSeg()
{
   String scriptfile=prefix;
   scriptfile.Add("_ibd1vsibd2.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_ibd1vsibd2.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix);
   fprintf(fp, "data<-c()\n");
   fprintf(fp, "if(file.exists(\"%s.seg\")) data <- read.table(file=\"%s.seg\", header=T)\n",
      (const char*)prefix,(const char*)prefix);
   plotIBD1vsIBD2(fp);
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   sprintf(command, "ps2pdf %s_ibd1vsibd2.ps", (const char*)prefix);
   system(command);
   printf("  IBD1 vs IBD2 plot is generated in %s_ibd1vsibd2.pdf\n", (const char*)prefix);
}

void Engine::plotIBD1vsIBD2(FILE *fp)
{
   // Page 1 plot: IBD1 vs IBD2
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
   fprintf(fp, "u <- par(\"usr\")\n");
   fprintf(fp, "legend(u[2], u[4], xjust=1, yjust=1,\n");
   fprintf(fp, "legend=c(\"Inferred MZ\", \"Inferred PO\", \"Inferred FS\", \"Inferred 2nd\", \"Inferred 3rd\", \"Inferred 4th\", \"Inferred UN\"),\n");
   fprintf(fp, "col=c(\"purple\", \"red\", \"green\", \"blue\", \"magenta\", \"gold\", \"black\"),\n");
   fprintf(fp, "text.col = c(\"purple\", \"red\", \"green\", \"blue\", \"magenta\", \"gold\", \"black\"), pch = 19, cex = 1.2)\n");
   fprintf(fp, "}\n");
}

void Engine::plotIBD2()
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
   fprintf(fp, "u <- par(\"usr\")\n");
   fprintf(fp, "legend(u[2], u[3], xjust=1, yjust=0,\n");
   fprintf(fp, "legend=c(\"MZ Twin\", \"Parent-Offspring\", \"Full Siblings\", \"2nd-Degree\", \"More Distant\", \"Unrelated\"),\n");
   fprintf(fp, "col=c(\"purple\", \"red\", \"green\", \"blue\", \"gold\", \"black\"),\n");
   fprintf(fp, "text.col = c(\"purple\", \"red\", \"green\", \"blue\", \"gold\", \"black\"), pch = 19, cex = 1.2)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   sprintf(command, "ps2pdf %s_ibd2plot.ps", (const char*)prefix);
   system(command);
   printf("  Relationship plot is generated in %s_ibd2plot.pdf\n", (const char*)prefix);
}

void Engine::plotHetConcvsIBD2()
{
   String scriptfile=prefix;
   scriptfile.Add("_ibd2plot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
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
   sprintf(command, "ps2pdf %s_ibd2plot.ps", (const char*)prefix);
   system(command);
   printf("  Relationship plot is generated in %s_ibd2plot.pdf\n", (const char*)prefix);
}

void Engine::plotGenderError(IntArray & plotx, Vector & ploty, IntArray & plotz, double xHeterozygosity, int gendererrorCount)
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
   fprintf(fp, "data <- read.table(file=\"%s\", header=T)\n",
      (const char*)genderplotdata);
   fprintf(fp, "postscript(\"%s\", paper=\"letter\", horizontal=T)\n",
      (const char*)genderplot);
   fprintf(fp, "isFemale<-data$SEX==2\n");
   fprintf(fp, "isMale<-data$SEX==1\n");
   fprintf(fp, "isUnknown<-data$SEX==0\n");
   fprintf(fp, "cutoff<-max(data$YCount)/2\n");
   fprintf(fp, "plot(data$YCount[isFemale],data$xHeterozygosity[isFemale], type=\"p\",\n");
   fprintf(fp, "     col=\"red\", xlim=c(0,max(data$YCount)), ylim=c(0,max(data$xHeterozygosity)),\n");
   if(gendererrorCount){
      fprintf(fp, "     main=\"Gender Checking in %s Samples (%d Samples Mislabeled)\", \n",
         (const char*)prefix, gendererrorCount);
      fprintf(fp, "     xlab=\"# Y-Chr SNPs\", ylab=\"X-Chr Heterozygosity\")\n");
   }else{
      fprintf(fp, "     main=\"Gender Checking in %s Samples\", \n", (const char*)prefix);
      fprintf(fp, "     xlab=\"# Y-Chr SNPs\", ylab=\"X-Chr Heterozygosity\")\n");
   }
   fprintf(fp, "points(data$YCount[isMale], data$xHeterozygosity[isMale], col=\"blue\")\n");
   fprintf(fp, "points(data$YCount[isUnknown], data$xHeterozygosity[isUnknown], col=\"black\")\n");
   fprintf(fp, "points(data$YCount[isFemale&data$YCount>cutoff], data$xHeterozygosity[isFemale&data$YCount>cutoff], col=\"red\")\n");
   fprintf(fp, "points(data$YCount[isMale&data$YCount<cutoff], data$xHeterozygosity[isMale&data$YCount<cutoff], col=\"blue\")\n");
   fprintf(fp, "abline(v=cutoff, col=\"purple\", lty=2)\n");
   fprintf(fp, "abline(v=cutoff*2/3, col=\"purple\")\n");
   fprintf(fp, "abline(v=cutoff*4/3, col=\"purple\")\n");
   fprintf(fp, "abline(h=%.4lf, col=\"purple\")\n", xHeterozygosity);
   fprintf(fp, "u <- par(\"usr\")\n");
   fprintf(fp, "legend(u[2], u[4], xjust=1, yjust=1,\n");
   fprintf(fp, "    legend=c(\"Female\", \"Male\", \"Unknown\"), col=c(\"red\", \"blue\", \"black\"),\n");
   fprintf(fp, "    text.col = c(\"red\", \"blue\", \"black\"), pch = 19, cex = 1.2)\n");
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   sprintf(command, "ps2pdf %s", (const char*)genderplot);
   system(command);
   printf("  Gender plot are generated in %s_gender_autoplot.pdf\n", (const char*)prefix);
}


void Engine::plotPopStructure()
{
   String scriptfile=prefix;
   scriptfile.Add("_pcplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_pcplot.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix);
   fprintf(fp, "data <- read.table(file=\"%spc.ped\", header=F)\n", (const char*)prefix);
   fprintf(fp, "plot(data[,7], data[,8], type=\"p\", xlab=\"PC1\", ylab=\"PC2\", main = \"Population Structure in %s\")\n",
      (const char*)prefix);
   if(projectFlag){
      fprintf(fp, "points(data[data[,6]==2,7], data[data[,6]==2,8], col = \"red\")\n");
      fprintf(fp, "u <- par(\"usr\")\n");
      fprintf(fp, "legend(u[2], u[4], xjust=1, yjust=1,\n");
      fprintf(fp, "legend=c(\"Reference Sample Used for PCA\", \"Sample Projected to PC Space of Reference Samples\"),\n");
      fprintf(fp, "col=c(\"black\",\"red\"),text.col = c(\"black\", \"red\"), pch = 19, cex = 1.2)\n");
   }
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   sprintf(command, "ps2pdf %s_pcplot.ps", (const char*)prefix);
   system(command);
   printf("  Population structure plot is generated in %s_pcplot.pdf\n", (const char*)prefix);
}
