#ifndef __rplot_h__
#define __rplot_h__

#include "IntArray.h"
#include "MathVector.h"

void plotMIerror(const char *prefix);
void plotUniqueFamily(const char *prefix, int degree, const char *analysis);
void plotDuplicate(const char *prefix);
void plotBuild(const char *prefix);
void plotSplitped(const char *prefix);
void plotCluster(const char *prefix);
void plotGenderError(const char *prefix, IntArray & plotx, Vector & ploty, IntArray & plotz, double xHeterozygosity, int gendererrorCount);
void plotRelationship(const char *prefix);
void plotIBDSeg(const char *prefix);
void plotPopStructure(const char *prefix, int projectFlag);

// not released yet
void plotAUCmapping(const char *prefix, int SEXCHR);
void plotNPL(const char *prefix, int SEXCHR);
void plotHEreg(const char *prefix, int SEXCHR);
void plotIBDmapping(const char *prefix, int SEXCHR);
void plotROHmapping(const char *prefix, const char *stratName, int SEXCHR);
void plotROHforQT(const char *prefix, int SEXCHR);
void plotPopROH(const char *prefix, int SEXCHR);
void plotPopDist(const char *prefix);
void plotAncestry(const char *prefix);

#endif
