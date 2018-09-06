////////////////////////////////////////////////////////////////////// 
// merlin/MerlinSort.h 
// (c) 2000-2007 Goncalo Abecasis
// 
// This file is distributed as part of the MERLIN source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile MERLIN.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Tuesday December 18, 2007
// 
 
#ifndef __MERLINSORT_H__
#define __MERLINSORT_H__

#include "Pedigree.h"

// This routine sorts families so that densely genotyped individuals
// appear before those with more missing data, which empirically appears
// to reduce the average size of gene flow trees
//

void SortFamilies(Pedigree & ped);

#endif
  
