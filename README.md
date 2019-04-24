TOPMed Variant Calling Pipeline (latest: Freeze 8)
==================================================

Overview of this repository
----------------------------

This repository is intended to provide a copy of software tools used for producing TOPMed Year 1 Freeze 5 variant calls and genotypes with a comprehensive documentation that allows investigators to understand the methods and reproduce the variant calls from the same set of aligned sequence reads.

This repository reflects specific versions of software tools that are
under active development in the Center for Statistical Genetics
(CSG). Most of the latest version of these software tools can be
accessed through multiple repositories linked as submodules, this repository is focused on a freeze of software tools that can reproduce a variant calls compatible to the latest TOPMed Freeze (Freeze 6 currently).


Outline of the variant calling procedure
----------------------------------------

Our ``GotCloud vt`` pipeline detects and genotype variants from a list of aligned sequence reads. Specifically, the pipeline consist of the following six key steps. Most of these procedure will be integrated into ``GotCloud`` software package later this year. 

Our ``GotCloud vt`` pipeline detects and genotype variants from a list of aligned sequence reads. Specifically, the pipeline consist of the following six key steps. Most of these procedure will be integrated into ``GotCloud`` software package later this year. 

1. **Sample quality control** : For each sequenced genome (in BAM/CRAMs), the genetic ancestry, sequence contamination, and biological sex are inferred using ``cramore cram-verify-bam`` and ``cramore vcf-normalize-depth``. 
2. **Variant detection** : For each sequenced genome (in BAM/CRAMs), candidate variants are detected by ``vt discover2`` software tools, separated by each chromosome. The candidate variants are normalized by ``vt normalize`` algorithm. 
3. **Variant consolidation** : For each chromosome, the called variant sites are merged across the genomes, accounting for overlap of variants between genomes, using ``cramore vcf-merge-candidate-variants``, ``vt annotate_indels`` software tool.
4. **Genotype and feature collection** : For each 100kb chunk of genome, the genotyping module implemented in ``cramore dense-genotypes`` collects individual genotypes and variant features across the merged sites by iterating each sequence genome focusing on the selected region.  
5. **Variant filtering** : We use the inferred pedigree of related and duplicated samples to calculate the Mendlian consistency statistics using ``king``, ``vcf-infer-ped``, ``vt milk-filter``, and train variant classifier using Support Vector Machine (SVM) implemented in the ``libsvm`` software package and ``run-svm-filter`` software tool.


![TOPMed Variant Calling Overview](topmed_variant_calling_overview.png)


Installation
-------------
First, clone the repository by recursively cloning each submodule.
```
   $ git clone --recurse-submodules -j8 https://github.com/statgen/topmed_variant_calling
   (Use --recursive instead of --recurse-submodules for git version 2.12 or lower) 
```

Next, build each submodule using the following sets of commands 
```
   $ cd libsvm/; make; cd ..
   $ cd apigenome; autoreconf -vfi; ./configure --prefix $PWD; make; make install; cd ..
   $ cd cramore; autoreconf -vcf; ./configure; make; cd ..
   $ cd libStatGen; make; cd ..
   $ cd bamUtil; make; cd ..
   $ cd invNorm; make; cd ..
   $ cd htslib; autoheader; autoconf; ./configure; make; cd ..
   $ cd samtools; authoheader; autoconf -Wno-syntax; ./configure; make; cd ..
   $ cd bcftools; make; cd ..
   $ cd king; g++ -O3 -c *.cpp; g++ -O3 -o king *.o -lz; cd ..
```

Performing variant calling with example data
----------------------------------------------
To produce variant calls using this pipeline, the following input
files are neded:

 1. Aligned sequenced reads in BAM or CRAM format. Each BAM and CRAM file should contain one sample per subject. It also must be indexed using ``samtools index`` or equivalent software tools.
 2. A sequence index file. Each line should contain [Sample ID] [Full Path to the BAM/CRAM file] [Contamination Estimates -- put zero if unknown]. See ``data/trio_data.index`` for example.
 3. Genomic resource files, such as FASTA, dbSNP, HapMap files.
 
Here we use 107 public samples from the TOPMed project to document a
reproducible variant calling pipeline that resembles the latest
TOPMed variant calling pipeline. In order to do so, you need do
download the following two sets of files.

1. Download 107 CRAMs from the public GCS bucket. (The location will be
   posted here asap, and contact hmkang@umich.edu in the meantime.) The
   ``examples/index/list.107.local.crams.index`` file contains the
   sample ID and CRAM file path. The CRAM file and the corresponding
   indices (.cram.crai) must present before running the examples.
2. Download the resource files for the variant calling. The tarball
   package will be available at
   ftp://share.sph.umich.edu/1000genomes/fullProject/hg38_resources
   . The ``resources/`` directory must be under ``examples/``directory
   to conform the examples below. 
   
The TOPMed variant calling was performed on the Google Cloud. The
software tool `cloudify` supports running jobs on the Google Cloud,
but these examples are configured to run in a local
computer. Cloud-based versions will also be available later.

### Step 1. Sample QC and Variant Detection 

First, change your current directory to
``topmed_variant_calling/examples``, and run the following command.

```
 $ ../apigenome/bin/cloudify --cmd ../scripts/run-discovery-local.cmd 
```

Sample QC and Variant Detection
-------------------------------
Variant detection from each sequence (ang aligned) genome is performed by ``vt discover2`` software tool. The script ``step-1-detect-variants.pl`` provide a mean to automate the variant detection across a large number of sequence genome.

The variant detection algorithm consider a variant as a potential candidate variant if there exists a mismatch between the aligned sequence reads and the reference genome. Because such a mismatch can easily occur by random errors, only potential candidate variants passing the following criteria are considered to be ***candidate variants*** in the next steps.

1. At least two identical evidence of variants must be observed from aligned sequence reads. 
  1. Each individual evidence will be normalized using the normalization algorithm implemented in ``vt normalize`` software tools.
  1. Only evidence on the reads with mapping quality 20 or greater will be considered.
  1. Duplicate reads, QC-passed reads, supplementary reads, secondary reads will be ignored. 
  1. Evidence of variant within overlapping fragments of read pairs will not be double counted. Either end of the overlapping read pair will be soft-clipped using ``bam clipOverlap`` software tool.  
1. Assuming per-sample heterozygosity of 0.1%, the posterior probability of having variant at the position should be greater than 50%. The method is equivalent to the `glfSingle` model described in http://www.ncbi.nlm.nih.gov/pubmed/25884587

The variant detection step is required only once per sequenced genome, when multiple freezes of variant calls are produced over the course of time.

 
Variant Consolidation
---------------------
Variants detected from the discovery step will be merged across all samples. This step is implemented in the ``step-2-detect-variants.pl`` scripts.

1. Each non-reference allele normalized by ``vt normalize`` algorithm is merged across the samples, and unique alleles are printed as biallelic candidate variants. The algorithm is published at http://www.ncbi.nlm.nih.gov/pubmed/25701572
2. If there are alleles overlapping with other SNPs and Indels, ``overlap_snp`` and ``overlap_indel`` filters are added in the ``FILTER`` column of the corresponding variant.
3. If there are tandem repeats with 2 or more repeats with total repeat length of 6bp or longer, the variant is annotated as potential VNTR (Variant Number Tandem Repeat), and ``overlap_vntr`` filters are added to the variant overlapping with the repeat track of the putative VNTR.     


Variant Genotyping and Feature Collection
-----------------------------------------
The genotyping step iterate all the merged variant site across the sample. It iterates each BAM/CRAM files one at a time sequentially for each 1Mb chunk to perform contamination-adjusted genotyping and annotation of variant features for filtering. The following variant features are calculated during the genotyping procedure. 

 * AVGDP : Average read depth per sample
 * AC : Non-reference allele count
 * AN : Total number of alleles
 * GC : Genotype count
 * GN : Total genotype counts
 * HWEAF : Allele frequency estimated from PL under HWE
 * HWDAF : Genotype frequency estimated from PL under HWD
 * IBC : [ Obs(Het) – Exp(Het) ] / Exp[Het]
 * HWE_SLP : -log(HWE likelihood-ratio test p-value) ⨉ sign(IBC)
 * ABE : Average fraction [#Ref Allele] across all heterozygotes
 * ABZ : Z-score for tesing deviation of ABE from expected value (0.5)
 * BQZ: Z-score testing association between allele and base qualities
 * CYZ: Z-score testing association between allele and the sequencing cycle
 * STZ : Z-score testing association between allele and strand
 * NMZ : Z-score testing association between allele and per-read mismatches
 * IOR : log [ Obs(non-ref, non-alt alleles) / Exp(non-ref, non-alt alleles) ]
 * NM1 : Average per-read mismatches for non-reference alleles
 * NM0 : Average per-read mismatches for reference alleles

The genotyping was done by adjusting for potential contamination. It uses adjusted genotype likelihood similar to the published method https://github.com/hyunminkang/cleancall, but does not use estimated population allele frequency for the sake of computational efficiency. It conservatively assumes that probability of observing non-reference read given homozygous reference genotype is equal to the half of the estimated contamination level, (or 1%, whichever is greater). The probability of observing reference reads given homozygous non-reference genotype is calculated in a similar way. This adjustment makes the heterozygous call more conservatively when the reference and non-reference allele is strongly imbalanced. For example, if 45 reference alleles and 5 non-reference alleles are observed at Q40, the new method calls it as homozygous reference genotype while the original method ignoring potential contamination calls it as heterozygous genotype. This adjustment improves the genotype quality of contaminated samples by reducing genotype errors by several folds.

Variant Filtering
-----------------
The variant filtering in TOPMed Freeze 6 were performed by (1) first calculating Mendelian consistency scores using known familial relatedness and duplicates, and (2) training SVM classifier between the known variant site (positive labels) and the Mendelian inconsistent variants (negative labels). 

The negative labels are defined if the Bayes Factor for Mendelian consistency quantified as ``Pr(Reads | HWE, Pedigree) / Pr(Reads | HWD, no Pedigree )`` less than 0.001. Also variant is marked as negative labels if 3 or more samples show 20% of non-reference Mendelian discordance within families or genotype discordance between duplicated samples.

The positive labels are the SNPs found polymorphic either in the 1000 Genomes Omni2.5 array or in HapMap 3.3, with additional evidence of being polymorphic from the sequenced samples. Variants eligible to be marked both positive and negative labels are discarded from the labels. The SVM scores trained and predicgted by ``libSVM`` software tool will be annotated in the VCF file. 

Two additional hard filtering was applied additionally. First is excessive heterozygosity filter ``(EXHET)``, if the Hardy-Weinberg disequilbrium p-value was less than 1e-6 in the direction of excessive heterozygosity. Another filter is Mendelian discordance filter ``(DISC)``, with 3 or more Mendelian discordance or duplicate discordance observed from the samples.

Questions
---------
For further questions, pleast contact Hyun Min Kang (hmkang@umich.edu) and Jonathon LeFaive (lefaivej@umich.edu).
