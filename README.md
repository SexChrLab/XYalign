# XYalign

# Inferring sex chromosome ploidy in NGS data
Slide show here: https://docs.google.com/presentation/d/1OB2d_mu5zC742N_NKfzHjVpUm4BFtm5lUzniLLI--OQ/edit?usp=sharing

### Publication  Links
- chromosome X: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2665286/
- chromosome Y: http://www.nature.com/nature/journal/v423/n6942/full/nature01722.html
- Repetitive mapping solutions: https://www.ncbi.nlm.nih.gov/pubmed/22124482


### List of Goals: Assess X/Y ploidy and correct for misalignment
1. Extract input chromosomes - recommend chrX, chrY, chr19 - from BAM (can input any autosome)

2. Infer sex chromosome ploidy from WGS data relative to autosomal ploidy 
  + XX 
  + XY
  + XXY
  + X0
  + And all other combinations
  Use 
  A. Quality
  B. Read Depth
  C. Allele Balance
  D. Ampliconic/Palindromic/CNV filter
  
 Typical expectations for heterozygous calls under different sex chromosome complements: 
  
  Genotype | X_call | Y_call
  --- | --- |  --- 
  XX | het | none
  XY | hap | hap
  X0 | hap | none or partial_hap
  XXY | het or hap | hap 
  XYY | hap | hap
  XXX | het | none

  Note: Half of 47,XXY are paternal in origin -> do not expect het sites: http://humupd.oxfordjournals.org/content/9/4/309.full.pdf

  Expectations for depth under different sex chromosome complements: 

  Genotype | X_depth | Y_depth
  --- | --- |  --- 
  XX | 2x | 0x
  XY | 1x | 1x
  X0 | 1x | 0x (or partial)
  XXY | 2x | 1x 
  XYY | 1x | 2x
  XXX | 3x | 0x

  
3. IF - If we infer there are no Y chromosomes in the sample, conduct re-mapping to increase confidence in X-linked alleles.
  Strip reads from X and Y
  Remap all X & Y reads to the X chromosome only 
  Remove X and Y from the input BAM file
  Merge the empty Y and the remapped X chromosome into the BAM

4. Assessment of 1000 genomes high coverage data
  Compare SNV and CNV variant calling in 1000 genomes high coverage before/after running this pipeline
  Test how different alignment algorithms, parameters, and reference sequences affect variant calling in different regions of the X and Y
  Compare variant calling with the "Gold Standard" reference individual

Other goals: Because I think we have to address this if we want to get a really good handle on #2 given the extremely high copy number variable regions on X and Y - the ampliconic regions. Likely we will masking them out to infer #2, which will be easiest, but then we can have an extended goal to see characterize variations in these regions.



### Known problems and complications
- High sequence identity between X and Y
- Higher amount of sequence repeats




### Group Members
Name | email | github ID
--- | --- |  --- 
Madeline Couse| mhcouse@gmail.com | @Madelinehazel
Bruno Grande | bgrande@sfu.ca | @brunogrande
Eric Karlins | karlinser@mail.nih.gov | @ekarlins
Tanya Phung | tnphung@ucla.edu | @tnphung
Phillip Richmond | phillip.a.richmond@gmail.com | @Phillip-a-Richmond 
Tim Webster | timothy.h.webster@asu.edu | @thw17
Whitney Whitford | whitney.whitford@auckland.ac.nz | @whitneywhitford
Melissa A. Wilson Sayres | melissa.wilsonsayres@asu.edu | @mwilsonsayres

