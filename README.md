# XYalign

   ```
    Copyright 2016 by Melissa A. Wilson Sayres

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    ```



## Inferring sex chromosome ploidy in NGS data
###XYalign: Hacking sex chromosome variation

Madeline Couse, Bruno Grande, Eric Karlins, Tanya Phung, Phillip Richmond, Timothy H. Webster, Whitney Whitford, Melissa A. Wilson Sayres

Sex chromosome aneuploidies are currently estimated to be as common as 1/400 in humans. Atypical ploidy will affect variant calling and measures of genomic variation that are central to most clinical genomic studies. Further, the high degree of similarity between gametologous sequences on the X and Y chromosomes can lead to the misalignment of sequencing reads and substantially affect variant calling. Here we present XYalign, a new tool that (1) quickly infers sex chromosome ploidy in NGS data (DNA and RNA), (2) remaps reads based on the inferred sex chromosome complement of the individual, and (3) outputs quality, depth, and allele-balance metrics across the sex chromosomes.

October 17, 2016 slide show here: https://docs.google.com/presentation/d/1OB2d_mu5zC742N_NKfzHjVpUm4BFtm5lUzniLLI--OQ/edit?usp=sharing


### List of Goals: Assess X/Y ploidy and correct for misalignment
1. Extract input chromosomes - recommend chrX, chrY, autosome (e.g., chr19) - from BAM

2. Infer sex chromosome ploidy from WGS data relative to autosomal ploidy
  + E.g., XX, XY, XXY, X0
  + And all other combinations

  Use
  + Quality
  + Read Depth
  + Allele Balance
  + Ampliconic/Palindromic/CNV filter

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
  + Strip reads from X and Y
  + Remap all X & Y reads to the X chromosome only
  + Remove X and Y from the input BAM file
  + Merge the empty Y and the remapped X chromosome into the BAM

4. IF - If we infer there are Y chromosomes in the sample, conduct re-mapping to increase confidence in X-linked alleles.
  + Strip reads from X and Y
  + Remap all X & Y reads to the X chromosome, and the accompanying Y chromosmoe with PARs masked out.
  + Remove X and Y from the input BAM file
  + Merge the empty Y and the remapped X chromosome into the BAM

5. Assessment:
  + Variant calling in 1000 genomes high coverage data before/after XYalign
  + Test how different alignment algorithms, parameters, and reference sequences affect X & Y variant calling
  + Compare variant calling with the "Gold Standard" reference individual
  + Assess XYalign on simulated sex chromosome aneuploidy data

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
