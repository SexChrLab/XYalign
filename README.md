<p align="center">
  <img src="https://github.com/WilsonSayresLab/XYalign/blob/master/Files/XYlogo.png" width="500"/>
</p>

# XYalign: Inferring Sex Chromosome Ploidy in NGS Data
Timothy H. Webster, Tanya Phung, Madeline Couse, Bruno Grande, Eric Karlins, Phillip Richmond, Whitney Whitford, Melissa A. Wilson Sayres

Sex chromosome aneuploidies are currently estimated to be as common as 1/400 in humans. Atypical ploidy will affect variant calling and measures of genomic variation that are central to most clinical genomic studies. Further, the high degree of similarity between gametologous sequences on the X and Y chromosomes can lead to the misalignment of sequencing reads and substantially affect variant calling. Here we present XYalign, a new tool that (1) quickly infers sex chromosome ploidy in NGS data (DNA and RNA), (2) remaps reads based on the inferred sex chromosome complement of the individual, and (3) outputs quality, depth, and allele-balance metrics across the sex chromosomes.

October 17, 2016 slide show of concept and inital results for Hackseq: https://docs.google.com/presentation/d/1OB2d_mu5zC742N_NKfzHjVpUm4BFtm5lUzniLLI--OQ/edit?usp=sharing

## Please note the XYalign is still in pre-release stages
Changes are made to the master branch quite frequently, so there still might be a number of known or unknown bugs at any given time.  If you'd like to try things before release (which will be early December), please feel free to contact us (either via an Issue or by email to Tim Webster - address at the bottom of this page) to inquire about the current status of the program.


## To-Do Items Before Release

Program
- [ ] Write commands, time, etc. to logfile
- [ ] Add support for arbitrary bwa flags
- [ ] Add support for additional mappers
- [ ] Handle duplicate reads
- [x] Write @PG line to final bam
- [ ] Accept and integrate external bed files for excluding and including regions
- [ ] Incorporate likelihood model for ploidy estimation
- [ ] Second round of bam window analysis and plotting after remapping and filtering low-quality regions

Testing
- [ ] Simulate and assemble XX, XY, XXY, and XO genomes 
  - [ ] Test XYalign on these genomes
- [ ] Download data for XX and XY:
  - [ ] high coverage whole genome
    - [ ] Test XYalign
  - [x] low coverage whole genome
    - [ ] Test XYalign
  - [x] exome 
    - [ ] Test XYalign
- [ ] In an example case, plot runtimes with increasing numbers of threads (e.g., 1, 2, 4, 8, 16)


## Quick Start
### What you need
Minimally, you'll need:

1. A BAM or CRAM file

2. The reference genome against which reads were mapped to create the BAM/CRAM file in (1)

3. A .fai index of the reference genome in (2) located in the same directory as the reference.  This can be generated using the command ```samtools faidx <reference_fasta>```.  See the [samtools documentation](http://www.htslib.org/doc/samtools.html) for more information 

4. An environment with a host of python packages (numpy, pandas, matplotlib, seaborn, pysam, and pybedtools) and external programs (platypus, bwa, samtools, and bbmap) installed.  Note that we only support Python 2.7 right now (because of Platypus). Probably the easiest way to do this is to download [miniconda](http://conda.pydata.org/miniconda.html) and let it append its path to your .bashrc file.  You should then be able to set up and environment with the following commands:
```
conda config --add channels r

conda config --add channels conda-forge

conda config --add channels bioconda

conda create -n xyalign_env python=2.7 pysam pybedtools numpy pandas matplotlib seaborn platypus-variant bwa bbmap samtools sambamba
```
This will work in within a Linux operating system.  As of right now, bioconda won't install platypus on Macs, so Mac users will have to replace the final command with:
```
conda create -n xyalign_env python=2.7 pysam pybedtools numpy pandas matplotlib seaborn bwa bbmap samtools sambamba
```
and [install platypus on their own](http://www.well.ox.ac.uk/platypus).

After setup, this environment can be loaded with the command (same on Linux and Macs):
```
source activate xyalign_env
```
See [Anaconda's documentation](http://conda.pydata.org/docs/using/envs.html) for details on working with environments.

### Running XYalign
XYalign is located in ```XYalign/xyalign/xyalign.py``` and can be run with a command along the lines of:
```
python <path/to/xyalign.py> --ref <path/to/reference.fasta> --bam </path/to/bam --sample_id <name_of_sample> --output_dir <path/to/outputdirectory - will be created if doesn't already exist> --cpus <number of cores/threads to use>
```
You can see a full list of options with the command ```python /path/to/xyalign.py -h```

## What XYalign does

## Goals, Problems, and Contributors
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
