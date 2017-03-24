<p align="center">
  <img src="https://github.com/WilsonSayresLab/XYalign/blob/master/Files/XYlogo.png" width="500"/>
</p>

# XYalign: Inferring Sex Chromosome Ploidy in NGS Data
Timothy H. Webster, Tanya Phung, Madeline Couse, Bruno Grande, Eric Karlins, Phillip Richmond, Whitney Whitford, Melissa A. Wilson Sayres

The high degree of similarity between gametologous sequences on the sex chromosomes can lead to the misalignment of sequencing reads and substantially affect variant calling. Here we present XYalign, a new tool that (1) quickly infers sex chromosome ploidy in NGS data, (2) remaps reads based on the inferred sex chromosome complement of the individual, and (3) outputs quality, depth, and allele-balance metrics across chromosomes.

[October 17, 2016 slide show of concept and inital results for Hackseq](https://docs.google.com/presentation/d/1OB2d_mu5zC742N_NKfzHjVpUm4BFtm5lUzniLLI--OQ/edit?usp=sharing)

[November 18, 2016 poster presentation for the Evolutionary Genomics of Sex 2016 Meeting](https://figshare.com/articles/XYalign_Inferring_and_Correcting_for_Sex_Chromosome_Ploidy_in_Next-Generation_Sequencing_Data/4292924/1)

## Please note the XYalign is still in pre-release stages
Changes are made to the master branch quite frequently, so there still might be a number of known or unknown bugs at any given time.  If you'd like to try things before release (which will be mid-April), please feel free to contact us (either via an Issue or by email to Tim Webster - address at the bottom of this page) to inquire about the current status of the program.

## Using XYalign

See full documentation at [Read The Docs](http://xyalign.readthedocs.io/en/latest/index.html)

## To-Do Items

*Before release*
- [x] Incorporate ploidy caller
- [x] Default to calculating mean depth per chromosome
	- [x] add flag to allow for whole-genome thresholds
- [ ] Make bwa indexing conditional for withY if no reference_mask provided
- [ ] Clean up parse_args() and update checks/validation
- [x] Implement check for dependencies upon loading

*Down the road*
- [ ] Add thorough module testing
- [ ] Test exact depth (vs. current approximation) - time and results
	- [ ] Implement exact depth if not too slow or if results differ substantially
	- [ ] Implement bam analyses in cython
- [ ] Add "batch" option (for CHARACTERIZE_SEX_CHROMS and/or ANALYZE_BAM)
	- not necessary for now, as can be easily done with a bash command
- [ ] Add support for additional mappers
	- not clear if this is necessary or valuable
- [ ] Handle duplicate reads
	- not clear if this is necessary or valuable
- [ ] Output high and low-quality bed files after second round of bam analyses
- [ ] Add histogram plotting for depth and mapq
- [ ] Add cram and sam support
- [ ] Implement some kind of checkpointing functionality
- [ ] Run bwa indexing in parallel
	- this might be a bit more difficult that it seems due to memory requirements for large genomes
- [ ] Implement subparser

## Group Members
Name | email | github ID
--- | --- |  ---
Tim Webster | timothy.h.webster@asu.edu | @thw17
Madeline Couse| mhcouse@gmail.com | @Madelinehazel
Bruno Grande | bgrande@sfu.ca | @brunogrande
Eric Karlins | karlinser@mail.nih.gov | @ekarlins
Tanya Phung | tnphung@ucla.edu | @tnphung
Phillip Richmond | phillip.a.richmond@gmail.com | @Phillip-a-Richmond
Whitney Whitford | whitney.whitford@auckland.ac.nz | @whitneywhitford
Melissa A. Wilson Sayres | melissa.wilsonsayres@asu.edu | @mwilsonsayres
