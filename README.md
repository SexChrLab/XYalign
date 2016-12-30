<p align="center">
  <img src="https://github.com/WilsonSayresLab/XYalign/blob/master/Files/XYlogo.png" width="500"/>
</p>

# XYalign: Inferring Sex Chromosome Ploidy in NGS Data
Timothy H. Webster, Tanya Phung, Madeline Couse, Bruno Grande, Eric Karlins, Phillip Richmond, Whitney Whitford, Melissa A. Wilson Sayres

Sex chromosome aneuploidies are currently estimated to be as common as 1/400 in humans. Atypical ploidy will affect variant calling and measures of genomic variation that are central to most clinical genomic studies. Further, the high degree of similarity between gametologous sequences on the X and Y chromosomes can lead to the misalignment of sequencing reads and substantially affect variant calling. Here we present XYalign, a new tool that (1) quickly infers sex chromosome ploidy in NGS data (DNA and RNA), (2) remaps reads based on the inferred sex chromosome complement of the individual, and (3) outputs quality, depth, and allele-balance metrics across the sex chromosomes.

[October 17, 2016 slide show of concept and inital results for Hackseq](https://docs.google.com/presentation/d/1OB2d_mu5zC742N_NKfzHjVpUm4BFtm5lUzniLLI--OQ/edit?usp=sharing)

[November 18, 2016 poster presentation for the Evolutionary Genomics of Sex 2016 Meeting](https://figshare.com/articles/XYalign_Inferring_and_Correcting_for_Sex_Chromosome_Ploidy_in_Next-Generation_Sequencing_Data/4292924/1)

## Please note the XYalign is still in pre-release stages
Changes are made to the master branch quite frequently, so there still might be a number of known or unknown bugs at any given time.  If you'd like to try things before release (which will be mid-to-late December), please feel free to contact us (either via an Issue or by email to Tim Webster - address at the bottom of this page) to inquire about the current status of the program.

## Using XYalign

See full documentation at [Read The Docs](http://xyalign.readthedocs.io/en/latest/index.html)

## To-Do Items Before Release

Program
- [x] Replace current logging and printing with logging module
- [x] Fix permutation test
- [x] Add KS 2 Sample Test
- [x] Add bootstrap
- [ ] Incorporate likelihood model for ploidy estimation
- [x] Add option to create new references up front (and not have to run process for every sample)
	- [x] Add option to directly read in preformatted references
- [ ] Test exact depth (vs. current approximation) - time and results
	- [ ] Implement exact depth if not too slow or if results differ substantially
- [x] Add ability to take targets instead of fixed window sizes for traverse_bam
- [x] Ensure all files have sample id in prefix to allow writing multiple samples to same directory
- [ ] Add "batch" option (for CHARACTERIZE_SEX_CHROMS and/or ANALYZE_BAM)
- [x] Create classes for reference fastas and bam/cram files
- [x] Add support for arbitrary bwa flags
- [ ] Add support for additional mappers
	- not clear if this is necessary or valuable
- [ ] Handle duplicate reads
	- not clear if this is necessary or valuable
- [x] Write @PG line to final bam
- [x] Accept and integrate external bed files for excluding regions
- [x] Second round of bam window analysis and plotting after remapping and filtering low-quality regions
- [ ] Output high and low-quality bed files after second round of bam analyses
- [ ] Add histogram plotting for depth and mapq
- [ ] Add cram and sam support
	- probably not necessary for initial release
- [ ] Clean up parse_args() and update checks/validation
- [x] Allow PREPARE_REFERENCE and REMAPPING to run without window_size or target_bed
- [x] Move subsections of pipeline to individual functions to clean up main()
- [ ] Implement check for dependencies upon loading
- [ ] Add thorough module testing (e.g. unittest, doctest, etc.)
- [x] Write documentation
- [ ] Implement some kind of checkpointing functionality
	- probably not necessary for initial release

Testing
- [x] Simulate and assemble XX, XY, XXY, and XO genomes
	- [ ] Test XYalign on these genomes
- [x] Download data for XX and XY:
	- [x] high coverage whole genome
		- [x] Test XYalign
	- [x] low coverage whole genome
		- [x] Test XYalign
	- [x] exome
		- [x] Test XYalign
- [ ] In an example case, plot runtimes with increasing numbers of threads (e.g., 1, 2, 4, 8, 16)


### Group Members
Name | email | github ID
--- | --- |  ---
Tim Webster | timothy.h.webster@asu.edu | @thw17
Tanya Phung | tnphung@ucla.edu | @tnphung
Madeline Couse| mhcouse@gmail.com | @Madelinehazel
Bruno Grande | bgrande@sfu.ca | @brunogrande
Eric Karlins | karlinser@mail.nih.gov | @ekarlins
Phillip Richmond | phillip.a.richmond@gmail.com | @Phillip-a-Richmond
Whitney Whitford | whitney.whitford@auckland.ac.nz | @whitneywhitford
Melissa A. Wilson Sayres | melissa.wilsonsayres@asu.edu | @mwilsonsayres
