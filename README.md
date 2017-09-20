<p align="center">
  <img src="https://github.com/WilsonSayresLab/XYalign/blob/master/Files/XYlogo.png" width="500"/>
</p>

# XYalign: Inferring Sex Chromosome Ploidy in NGS Data
Timothy H. Webster, Madeline Couse, Bruno Grande, Eric Karlins, Tanya Phung, Phillip Richmond, Whitney Whitford, Melissa A. Wilson Sayres

## Background
The high degree of similarity between gametologous sequences on the sex chromosomes can lead to the misalignment of sequencing reads and substantially affect variant calling. Here we present XYalign, a new tool that (1) quickly infers sex chromosome ploidy in NGS data, (2) remaps reads based on the inferred sex chromosome complement of the individual, and (3) outputs quality, depth, and allele-balance metrics across chromosomes.

[October 17, 2016 slide show of concept and inital results for Hackseq](https://docs.google.com/presentation/d/1OB2d_mu5zC742N_NKfzHjVpUm4BFtm5lUzniLLI--OQ/edit?usp=sharing)

[November 18, 2016 poster presentation for the Evolutionary Genomics of Sex 2016 Meeting](https://figshare.com/articles/XYalign_Inferring_and_Correcting_for_Sex_Chromosome_Ploidy_in_Next-Generation_Sequencing_Data/4292924/1)

## Using XYalign

See full documentation at [Read The Docs](http://xyalign.readthedocs.io/en/latest/index.html)

Post any questions you have at the [XYalign Google Group](https://groups.google.com/forum/#!forum/xyalign)

## If you cloned XYalign before April 5, 2017 please clone a fresh version
On April 5th, we removed a host of large files using BFG Repo Cleaner, which changed
the commit history of the repo.  Because of this, ```git pull``` will no longer work on clones
prior to this cleaning.  Therefore, if you cloned before this date (or can't remember when you
cloned), it's probably a good idea to clone a fresh version.

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
