<p align="center">
  <img src="https://github.com/WilsonSayresLab/XYalign/blob/master/Files/XYlogo.png" width="500"/>
</p>

## Background
The high degree of similarity between gametologous sequences on the sex chromosomes can lead to the misalignment of sequencing reads and substantially affect variant calling. Here we present XYalign, a new tool that (1) quickly infers sex chromosome ploidy in NGS data, (2) remaps reads based on the inferred sex chromosome complement of the individual, and (3) outputs quality, depth, and allele-balance metrics across chromosomes.

## Preprint

Please see our preprint for more information:

*Identifying, understanding, and correcting technical biases on the sex chromosomes in next-generation sequencing data.* 2018. Webster TH; Couse M; Grande BM; Karlins E;
Phung T; Richmond PA; Whitford W; Wilson Sayres MA. bioRxiv 346940; doi: https://doi.org/10.1101/346940

If you use XYalign or discuss/correct for bias in mapping on the sex chromosomes, please cite this preprint.

## Using XYalign

See full documentation at [Read The Docs](http://xyalign.readthedocs.io/en/latest/index.html) -- Under construction

Post any questions you have at the [XYalign Google Group](https://groups.google.com/forum/#!forum/xyalign)

Post any bugs/issues to [XYalign's issues page on Github](https://github.com/WilsonSayresLab/XYalign/issues)

## Quick start and examples

### Installing XYalign
XYalign has only been tested on Linux and Mac systems. We recommend users install and manage XYalign (and programming environments) using
Conda. To do this

1. First download and install either
[Anaconda](https://www.continuum.io/downloads)
or [Miniconda](http://conda.pydata.org/miniconda.html) (both work well,
Miniconda is a lightweight version of Anaconda).

2. Finish installation with the following commands to install XYalign and all
of its dependancies in an environment called "xyalign_env":

```
conda config --add channels r

conda config --add channels defaults

conda config --add channels conda-forge

conda config --add channels bioconda

conda create -n xyalign_env xyalign

```

3. Load your new environment (containing XYalign and all related programs) with:

```
source activate xyalign_env
```

See [Bioconda](https://bioconda.github.io/) and [Conda](https://conda.io/docs/user-guide/tasks/manage-environments.html) documentation
for more information.

### Prepare a sex-specific reference genome
Assuming XYalign is installed correctly with all associated programs and is available
in your ``PATH`` (see "Installing XYalign above"), you can use the command
(assume the following is on one line):

```
xyalign --PREPARE_REFERENCE --ref reference.fasta
--xx_ref_out /path/to/reference.XXonly.fasta
--xy_ref_out /path/to/reference.XY.fasta
--x_chromosome chrX
--y_chromosome chrY
--reference_mask mask.bed
```

In the above command, ``reference.fasta`` is the original reference genome,
``/path/to/reference.XXonly.fasta`` and ``/path/to/reference.XY.fasta`` are the
full paths to and names of the desired output references for XX and XY samples,
respectively. ``chrX`` and ``chrY`` are the *exact* names of the X and Y chromosome
scaffolds in the assembly. ``mask.bed`` is some bed file containing regions that
should be masked in *both* output fastas.

### Analyze a single bam file to explore sex chromosome content, etc.
You can use the command (assume the following is on one line):

```
xyalign --CHARACTERIZE_SEX_CHROMS
--ref reference.fasta
--bam sample1.bam
--output_dir sample1_results
--sample_id sample1
--cpus 4
--window_size 5000
--chromosomes chr19 chrX chrY
--x_chromosome chrX
--y_chromosome chrY
```

In the above command, ``reference.fasta`` is the full path to the reference genome
*used to generate the bam file*, ``sample1.bam`` is the full path to the bam file
``sample1_results`` is our desired output directory, and ``sample1`` is the name of
our sample. we're using four cores (``--cpus 4``) and 5kb nonoverlapping
windows for analysis. We're analyzing three chromosomes named ``chr19``,
``chrX``, and ``chrY``, and our X and Y scaffolds in the reference are named
``chrX`` and ``chrY``.

Our output of interest will be in ``sample1_results/plots``
and ``sample1_results/results``. Tables (.csv) of depth and mapq measurements per window
will in ``sample1_results/bed`` with "full_dataframe" in their file names. BED files containing windows passing ("highquality") and failing ("lowquality") filtering
thresholds will also be in ``sample1_results/bed``.

Relevant flags for filtering variants include:

```
	--variant_site_quality
	--variant_genotype_quality
	--variant_depth
```

Relevant flags for filtering windows include:

```
	--mapq_cutoff
	--min_depth_filter
	--max_depth_filter
	--min_variant_count
```

You can get details about these (and more) flags with the command:

```
	xyalign -h
```

### Analyze multiple bam files to determine sex chromosome complement, identify sex chromosome scaffolds, etc.

```
xyalign --CHROM_STATS
--chromosomes chr1 chr8 chr19 chrX chrY
--bam sample1.bam sample2.bam sample3.bam
--ref null
--sample_id bam_comparison1
--output_dir bam_comparison1_results
```

In the above command, we're analyzing five chromosomes in three different bam files.
We provide ``null`` as our reference because it's not used in these analyses.
``--sample_id`` now becomes the name of our comparison (it's used in file names, etc.)
and our output will be located in ``bam_comparison1_results/results``. We could also use
``--use_counts`` to force XYalign to simply use counts of reads on each chromosome in
comparisons.

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
