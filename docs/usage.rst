Usage
=====

The Basics
-----------

Requirements
~~~~~~~~~~~~

Currently, XYalign requires two files to run: a **bam file** and the **reference fasta file**
used to generate it.  It also requires a *list of chromosomes* to analyze (at least 3, including X and Y),
the *name of the X chromosome*, and the *name of the Y chromosome*. The chromosome names must *exactly* match
those in the bam header and reference fasta - 'chr19' is not equivalent to '19', for example.

As we are still finalizing methods for estimating sex chromosome complement (e.g., XX, XY, XXY, X0), for the
time being, **you also need to provide** ``--y_present`` or ``--y_absent`` to indicate whether XYalign should
treat the sample as possessing a Y chromosome or not.

You also need a variety of python packages and external programs installed.  See
:doc:`installation` for more information.

The Pipeline
~~~~~~~~~~~~

Xyalign is composed of the following modules that can be thought of as steps in the pipeline::

	PREPARE_REFERENCE
	ANALYZE_BAM
	CHARACTERIZE_SEX_CHROMS
	REMAPPING

Each of these modules can be invoked as a command line flag with no arguments
(e.g., ``--PREPARE_REFERENCE``), and XYalign will execute *only that module*.  If no flags
are provided, XYalign will run the full pipeline in the following order: PREPARE_REFERENCE ->
ANALYZE_BAM -> CHARACTERIZE_SEX_CHROMS -> REMAPPING -> ANALYZE_BAM.  This will:

	1. Prepare two reference genomes - one with the Y chromosome masked, the other with both X and Y
	unmasked.  In both cases, XYalign will optionally mask other regions of the genome provided in an
	input bed file (using the flag ``--reference_mask <file1.bed> <file2.bed> ...``).

	2. Analyze the bam file to calculate metrics related to read balance, read depth, and mapping quality.
	Read depth and mapping quality are calculated in windows, and either ``--window_size <integer window size>``
	or ``--target_bed <path to target bed file>`` must be provided.  ``--window_size`` is the fixed size
	of windows to use in a sliding window analysis in bases (e.g., ``10000`` for 10 kb windows).  ``--target_bed``
	is a bed file of targets to use as windows, e.g. exome capture targets.

	3. Plot read balance, depth, and mapq for each chromosome, and output bed files of high
	and low quality regions.

	4. Run a series of tests comparing ANALYZE_BAM metrics for each chromosome. If the flag
	``--CHARACTERIZE_SEX_CHROMS`` is invoked, XYalign will carry out the bam analysis steps above
	and then proceed to these tests.

	5. Strip and sort reads mapping to the sex chromosomes, map to the reference with
	the appropriate masking (step 1) based on the results of step 4, and replace the sex
	chromosome alignments in the original bam file with these new ones.

	6. Analyze the new bam file as in steps 4 and 5.

Suggested Command Lines
~~~~~~~~~~~~~~~~~~~~~~~

Below we highlight example command lines, as well as useful optional flags for
each module (PREPARE_REFERENCE, ANALYZE_BAM, CHARACTERIZE_SEX_CHROMS, REMAPPING)
as well as the full pipeline.  You can find a complete list of command line flags,
their descriptions, and their defaults listed at the bottom of the page.  You can
also get this list from the command line::

	python xyalign.py -h

In all examples, ``reference.fasta`` is our input reference in fasta format, ``input.bam``
is our input bam file (created using ``reference.fasta``), ``sample1`` is the ID of our
sample, and ``sample1_output`` is the name of our desired output directory.  We'll
analyze chromosomes named 'chr19', 'chrX', and 'chrY', with chrX representing the X chromosome
and chrY representing the Y chromosome (these are XYalign defaults).  We'll assume that all programs are in
our ``PATH`` and can be invoked by typing the program name from the command line
without any associated path (e.g., ``samtools``).  We'll also assume that we're
working on a cluster with 4 cores available to XYalign.

1. PREPARE_REFERENCE
::

	python xyalign.py --PREPARE_REFERENCE --ref reference.fasta --bam input.bam \
	--output_dir sample1_output --sample_id sample1 --cpus 4 --reference_mask mask.bed

Here, ``mask.bed`` is a bed file containing regions to mask in *both* output reference
genomes (e.g., coordinates for the pseudoautosomal regions on the Y chromosome).  More
than one can be included as well (e.g., ``--reference_mask mask.bed mask2.bed``).

This will output two reference genomes, one with the Y chromosome completely masked
(defaults to ``sample1_output/reference/xyalign_noY.masked.fa``) and one with
an unmasked Y (defaults to ``sample1_output/reference/xyalign_withY.masked.fa``). These
defaults can be changed with the ``--xx_ref_out`` and ``--xy_ref_out`` flags.

2. ANALYZE_BAM
::

	python xyalign.py --ANALYZE_BAM --ref reference.fasta --bam input.bam \
	--output_dir sample1_output --sample_id sample1 --cpus 4 --window_size 10000

Here, 10000 is the fixed window size to use in (nonoverlapping) sliding window
analyses of the bam file.  If you're working with targeted sequencing data (e.g. exome),
you can provide a list of regions to use instead of windows.  For example, if your
regions are in ``targets.bed`` you would add the flag: ``--targed_bed targets.bed``.

This command line will default to a minimum quality of 20 (SNP), minimum
mapping quality of 20 (bam window), and a depth filter of 4 (bam window).  These
can be set with the flags ``--variant_quality_cutoff``, ``--mapq_cutoff``, and
``--depth_filter``, respectively. A quick note about the depth filter: the formula
we use is based off of Li's (2004, Bioinformatics 30: 2843-2851) analysis of artifacts
in sequencing data.  For SNPs, he recommends a depth_filter of 3 or 4 using the equation
mean_depth +- (depth_filter * square_root(mean_depth)).  We use this formula for our filter
as well.

This will output a series of plots in ``sample1_output/plots`` and bed files containing
high and low quality windows in ``sample1_output/plots``.

3. CHARACTERIZE_SEX_CHROMS
::

	python xyalign.py --CHARACTERIZE_SEX_CHROMS --ref reference.fasta --bam input.bam \
	--output_dir sample1_output --sample_id sample1 --cpus 4 --window_size 10000

Settings here are identical to 3 because the first step of CHARACTERIZE_SEX_CHROMS
involves running ANALYZE_BAM.

In addition to everything in ANALYZE_BAM, CHARACTERIZE_SEX_CHROMS will output the
results of a series of statistical tests in ``sample1_output/results``.

4. REMAPPING
::

	python xyalign.py --REMAPPING --ref reference.fasta --bam input.bam \
	--output_dir sample1_output --sample_id sample1 --cpus 4 \
	--xx_ref_in sample1_output/reference/xyalign_noY.masked.fa \
	--xy_ref_in sample1_output/reference/xyalign_withY.masked.fa \
	--y_absent

Here, we've input our reference genomes generated in step 1 (if we don't, XYalign
will repeat that step).  We've also used the flag ``--y_absent`` to indicate that
there is no Y chromosome in our sample (perhaps as the result of step 3, or outside
knowledge).  If a Y is present, we would have used ``--y_present`` instead.  REMAPPING
requires one of those two flags, as it does not involve any steps to estimate
sex chromosome content (those are carried out in CHARACTERIZE_SEX_CHROMS).

5. Full pipeline

And if we want to run the full XYalign pipeline on a sample, we'd use a command line
along the lines of::

	python xyalign.py --PREPARE_REFERENCE --ref reference.fasta --bam input.bam \
	--output_dir sample1_output --sample_id sample1 --cpus 4 --reference_mask mask.bed \
	--window_size 10000

We could optionally provided preprocessed reference genomes with ``--xx_ref_in``
and ``--xx_ref_in``, as in 4.  We could also use ``--y_absent`` or ``--y_present``
to force mapping to a certain reference.

.. note::

	We are currently experimenting with methods for determining the presence or
	absence of a Y chromosome, so either ``--y_absent`` or ``--y_present`` is
	required for the time being until we've finalized the implementation of ploidy
	estimation.

Recommendations for Incorporating XYalign into Pipelines
--------------------------------------------------------

While the full XYalign pipeline will be useful in certain situations, we feel that
the following pipeline is better suited to most users' needs and will save time and space.

1. Use XYalign PREPARE_REFERENCE to prepare Y present and Y absent genomes.

2. Preliminarily map reads to the standard reference (or Y present) and sort the bam file
using any mapper and sorting algorithm.

3. Run CHARACTERIZE_SEX_CHROMS, to analyze the bam file, output plots, and estimate
ploidy.

4. Remap reads to the fasta produced in 1 corresponding to the sex chromosome
complement characterized in 3.  E.g., if Y is not detected, map to Y absent.  This time
run full pipeline of mapping, sorting, removing duplicates, etc., using users' preferred
tools/pipeline.

5. Optionally run ANALYZE_BAM on bam file produced in 4.

6. Call variants using user-preferred caller.

7. Analyze variants taking into account ploidy estimated in 3, and consider masking
low quality regions using bed files output in 5.

Exome data
----------

XYalign handles exome data, with a few minor considerations.  In particular, either setting
``--window_size`` to a smaller value, perhaps 5000 or less, or inputting
targets to use in lieu of windows (``--target_bed targets.bed``) will be critical
for getting more accurate window measures.  In addition, users should manually
check the results of CHARACTERIZE_SEX_CHROMS for a number of samples to get a feel
for expected values on the sex chromosomes, as these values are likely to vary among
experimental design (especially among different capture kits).

Nonhuman genomes
----------------

Analyzing arbitrary chromosomes
-------------------------------

Using XYalign as a Python library
---------------------------------

Full List of Command-Line Flags
-------------------------------
