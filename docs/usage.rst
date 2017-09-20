Usage Overview
==============

The Basics
-----------

Requirements
~~~~~~~~~~~~

The different modules of XYalign have slightly different requirements, but in
general you'll need: a **bam file** and the **reference fasta file**
used to generate it.  XYalign also requires a *list of chromosomes* to analyze,
the *name of the X chromosome*, and the *name of the Y chromosome* (if in the assembly). The chromosome names must *exactly* match
those in the bam header and reference fasta - 'chr19' is not equivalent to '19', for example.

You also need a variety of python packages and external programs installed.  See
:doc:`installation` for more information.

The Pipeline
~~~~~~~~~~~~

Xyalign is composed of the following modules that can be thought of as steps in the pipeline (with the exception of CHROM_STATS)::

	PREPARE_REFERENCE
	ANALYZE_BAM
	CHARACTERIZE_SEX_CHROMS
	STRIP_READS
	REMAPPING
	CHROM_STATS

Each of these modules can be invoked as a command line flag with no arguments
(e.g., ``--PREPARE_REFERENCE``), and XYalign will execute *only that module*.  If no flags
are provided, XYalign will run the full pipeline in the following order: PREPARE_REFERENCE ->
ANALYZE_BAM -> CHARACTERIZE_SEX_CHROMS -> STRIP_READS -> REMAPPING -> ANALYZE_BAM.  This will:

	1. Prepare two reference genomes - one with the Y chromosome masked, the other with both X and Y
	unmasked.  In both cases, XYalign will optionally mask other regions of the genome provided in an
	input bed file (using the flag ``--reference_mask <file1.bed> <file2.bed> ...``).

	2. Analyze the bam file to calculate metrics related to read balance, read depth, and mapping quality.
	Read depth and mapping quality are calculated in windows, and either ``--window_size <integer window size>``
	or ``--target_bed <path to target bed file>`` must be provided.  ``--window_size`` is the fixed size
	of windows to use in a nonoverlapping sliding window analysis in bases (e.g., ``10000`` for 10 kb windows).  ``--target_bed``
	is a bed file of targets to use as windows, e.g. exome capture targets.

	3. Plot read balance, depth, and mapq for each chromosome, and output bed files of high
	and low quality regions, based on either default or user-defined thresholds.

	4. Run a series of tests comparing ANALYZE_BAM metrics for each chromosome. If the flag
	``--CHARACTERIZE_SEX_CHROMS`` is invoked, XYalign will carry out the bam analysis steps above
	and then proceed to these tests.

	5. Strip and sort reads mapping to the sex chromosomes, map to the reference with
	the appropriate masking (step 1) based on the results of step 4, and replace the sex
	chromosome alignments in the original bam file with these new ones.

	6. Analyze the new bam file as in steps 4 and 5.

``CHROM_STATS`` provides quicker, coarser statistics and is designed for cases in which a reference genome is well-understood
and when multiple samples are available.

Suggested Command Lines
~~~~~~~~~~~~~~~~~~~~~~~

Below we highlight example command lines, as well as useful optional flags for
each module (PREPARE_REFERENCE, ANALYZE_BAM, CHARACTERIZE_SEX_CHROMS, STRIP_READS, REMAPPING)
as well as the full pipeline.  You can find a complete list of command line flags,
their descriptions, and their defaults listed at the bottom of the page.  You can
also get this list from the command line::

	xyalign -h

In all examples, ``reference.fasta`` is our input reference in fasta format, ``input.bam``
is our input bam file (created using ``reference.fasta``), ``sample1`` is the ID of our
sample, and ``sample1_output`` is the name of our desired output directory.  We'll
analyze chromosomes named 'chr19', 'chrX', and 'chrY', with chrX representing the X chromosome
and chrY representing the Y chromosome.  We'll assume that all programs are in
our ``PATH`` and can be invoked by typing the program name from the command line
without any associated path (e.g., ``samtools``).  We'll also assume that we're
working on a cluster with 4 cores available to XYalign.

1. PREPARE_REFERENCE
::

	xyalign --PREPARE_REFERENCE --ref reference.fasta \
	--output_dir sample1_output --sample_id sample1 --cpus 4 --reference_mask mask.bed \
	--x_chromosome chrX --y_chromosome chrY

Here, ``mask.bed`` is a bed file containing regions to mask in *both* output reference
genomes (e.g., coordinates for the pseudoautosomal regions on the Y chromosome).  More
than one can be included as well (e.g., ``--reference_mask mask.bed mask2.bed``).

This will output two reference genomes, one with the Y chromosome completely masked
(defaults to ``sample1_output/reference/xyalign_noY.masked.fa``) and one with
an unmasked Y (defaults to ``sample1_output/reference/xyalign_withY.masked.fa``). These
default *names* can be changed with the ``--xx_ref_out`` and ``--xy_ref_out`` flags.
Note, however, that files will still be deposited in ``sample1_output/reference``.

2. ANALYZE_BAM
::

	xyalign --ANALYZE_BAM --ref reference.fasta --bam input.bam \
	--output_dir sample1_output --sample_id sample1 --cpus 4 --window_size 10000 \
	--chromosomes chr19 chrX chrY --x_chromosome chrX --y_chromosome chrY

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

This will output a series of plots in ``sample1_output/plots``, bed files containing
high and low quality windows in ``sample1_output/bed``, and the entire dataframe
with values for each measure in each window in ``sample1_output/bed``.

3. CHARACTERIZE_SEX_CHROMS
::

	xyalign --CHARACTERIZE_SEX_CHROMS --ref reference.fasta --bam input.bam \
	--output_dir sample1_output --sample_id sample1 --cpus 4 --window_size 10000 \
	--chromosomes chr19 chrX chrY --x_chromosome chrX --y_chromosome chrY

Settings here are identical to 3 because the first step of CHARACTERIZE_SEX_CHROMS
involves running ANALYZE_BAM.

In addition to everything in ANALYZE_BAM, CHARACTERIZE_SEX_CHROMS will output the
results of a series of statistical tests in ``sample1_output/results``.

4. STRIP_READS
::
	xyalign --STRIP_READS --ref reference.fasta --bam input.bam \
	--output_dir sample1_output --sample_id sample1 --cpus 4 \
	--chromosomes chr1 chr2 chr3 chr4 chr5

This will strip the reads, by read group, from chromosomes 1-5 and output
a pair of fastqs per read group, as well as the read groups themselves, and a
text file connecting fastqs with their respective read groups in the directory
`` sample1_output/fastq ``.  If we were working with single-end reads, we would
have had to include the flag `` --single_end ``.  Here, the reference file isn't
used at all (it's a general requirement of XYalign), so a dummy file can be used
in its place.  To strip reads from the entire genome (including unmapped), use
`` --chromosomes ALL``

5. REMAPPING
::

	xyalign --REMAPPING --ref reference.fasta --bam input.bam \
	--output_dir sample1_output --sample_id sample1 --cpus 4 \
	--chromosomes chr19 chrX chrY --x_chromosome chrX --y_chromosome chrY \
	--xx_ref_in sample1_output/reference/xyalign_noY.masked.fa \
	--xy_ref_in sample1_output/reference/xyalign_withY.masked.fa \
	--y_absent

Here, we've input our reference genomes generated in step 1 (if we don't, XYalign
will repeat that step).  We've also used the flag ``--y_absent`` to indicate that
there is no Y chromosome in our sample (perhaps as the result of step 3, or outside
knowledge).  If a Y is present, we would have used ``--y_present`` instead.  REMAPPING
requires one of those two flags, as it does not involve any steps to estimate
sex chromosome content (those are carried out in CHARACTERIZE_SEX_CHROMS). Note that
REMAPPING will run STRIP_READS first.

5. Full pipeline

And if we want to run the full XYalign pipeline on a sample, we'd use a command line
along the lines of::

	xyalign --ref reference.fasta --bam input.bam \
	--output_dir sample1_output --sample_id sample1 --cpus 4 --reference_mask mask.bed \
	--window_size 10000 \ --chromosomes chr19 chrX chrY --x_chromosome chrX --y_chromosome chrY

We could have optionally provided preprocessed reference genomes with ``--xx_ref_in``
and ``--xx_ref_in``, as in 4.  We could have also used ``--y_absent`` or ``--y_present``
to force mapping to a certain reference.  Because we didn't include either of these
two flags, XYalign will use ``--sex_chrom_calling_threshold`` to determine the
sex chromosome complement (default is 2.0).

6. CHROM_STATS
::
	xyalign --CHROM_STATS --use_counts --bam input1.bam input2.bam input3.bam --ref null \
	--output_dir directory_name --sample_id analysis_name --chromosomes chr19 chrX chrY

Here, ``--use_counts`` simply grabs the number of reads mapped to each chromosome from the
bam index. It's by far the fastest, yet coarsest option. Running without this flag
will calculated depth and mapq along each chromosome for more detail, but this will take longer.


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

XYalign - Speed and Memory
--------------------------
The minimum memory requirements for XYalign are determined by external programs,
rather than any internal code.  Right now, the major limiting step is bwa indexing
of reference genomes which requires 5-6 GB of memory to index a human-sized genome.  In addition,
in certain situations (e.g., removing all reads from deep coverage genome data with
a single - or no - read group) the STRIP_READS module will require a great deal
of memory to sort and match paired reads (the memory requirement is that of the
external program repair.sh).

The slowest parts of the pipeline also all involve steps relying on external programs, such as
genome preparation, variant calling, read mapping, swapping sex chromosome alignments, etc.
In almost all cases, you'll see substantial increases in the speed of the pipeline by increasing the
number of threads/cores.  You must provide information about the number of cores available
to XYalign with the ``--cpus`` flag (XYalign will assume only a single thread is
available unless this flag is set).

Exome data
----------

XYalign handles exome data, with a few minor considerations.  In particular, either setting
``--window_size`` to a smaller value, perhaps 5000 or less, or inputting
targets instead of a window size (``--target_bed targets.bed``) will be critical
for getting more accurate window measures.  In addition, users should manually
check the results of CHARACTERIZE_SEX_CHROMS for a number of samples to get a feel
for expected values on the sex chromosomes, as these values are likely to vary among
experimental design (especially among different capture kits).

Nonhuman genomes
----------------

XYalign will theoretically work with any genome, and on any combination of chromosomes
or scaffolds (see more on the latter below).  Simply provide the names of the
chromosomes/scaffolds to analyze and the names of the sex chromosomes (e.g.,
``--chromosomes chr1a chr1b chr2 lga lgb --x_chromosome lga --y_chromosome lgb``
if our x_linked scaffold was lga and y_linked scaffold was lgb, and we wanted
to compare these scaffolds to chromosomes: chr1a chr1b and chr2). However,
please note that, as of right now, XYalign does not support multiple X or Y
chromosomes/scaffolds (we are planning on supporting this soon though).

Keep in mind, however, that read balance, mapq, and depth ratios might differ
among organisms, so default XYalign settings will likely not be appropriate in
most cases.  Instead, if multiple samples are available, we recommend running
XYalign's CHARACTERIZE_SEX_CHROMS  on each sample (steps 2-3 in
"Recommendations for Incorporating XYalign into pipelines" above)
using the same output directory for all samples.  One can then quickly concatenate
results (we recommend starting with bootstrap results) and plot them to look
for clustering of samples (see the XYalign publication for examples of this).

Analyzing arbitrary chromosomes
-------------------------------

Currently, XYalign requires a minimum of two chromosomes (an "autosome and an "x chromosome")
for analyses in ANALYZE_BAM and CHARACTERIZE_SEX_CHROMS (and therefore, the whole pipeline)
These chromosomes, however, can be arbitrary. Below, we highlight two example cases:
looking for evidence of Trisomy 21 in human samples,
and running the full XYalign pipeline on a ZW sample (perhaps a bird, squamate reptile, or moth).

If one wanted to look for evidence of Trisomy 21 in human data mapped to hg19 (which uses
"chr" in chromosome names), s/he could use a command along the lines of::

	xyalign --CHARACTERIZE_SEX_CHROMS --ref reference.fasta --bam input.bam \
	--output_dir sample1_output --sample_id sample1 --cpus 4 --window_size 10000 \
	--chromosomes chr1 chr10 chr19 chr21 --x_chromosome chr21

This would run the CHARACTERIZE_SEX_CHROMS module, systematically comparing
``chr21`` to ``chr1``, ``chr10``, and ``chr19``.

To run the full pipeline on a ZW sample (in ZZ/ZW systems, males are ZZ and females
are ZW), one could simply run a command like (assuming the Z scaffold was named
"scaffoldz" and the W scaffold was named "scaffoldw")::

	xyalign --ref reference.fasta --bam input.bam \
	--output_dir sample1_output --sample_id sample1 --cpus 4 --reference_mask mask.bed \
	--window_size 10000 --chromosomes scaffold1 scaffoldz scaffoldw --x_chromosome scaffoldz \
	--y_chromosome scaffoldw

In this example, it's important that the the "X" and "Y" chromosomes are assigned in this way
because PREPARE_REFERENCE (the first step in the full pipeline) will create two
reference genomes: one with the "Y" completely masked, and one with both "X" and "Y"
unmasked.  This command will therefore create the appropriate references (a ZW and
a Z only).  Other organisms or uses might not require this consideration.

Using XYalign as a Python library
---------------------------------
All modules in the XYalign/xyalign directory are designed to support the command
line program XYalign.  However, some classes and functions might be of use in other
circumstances. If you've installed XYalign as described in :doc:`installation`, then you
should be able to import XYalign libraries just like you would for any other Python package. E.g.::

	from xyalign import bam

Or::

	import xyalign.bam


Full List of Command-Line Flags
-------------------------------
This list can also be produced with the command::
	xyalign -h

::
-h, --help            show this help message and exit
--bam [BAM [BAM ...]]
					  Full path to input bam files. If more than one
					  provided, only the first will be used for modules
					  other than --CHROM_STATS
--cram [CRAM [CRAM ...]]
					  Full path to input cram files. If more than one
					  provided, only the first will be used for modules
					  other than --CHROM_STATS. Not currently supported.
--sam [SAM [SAM ...]]
					  Full path to input sam files. If more than one
					  provided, only the first will be used for modules
					  other than --CHROM_STATS. Not currently supported.
--ref REF             REQUIRED. Path to reference sequence (including file
					  name).
--output_dir OUTPUT_DIR, -o OUTPUT_DIR
					  REQUIRED. Output directory. XYalign will create a
					  directory structure within this directory
--chromosomes [CHROMOSOMES [CHROMOSOMES ...]], -c [CHROMOSOMES [CHROMOSOMES ...]]
					  Chromosomes to analyze (names must match reference
					  exactly). For humans, we recommend at least chr19,
					  chrX, chrY. Generally, we suggest including the sex
					  chromosomes and at least one autosome. To analyze all
					  chromosomes use '--chromosomes ALL' or '--chromosomes
					  all'.
--x_chromosome [X_CHROMOSOME [X_CHROMOSOME ...]], -x [X_CHROMOSOME [X_CHROMOSOME ...]]
					  Names of x-linked scaffolds in reference fasta (must
					  match reference exactly).
--y_chromosome [Y_CHROMOSOME [Y_CHROMOSOME ...]], -y [Y_CHROMOSOME [Y_CHROMOSOME ...]]
					  Names of y-linked scaffolds in reference fasta (must
					  match reference exactly). Defaults to chrY. Give None
					  if using an assembly without a Y chromosome
--sample_id SAMPLE_ID, -id SAMPLE_ID
					  Name/ID of sample - for use in plot titles and file
					  naming. Default is sample
--cpus CPUS           Number of cores/threads to use. Default is 1
--xmx XMX             Memory to be provided to java programs via -Xmx. E.g.,
					  use the flag '--xmx 4g' to pass '-Xmx4g' as a flag
					  when running java programs (currently just repair.sh).
					  Default is 'None' (i.e., nothing provided on the
					  command line), which will allow repair.sh to
					  automatically allocate memory. Note that if you're
					  using --STRIP_READS on deep coverage whole genome
					  data, you might need quite a bit of memory, e.g. '--
					  xmx 16g', '--xmx 32g', or more depending on how many
					  reads are present per read group.
--fastq_compression {0,1,2,3,4,5,6,7,8,9}
					  Compression level for fastqs output from repair.sh.
					  Between (inclusive) 0 and 9. Default is 3. 1 through 9
					  indicate compression levels. If 0, fastqs will be
					  uncompressed.
--single_end          Include flag if reads are single-end and NOT paired-
					  end.
--version, -V         Print version and exit.
--no_cleanup          Include flag to preserve temporary files.
--PREPARE_REFERENCE   This flag will limit XYalign to only preparing
					  reference fastas for individuals with and without Y
					  chromosomes. These fastas can then be passed with each
					  sample to save subsequent processing time.
--CHROM_STATS         This flag will limit XYalign to only analyzing
					  provided bam files for depth and mapq across entire
					  chromosomes.
--ANALYZE_BAM         This flag will limit XYalign to only analyzing the bam
					  file for depth, mapq, and (optionally) read balance
					  and outputting plots.
--CHARACTERIZE_SEX_CHROMS
					  This flag will limit XYalign to the steps required to
					  characterize sex chromosome content (i.e., analyzing
					  the bam for depth, mapq, and read balance and running
					  statistical tests to help infer ploidy)
--REMAPPING           This flag will limit XYalign to only the steps
					  required to strip reads and remap to masked
					  references. If masked references are not provided,
					  they will be created.
--STRIP_READS         This flag will limit XYalign to only the steps
					  required to strip reads from a provided bam file.
--logfile LOGFILE     Name of logfile. Will overwrite if exists. Default is
					  sample_xyalign.log
--reporting_level {DEBUG,INFO,ERROR,CRITICAL}
					  Set level of messages printed to console. Default is
					  'INFO'. Choose from (in decreasing amount of
					  reporting) DEBUG, INFO, ERROR or CRITICAL
--platypus_path PLATYPUS_PATH
					  Path to platypus. Default is 'platypus'. If platypus
					  is not directly callable (e.g., '/path/to/platypus' or
					  '/path/to/Playpus.py'), then provide path to python as
					  well (e.g., '/path/to/python /path/to/platypus'). In
					  addition, be sure provided python is version 2. See
					  the documentation for more information about setting
					  up an anaconda environment.
--bwa_path BWA_PATH   Path to bwa. Default is 'bwa'
--samtools_path SAMTOOLS_PATH
					  Path to samtools. Default is 'samtools'
--repairsh_path REPAIRSH_PATH
					  Path to bbmap's repair.sh script. Default is
					  'repair.sh'
--shufflesh_path SHUFFLESH_PATH
					  Path to bbmap's shuffle.sh script. Default is
					  'shuffle.sh'
--sambamba_path SAMBAMBA_PATH
					  Path to sambamba. Default is 'sambamba'
--bedtools_path BEDTOOLS_PATH
					  Path to bedtools. Default is 'bedtools'
--platypus_calling {both,none,before,after}
					  Platypus calling withing the pipeline (before
					  processing, after processing, both, or neither).
					  Options: both, none, before, after.
--no_variant_plots    Include flag to prevent plotting read balance from VCF
					  files.
--no_bam_analysis     Include flag to prevent depth/mapq analysis of bam
					  file. Used to isolate platypus_calling.
--skip_compatibility_check
					  Include flag to prevent check of compatibility between
					  input bam and reference fasta
--no_perm_test        Include flag to turn off permutation tests.
--no_ks_test          Include flag to turn off KS Two Sample tests.
--no_bootstrap        Include flag to turn off bootstrap analyses. Requires
					  either --y_present, --y_absent, or
					  --sex_chrom_calling_threshold if running full
					  pipeline.
--variant_site_quality VARIANT_SITE_QUALITY, -vsq VARIANT_SITE_QUALITY
					  Consider all SNPs with a site quality (QUAL) greater
					  than or equal to this value. Default is 30.
--variant_genotype_quality VARIANT_GENOTYPE_QUALITY, -vgq VARIANT_GENOTYPE_QUALITY
					  Consider all SNPs with a sample genotype quality
					  greater than or equal to this value. Default is 30.
--variant_depth VARIANT_DEPTH, -vd VARIANT_DEPTH
					  Consider all SNPs with a sample depth greater than or
					  equal to this value. Default is 4.
--platypus_logfile PLATYPUS_LOGFILE
					  Prefix to use for Platypus log files. Will default to
					  the sample_id argument provided
--homogenize_read_balance HOMOGENIZE_READ_BALANCE
					  If True, read balance values will be transformed by
					  subtracting each value from 1. For example, 0.25 and
					  0.75 would be treated equivalently. Default is False.
--min_variant_count MIN_VARIANT_COUNT
					  Minimum number of variants in a window for the read
					  balance of that window to be plotted. Note that this
					  does not affect plotting of variant counts. Default is
					  1, though we note that many window averages will be
					  meaningless at this setting.
--reference_mask [REFERENCE_MASK [REFERENCE_MASK ...]]
					  Bed file containing regions to replace with Ns in the
					  sex chromosome reference. Examples might include the
					  pseudoautosomal regions on the Y to force all
					  mapping/calling on those regions of the X chromosome.
					  Default is None.
--xx_ref_out XX_REF_OUT
					  Desired name for masked output fasta for samples
					  WITHOUT a Y chromosome (e.g., XX, XXX, XO, etc.).
					  Defaults to 'xyalign_noY.masked.fa'. Will be output in
					  the XYalign reference directory.
--xy_ref_out XY_REF_OUT
					  Desired name for masked output fasta for samples WITH
					  a Y chromosome (e.g., XY, XXY, etc.). Defaults to
					  'xyalign_withY.masked.fa'. Will be output in the
					  XYalign reference directory.
--xx_ref_in XX_REF_IN
					  Path to preprocessed reference fasta to be used for
					  remapping in X0 or XX samples. Default is None. If
					  none, will produce a sample-specific reference for
					  remapping.
--xy_ref_in XY_REF_IN
					  Path to preprocessed reference fasta to be used for
					  remapping in samples containing Y chromosome. Default
					  is None. If none, will produce a sample-specific
					  reference for remapping.
--read_group_id READ_GROUP_ID
					  If read groups are present in a bam file, they are
					  used by default in remapping steps. However, if read
					  groups are not present in a file, there are two
					  options for proceeding. If '--read_group_id None' is
					  provided (case sensitive), then no read groups will be
					  used in subsequent mapping steps. Otherwise, any other
					  string provided to this flag will be used as a read
					  group ID. Default is '--read_group_id xyalign'
--bwa_flags BWA_FLAGS
					  Provide a string (in quotes, with spaces between
					  arguments) for additional flags desired for BWA
					  mapping (other than -R and -t). Example: '-M -T 20 -v
					  4'. Note that those are spaces between arguments.
--sex_chrom_bam_only  This flag skips merging the new sex chromosome bam
					  file back into the original bam file (i.e., sex chrom
					  swapping). This will output a bam file containing only
					  the newly remapped sex chromosomes.
--sex_chrom_calling_threshold SEX_CHROM_CALLING_THRESHOLD
					  This is the *maximum* filtered X/Y depth ratio for an
					  individual to be considered as having heterogametic
					  sex chromsomes (e.g., XY) for the REMAPPING module of
					  XYalign. Note here that X and Y chromosomes are simply
					  the chromosomes that have been designated as X and Y
					  via --x_chromosome and --y_chromosome. Keep in mind
					  that the ideal threshold will vary according to sex
					  determination mechanism, sequence homology between the
					  sex chromosomes, reference genome, sequencing methods,
					  etc. See documentation for more detail. Default is
					  2.0, which we found to be reasonable for exome, low-
					  coverage whole-genome, and high-coverage whole-genome
					  human data.
--y_present           Overrides sex chr estimation by XYalign and remaps
					  with Y present.
--y_absent            Overrides sex chr estimation by XY align and remaps
					  with Y absent.
--window_size WINDOW_SIZE, -w WINDOW_SIZE
					  Window size (integer) for sliding window calculations.
					  Default is 50000. Default is None. If set to None,
					  will use targets provided using --target_bed.
--target_bed TARGET_BED
					  Bed file containing targets to use in sliding window
					  analyses instead of a fixed window width. Either this
					  or --window_size needs to be set. Default is None,
					  which will use window size provided with
					  --window_size. If not None, and --window_size is None,
					  analyses will use targets in provided file. Must be
					  typical bed format, 0-based indexing, with the first
					  three columns containing the chromosome name, start
					  coordinate, stop coordinate.
--exact_depth         Calculate exact depth within windows, else use much
					  faster approximation. *Currently exact is not
					  implemented*. Default is False.
--whole_genome_threshold
					  This flag will calculate the depth filter threshold
					  based on all values from across the genome. By
					  default, thresholds are calculated per chromosome.
--mapq_cutoff MAPQ_CUTOFF, -mq MAPQ_CUTOFF
					  Minimum mean mapq threshold for a window to be
					  considered high quality. Default is 20.
--min_depth_filter MIN_DEPTH_FILTER
					  Minimum depth threshold for a window to be considered
					  high quality. Calculated as mean depth *
					  min_depth_filter. So, a min_depth_filter of 0.2 would
					  require at least a minimum depth of 2 if the mean
					  depth was 10. Default is 0.0 to consider all windows.
--max_depth_filter MAX_DEPTH_FILTER
					  Maximum depth threshold for a window to be considered
					  high quality. Calculated as mean depth *
					  max_depth_filter. So, a max_depth_filter of 4 would
					  require depths to be less than or equal to 40 if the
					  mean depth was 10. Default is 10000.0 to consider all
					  windows.
--num_permutations NUM_PERMUTATIONS
					  Number of permutations to use for permutation
					  analyses. Default is 10000
--num_bootstraps NUM_BOOTSTRAPS
					  Number of bootstrap replicates to use when
					  bootstrapping mean depth ratios among chromosomes.
					  Default is 10000
--ignore_duplicates   Ignore duplicate reads in bam analyses. Default is to
					  include duplicates.
--marker_size MARKER_SIZE
					  Marker size for genome-wide plots in matplotlib.
					  Default is 10.
--marker_transparency MARKER_TRANSPARENCY, -mt MARKER_TRANSPARENCY
					  Transparency of markers in genome-wide plots. Alpha in
					  matplotlib. Default is 0.5
--coordinate_scale COORDINATE_SCALE
					  For genome-wide scatter plots, divide all coordinates
					  by this value.Default is 1000000, which will plot
					  everything in megabases.
--use_counts          If True, get counts of reads per chromosome for
					  CHROM_STATS, rather than calculating mean depth and
					  mapq. Much faster, but provides less information.
					  Default is False
