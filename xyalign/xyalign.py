# XYalign main program
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import argparse
import collections
import csv
import logging
import os
import subprocess
import sys
import time
import pandas as pd
import xyalign.assemble as assemble
import xyalign.bam as bam
import xyalign.ploidy as ploidy
import xyalign.reftools as reftools
import xyalign.utils as utils
import xyalign.variants as variants
# from xyalign import assemble, bam, ploidy, reftools, utils, variants

# Grab XYalign version from _version.py in the xyalign directory
dir = os.path.dirname(__file__)
version_py = os.path.join(dir, "_version.py")
exec(open(version_py).read())


def parse_args():
	"""
	Parse command-line arguments

	Returns
	-------

	Parser argument namespace
	"""
	parser = argparse.ArgumentParser(description="XYalign")

	parser.add_argument(
		"--bam", nargs="*", help="Full path to input bam files. If more "
		"than one provided, only the first will be used for modules other "
		"than --CHROM_STATS")

	parser.add_argument(
		"--cram", nargs="*", help="Full path to input cram files. If more "
		"than one provided, only the first will be used for modules other "
		"than --CHROM_STATS. Not currently supported.")

	parser.add_argument(
		"--sam", nargs="*", help="Full path to input sam files. If more "
		"than one provided, only the first will be used for modules other "
		"than --CHROM_STATS. Not currently supported.")

	parser.add_argument(
		"--ref", required=True,
		help="REQUIRED. Path to reference sequence (including file name).")

	parser.add_argument(
		"--output_dir", "-o", required=True,
		help="REQUIRED. Output directory. XYalign will create a directory "
		"structure within this directory")

	parser.add_argument(
		"--chromosomes", "-c", nargs="*", default=[None],
		help="Chromosomes to analyze (names must match reference exactly). "
		"For humans, we recommend at least chr19, chrX, chrY.  Generally, we "
		"suggest including the sex chromosomes and at least one autosome. "
		"To analyze all chromosomes use '--chromosomes ALL' or "
		"'--chromosomes all'.")

	parser.add_argument(
		"--x_chromosome", "-x", nargs="*", default=[None],
		help="Names of x-linked scaffolds in reference fasta (must match "
		"reference exactly).  ")

	parser.add_argument(
		"--y_chromosome", "-y", nargs="*", default=[None],
		help="Names of y-linked scaffolds in reference fasta (must match "
		"reference exactly). Defaults to chrY. Give None if using an assembly "
		"without a Y chromosome")

	parser.add_argument(
		"--sample_id", "-id", default="sample",
		help="Name/ID of sample - for use in plot titles and file naming. "
		"Default is sample")

	parser.add_argument(
		"--cpus", type=int, default=1,
		help="Number of cores/threads to use. Default is 1")

	parser.add_argument(
		"--xmx", default="None",
		help="Memory to be provided to java programs via -Xmx.  E.g., use the "
		"flag '--xmx 4g' to pass '-Xmx4g' as a flag when running java "
		"programs (currently just repair.sh). Default is 'None' (i.e., nothing "
		"provided on the command line), which will allow repair.sh to "
		"automatically allocate memory. Note that if you're using --STRIP_READS "
		"on deep coverage whole genome data, you might need quite a bit of memory, "
		"e.g. '--xmx 16g', '--xmx 32g', or more depending on how many reads "
		"are present per read group.")

	parser.add_argument(
		"--fastq_compression", default=3, type=int, choices=range(0, 10),
		help="Compression level for fastqs output from repair.sh. Between "
		"(inclusive) 0 and 9.  Default is 3.  1 through 9 indicate "
		"compression levels. If 0, fastqs will be uncompressed.")

	parser.add_argument(
		"--single_end", action="store_true", default=False,
		help="Include flag if reads are single-end and NOT paired-end.")

	parser.add_argument(
		"--version", "-V", action="version",
		version="XYalign {}".format(__version__),
		help="Print version and exit.")

	parser.add_argument(
		"--no_cleanup", action="store_true", default=False,
		help="Include flag to preserve temporary files.")

	# Options to run specific parts of the pipeline
	pipeline_group = parser.add_mutually_exclusive_group(required=False)
	pipeline_group.add_argument(
		"--PREPARE_REFERENCE", action="store_true", default=False,
		help="This flag will limit XYalign to only preparing reference fastas "
		"for individuals with and without Y chromosomes.  These fastas can "
		"then be passed with each sample to save subsequent processing time.")

	pipeline_group.add_argument(
		"--CHROM_STATS", action="store_true", default=False,
		help="This flag will limit XYalign to only analyzing provided bam files "
		"for depth and mapq across entire chromosomes.")

	pipeline_group.add_argument(
		"--ANALYZE_BAM", action="store_true", default=False,
		help="This flag will limit XYalign to only analyzing the bam file for "
		"depth, mapq, and (optionally) read balance and outputting plots.")

	pipeline_group.add_argument(
		"--CHARACTERIZE_SEX_CHROMS", action="store_true", default=False,
		help="This flag will limit XYalign to the steps required to "
		"characterize sex chromosome content (i.e., analyzing the bam "
		"for depth, mapq, and read balance and running statistical tests "
		"to help infer ploidy)")

	pipeline_group.add_argument(
		"--REMAPPING", action="store_true", default=False,
		help="This flag will limit XYalign to only the steps required to "
		"strip reads and remap to masked references. If masked references "
		"are not provided, they will be created.")

	pipeline_group.add_argument(
		"--STRIP_READS", action="store_true", default=False,
		help="This flag will limit XYalign to only the steps required to "
		"strip reads from a provided bam file.")

	# Logging options
	parser.add_argument(
		"--logfile", default=None,
		help="Name of logfile.  Will overwrite if exists.  Default is "
		"sample_xyalign.log")

	parser.add_argument(
		"--reporting_level", default='INFO',
		choices=["DEBUG", "INFO", "ERROR", "CRITICAL"],
		help="Set level of messages printed to console. Default is 'INFO'. "
		"Choose from (in decreasing amount of reporting) DEBUG, INFO, ERROR "
		"or CRITICAL")

	# Program paths
	parser.add_argument(
		"--platypus_path", default="platypus",
		help="Path to platypus.  Default is 'platypus'.  If platypus is not "
		"directly callable (e.g., '/path/to/platypus' or "
		"'/path/to/Playpus.py'), then provide path to python as well (e.g., "
		"'/path/to/python /path/to/platypus').  In addition, be sure provided "
		"python is version 2.  See the documentation for more information "
		"about setting up an anaconda environment.")

	parser.add_argument(
		"--bwa_path", default="bwa",
		help="Path to bwa. Default is 'bwa'")

	parser.add_argument(
		"--samtools_path", default="samtools",
		help="Path to samtools. Default is 'samtools'")

	parser.add_argument(
		"--repairsh_path", default="repair.sh",
		help="Path to bbmap's repair.sh script. Default is 'repair.sh'")

	parser.add_argument(
		"--shufflesh_path", default="repair.sh",
		help="Path to bbmap's shuffle.sh script. Default is 'shuffle.sh'")

	parser.add_argument(
		"--sambamba_path", default="sambamba",
		help="Path to sambamba. Default is 'sambamba'")

	parser.add_argument(
		"--bedtools_path", default="bedtools",
		help="Path to bedtools. Default is 'bedtools'")

	# Options for turning on/off parts of the pipeline
	parser.add_argument(
		"--platypus_calling", default="both",
		choices=["both", "none", "before", "after"],
		help="Platypus calling withing the pipeline "
		"(before processing, after processing, both, "
		"or neither). Options: both, none, before, after.")

	parser.add_argument(
		"--no_variant_plots", action="store_true", default=False,
		help="Include flag to prevent plotting read balance from VCF files.")

	parser.add_argument(
		"--no_bam_analysis", action="store_true", default=False,
		help="Include flag to prevent depth/mapq analysis of bam file. "
		"Used to isolate platypus_calling.")

	parser.add_argument(
		"--skip_compatibility_check", action="store_true", default=False,
		help="Include flag to prevent check of compatibility between input bam "
		"and reference fasta")

	parser.add_argument(
		"--no_perm_test", action="store_true", default=False,
		help="Include flag to turn off permutation tests.")

	parser.add_argument(
		"--no_ks_test", action="store_true", default=False,
		help="Include flag to turn off KS Two Sample tests.")

	parser.add_argument(
		"--no_bootstrap", action="store_true", default=False,
		help="Include flag to turn off bootstrap analyses. Requires either "
		"--y_present, --y_absent, or --sex_chrom_calling_threshold "
		"if running full pipeline.")

	# Variant Flags
	parser.add_argument(
		"--variant_site_quality", "-vsq", type=int, default=30,
		help="Consider all SNPs with a site quality (QUAL) greater than or "
		"equal to this value. Default is 30.")

	parser.add_argument(
		"--variant_genotype_quality", "-vgq", type=int, default=30,
		help="Consider all SNPs with a sample genotype quality greater than or "
		"equal to this value. Default is 30.")

	parser.add_argument(
		"--variant_depth", "-vd", type=int, default=4,
		help="Consider all SNPs with a sample depth greater than or "
		"equal to this value. Default is 4.")

	parser.add_argument(
		"--platypus_logfile", default=None,
		help="Prefix to use for Platypus log files.  Will default to the "
		"sample_id argument provided")

	parser.add_argument(
		"--homogenize_read_balance", type=bool, default=False,
		help="If True, read balance values will be transformed by subtracting "
		"each value from 1. For example, 0.25 and 0.75 would be treated "
		"equivalently. Default is False.")

	parser.add_argument(
		"--min_variant_count", type=int, default=1,
		help="Minimum number of variants in a window for the read balance of "
		"that window to be plotted. Note that this does not affect plotting "
		"of variant counts. Default is 1, though we note that many window "
		"averages will be meaningless at this setting.")

	# Reference-related Flags
	parser.add_argument(
		"--reference_mask", nargs="*",
		help="Bed file containing regions to replace with Ns in the sex "
		"chromosome reference.  Examples might include the pseudoautosomal "
		"regions on the Y to force all mapping/calling on those regions of the "
		"X chromosome.  Default is None.")

	parser.add_argument(
		"--xx_ref_out_name", default=None,
		help="Desired name for masked output fasta for "
		"samples WITHOUT a Y chromosome (e.g., XX, XXX, XO, etc.). "
		"Defaults to 'xyalign_noY.masked.fa'. Will be output "
		"in the XYalign reference directory.")

	parser.add_argument(
		"--xy_ref_out_name", default=None,
		help="Desired name for masked output fasta for "
		"samples WITH a Y chromosome (e.g., XY, XXY, etc.). "
		"Defaults to 'xyalign_withY.masked.fa'. Will be output "
		"in the XYalign reference directory.")

	parser.add_argument(
		"--xx_ref_out", default=None,
		help="Desired path to and name of masked output fasta for "
		"samples WITHOUT a Y chromosome (e.g., XX, XXX, XO, etc.). "
		"Overwrites if exists. "
		"Use if you would like output somewhere other than XYalign reference "
		"directory. Otherwise, use --xx_ref_name.")

	parser.add_argument(
		"--xy_ref_out", default=None,
		help="Desired path to and name of masked output fasta for "
		"samples WITH a Y chromosome (e.g., XY, XXY, etc.). Overwrites if exists. "
		"Use if you would like output somewhere other than XYalign reference "
		"directory. Otherwise, use --xy_ref_name.")

	parser.add_argument(
		"--xx_ref_in", default=None,
		help="Path to preprocessed reference fasta to be used for remapping in "
		"X0 or XX samples.  Default is None.  If none, will produce a "
		"sample-specific reference for remapping.")

	parser.add_argument(
		"--xy_ref_in", default=None,
		help="Path to preprocessed reference fasta to be used for remapping in "
		"samples containing Y chromosome.  Default is None.  If none, will "
		"produce a sample-specific reference for remapping.")

	parser.add_argument(
		"--bwa_index", type=bool, default=False,
		help="If True, index with BWA during PREPARE_REFERENCE. Only relevant"
		"when running the PREPARE_REFERENCE module by itself. Default is False.")

	# Mapping/remapping arguments
	parser.add_argument(
		"--read_group_id", default="xyalign", type=str,
		help="If read groups are present in a bam file, they are used by default in "
		"remapping steps.  However, if read groups are not present in a file, "
		"there are two options for proceeding. If '--read_group_id None' is "
		"provided (case sensitive), then no read groups will be used in "
		"subsequent mapping steps. Otherwise, any other string provided to "
		"this flag will be used as a read group ID.  "
		"Default is '--read_group_id xyalign'")

	parser.add_argument(
		"--bwa_flags", type=str, default="",
		help="Provide a string (in quotes, with spaces between arguments) "
		"for additional flags desired for BWA mapping (other than -R and -t). "
		"Example: '-M -T 20 -v 4'.  Note that those are spaces between "
		"arguments.")

	parser.add_argument(
		"--sex_chrom_bam_only", action="store_true", default=False,
		help="This flag skips merging the new sex chromosome bam file back "
		"into the original bam file (i.e., sex chrom swapping). This will "
		"output a bam file containing only the newly remapped sex chromosomes.")

	parser.add_argument(
		"--sex_chrom_calling_threshold", type=float, default=2.0,
		help="This is the *maximum* filtered X/Y depth ratio for an individual "
		"to be considered as having heterogametic sex chromsomes (e.g., XY) for "
		"the REMAPPING module of XYalign. "
		"Note here that X and Y chromosomes are simply the chromosomes that "
		"have been designated as X and Y via --x_chromosome and --y_chromosome. "
		"Keep in mind that the ideal threshold will vary according to sex "
		"determination mechanism, sequence homology between the sex chromosomes, "
		"reference genome, sequencing methods, etc.  See documentation for more "
		"detail. Default is 2.0, which we found to be reasonable for exome, "
		"low-coverage whole-genome, and high-coverage whole-genome human data.")

	group_y = parser.add_mutually_exclusive_group(required=False)

	group_y.add_argument(
		"--y_present", action="store_true", default=False,
		help="Overrides sex chr estimation by XYalign and remaps with Y present.")

	group_y.add_argument(
		"--y_absent", action="store_true", default=False,
		help="Overrides sex chr estimation by XY align and remaps with Y absent.")

	# Bam Analysis and Characterize Sex Chrom Flags
	parser.add_argument(
		"--window_size", "-w", default=None,
		help="Window size (integer) for sliding window calculations. Default "
		"is 50000.  Default is None. If set to None, will use targets "
		"provided using --target_bed.")

	parser.add_argument(
		"--target_bed", default=None,
		help="Bed file containing targets to use in sliding window analyses "
		"instead of a fixed window width. Either this or --window_size needs "
		"to be set.  Default is None, which will use window size provided "
		"with --window_size.  If not None, and --window_size is None, analyses "
		"will use targets in provided file.  Must be typical bed format, "
		"0-based indexing, with the first three columns containing "
		"the chromosome name, start coordinate, stop coordinate.")

	parser.add_argument(
		"--exact_depth", action="store_true", default=False,
		help="Calculate exact depth within windows, else use much faster "
		"approximation. *Currently exact is not implemented*. Default is False.")

	parser.add_argument(
		"--whole_genome_threshold", action="store_true", default=False,
		help="This flag will calculate the depth filter threshold based on "
		"all values from across the genome.  By default, thresholds are "
		"calculated per chromosome.")

	parser.add_argument(
		"--mapq_cutoff", "-mq", type=int, default=20,
		help="Minimum mean mapq threshold for a window to be "
		"considered high quality. Default is 20.")

	parser.add_argument(
		"--min_depth_filter", type=float, default=0.0,
		help="Minimum depth threshold for a window to be considered high "
		"quality. Calculated as mean depth * min_depth_filter. So, a "
		"min_depth_filter of 0.2 would require at least a minimum depth "
		"of 2 if the mean depth was 10. Default is 0.0 to consider all windows.")

	parser.add_argument(
		"--max_depth_filter", type=float, default=10000.0,
		help="Maximum depth threshold for a window to be considered high "
		"quality. Calculated as mean depth * max_depth_filter. So, a "
		"max_depth_filter of 4 would require depths to be less than or "
		"equal to 40 if the mean depth was 10. "
		"Default is 10000.0 to consider all windows.")

	parser.add_argument(
		"--num_permutations", type=int, default=10000,
		help="Number of permutations to use for permutation analyses. "
		"Default is 10000")

	parser.add_argument(
		"--num_bootstraps", type=int, default=10000,
		help="Number of bootstrap replicates to use when bootstrapping mean "
		"depth ratios among chromosomes. Default is 10000")

	parser.add_argument(
		"--ignore_duplicates", action="store_true", default=False,
		help="Ignore duplicate reads in bam analyses. Default is to include "
		"duplicates.")

	# Plotting flags
	parser.add_argument(
		"--marker_size", type=float, default=10.0,
		help="Marker size for genome-wide plots in matplotlib. Default is 10.")

	parser.add_argument(
		"--marker_transparency", "-mt", type=float, default=0.5,
		help="Transparency of markers in genome-wide plots.  "
		"Alpha in matplotlib.  Default is 0.5")

	parser.add_argument(
		"--coordinate_scale", type=int, default=1000000,
		help="For genome-wide scatter plots, divide all coordinates by this value."
		"Default is 1000000, which will plot everything in megabases.")

	parser.add_argument(
		"--include_fixed", type=bool, default=False,
		help="Default is False, which removes read balances less than or equal to "
		"0.05 and equal to 1.0 for histogram plotting. True will include "
		"all values. Extreme values removed by default because they often swamp "
		"out the signal of the rest of the distribution.")

	# CHROM_STATS flags
	parser.add_argument(
		"--use_counts", action="store_true", default=False,
		help="If True, get counts of reads per chromosome for CHROM_STATS, rather "
		"than calculating mean depth and mapq. Much faster, but provides less "
		"information. Default is False")

	args = parser.parse_args()

	# Validate Arguments

	mods = [
		args.PREPARE_REFERENCE, args.ANALYZE_BAM, args.CHARACTERIZE_SEX_CHROMS,
		args.REMAPPING, args.STRIP_READS, args.CHROM_STATS]
	if not any(mods) is True:
		full_pipeline = True
	else:
		full_pipeline = False

	# Ensure only one of bam, sam, or cram is set

	bam_flags = [args.bam, args.cram, args.sam]
	bam_flags = [x for x in bam_flags if x is not None]
	if len(bam_flags) > 1:
		sys.exit(
			"Error. --bam, --cram, and --sam are mutally exclusive, please "
			"provide no more than one. Submitted values are: --bam {} "
			"--cram {} --sam {}".format(args.bam, args.cram, args.sam))
	elif len(bam_flags) == 0:
		if args.PREPARE_REFERENCE is not True:
			sys.exit(
				"Error. All modules other than PREPARE_REFERENCE require an input "
				"bam file.")

	# Cram and SAM files are currently unsupported
	if args.cram is not None:
		sys.exit(
			"Error. XYalign does not currently support cram files. "
			"Please provide a bam file via --bam instead.")

	if args.sam is not None:
		sys.exit(
			"Error. XYalign does not currently support sam files. "
			"Please provide a bam file via --bam instead.")

	# Exact depth is not currently implemented
	if args.exact_depth is True:
		sys.exit(
			"Error. Exact depth is not currently implemented. Exiting.")

	# Validate ploidy test arguments
	test_flags = [args.no_perm_test, args.no_bootstrap, args.no_ks_test]
	if full_pipeline is True:
		if args.no_bootstrap is True:
			if args.y_present is False and args.y_absent is False:
				sys.exit(
					"Error. Either --y_present or --y_absent needs to be "
					"included with if running the full pipeline "
					"with --no_bootstrap. XYalign uses the results of the "
					"bootstrap for sex chromosome complement inference.")
	if args.REMAPPING is True:
		if any([args.y_present, args.y_absent]) is False:
				sys.exit(
					"Error. Either --y_present or --y_absent needs to be "
					"included with if runnning --REMAPPING.  If you are "
					"interested in both sex chromosome complement inference "
					"and remapping based on the results of said inference, "
					"consider running the full pipeline.")

	# Validate chromosome arguments
	if any(
		[full_pipeline, args.ANALYZE_BAM, args.CHARACTERIZE_SEX_CHROMS,
			args.STRIP_READS, args.CHROM_STATS]) is True:
		if args.chromosomes == [None]:
			sys.exit("Please provide chromosome names to analyze (--chromosomes)")
	if any(
		[full_pipeline, args.CHARACTERIZE_SEX_CHROMS]) is True:
		if len(args.chromosomes) < 2:
			sys.exit(
				"You only provided a single chromosome to analyze with "
				"--chromosomes.  CHARACTERIZE_SEX_CHROMS requires at least two "
				"chromosomes, including one designated as an X chromosome.")
	if any(
		[full_pipeline, args.CHARACTERIZE_SEX_CHROMS, args.REMAPPING]) is True:
		if args.x_chromosome == [None]:
			sys.exit(
				"Please provide an 'X chromosome' to analyze (--x_chromosome). "
				"Note that this does not have to actually be an x-linked sequence, "
				"but this designation is required for the tests in "
				"CHARACTERIZE_SEX_CHROMS, which involve pairwise comparisons "
				"between the 'X chromosome', the 'Y chromosome', and all other "
				"'autosomes' (see documentation for more details).")
	if any(
		[full_pipeline, args.PREPARE_REFERENCE, args.REMAPPING]) is True:
		if args.y_chromosome == [None]:
			sys.exit(
				"Please provide an 'Y chromosome' for the full pipeline, "
				"--PREPARE_REFERENCE, or --REMAPPING.  It is required for the "
				"creation and use of an alternate reference (e.g., "
				"masking the Y chromosome for the XX reference)")

	# Validate bwa arguments
	bwa_args = [str(x).strip() for x in args.bwa_flags.split()]
	red_list = ["-rm", "rm", "-rf", "rf", "-RM", "RM", "-RF", "RF"]
	if any(x in bwa_args for x in red_list):
		sys.exit(
			"Found either rm or rf in your bwa flags. Exiting to prevent "
			"unintended shell consequences")
	yellow_list = ["-R", "-t"]
	if any(x in bwa_args for x in yellow_list):
		sys.exit(
			"Found either -R or -t in bwa flags.  These flags are already used "
			"in XYalign.  Please remove.")

	# Validate sliding window options
	if any(
		[full_pipeline, args.ANALYZE_BAM, args.CHARACTERIZE_SEX_CHROMS]) is True:
		if args.window_size is not None and args.window_size != "None":
			if args.window_size.isdigit() is False:
				sys.exit(
					"--window_size needs to be either None or a positive integer. "
					"Exiting.")
		else:
			if args.target_bed is None:
				sys.exit(
					"If --window_size is None, --target_bed needs to be used. "
					"Please set one of these two flags if running ANALYZE_BAM, "
					"CHARACTERIZE_SEX_CHROMS, or the full pipeline. Exiting.")
			elif os.path.exists(args.target_bed) is False:
				sys.exit(
					"Invalid file provided with --target_bed. Check path. Exiting.")

	# Return arguments namespace
	return args


def ref_prep(
	ref_obj, ref_mask, ref_dir, xx, xy, y_chromosome,
	samtools_path, bwa_path, bwa_index):
	"""
	Reference prep part of XYalign pipeline.

	* Creates two reference fasta files.  Both will include masks provied with
	ref_mask.  One will additionally have the entire Y chromosome
	hard masked.

	* Indexes (.fai, .dict, and optionally bwa indices) both new references

	Parameters
	----------
	ref_obj : RefFasta() object
		A reftools.RefFasta() object of a fasta reference file to be processed
	ref_mask : list or None
		List of files to use to *hard-mask* references. None will ignore masking.
	ref_dir : str
		Path to output directory
	xx : str
		Path to XX output reference
	xy : str
		Path to XY output reference
	y_chromosome : str
		Name of Y chromosome in fasta
	samtools_path : str
		The path to samtools (i.e, "samtools" if in path)
	bwa_path : str
		The path to bwa (i.e, "bwa" if in path)
	bwa_index : bool
		If True, create bwa indices. Don't if False.

	Returns
	-------

	tuple
		Paths to two masked references (y_masked, y_unmasked)
	"""
	logger = logging.getLogger("xyalign")
	# Combine masks, if more than one present
	if ref_mask is not None:
		if len(ref_mask) > 1:
			reference_mask = utils.merge_bed_files(
				"{}/reference_mask.merged.bed".format(
					ref_dir), *ref_mask)
		else:
			reference_mask = ref_mask[0]
	else:
		reference_mask = None
	# Create masked noY reference
	y_mask = ref_obj.chromosome_bed("{}/Y.bed".format(
		ref_dir), y_chromosome)
	if reference_mask is not None:
		noy_out = ref_obj.mask_reference(
			utils.merge_bed_files(
				"{}/reference_mask.maskY.merged.bed".format(
					ref_dir), reference_mask, y_mask), xx)
	else:
		noy_out = ref_obj.mask_reference(y_mask, xx)
	noy_ref = reftools.RefFasta(noy_out, samtools_path, bwa_path)
	noy_ref.seq_dict()
	# Create masked withY reference
	if reference_mask is not None:
		withy_out = ref_obj.mask_reference(reference_mask, xy)
	else:
		logger.info(
			"No reference mask provided, copying full reference to {} as "
			"XY reference to prevent damage, modification, etc. to original "
			"reference.".format(xy))
		subprocess.call(["cp", ref_obj.filepath, xy])
		withy_out = xy
	withy_ref = reftools.RefFasta(withy_out, samtools_path, bwa_path)
	withy_ref.seq_dict()
	if bwa_index is True:
		noy_ref.index_bwa()
		withy_ref.index_bwa()
	return (noy_ref, withy_ref)


def chrom_stats(bam_obj_list, chrom_list, use_counts):
	"""
	Runs chrom stats module.

	Calculates mean depth and mapq across entire scaffolds for a list of bam
	files

	Input
	-----

	bam_obj_list : list
		List of bam.BamFile() objects. Need to have been created using same
		reference (i.e., seq names and lengths are the same)
	chrom_list : list
		List of chromosome names to analyze
	use_counts : bool
		If True, will just grab counts from bam index INSTEAD of traversing
		for depth and mapq

	Returns
	-------

	tuple
		Tuple containing two dictionaries with results for
		depth and mapq, respectively. Or, if use_counts is True, returns a
		tuple containing the count dictionary and None.

	"""
	logger = logging.getLogger("xyalign")
	logger.info(
		"Running CHROM_STATS on following bam files: {}".format(
			", ".join([x.filepath for x in bam_obj_list])))

	comp_check = utils.check_compatibility_bam_list(bam_obj_list)
	if comp_check is False:
		logger.error(
			"Bam files incompatible - ensure same reference used for all files "
			"and that sequence portion of all bam headers match. Exiting")
		logging.shutdown()
		sys.exit(1)

	if chrom_list == ["ALL"] or chrom_list == ["all"]:
		chrom_list = list(bam_obj_list[0].chromosome_names())

	missing_chroms = bam_obj_list[0].check_chrom_in_bam(chrom_list)
	if len(missing_chroms) != 0:
		logger.error(
			"One or more chromosomes provided via --chromosomes not "
			"present in bam files. Exiting.")
		logging.shutdown()
		sys.exit(1)

	if use_counts is False:
		chrom_depth_dict = collections.OrderedDict()
		chrom_mapq_dict = collections.OrderedDict()

		chrom_depth_dict["header"] = [
			os.path.basename(x.filepath) for x in bam_obj_list]
		chrom_depth_dict["header"].insert(0, "chrom")
		chrom_mapq_dict["header"] = [
			os.path.basename(x.filepath) for x in bam_obj_list]
		chrom_mapq_dict["header"].insert(0, "chrom")
		for chromosome in chrom_list:
			chrom_results = [
				x.chrom_stats(chromosome, False) for x in bam_obj_list]
			chrom_depth_dict[chromosome], chrom_mapq_dict[chromosome] = zip(
				*chrom_results)

			chrom_depth_dict[chromosome] = (chromosome,) + chrom_depth_dict[chromosome]
			chrom_mapq_dict[chromosome] = (chromosome,) + chrom_mapq_dict[chromosome]

		return (chrom_depth_dict, chrom_mapq_dict)

	else:
		chrom_count_dict = collections.OrderedDict()

		chrom_count_dict["header"] = [
			os.path.basename(x.filepath) for x in bam_obj_list]
		chrom_count_dict["header"].insert(0, "chrom")

		for chromosome in chrom_list:
			chrom_count_dict[chromosome] = [chromosome]

		for i in bam_obj_list:
			idx_stats = i.chrom_counts()
			for k in idx_stats:
				if k.contig in chrom_count_dict:
					chrom_count_dict[k.contig].append(int(k.mapped))

		return (chrom_count_dict, None)


def bam_analysis(
	input_bam_obj, platypus_calling, platypus_path, vcf_log, ref_obj,
	input_chroms, cpus, out_vcf, no_variant_plots, window_size, target_bed,
	sample_id, readbalance_prefix, variant_site_quality,
	variant_genotype_quality, variant_depth, marker_size, marker_transparency,
	homogenize_read_balance, data_frame_readbalance, min_variant_count,
	no_bam_analysis, ignore_duplicates, exact_depth, whole_genome_threshold,
	mapq_cutoff, min_depth_filter, max_depth_filter, depth_mapq_prefix,
	bam_data_frame, output_bed_high, output_bed_low, use_bed_for_platypus,
	coordinate_scale, fixed):
	"""
	Runs bam analyis part of XYalign pipeline on bam file.

	* (Optionally) calls variants using Platypus

	* (Optionally) parses and filters Platypus vcf, and plots read balance

	* (Optionally) Calculates window based metrics from the bam file:
		depth and mapq

	* (optionally) Plots window-based metrics

	* Outputs two bed files: high quality windows, and low quality windows.

	Parameters
	----------
	input_bam_obj : bam.BamFile() object
	platypus_calling : bool
		If True, will call and analyze variants
	platypus_path : str
		Command to call platypus (e.g, "platypus")
	vcf_log : str
		Path to file for platypus log
	ref_obj : reftools.RefFasta() object
	input_chroms : list
		Chromosomes to analyze
	cpus : int
		Number of threads/cpus
	out_vcf : str
		Output vcf path/name
	no_variant_plots : bool
		If True, will not plot read balance
	window_size : int or None
		Window size for sliding window analyses (both bam and vcf). If None,
		will use regions in target_bed
	target_bed : str or None
		Path to bed file containing targets to use in sliding window analyses
	sample_id : str
	readbalance_prefix : str
		Prefix, including full path, to use for output files for
		readbalance analyses
	variant_site_quality : int
		Minimum site quality (PHRED) for a site to be included in
		readbalance analyses
	variant_genotype_quality : int
		Minimum genotype quality for a site to be included in read balance analyses
	variant_depth : int
		Minimum depth for a site to be included in read balance analyses
	marker_size : float
		Marker size for plotting genome scatter plots
	marker_transparency: float
		Value to use for marker transparency in genome scatter plots
	homogenize_read_balance : bool
		If true, will subtract values less than 0.5 from 1. I.e., 0.25 and 0.75
		would be treated equivalently
	data_frame_readbalance: str
		Path of output file for full read balance dataframe
	min_variant_count : int
		Minimum number of variants in a given window for the window to be
		plotted in window-based read balance analyses
	no_bam_analysis : bool
		If True, no bam analyses will take place
	ignore_duplicates : bool
		If True, duplicates excluded from bam analyses
	exact_depth : bool
		If True, exact depth calculated in each window. Else, a much faster
		approximation will be used
	whole_genome_threshold : bool
		If True, values for depth filters will be calculated using mean from
		across all chromosomes included in analyses. Else, mean will be taken
		per chromosome
	min_depth_filter : float
		Minimum depth threshold for a window to be considered high. Calculated
		as mean depth * min_depth_filter.
	max_depth_filter : float
		Maximum depth threshold for a window to be considered high. Calculated
		as mean depth * min_depth_filter.
	depth_mapq_prefix : str
		Prefix, including full path, to be used for files output from depth and
		mapq analyses
	bam_data_frame : str
		Full path to output file for dataframe containing all data from bam
		analyses
	output_bed_high : str
		Full path to output bed containing high quality (i.e., passing
		filters) windows
	output_bed_low : str
		Full path to output bed containing low quality (i.e., failing
		filters) windows
	use_bed_for_platypus : bool
		If True, use output_bed_high as regions for Platypus calling
	coordinate_scale : int
		Divide all coordinates by this value for plotting. In most cases, 1000000
		will be ideal for eukaryotic genomes.
	fixed : bool
		If False, only plots histogram for values between 0.05 and 1.0
		(non-inclusive). If True, plots histogram of all variants.

	Returns
	-------
	tuple
		(list of pandas dataframes with passing windows,
		list of pandas dataframes with failing windows)
	"""
	logger = logging.getLogger("xyalign")
	# Depth/MAPQ
	if no_bam_analysis is not True:
		all_df = []
		pass_df = []
		fail_df = []
		for chromosome in input_chroms:
			if window_size is not None and window_size != "None":
				data = input_bam_obj.analyze_bam(
					chromosome, ignore_duplicates,
					exact_depth, int(window_size))
			else:
				data = input_bam_obj.analyze_bam(
					chromosome, ignore_duplicates,
					exact_depth, None, target_bed)
			if whole_genome_threshold is True:
				tup = utils.make_region_lists_genome_filters(
					data, mapq_cutoff,
					min_depth_filter, max_depth_filter)
			else:
				tup = utils.make_region_lists_chromosome_filters(
					data, mapq_cutoff,
					min_depth_filter, max_depth_filter)
			all_df.append(data)
			pass_df.append(tup[0])
			fail_df.append(tup[1])
			utils.plot_depth_mapq(
				data, depth_mapq_prefix, sample_id,
				input_bam_obj.get_chrom_length(chromosome), marker_size,
				marker_transparency, coordinate_scale)
		all_concat = pd.concat(all_df)
		all_concat.to_csv(
			bam_data_frame, index=False, sep="\t", quoting=csv.QUOTE_NONE)
		utils.output_bed(output_bed_high, *pass_df)
		utils.output_bed(output_bed_low, *fail_df)

	# Platypus
	if platypus_calling is True:
		if use_bed_for_platypus is True:
			a = input_bam_obj.platypus_caller(
				platypus_path, vcf_log,
				ref_obj.filepath, input_chroms, cpus, out_vcf,
				output_bed_high)
		else:
			a = input_bam_obj.platypus_caller(
				platypus_path, vcf_log,
				ref_obj.filepath, input_chroms, cpus, out_vcf,
				None)
		if a != 0:
			logger.error("Error in platypus calling on {}".format(
				input_bam_obj.filepath))
			logging.shutdown()
			sys.exit(1)
		noprocess_vcf_object = variants.VCFFile(out_vcf)
		if no_variant_plots is not True:
			if window_size is not None and window_size != "None":
				noprocess_vcf_object.plot_variants_per_chrom(
					chrom_list=input_chroms,
					sampleID=sample_id,
					output_prefix=readbalance_prefix,
					site_qual=variant_site_quality,
					genotype_qual=variant_genotype_quality,
					depth=variant_depth,
					MarkerSize=marker_size,
					MarkerAlpha=marker_transparency,
					bamfile_obj=input_bam_obj,
					variant_caller="platypus",
					homogenize=homogenize_read_balance,
					dataframe_out=data_frame_readbalance,
					min_count=min_variant_count,
					window_size=int(window_size),
					x_scale=coordinate_scale,
					target_file=None,
					include_fixed=fixed)
			else:
				noprocess_vcf_object.plot_variants_per_chrom(
					chrom_list=input_chroms,
					sampleID=sample_id,
					output_prefix=readbalance_prefix,
					site_qual=variant_site_quality,
					genotype_qual=variant_genotype_quality,
					depth=variant_depth,
					MarkerSize=marker_size,
					MarkerAlpha=marker_transparency,
					bamfile_obj=input_bam_obj,
					variant_caller="platypus",
					homogenize=homogenize_read_balance,
					dataframe_out=data_frame_readbalance,
					min_count=min_variant_count,
					window_size=None,
					x_scale=coordinate_scale,
					target_file=target_bed,
					include_fixed=fixed)
	return(pass_df, fail_df)


def ploidy_analysis(
	passing_df, failing_df, no_perm_test, no_ks_test, no_bootstrap,
	input_chroms, x_chromosome, y_chromosome, results_dir,
	num_permutations, num_bootstraps, sample_id):
	"""
	Runs the ploidy analysis part of XYalign.

	* Runs permutation test to systematically compare means between
	every possible pair of chromosomes

	* Runs K-S two sample test to systematically compare distributions between
	every possible pair of chromosomes

	* Bootstraps the mean depth ratio for every possible pair of chromosomes

	Parameters
	----------

	passing_df : list
		Passing pandas dataframes, one per chromosome
	failing_df : list
		Failing pandas dataframes, one per chromosome
	no_perm_test : bool
		If False, permutation test will be run
	no_ks_test : bool
		If False, KS test will be run
	no_bootstrap : bool
		If False, bootstrap analysis will be run
	input_chroms : list
		Chromosomes/scaffolds to analyze
	x_chromosome : list
		X-linked scaffolds
	y_chromosome : list
		Y-likned scaffolds
	results_dir : str
		Full path to directory to output results
	num_permutations : int
		Number of permutations
	num_bootstraps : int
		Number of bootstrap replicates
	sample_id : str

	Returns
	-------

	dictionary
		Results for each test. Keys: perm, ks, boot.
	"""
	# Permutations
	if no_perm_test is False:
		if y_chromosome != [None]:
			sex_chromosomes = x_chromosome + y_chromosome
			perm_res_x = []
			perm_res_y = []
		else:
			sex_chromosomes = x_chromosome
			perm_res_x = []
			perm_res_y = None
		autosomes = [
			x for x in input_chroms if x not in sex_chromosomes]
		for c in autosomes:
			perm_res_x.append(ploidy.permutation_test_chromosomes(
				pd.concat(passing_df), c,
				str(x_chromosome[0]), "chrom",
				"depth", num_permutations,
				results_dir + "/{}.{}_{}_permutation_results.txt".format(
					sample_id, c, str(x_chromosome[0]))))
			if perm_res_y is not None:
				perm_res_y.append(ploidy.permutation_test_chromosomes(
					pd.concat(passing_df), c,
					str(y_chromosome[0]), "chrom",
					"depth", num_permutations,
					results_dir + "/{}.{}_{}_permutation_results.txt".format(
						sample_id, c, str(y_chromosome[0]))))
		if perm_res_y is not None:
			sex_perm_res = ploidy.permutation_test_chromosomes(
				pd.concat(passing_df), str(x_chromosome[0]),
				str(y_chromosome[0]),
				"chrom", "depth", num_permutations,
				results_dir + "/{}.{}_{}_permutation_results.txt".format(
					sample_id, str(x_chromosome[0]), str(y_chromosome[0])))
		else:
			sex_perm_res = None

	# K-S Two Sample
	if no_ks_test is False:
		if y_chromosome != [None]:
			sex_chromosomes = x_chromosome + y_chromosome
			ks_res_x = []
			ks_res_y = []
		else:
			sex_chromosomes = x_chromosome
			ks_res_x = []
			ks_res_y = None
		autosomes = [
			x for x in input_chroms if x not in sex_chromosomes]
		for c in autosomes:
			ks_res_x.append(ploidy.ks_two_sample(
				pd.concat(passing_df), c,
				str(x_chromosome[0]), "chrom", "depth",
				results_dir + "/{}.{}_{}_ks_results.txt".format(
					sample_id, c, str(x_chromosome[0]))))
			if ks_res_y is not None:
				ks_res_y.append(ploidy.ks_two_sample(
					pd.concat(passing_df), c,
					str(y_chromosome[0]), "chrom", "depth",
					results_dir + "/{}.{}_{}_ks_results.txt".format(
						sample_id, c, str(y_chromosome[0]))))
		if ks_res_y is not None:
			sex_ks_res = ploidy.ks_two_sample(
				pd.concat(passing_df), str(x_chromosome[0]),
				str(y_chromosome[0]), "chrom", "depth",
				results_dir + "/{}.{}_{}_ks_results.txt".format(
					sample_id, str(x_chromosome[0]), str(y_chromosome[0])))
		else:
			sex_ks_res = None

	# Bootstrap
	if no_bootstrap is False:
		if y_chromosome != [None]:
			sex_chromosomes = x_chromosome + y_chromosome
			boot_res_x = []
			boot_res_y = []
		else:
			sex_chromosomes = x_chromosome
			boot_res_x = []
			boot_res_y = None
		autosomes = [
			x for x in input_chroms if x not in sex_chromosomes]
		for c in autosomes:
			boot_res_x.append(ploidy.bootstrap(
				pd.concat(passing_df), c,
				str(x_chromosome[0]), "chrom",
				"depth", num_bootstraps,
				results_dir + "/{}.{}_{}_bootstrap_results.txt".format(
					sample_id, c, str(x_chromosome[0]))))
			if boot_res_y is not None:
				boot_res_y.append(ploidy.bootstrap(
					pd.concat(passing_df), c,
					str(y_chromosome[0]), "chrom",
					"depth", num_bootstraps,
					results_dir + "/{}.{}_{}_bootstrap_results.txt".format(
						sample_id, c, str(y_chromosome[0]))))
		if boot_res_y is not None:
			sex_boot_res = ploidy.bootstrap(
				pd.concat(passing_df), str(x_chromosome[0]),
				str(y_chromosome[0]),
				"chrom", "depth", num_bootstraps,
				results_dir + "/{}.{}_{}_bootstrap_results.txt".format(
					sample_id, str(x_chromosome[0]), str(y_chromosome[0])))
		else:
			sex_boot_res = None
	return {
		"perm": [perm_res_x, perm_res_y, sex_perm_res],
		"ks": [ks_res_x, ks_res_y, sex_ks_res],
		"boot": [boot_res_x, boot_res_y, sex_boot_res]}


def remapping(
	input_bam_obj, y_pres, masked_references, samtools_path, sambamba_path,
	repairsh_path, shufflesh_path, bwa_path, bwa_flags, single_end, bam_dir,
	fastq_dir, sample_id, x_chromosome, y_chromosome, cpus,
	xmx, fastq_compression, cleanup, read_group_id):
	"""
	Runs remapping steps of XYalign.

	* Strips, sorts, and re-pair reads from the sex chromosomes (collecting read
	group information)

	* Maps (with sorting) reads (with read group information) to appropriate
	reference based on presence (or not) of Y chromosome

	* Merge bam files (if more than one read group)

	Parameters
	----------

	input_bam_obj : bam.BamFile() object
	y_pres : bool
		True if Y chromosome present in individual
	masked_references : tuple
		Masked reference objects (xx, xy)
	samtools_path : str
		Path/command to call samtools
	sambamba_path : str
		Path/command to call sambamba
	repairsh_path : str
		Path/command to call repair.sh
	shufflesh_path : str
		Path/command to call shuffle.sh
	bwa_path : str
		Path/command to call bwa
	bwa_flags : str
		Flags to use for bwa mapping
	single_end : bool
		If True, reads treated as single end
	bam_dir : str
		Path to output directory for bam files
	fastq_dir : str
		Path to output directory for fastq files
	sample_id : str
	x_chromosome : list
		X-linked scaffolds
	y_chromosome : list
		Y-linked scaffolds
	cpus : int
		Number of threads/cpus
	xmx : str
		Value to be combined with -Xmx for java programs (i.e., 4g would
		result in -Xmx4g)
	fastq_compression : int
		Compression level for fastq files. 0 leaves fastq files uncompressed.
		Otherwise values should be between 1 and 9 (inclusive), with
		larger values indicating more compression
	cleanup : bool
		If True, will delete temporary files
	read_group_id : str
		ID to use to add read group information

	Returns
	-------

	str
		Path to bam containing remapped sex chromsomes
	"""
	if y_pres is True:
		new_reference = masked_references[1]
	else:
		new_reference = masked_references[0]
	new_fastqs = input_bam_obj.strip_reads(
		repairsh_path, shufflesh_path, single_end, fastq_dir,
		sample_id, x_chromosome + y_chromosome, xmx,
		fastq_compression, cleanup, read_group_id)
	with open(new_fastqs[0]) as f:
		read_group_and_fastqs = [line.strip() for line in f]
		read_group_and_fastqs = [x.split() for x in read_group_and_fastqs]
	if new_fastqs[1] is not None:
		with open(new_fastqs[1]) as f:
			read_group_headers = [line.split() for line in f]
		temp_bam_list = []
		for i in read_group_and_fastqs:
			if i != [""]:
				rg_id = i[0]
				fastq_files = i[1:]
				for j in read_group_headers:
					for k in j:
						if k[0:2] == 'ID':
							if k[3:] == rg_id:
								rg_tag = "\t".join(j)
							break
				temp_bam = assemble.bwa_mem_mapping_sambamba(
					bwa_path, samtools_path, sambamba_path,
					new_reference, "{}/{}.sex_chroms.{}".format(
						bam_dir, sample_id, rg_id),
					fastq_files, cpus, rg_tag,
					[str(x).strip() for x in bwa_flags.split()])
				temp_bam_list.append(temp_bam)
		if len(temp_bam_list) < 2:
			new_bam = temp_bam_list[0]
		else:
			new_bam = bam.samtools_merge(
				samtools_path, temp_bam_list, "{}/{}.sex_chroms".format(
					bam_dir, sample_id), cpus)
	else:
		fastq_files = read_group_and_fastqs[0][1:]
		new_bam = assemble.bwa_mem_mapping_sambamba(
			bwa_path, samtools_path, sambamba_path,
			new_reference, "{}/{}.sex_chroms.{}".format(
				bam_dir, sample_id, rg_id),
			fastq_files, cpus, "None",
			[str(x).strip() for x in bwa_flags.split()])
	return new_bam


def swap_sex_chroms(
	input_bam_obj, new_bam_obj, samtools_path, sambamba_path, x_chromosome,
	y_chromosome, bam_dir, sample_id, cpus, xyalign_params):
	"""
	Switches sex chromosmes from new_bam_file with those in original bam file

	Parameters
	----------

	input_bam_obj : bam.BamFile() object
		Original input bam file object
	new_bam_obj : bam.BamFile() object
		Bam file object containing newly mapped sex chromosomes (to insert)
	samtools_path : str
		Path/command to call samtools
	sambamba_path : str
		Path/command to call sambamba
	x_chromosome : list
		X-linked scaffolds
	y_chromosome : str
		Y-linked scaffolds
	bam_dir : str
		Path to bam output directory
	sample_id : str
	cpus : int
		Number of threads/cpus
	xyalign_params : dict
		Dictionary of xyalign_params to add to bam header

	Returns
	-------

	str
		Path to new bam file containing original autosomes and new sex chromosomes
	"""
	merged_bam = bam.switch_sex_chromosomes_sambamba(
		samtools_path, sambamba_path, input_bam_obj.filepath,
		new_bam_obj.filepath, x_chromosome + y_chromosome,
		bam_dir, sample_id, cpus, xyalign_params)
	return merged_bam


def main():
	# Version - placeholder for now - need to incorporate it into __init__.py
	citation = """
	XYalign: Inferring Sex Chromosome Ploidy in NGS Data

	Timothy H Webster, Tanya Phung, Madeline Couse, Bruno Grande, Eric Karlins,
	Phillip Richmond, Whitney Whitford, Melissa A. Wilson Sayres

	2017

	Version: {}
	""".format(__version__)

	# Start timer
	xyalign_start = time.time()

	# Grab arguments
	args = parse_args()

	# Create directory structure if not already in place
	utils.validate_dir(args.output_dir, "fastq")
	utils.validate_dir(args.output_dir, "bam")
	utils.validate_dir(args.output_dir, "reference")
	utils.validate_dir(args.output_dir, "bed")
	utils.validate_dir(args.output_dir, "vcf")
	utils.validate_dir(args.output_dir, "plots")
	utils.validate_dir(args.output_dir, "results")
	utils.validate_dir(args.output_dir, "logfiles")

	# Set up logfile
	logfile_path = os.path.join(args.output_dir, "logfiles")
	if args.logfile is not None:
		logfile = os.path.join(
			logfile_path, args.logfile)
	else:
		logfile = os.path.join(
			logfile_path, "{}_xyalign.log".format(
				args.sample_id))

	# Initiate logging
	logger = logging.getLogger("xyalign")
	logger.setLevel(logging.DEBUG)
	# Create File Handler
	fh = logging.FileHandler(logfile, mode="w")
	fh.setLevel(logging.DEBUG)
	# Create Console Handler
	ch = logging.StreamHandler()
	if args.reporting_level == "DEBUG":
		ch.setLevel(logging.DEBUG)
	elif args.reporting_level == "INFO":
		ch.setLevel(logging.INFO)
	elif args.reporting_level == "ERROR":
		ch.setLevel(logging.ERROR)
	elif args.reporting_level == "CRITICAL":
		ch.setLevel(logging.CRITICAL)
	else:
		ch.setLevel(logging.INFO)
	# Set Formatter
	formatter = logging.Formatter(
		'%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	fh.setFormatter(formatter)
	ch.setFormatter(formatter)
	# Add handlers to logger
	logger.addHandler(fh)
	logger.addHandler(ch)

	# Log citation
	logger.info("{}".format(citation))

	# Set up param dictionary
	xyalign_params_dict = {'ID': 'XYalign', 'VN': __version__, 'CL': []}
	p = ""
	for arg in args.__dict__:
		p = p + "{}={}, ".format(arg, args.__dict__[arg])
		xyalign_params_dict['CL'].append("{}={}".format(arg, args.__dict__[arg]))

	# Log parameters and pipeline start
	logger.info("Parameters: {}".format(p))
	logger.info("Beginning XYalign")

	# Check for external programs
	logger.info("Checking external programs")
	utils.validate_external_prog(args.repairsh_path, "bbmap's repair.sh")
	utils.validate_external_prog(args.bedtools_path, "bedtools")
	utils.validate_external_prog(args.bwa_path, "bwa")
	utils.validate_external_prog(args.platypus_path, "platypus")
	utils.validate_external_prog(args.samtools_path, "samtools")
	utils.validate_external_prog(args.sambamba_path, "sambamba")
	logger.info("All external program paths successfully checked")

	# Setup output paths
	fastq_path = os.path.join(args.output_dir, "fastq")
	bam_path = os.path.join(args.output_dir, "bam")
	reference_path = os.path.join(args.output_dir, "reference")
	bed_path = os.path.join(args.output_dir, "bed")
	vcf_path = os.path.join(args.output_dir, "vcf")
	plots_path = os.path.join(args.output_dir, "plots")
	results_path = os.path.join(args.output_dir, "results")

	# Create paths for output files
	# reference-related
	if args.xx_ref_out_name is not None:
		xx_out = os.path.join(reference_path, args.xx_ref_out_name)
	elif args.xx_ref_out is not None:
		xx_out = args.xx_ref_out
	else:
		xx_out = os.path.join(reference_path, "xyalign_noY.masked.fa")
	if args.xy_ref_out_name is not None:
		xy_out = os.path.join(reference_path, args.xy_ref_out_name)
	elif args.xy_ref_out is not None:
		xy_out = args.xy_ref_out
	else:
		xy_out = os.path.join(reference_path, "xyalign_withY.masked.fa")
	# variant/vcf related
	noprocessing_vcf = os.path.join(
		vcf_path, "{}.noprocessing.vcf".format(
			args.sample_id))
	postprocessing_vcf = os.path.join(
		vcf_path, "{}.postprocessing.vcf".format(
			args.sample_id))
	if args.platypus_logfile is not None:
		plat_log = args.platypus_logfile
	else:
		plat_log = args.sample_id
	noprocessing_vcf_log = os.path.join(
		logfile_path, "{}_noprocessing_platypus.log".format(
			plat_log))
	postprocessing_vcf_log = os.path.join(
		logfile_path, "{}_postprocessing_platypus.log".format(
			plat_log))
	readbalance_prefix_noprocessing = os.path.join(
		plots_path, "{}_noprocessing".format(args.sample_id))
	readbalance_prefix_postprocessing = os.path.join(
		plots_path, "{}_postprocessing".format(args.sample_id))
	# Depth/mapq related
	chrom_stats_mapq = os.path.join(
		results_path, "{}_chrom_stats_mapq.txt".format(args.sample_id))
	chrom_stats_depth = os.path.join(
		results_path, "{}_chrom_stats_depth.txt".format(args.sample_id))
	chrom_stats_count = os.path.join(
		results_path, "{}_chrom_stats_count.txt".format(args.sample_id))
	depth_mapq_prefix_noprocessing = os.path.join(
		plots_path, "{}_noprocessing".format(args.sample_id))
	depth_mapq_prefix_postprocessing = os.path.join(
		plots_path, "{}_postprocessing".format(args.sample_id))
	data_frame_preprocessing = os.path.join(
		bed_path, "{}_full_dataframe_depth_mapq_preprocessing.csv".format(
			args.sample_id))
	data_frame_postprocessing = os.path.join(
		bed_path, "{}_full_dataframe_depth_mapq_postprocessing.csv".format(
			args.sample_id))
	data_frame_readbalance_preprocessing = os.path.join(
		bed_path, "{}_full_dataframe_readbalance_preprocessing.csv".format(
			args.sample_id))
	data_frame_readbalance_postprocessing = os.path.join(
		bed_path, "{}_full_dataframe_readbalance_postprocessing.csv".format(
			args.sample_id))
	high_prefix = "{}_highquality_preprocessing".format(args.sample_id)
	output_bed_high_preprocessing = os.path.join(
		bed_path, "{}.bed".format(high_prefix))
	low_prefix = "{}_lowquality_preprocessing".format(args.sample_id)
	output_bed_low_preprocessing = os.path.join(
		bed_path, "{}.bed".format(low_prefix))
	high_prefix_postprocessing = "{}_highquality_postprocessing".format(
		args.sample_id)
	output_bed_high_postprocessing = os.path.join(
		bed_path, "{}.bed".format(high_prefix_postprocessing))
	low_prefix_postprocessing = "{}_lowquality_postprocessing".format(
		args.sample_id)
	output_bed_low_postprocessing = os.path.join(
		bed_path, "{}.bed".format(low_prefix_postprocessing))

	# Set cleanup
	if args.no_cleanup is True:
		cleanup = False
	else:
		cleanup = True

	######################################
	#            Run XYalign             #
	######################################
	# Reference Prep Only
	if args.PREPARE_REFERENCE is True:
		logger.info(
			"PREPARE_REFERENCE set, so only preparing reference fastas.")
		ref = reftools.RefFasta(args.ref, args.samtools_path, args.bwa_path)

		ref_prep(
			ref_obj=ref,
			ref_mask=args.reference_mask,
			ref_dir=reference_path,
			xx=xx_out,
			xy=xy_out,
			y_chromosome=args.y_chromosome,
			samtools_path=args.samtools_path,
			bwa_path=args.bwa_path,
			bwa_index=args.bwa_index)

		logger.info("PREPARE_REFERENCE complete.")
		logger.info("XYalign complete. Elapsed time: {} seconds".format(
			time.time() - xyalign_start))
		logging.shutdown()
		sys.exit(0)

	if args.CHROM_STATS is True:
		logger.info(
			"CHROM_STATS set, will iterate through bam files to calculate "
			"chromosome-wide averages.")
		bam_list = [bam.BamFile(x, args.samtools_path) for x in args.bam]

		chrom_stats_results = chrom_stats(
			bam_obj_list=bam_list,
			chrom_list=args.chromosomes,
			use_counts=args.use_counts)

		if args.use_counts is False:
			cs_depth_dict = chrom_stats_results[0]
			cs_mapq_dict = chrom_stats_results[1]

			with open(chrom_stats_depth, "w") as f:
				for i in cs_depth_dict:
					out_line = "\t".join([str(x) for x in cs_depth_dict[i]])
					f.write("{}\n".format(out_line))

			with open(chrom_stats_mapq, "w") as f:
				for i in cs_mapq_dict:
					out_line = "\t".join([str(x) for x in cs_mapq_dict[i]])
					f.write("{}\n".format(out_line))
		else:
			cs_count_dict = chrom_stats_results[0]

			with open(chrom_stats_count, "w") as f:
				for i in cs_count_dict:
					out_line = "\t".join([str(x) for x in cs_count_dict[i]])
					f.write("{}\n".format(out_line))

		logger.info("CHROM_STATS complete.")
		logger.info("XYalign complete. Elapsed time: {} seconds".format(
			time.time() - xyalign_start))
		logging.shutdown()
		sys.exit(0)

	input_bam = bam.BamFile(args.bam[0], args.samtools_path)

	if args.chromosomes == ["ALL"] or args.chromosomes == ["all"]:
		input_chromosomes = list(input_bam.chromosome_names())
	else:
		input_chromosomes = args.chromosomes

	if any(
		[args.ANALYZE_BAM, args.CHARACTERIZE_SEX_CHROMS, args.STRIP_READS]) is True:
		missing_chroms = input_bam.check_chrom_in_bam(input_chromosomes)
		if len(missing_chroms) != 0:
			logger.error(
				"One or more chromosomes provided via --chromosomes not "
				"present in bam file. Exiting.")
			logging.shutdown()
			sys.exit(1)

	# Strip reads only
	# This module is isolated first because it does not require a reference fasta
	if args.STRIP_READS is True:
		logger.info(
			"STRIP_READS set, so only stripping reads from {}.".format(
				input_bam.filepath))
		if args.chromosomes == ["ALL"] or args.chromosomes == ["all"]:
			stripped_fastqs = input_bam.strip_reads(
				args.repairsh_path, args.shufflesh_path, args.single_end,
				fastq_path, args.sample_id,
				[], args.xmx, args.fastq_compression, cleanup, args.read_group_id)
		else:
			stripped_fastqs = input_bam.strip_reads(
				args.repairsh_path, args.shufflesh_path, args.single_end,
				fastq_path, args.sample_id,
				input_chromosomes, args.xmx, args.fastq_compression, cleanup,
				args.read_group_id)
		logger.info("STRIP_READS complete. Output in {}".format(fastq_path))
		logger.info("XYalign complete. Elapsed time: {} seconds".format(
			time.time() - xyalign_start))
		logging.shutdown()
		sys.exit(0)

	# Other modules
	ref = reftools.RefFasta(args.ref, args.samtools_path, args.bwa_path)

	# Check to ensure bam and fasta are compatible and imports work
	if args.skip_compatibility_check is False:
		compatible = utils.check_bam_fasta_compatibility(input_bam, ref)
		if compatible is False:
			logger.error(
				"Exiting due to compatibility issues between {} and {}. Check: "
				"1) that this fasta was used when generating this bam, 2) that "
				"sequence lengths and names are identical between the two files. "
				"You can check 2) by comparing the bam header (e.g., "
				"samtools -H {}) with a sequence dictionary (.dict) created for "
				"{}.  If you think this is an error, you can use the flag "
				"--skip_compatibility_check.".format(
					input_bam.filepath, ref.filepath, input_bam.filepath,
					ref.filepath))
			logging.shutdown()
			sys.exit(1)

	# Stats Only
	if args.ANALYZE_BAM is True:
		logger.info(
			"ANALYZE_BAM set, so only running steps required "
			"for bam analysis.")
		if args.platypus_calling == "both" or args.platypus_calling == "before":
			call_variants = True
		else:
			call_variants = False
		bam_analysis(
			input_bam_obj=input_bam,
			platypus_calling=call_variants,
			platypus_path=args.platypus_path,
			vcf_log=noprocessing_vcf_log,
			ref_obj=ref,
			input_chroms=input_chromosomes,
			cpus=args.cpus,
			out_vcf=noprocessing_vcf,
			no_variant_plots=args.no_variant_plots,
			window_size=args.window_size,
			target_bed=args.target_bed,
			sample_id=args.sample_id,
			readbalance_prefix=readbalance_prefix_noprocessing,
			variant_site_quality=args.variant_site_quality,
			variant_genotype_quality=args.variant_genotype_quality,
			variant_depth=args.variant_depth,
			marker_size=args.marker_size,
			marker_transparency=args.marker_transparency,
			homogenize_read_balance=args.homogenize_read_balance,
			data_frame_readbalance=data_frame_readbalance_preprocessing,
			min_variant_count=args.min_variant_count,
			no_bam_analysis=args.no_bam_analysis,
			ignore_duplicates=args.ignore_duplicates,
			exact_depth=args.exact_depth,
			whole_genome_threshold=args.whole_genome_threshold,
			mapq_cutoff=args.mapq_cutoff,
			min_depth_filter=args.min_depth_filter,
			max_depth_filter=args.max_depth_filter,
			depth_mapq_prefix=depth_mapq_prefix_noprocessing,
			bam_data_frame=data_frame_preprocessing,
			output_bed_high=output_bed_high_preprocessing,
			output_bed_low=output_bed_low_preprocessing,
			use_bed_for_platypus=False,
			coordinate_scale=args.coordinate_scale,
			fixed=args.include_fixed)

		logger.info("ANALYZE_BAM complete.")
		logger.info("XYalign complete. Elapsed time: {} seconds".format(
			time.time() - xyalign_start))
		logging.shutdown()
		sys.exit(0)

	# Ploidy Estimation Only
	elif args.CHARACTERIZE_SEX_CHROMS is True:
		logger.info(
			"CHARACTERIZE_SEX_CHROMS set, so only running steps required "
			"for to characterize sex chromosome complement. Note that "
			"this involve both playtpus calling and bam analysis too.")

		if args.platypus_calling == "both" or args.platypus_calling == "before":
			call_variants = True
		else:
			call_variants = False
		bam_dfs = bam_analysis(
			input_bam_obj=input_bam,
			platypus_calling=call_variants,
			platypus_path=args.platypus_path,
			vcf_log=noprocessing_vcf_log,
			ref_obj=ref,
			input_chroms=input_chromosomes,
			cpus=args.cpus,
			out_vcf=noprocessing_vcf,
			no_variant_plots=args.no_variant_plots,
			window_size=args.window_size,
			target_bed=args.target_bed,
			sample_id=args.sample_id,
			readbalance_prefix=readbalance_prefix_noprocessing,
			variant_site_quality=args.variant_site_quality,
			variant_genotype_quality=args.variant_genotype_quality,
			variant_depth=args.variant_depth,
			marker_size=args.marker_size,
			marker_transparency=args.marker_transparency,
			homogenize_read_balance=args.homogenize_read_balance,
			data_frame_readbalance=data_frame_readbalance_preprocessing,
			min_variant_count=args.min_variant_count,
			no_bam_analysis=args.no_bam_analysis,
			ignore_duplicates=args.ignore_duplicates,
			exact_depth=args.exact_depth,
			whole_genome_threshold=args.whole_genome_threshold,
			mapq_cutoff=args.mapq_cutoff,
			min_depth_filter=args.min_depth_filter,
			max_depth_filter=args.max_depth_filter,
			depth_mapq_prefix=depth_mapq_prefix_noprocessing,
			bam_data_frame=data_frame_preprocessing,
			output_bed_high=output_bed_high_preprocessing,
			output_bed_low=output_bed_low_preprocessing,
			use_bed_for_platypus=False,
			coordinate_scale=args.coordinate_scale,
			fixed=args.include_fixed)

		ploidy_results = ploidy_analysis(
			passing_df=bam_dfs[0],
			failing_df=bam_dfs[1],
			no_perm_test=args.no_perm_test,
			no_ks_test=args.no_ks_test,
			no_bootstrap=args.no_bootstrap,
			input_chroms=input_chromosomes,
			x_chromosome=args.x_chromosome,
			y_chromosome=args.y_chromosome,
			results_dir=results_path,
			num_permutations=args.num_permutations,
			num_bootstraps=args.num_bootstraps,
			sample_id=args.sample_id)

		logger.info("CHARACTERIZE_SEX_CHROMS complete.")
		logger.info("XYalign complete. Elapsed time: {} seconds".format(
			time.time() - xyalign_start))
		logging.shutdown()
		sys.exit(0)

	# Remapping Only
	elif args.REMAPPING is True:
		logger.info(
			"REMAPPING set, so only running steps required for remapping. "
			"--y_present or --y_absent must be set")
		if args.y_present is True:
			y_present = True
			logger.info("--y_present provided, so remapping for XY individual")
		elif args.y_absent is True:
			logger.info("--y_absent provided, so remapping for XX individual")
			y_present = False
		else:
			logger.error(
				"--y_present or --y_absent required for --REMAPPING. "
				"Exiting.")
			logging.shutdown()
			sys.exit(1)
		if args.xx_ref_in is None or args.xy_ref_in is None:
			logger.info(
				"Input masked reference not provided for both "
				"XX and XY mapping, so creating both")

			masked_refs = ref_prep(
				ref_obj=ref,
				ref_mask=args.reference_mask,
				ref_dir=reference_path,
				xx=xx_out,
				xy=xy_out,
				y_chromosome=args.y_chromosome,
				samtools_path=args.samtools_path,
				bwa_path=args.bwa_path,
				bwa_index=True)

		else:
			xx = reftools.RefFasta(
				args.xx_ref_in, args.samtools_path, args.bwa_path)
			xy = reftools.RefFasta(
				args.xy_ref_in, args.samtools_path, args.bwa_path)
			xx.conditional_index_bwa()
			xx.conditional_seq_dict()
			xy.conditional_index_bwa()
			xy.conditional_seq_dict()
			masked_refs = (xx, xy)

		sex_chrom_bam = bam.BamFile(
			remapping(
				input_bam_obj=input_bam,
				y_pres=y_present,
				masked_references=masked_refs,
				samtools_path=args.samtools_path,
				sambamba_path=args.sambamba_path,
				repairsh_path=args.repairsh_path,
				shufflesh_path=args.shufflesh_path,
				bwa_path=args.bwa_path,
				bwa_flags=args.bwa_flags,
				single_end=args.single_end,
				bam_dir=bam_path,
				fastq_dir=fastq_path,
				sample_id=args.sample_id,
				x_chromosome=args.x_chromosome,
				y_chromosome=args.y_chromosome,
				cpus=args.cpus,
				xmx=args.xmx,
				fastq_compression=args.fastq_compression,
				cleanup=cleanup,
				read_group_id=args.read_group_id),
			args.samtools_path)

		if args.sex_chrom_bam_only is True:
			final_bam = sex_chrom_bam
		else:
			final_bam = bam.BamFile(
				swap_sex_chroms(
					input_bam_obj=input_bam,
					new_bam_obj=sex_chrom_bam,
					samtools_path=args.samtools_path,
					sambamba_path=args.sambamba_path,
					x_chromosome=args.x_chromosome,
					y_chromosome=args.y_chromosome,
					bam_dir=bam_path,
					sample_id=args.sample_id,
					cpus=args.cpus,
					xyalign_params=xyalign_params_dict),
				args.samtools_path)

		logger.info("REMAPPING complete.")
		logger.info("XYalign complete. Elapsed time: {} seconds".format(
			time.time() - xyalign_start))
		logging.shutdown()
		sys.exit(0)

	# Full Pipeline
	else:
		logger.info(
			"Running entire XYalign pipeline")
		missing_chroms = input_bam.check_chrom_in_bam(input_chromosomes)
		if len(missing_chroms) != 0:
			logger.error(
				"One or more chromosomes provided via --chromosomes not "
				"present in bam file. Exiting.")
			logging.shutdown()
			sys.exit(1)
		if args.xx_ref_in is None or args.xy_ref_in is None:
			logger.info(
				"Input masked reference not provided for both "
				"XX and XY mapping, so creating both")

			masked_refs = ref_prep(
				ref_obj=ref,
				ref_mask=args.reference_mask,
				ref_dir=reference_path,
				xx=xx_out,
				xy=xy_out,
				y_chromosome=args.y_chromosome,
				samtools_path=args.samtools_path,
				bwa_path=args.bwa_path,
				bwa_index=True)

		else:
			xx = reftools.RefFasta(
				args.xx_ref_in, args.samtools_path, args.bwa_path)
			xy = reftools.RefFasta(
				args.xy_ref_in, args.samtools_path, args.bwa_path)
			xx.conditional_index_bwa()
			xx.conditional_seq_dict()
			xy.conditional_index_bwa()
			xy.conditional_seq_dict()
			masked_refs = (xx, xy)

		if args.platypus_calling == "both" or args.platypus_calling == "before":
			call_variants = True
		else:
			call_variants = False
		bam_dfs = bam_analysis(
			input_bam_obj=input_bam,
			platypus_calling=call_variants,
			platypus_path=args.platypus_path,
			vcf_log=noprocessing_vcf_log,
			ref_obj=ref,
			input_chroms=input_chromosomes,
			cpus=args.cpus,
			out_vcf=noprocessing_vcf,
			no_variant_plots=args.no_variant_plots,
			window_size=args.window_size,
			target_bed=args.target_bed,
			sample_id=args.sample_id,
			readbalance_prefix=readbalance_prefix_noprocessing,
			variant_site_quality=args.variant_site_quality,
			variant_genotype_quality=args.variant_genotype_quality,
			variant_depth=args.variant_depth,
			marker_size=args.marker_size,
			marker_transparency=args.marker_transparency,
			homogenize_read_balance=args.homogenize_read_balance,
			data_frame_readbalance=data_frame_readbalance_preprocessing,
			min_variant_count=args.min_variant_count,
			no_bam_analysis=args.no_bam_analysis,
			ignore_duplicates=args.ignore_duplicates,
			exact_depth=args.exact_depth,
			whole_genome_threshold=args.whole_genome_threshold,
			mapq_cutoff=args.mapq_cutoff,
			min_depth_filter=args.min_depth_filter,
			max_depth_filter=args.max_depth_filter,
			depth_mapq_prefix=depth_mapq_prefix_noprocessing,
			bam_data_frame=data_frame_preprocessing,
			output_bed_high=output_bed_high_preprocessing,
			output_bed_low=output_bed_low_preprocessing,
			use_bed_for_platypus=False,
			coordinate_scale=args.coordinate_scale,
			fixed=args.include_fixed)

		ploidy_results = ploidy_analysis(
			passing_df=bam_dfs[0],
			failing_df=bam_dfs[1],
			no_perm_test=args.no_perm_test,
			no_ks_test=args.no_ks_test,
			no_bootstrap=args.no_bootstrap,
			input_chroms=input_chromosomes,
			x_chromosome=args.x_chromosome,
			y_chromosome=args.y_chromosome,
			results_dir=results_path,
			num_permutations=args.num_permutations,
			num_bootstraps=args.num_bootstraps,
			sample_id=args.sample_id)

		if args.y_present is True:
			y_present = True
			logger.info("--y_present provided, so remapping for XY individual")
		elif args.y_absent is True:
			y_present = False
			logger.info("--y_absent provided, so remapping for XX individual")
		else:
			if ploidy_results["boot"][2][2] > args.sex_chrom_calling_threshold:
				y_present = False
				logger.info(
					"X/Y depth ratio ({}) > {}. Y inferred to be absent.".format(
						ploidy_results["boot"][2][2], args.sex_chrom_calling_threshold))
			else:
				y_present = True
				logger.info(
					"X/Y depth ratio ({}) <= {}. Y inferred to be present.".format(
						ploidy_results["boot"][2][2], args.sex_chrom_calling_threshold))

		sex_chrom_bam = bam.BamFile(
			remapping(
				input_bam_obj=input_bam,
				y_pres=y_present,
				masked_references=masked_refs,
				samtools_path=args.samtools_path,
				sambamba_path=args.sambamba_path,
				repairsh_path=args.repairsh_path,
				shufflesh_path=args.shufflesh_path,
				bwa_path=args.bwa_path,
				bwa_flags=args.bwa_flags,
				single_end=args.single_end,
				bam_dir=bam_path,
				fastq_dir=fastq_path,
				sample_id=args.sample_id,
				x_chromosome=args.x_chromosome,
				y_chromosome=args.y_chromosome,
				cpus=args.cpus,
				xmx=args.xmx,
				fastq_compression=args.fastq_compression,
				cleanup=cleanup,
				read_group_id=args.read_group_id),
			args.samtools_path)

		if args.sex_chrom_bam_only is True:
			final_bam = sex_chrom_bam
		else:
			final_bam = bam.BamFile(
				swap_sex_chroms(
					input_bam_obj=input_bam,
					new_bam_obj=sex_chrom_bam,
					samtools_path=args.samtools_path,
					sambamba_path=args.sambamba_path,
					x_chromosome=args.x_chromosome,
					y_chromosome=args.y_chromosome,
					bam_dir=bam_path,
					sample_id=args.sample_id,
					cpus=args.cpus,
					xyalign_params=xyalign_params_dict),
				args.samtools_path)

		if args.platypus_calling == "both" or args.platypus_calling == "after":
			call_variants = True
		else:
			call_variants = False
		bam_analysis(
			input_bam_obj=final_bam,
			platypus_calling=call_variants,
			platypus_path=args.platypus_path,
			vcf_log=postprocessing_vcf_log,
			ref_obj=ref,
			input_chroms=input_chromosomes,
			cpus=args.cpus,
			out_vcf=postprocessing_vcf,
			no_variant_plots=args.no_variant_plots,
			window_size=args.window_size,
			target_bed=args.target_bed,
			sample_id=args.sample_id,
			readbalance_prefix=readbalance_prefix_postprocessing,
			variant_site_quality=args.variant_site_quality,
			variant_genotype_quality=args.variant_genotype_quality,
			variant_depth=args.variant_depth,
			marker_size=args.marker_size,
			marker_transparency=args.marker_transparency,
			homogenize_read_balance=args.homogenize_read_balance,
			data_frame_readbalance=data_frame_readbalance_postprocessing,
			min_variant_count=args.min_variant_count,
			no_bam_analysis=args.no_bam_analysis,
			ignore_duplicates=args.ignore_duplicates,
			exact_depth=args.exact_depth,
			whole_genome_threshold=args.whole_genome_threshold,
			mapq_cutoff=args.mapq_cutoff,
			min_depth_filter=args.min_depth_filter,
			max_depth_filter=args.max_depth_filter,
			depth_mapq_prefix=depth_mapq_prefix_postprocessing,
			bam_data_frame=data_frame_postprocessing,
			output_bed_high=output_bed_high_postprocessing,
			output_bed_low=output_bed_low_postprocessing,
			use_bed_for_platypus=True,
			coordinate_scale=args.coordinate_scale,
			fixed=args.include_fixed)

		logger.info("XYalign complete. Elapsed time: {} seconds".format(
			time.time() - xyalign_start))
		logging.shutdown()
		sys.exit(0)


if __name__ == "__main__":
	main()
