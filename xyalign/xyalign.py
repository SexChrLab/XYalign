# XYalign main program

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
import pysam
import assemble
import bam
import ploidy
import reftools
import utils
import variants

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
		"--cleanup", type=bool, default=True,
		help="Default is True. Flag determines whether temporary files are "
		"deleted or not.  '--cleanup False' will preserve all temporary "
		"files.")

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
		"--homogenize_read_balance", type=bool, default=True,
		help="If True, read balance values will be transformed by subtracting "
		"each value from 1. For example, 0.25 and 0.75 would be treated "
		"equivalently. Default is True.")

	parser.add_argument(
		"--min_variant_count", type=int, default=2,
		help="Minimum number of variants in a window for the read balance of "
		"that window to be plotted. Note that this does not affect plotting "
		"of variant counts. Default is 2.")

	# Reference-related Flags
	parser.add_argument(
		"--reference_mask", nargs="*",
		help="Bed file containing regions to replace with Ns in the sex "
		"chromosome reference.  Examples might include the pseudoautosomal "
		"regions on the Y to force all mapping/calling on those regions of the "
		"X chromosome.  Default is None.")

	parser.add_argument(
		"--xx_ref_out", default=None,
		help="Desired name for masked output fasta for "
		"samples WITHOUT a Y chromosome (e.g., XX, XXX, XO, etc.). "
		"Defaults to 'xyalign_noY.masked.fa'. Will be output "
		"in the XYalign reference directory.")

	parser.add_argument(
		"--xy_ref_out", default=None,
		help="Desired name for masked output fasta for "
		"samples WITH a Y chromosome (e.g., XY, XXY, etc.). "
		"Defaults to 'xyalign_withY.masked.fa'. Will be output "
		"in the XYalign reference directory.")

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
	ref_obj, ref_mask, ref_dir, xx, xy, y_chromosome, samtools_path, bwa_path):
	"""
	Reference prep part of XYalign pipeline.

	* Creates two reference fasta files.  Both will include masks provied with
	--reference_mask.  One will additionally have the entire Y chromosome
	hard masked.

	* Indexes (.fai, .dict, and bwa indices) both new references

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

	Returns
	-------

	tuple
		Paths to two masked references (y_masked, y_unmasked")
	"""
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
	return (noy_ref, withy_ref)


def chrom_stats(bam_obj_list, chrom_list):
	"""
	Runs chrom stats module.

	Calculates mean depth and mapq across entire scaffolds for a list of bam
	files

	Input
	-----

	bam_obj_list : list
		List of bam.BamFile() objects. Need to have been created using same
		reference (i.e., seq names and lengths are the same)

	Returns
	-------

	"""
	logger.info(
		"Running CHROM_STATS on following bam files: {}".format(
			", ".join([x.filepath for x in bam_obj_list])))

	comp_check = check_compatibility_bam_list(bam_obj_list)
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

	chrom_depth_dict = collections.OrderedDict()
	chrom_mapq_dict = collections.OrderedDict()

	chrom_depth_dict["header"] = [os.path.basename(x.filepath) for x in bam_obj_list]
	for chromosome in chrom_list:
		chrom_results = [x.chrom_stats(chromosome, False) for x in bam_obj_list]
		chrom_depth_dict[chromosome], chrom_mapq_dict[chromosome] = zip(*chrom_results)

		chrom_depth_dict[chromosome].insert(0, chromosome)
		chrom_mapq_dict[chromosome].insert(0, chromosome)

	return (chrom_depth_dict, chrom_mapq_dict)


def bam_analysis_noprocessing():
	"""
	Runs bam analyis part of XYalign pipeline on unprocessed bam file.

	* (Optionally) calls variants using Platypus

	* (Optionally) parses and filters Platypus vcf, and plots read balance

	* (Optionally) Calculates window based metrics from the bam file:
		depth and mapq

	* (optionally) Plots window-based metrics

	* Outputs two bed files: high quality windows, and low quality windows.

	Returns
	-------
	tuple
		(list of pandas dataframes with passing windows,
		list of pandas dataframes with failing windows)
	"""
	# Platypus
	if args.platypus_calling in ["both", "before"]:
		a = input_bam.platypus_caller(
			args.platypus_path, noprocessing_vcf_log,
			ref.filepath, input_chromosomes, args.cpus, noprocessing_vcf,
			None)
		if a != 0:
			logger.error("Error in platypus calling on {}".format(
				input_bam.filepath))
			logging.shutdown()
			sys.exit(1)
		noprocess_vcf_object = variants.VCFFile(noprocessing_vcf)
		if args.no_variant_plots is not True:
			if args.window_size is not None and args.window_size != "None":
				noprocess_vcf_object.plot_variants_per_chrom(
					input_chromosomes, args.sample_id, readbalance_prefix_noprocessing,
					args.variant_site_quality,
					args.variant_genotype_quality,
					args.variant_depth,
					args.marker_size,
					args.marker_transparency, input_bam, "platypus",
					args.homogenize_read_balance, data_frame_readbalance_preprocessing,
					args.min_variant_count, int(args.window_size))
			else:
				noprocess_vcf_object.plot_variants_per_chrom(
					input_chromosomes, args.sample_id, readbalance_prefix_noprocessing,
					args.variant_site_quality,
					args.variant_genotype_quality,
					args.variant_depth,
					args.marker_size,
					args.marker_transparency, input_bam, "platypus",
					args.homogenize_read_balance, data_frame_readbalance_preprocessing,
					args.min_variant_count, None, args.target_bed)
	# Depth/MAPQ
	if args.no_bam_analysis is not True:
		all_df = []
		pass_df = []
		fail_df = []
		for chromosome in input_chromosomes:
			if args.window_size is not None and args.window_size != "None":
				data = input_bam.analyze_bam(
					chromosome, args.ignore_duplicates,
					args.exact_depth, int(args.window_size))
			else:
				data = input_bam.analyze_bam(
					chromosome, args.ignore_duplicates,
					args.exact_depth, None, args.target_bed)
			if args.whole_genome_threshold is True:
				tup = utils.make_region_lists_genome_filters(
					data, args.mapq_cutoff,
					args.min_depth_filter, args.max_depth_filter)
			else:
				tup = utils.make_region_lists_chromosome_filters(
					data, args.mapq_cutoff,
					args.min_depth_filter, args.max_depth_filter)
			all_df.append(data)
			pass_df.append(tup[0])
			fail_df.append(tup[1])
			utils.plot_depth_mapq(
				data, depth_mapq_prefix_noprocessing, args.sample_id,
				input_bam.get_chrom_length(chromosome), args.marker_size,
				args.marker_transparency)
		all_concat = pd.concat(all_df)
		all_concat.to_csv(
			data_frame_preprocessing, index=False, sep="\t", quoting=csv.QUOTE_NONE)
		utils.output_bed(output_bed_high, *pass_df)
		utils.output_bed(output_bed_low, *fail_df)
	return(pass_df, fail_df)


def ploidy_analysis(passing_df, failing_df):
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

	Returns
	-------

	dictionary
		Results for each test. Keys: perm, ks, boot.
	"""
	# Permutations
	if args.no_perm_test is False:
		if args.y_chromosome != [None]:
			sex_chromosomes = args.x_chromosome + args.y_chromosome
			perm_res_x = []
			perm_res_y = []
		else:
			sex_chromosomes = args.x_chromosome
			perm_res_x = []
			perm_res_y = None
		autosomes = [
			x for x in input_chromosomes if x not in sex_chromosomes]
		for c in autosomes:
			perm_res_x.append(ploidy.permutation_test_chromosomes(
				pd.concat(passing_df), c,
				str(args.x_chromosome[0]), "chrom",
				"depth", args.num_permutations,
				results_path + "/{}.{}_{}_permutation_results.txt".format(
					args.sample_id, c, str(args.x_chromosome[0]))))
			if perm_res_y is not None:
				perm_res_y.append(ploidy.permutation_test_chromosomes(
					pd.concat(passing_df), c,
					str(args.y_chromosome[0]), "chrom",
					"depth", args.num_permutations,
					results_path + "/{}.{}_{}_permutation_results.txt".format(
						args.sample_id, c, str(args.y_chromosome[0]))))
		if perm_res_y is not None:
			sex_perm_res = ploidy.permutation_test_chromosomes(
				pd.concat(passing_df), str(args.x_chromosome[0]),
				str(args.y_chromosome[0]),
				"chrom", "depth", args.num_permutations,
				results_path + "/{}.{}_{}_permutation_results.txt".format(
					args.sample_id, str(args.x_chromosome[0]), str(args.y_chromosome[0])))
		else:
			sex_perm_res = None

	# K-S Two Sample
	if args.no_ks_test is False:
		if args.y_chromosome != [None]:
			sex_chromosomes = args.x_chromosome + args.y_chromosome
			ks_res_x = []
			ks_res_y = []
		else:
			sex_chromosomes = args.x_chromosome
			ks_res_x = []
			ks_res_y = None
		autosomes = [
			x for x in input_chromosomes if x not in sex_chromosomes]
		for c in autosomes:
			ks_res_x.append(ploidy.ks_two_sample(
				pd.concat(passing_df), c,
				str(args.x_chromosome[0]), "chrom", "depth",
				results_path + "/{}.{}_{}_ks_results.txt".format(
					args.sample_id, c, str(args.x_chromosome[0]))))
			if ks_res_y is not None:
				ks_res_y.append(ploidy.ks_two_sample(
					pd.concat(passing_df), c,
					str(args.y_chromosome[0]), "chrom", "depth",
					results_path + "/{}.{}_{}_ks_results.txt".format(
						args.sample_id, c, str(args.y_chromosome[0]))))
		if ks_res_y is not None:
			sex_ks_res = ploidy.ks_two_sample(
				pd.concat(passing_df), str(args.x_chromosome[0]),
				str(args.y_chromosome[0]), "chrom", "depth",
				results_path + "/{}.{}_{}_ks_results.txt".format(
					args.sample_id, str(args.x_chromosome[0]), str(args.y_chromosome[0])))
		else:
			sex_ks_res = None

	# Bootstrap
	if args.no_bootstrap is False:
		if args.y_chromosome != [None]:
			sex_chromosomes = args.x_chromosome + args.y_chromosome
			boot_res_x = []
			boot_res_y = []
		else:
			sex_chromosomes = args.x_chromosome
			boot_res_x = []
			boot_res_y = None
		autosomes = [
			x for x in input_chromosomes if x not in sex_chromosomes]
		for c in autosomes:
			boot_res_x.append(ploidy.bootstrap(
				pd.concat(passing_df), c,
				str(args.x_chromosome[0]), "chrom",
				"depth", args.num_bootstraps,
				results_path + "/{}.{}_{}_bootstrap_results.txt".format(
					args.sample_id, c, str(args.x_chromosome[0]))))
			if boot_res_y is not None:
				boot_res_y.append(ploidy.bootstrap(
					pd.concat(passing_df), c,
					str(args.y_chromosome[0]), "chrom",
					"depth", args.num_bootstraps,
					results_path + "/{}.{}_{}_bootstrap_results.txt".format(
						args.sample_id, c, str(args.y_chromosome[0]))))
		if boot_res_y is not None:
			sex_boot_res = ploidy.bootstrap(
				pd.concat(passing_df), str(args.x_chromosome[0]),
				str(args.y_chromosome[0]),
				"chrom", "depth", args.num_bootstraps,
				results_path + "/{}.{}_{}_bootstrap_results.txt".format(
					args.sample_id, str(args.x_chromosome[0]), str(args.y_chromosome[0])))
		else:
			sex_boot_res = None
	return {
		"perm": [perm_res_x, perm_res_y, sex_perm_res],
		"ks": [ks_res_x, ks_res_y, sex_ks_res],
		"boot": [boot_res_x, boot_res_y, sex_boot_res]}


def remapping():
	"""
	Runs remapping steps of XYalign.

	* Strips, sorts, and re-pair reads from the sex chromosomes (collecting read
	group information)

	* Maps (with sorting) reads (with read group information) to appropriate
	reference based on presence (or not) of Y chromosome

	* Merge bam files (if more than one read group)

	Returns
	-------

	str
		Path to bam containing remapped sex chromsomes
	"""
	if y_present is True:
		new_reference = masked_refs[1]
	else:
		new_reference = masked_refs[0]
	new_fastqs = input_bam.strip_reads(
		args.repairsh_path, args.shufflesh_path, args.single_end, fastq_path,
		args.sample_id, args.x_chromosome + args.y_chromosome, args.xmx,
		args.fastq_compression, args.cleanup, args.read_group_id)
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
					args.bwa_path, args.samtools_path, args.sambamba_path,
					new_reference, "{}/{}.sex_chroms.{}.".format(
						bam_path, args.sample_id, rg_id),
					fastq_files, args.cpus, rg_tag,
					[str(x).strip() for x in args.bwa_flags.split()])
				temp_bam_list.append(temp_bam)
		if len(temp_bam_list) < 2:
			new_bam = temp_bam_list[0]
		else:
			new_bam = bam.samtools_merge(
				args.samtools_path, temp_bam_list, "{}/{}.sex_chroms".format(
					bam_path, args.sample_id), args.cpus)
	else:
		fastq_files = read_group_and_fastqs[0][1:]
		new_bam = assemble.bwa_mem_mapping_sambamba(
			args.bwa_path, args.samtools_path, args.sambamba_path,
			new_reference, "{}/{}.sex_chroms.{}.".format(
				bam_path, args.sample_id, rg_id),
			fastq_files, args.cpus, "None",
			[str(x).strip() for x in args.bwa_flags.split()])
	return new_bam


def swap_sex_chroms(new_bam_file):
	"""
	Switches sex chromosmes from new_bam_file with those in original bam file

	Parameters
	----------

	new_bam_file : str
		Path to bam file containing newly mapped sex chromosomes (to insert)

	Returns
	-------

	str
		Path to new bam file containing original autosomes and new sex chromosomes
	"""
	merged_bam = bam.switch_sex_chromosomes_sambamba(
		args.samtools_path, args.sambamba_path, input_bam.filepath,
		new_bam_file.filepath, args.x_chromosome + args.y_chromosome,
		bam_path, args.sample_id, args.cpus, xyalign_params_dict)
	return merged_bam


def bam_analysis_postprocessing():
	"""
	Runs bam analyis part of XYalign pipeline on processed bam file.

	* (Optionally) calls variants using Platypus

	* (Optionally) parses and filters Platypus vcf, and plots read balance

	* Calculates window based metrics from the bam file: depth and mapq

	* Plots window-based metrics

	* Outputs two bed files: high quality windows, and low quality windows.

	Returns
	-------

	tuple
		(list of pandas dataframes with passing windows,
		list of pandas dataframes with failing windows)
	"""
	# Depth/MAPQ
	if args.no_bam_analysis is not True:
		all_df = []
		pass_df = []
		fail_df = []
		for chromosome in input_chromosomes:
			if args.window_size is not None and args.window_size != "None":
				data = final_bam.analyze_bam(
					chromosome, args.ignore_duplicates,
					args.exact_depth, int(args.window_size))
			else:
				data = final_bam.analyze_bam(
					chromosome, args.ignore_duplicates, args.exact_depth,
					None, args.target_bed)
			if args.whole_genome_threshold is True:
				tup = utils.make_region_lists_genome_filters(
					data, args.mapq_cutoff, args.depth_filter)
			else:
				tup = utils.make_region_lists_chromosome_filters(
					data, args.mapq_cutoff, args.depth_filter)
			all_df.append(data)
			pass_df.append(tup[0])
			fail_df.append(tup[1])
			utils.plot_depth_mapq(
				data, depth_mapq_prefix_postprocessing, args.sample_id,
				final_bam.get_chrom_length(chromosome), args.marker_size,
				args.marker_transparency)
		all_concat = pd.concat(all_df)
		all_concat.to_csv(
			data_frame_postprocessing, index=False, sep="\t", quoting=csv.QUOTE_NONE)
		utils.output_bed(output_bed_high_postprocessing, *pass_df)
		utils.output_bed(output_bed_low_postprocessing, *fail_df)

	# Platypus
	include_bed = output_bed_high_postprocessing
	if args.platypus_calling in ["both", "after"]:
		a = input_bam.platypus_caller(
			args.platypus_path, postprocessing_vcf_log,
			ref.filepath, input_chromosomes, args.cpus, postprocessing_vcf,
			include_bed)
		if a != 0:
			logger.error("Error in platypus calling on {}".format(
				final_bam.filepath))
			logging.shutdown()
			sys.exit(1)
		postprocess_vcf_object = variants.VCFFile(postprocessing_vcf)
		if args.no_variant_plots is not True:
			if args.window_size is not None and args.window_size != "None":
				postprocess_vcf_object.plot_variants_per_chrom(
					input_chromosomes, args.sample_id, readbalance_prefix_postprocessing,
					args.variant_site_quality,
					args.variant_genotype_quality,
					args.variant_depth,
					args.marker_size,
					args.marker_transparency, final_bam, "platypus",
					args.homogenize_read_balance, data_frame_readbalance_postprocessing,
					args.min_variant_count, int(args.window_size))
			else:
				postprocess_vcf_object.plot_variants_per_chrom(
					input_chromosomes, args.sample_id, readbalance_prefix_postprocessing,
					args.variant_site_quality,
					args.variant_genotype_quality,
					args.variant_depth,
					args.marker_size,
					args.marker_transparency, final_bam, "platypus",
					args.homogenize_read_balance, data_frame_readbalance_postprocessing,
					args.min_variant_count, None, args.target_bed)

	return(pass_df, fail_df)


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
	if args.xx_ref_out is not None:
		xx_out = os.path.join(reference_path, args.xx_ref_out)
	else:
		xx_out = os.path.join(reference_path, "xyalign_noY.masked.fa")
	if args.xy_ref_out is not None:
		xy_out = os.path.join(reference_path, args.xy_ref_out)
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
	depth_mapq_prefix_noprocessing = os.path.join(
		plots_path, "{}_noprocessing".format(args.sample_id))
	depth_mapq_prefix_postprocessing = os.path.join(
		plots_path, "{}_postprocessing".format(args.sample_id))
	data_frame_preprocessing = os.path.join(
		bed_path, "{}_full_dataframe_depth_mapq_preprocessing.csv".format(args.sample_id))
	data_frame_postprocessing = os.path.join(
		bed_path, "{}_full_dataframe_depth_mapq_postprocessing.csv".format(args.sample_id))
	data_frame_readbalance_preprocessing = os.path.join(
		bed_path, "{}_full_dataframe_readbalance_preprocessing.csv".format(args.sample_id))
	data_frame_readbalance_postprocessing = os.path.join(
		bed_path, "{}_full_dataframe_readbalance_postprocessing.csv".format(args.sample_id))
	high_prefix = "{}_highquality_preprocessing".format(args.sample_id)
	output_bed_high = os.path.join(
		bed_path, "{}.bed".format(high_prefix))
	low_prefix = "{}_lowquality_preprocessing".format(args.sample_id)
	output_bed_low = os.path.join(
		bed_path, "{}.bed".format(low_prefix))
	high_prefix_postprocessing = "{}_highquality_postprocessing".format(
		args.sample_id)
	output_bed_high_postprocessing = os.path.join(
		bed_path, "{}.bed".format(high_prefix_postprocessing))
	low_prefix_postprocessing = "{}_lowquality_postprocessing".format(
		args.sample_id)
	output_bed_low_postprocessing = os.path.join(
		bed_path, "{}.bed".format(low_prefix_postprocessing))

	######################################
	#            Run XYalign             #
	######################################
	# Reference Prep Only
	if args.PREPARE_REFERENCE is True:
		logger.info(
			"PREPARE_REFERENCE set, so only preparing reference fastas.")
		ref = reftools.RefFasta(args.ref, args.samtools_path, args.bwa_path)
		ref_prep(
			ref, args.reference_mask, reference_path,
			xx_out, xy_out, args.y_chromosome, args.samtools_path, args.bwa_path)
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
		chrom_stats_results = chrom_stats(bam_list, args.chromosomes)

		with open(chrom_stats_depth, "w") as f:
			for i in chrom_stats_results[0]:
				out_line = "\t".join(i)
				f.write("{}\n".format(out_line))

		with open(chrom_stats_mapq, "w") as f:
			for i in chrom_stats_results[1]:
				out_line = "\t".join(i)
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
				[], args.xmx, args.fastq_compression, args.cleanup, args.read_group_id)
		else:
			stripped_fastqs = input_bam.strip_reads(
				args.repairsh_path, args.shufflesh_path, args.single_end,
				fastq_path, args.sample_id,
				input_chromosomes, args.xmx, args.fastq_compression, args.cleanup,
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
		bam_analysis_noprocessing()
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
		bam_dfs = bam_analysis_noprocessing()
		ploidy_results = ploidy_analysis(bam_dfs[0], bam_dfs[1])
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
				ref, args.reference_mask, reference_path,
				xx_out, xy_out, args.y_chromosome, args.samtools_path, args.bwa_path)
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
		sex_chrom_bam = bam.BamFile(remapping(), args.samtools_path)
		if args.sex_chrom_bam_only is True:
			final_bam = sex_chrom_bam
		else:
			final_bam = bam.BamFile(
				swap_sex_chroms(sex_chrom_bam), args.samtools_path)
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
				ref, args.reference_mask, reference_path,
				xx_out, xy_out, args.y_chromosome, args.samtools_path, args.bwa_path)
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
		bam_dfs = bam_analysis_noprocessing()
		ploidy_results = ploidy_analysis(bam_dfs[0], bam_dfs[1])
		if args.y_present is True:
			y_present = True
			logger.info("--y_present provided, so remapping for XY individual")
		elif args.y_absent is True:
			y_present = False
			logger.info("--y_absent provided, so remapping for XX individual")
		else:
			if ploidy_results["boot"][2] > args.sex_chrom_calling_threshold:
				y_present = False
				logger.info(
					"X/Y depth ratio ({}) > {}. Y inferred to be absent.".format(
						ploidy_results["boot"][2], args.sex_chrom_calling_threshold))
			else:
				y_present = True
				logger.info(
					"X/Y depth ratio ({}) <= {}. Y inferred to be present.".format(
						ploidy_results["boot"][2], args.sex_chrom_calling_threshold))
		sex_chrom_bam = bam.BamFile(remapping(), args.samtools_path)
		if args.sex_chrom_bam_only is True:
			final_bam = sex_chrom_bam
		else:
			final_bam = bam.BamFile(
				swap_sex_chroms(sex_chrom_bam), args.samtools_path)
		bam_analysis_postprocessing()
		logger.info("XYalign complete. Elapsed time: {} seconds".format(
			time.time() - xyalign_start))
		logging.shutdown()
		sys.exit(0)

if __name__ == "__main__":
	main()
