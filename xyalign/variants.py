# Part of XYalign
# Functions for calling and processing variants

from __future__ import division
import logging
import subprocess
import time
import pysam
import bam
import cyvcf2
# Matplotlib needs to be called in this way to set the display variable
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt


# Create logger for variants submodule
variants_logger = logging.getLogger("xyalign.variants")


class VCFFile():
	"""
	A class for working with external vcf files.

	Attributes
	----------

	filepath : str
		Full path to external vcf file
	bgzip : str
		Full path to bgzip. Default = 'bgzip'
	tabix : str
		Full path to tabix. Default = "tabix"

	"""
	def __init__(
		self, filepath, bgzip="bgzip", tabix="tabix",
		no_initial_compress=False):
		self.filepath = filepath
		self.bgzip = bgzip
		self.tabix = tabix
		self.logger = logging.getLogger("xyalign.variants.VCFFile")
		self.logger.info("Creating a VCFFile instance for {}".format(
			self.filepath))
		if no_initial_compress is False:
			if self.is_bgzipped() is False:
				self.compress_vcf()
				self.index_vcf()

	def is_bgzipped(self):
		"""
		Checks to see if vcf file is gzipped, simply by looking for a .gz  or
		.bgz ending.
		If .gz or .bgz ending exists, assumes file is compressed using bgzip.

		Returns
		-------

		bool
			True if ends in .gz, False otherwise

		"""
		self.logger.info("Checking for .gz or .bgz ending in {}".format(
			self.filepath))
		if self.filepath[-3:] == ".gz":
			self.logger.info("Ends in .gz")
			return True
		elif self.filepath[-4:] == ".bgz":
			self.logger.info("Ends in .bgz")
			return True
		else:
			self.logger.info("Does not end in .gz or .bgz")
			return False

	def compress_vcf(self):
		"""
		Compresses vcf file using bgzip.

		Returns
		-------

		bool
			True if successful

		Raises
		-------

		RuntimeError
			If return code from external call is not 0

		"""
		self.logger.info("Compressing vcf file {}".format(self.filepath))
		bgzip_start = time.time()
		rc = subprocess.call([self.bgzip, self.filepath])
		if rc == 0:
			self.logger.info("Compression complete. Elapsed time: {} seconds".format(
				time.time() - bgzip_start))
			self.filepath = self.filepath + ".gz"
			return True
		else:
			self.logger.error("Unable to compress vcf file: {}. Exiting.".format(
				self.filepath))
			logging.shutdown()
			raise RuntimeError("Unable to compress vcf file. Exiting.")

	def index_vcf(self):
		"""
		Indexes vcf file using tabix.  If file does not end in .gz, will
		compress with bgzip (by calling self.compress_vcf).

		Note: Files MUST be compressed using bgzip.

		Returns
		-------

		bool
			True if successful.

		Raises
		------

		RuntimeError
			If return code from external call is not 0.
		"""
		self.logger.info("Indexing vcf file: {}".format(self.filepath))
		index_start = time.time()
		rc = subprocess.call([self.tabix, "-p", "vcf", self.filepath])
		if rc == 0:
			self.logger.info("Indexing complete. Elapsed time: {} seconds.".format(
				time.time() - index_start))
			return True
		else:
			self.logger.info("Unable to index vcf file: {}. Exiting".format(
				self.filepath))
			logging.shutdown()
			raise RuntimeError("unable to index vcf file. Exiting.")

	def read_balance_per_window(self, chrom, window_size, target_file=None):
		"""
		Calculates mean read balance per genomic window (defined by size or an
		external target bed file).  Currently only supports Platypus vcfs bgzipped
		and tabix indexed.

		Parameters
		----------

		chrom : str
			Name of the chromosome to analyze
		window_size
			If int, the window size to use for sliding window analyses, if None
			intervals from target_file
		target_file : str
			Path to bed_file containing regions to analyze instead of
			windows of a fixed size. Will only be engaged if window_size is None

		Returns
		-------

		pandas dataframe
			With columns: "chrom", "start", "stop", "balance", and "count"
		"""
		# vcf_start = time.time()
		# variants_logger.info(
		# 	"Traversing {} in {} to analyze depth and mapping quality".format(
		# 		chrom, self.filepath))
		# vcf = cyvcf2.VCF(vcf)
		pass


def platypus_caller(
	platypus_path, log_path, bam, ref, chroms, cpus, output_file,
	regions_file=None):
	"""
	Uses platypus to make variant calls on provided bam file

	Parameters
	----------

	platypus_path : str
		Path to platypus
	log_path : str
		Path to and name of desired log file for platypus
	bam : str
		Path to input bam (or cram) file
	ref : str
		Path to reference sequence
	chroms : list
		Chromosomes to call variants on, e.g., ["chrX", "chrY", "chr19"]
	cpus : int
		Number of threads/cores to use
	output_file : path
		Path to and name of the output vcf
	regions_file : {str, None}
		If not None, must be path to bed file containing regions to call variants
		in.  If None, calls in call regions of provided chromosomes. Default =
		None.

	Returns
	-------

	int
		Exit code of the platypus call

	"""
	platy_start = time.time()
	if regions_file is None:
		regions = ','.join(map(str, chroms))
	else:
		regions = regions_file
	command_line = [
		platypus_path, "callVariants", "--bamFiles", bam, "-o",
		output_file, "--refFile", ref, "--nCPU", str(cpus), "--regions", regions,
		"--assemble", "1", "--logFileName", log_path]
	variants_logger.info("Calling variants with command line: {}".format(
		" ".join(command_line)))
	return_code = subprocess.call(command_line)
	variants_logger.info(
		"Variant calling complete. Elapsed time: {} seconds".format(
			time.time() - platy_start))
	return return_code


def parse_platypus_VCF(filename, qual_cutoff, chrom):
	"""
	Parse vcf generated by Platypus to grab read balance. Note that this
	is hard-coded to Platypus (version 0.8.1) and will not generalize to vcfs
	generated with other programs (and, potentially, other versions of Platypus)

	Parameters
	----------

	filename : str
		Full path to the input vcf
	qual_cutoff : str
		Minimum (PHRED) quality at which sites should be included
	chrom : str
		Name of the chromosome to include

	Returns
	-------
	tuple
		three corresponding arrays of the same length:
			(position across the chromosome, site quality, read balance)

	"""
	parse_start = time.time()
	variants_logger.info("Parsing {} for read balance.".format(filename))
	infile = open("{}".format(filename), 'r')
	positions = []
	quality = []
	readBalance = []
	for line in infile:
		cols = line.strip('\n').split('\t')
		if cols[0] != chrom:
			continue
		pos = int(cols[1])
		qual = float(cols[5])
		if qual < qual_cutoff:
			continue
		try:
			TR = cols[7].split(';')[17].split('=')[1]
			TC = cols[7].split(';')[14].split('=')[1]
		except IndexError:
			variants_logger.error("Error parsing line: {}".format(line))
			continue
		if ',' in TR or ',' in TC:
			continue
		if (float(TR) == 0) or (float(TC) == 0):
			continue
		ReadRatio = float(TR) / float(TC)

		# Add to arrays
		readBalance.append(ReadRatio)
		positions.append(pos)
		quality.append(qual)
	variants_logger.info("Parsing complete. Elapsed time: {} seconds".format(
		time.time() - parse_start))
	return (positions, quality, readBalance)


def plot_variants_per_chrom(
	chrom_list, vcf_file, sampleID, output_prefix, qual_cutoff,
	MarkerSize, MarkerAlpha, bamfile_obj):
	"""
	Parses a vcf file and plots read balance in separate plots
	for each chromosome in the input list

	Parameters
	----------

	chrom_list : list
		Chromosomes to include
	vcf_file : str
		File (including path) of platypus vcf to analyze
	sampleID : str
		Sample ID (for plot titles)
	output_prefix : str
		Full path to and prefix of desired output plots
	qual_cutoff : int
		Minimum (Phred) quality to consider a site in the vcf
	MarkerSize : float
		Size of markers (matplotlib sizes) to use in the figure
	MarkerAlpha : float
		Transparency (matplotlib values, 0 to 1) of markers
	bamfile_obj : BamFile() object
		Used to get chromosome lengths only

	Returns
	-------
	int
		0

	"""
	plot_start = time.time()
	variants_logger.info("Plotting read balance from {} for chroms: {}".format(
		vcf_file, " ".join(chrom_list)))
	no_sites = []
	for i in chrom_list:
		parse_results = parse_platypus_VCF(vcf_file, qual_cutoff, i)
		if len(parse_results[0]) < 1:
			no_sites.append(i)
		else:
			plot_read_balance(
				i, parse_results[0], parse_results[2],
				sampleID, output_prefix, MarkerSize, MarkerAlpha, bamfile_obj)
			hist_read_balance(
				i, parse_results[2], sampleID, output_prefix)
	if len(no_sites) >= 1:
		variants_logger.info(
			"No variants passing filters on the following chromosomes: {}".format(
				" ".join(no_sites)))
	else:
		variants_logger.info(
			"All chromosomes had variant sites passing filters.")
	variants_logger.info(
		"Read balance plotting complete. Elapsed time: {} seconds".format(
			time.time() - plot_start))
	return 0


def plot_read_balance(
	chrom, positions, readBalance, sampleID, output_prefix, MarkerSize,
	MarkerAlpha, bamfile_obj):
	"""
	Plots read balance at each SNP along a chromosome

	Parameters
	----------

	chrom : str
		Name of the chromosome
	positions : numpy array
		Positions along the chromosome (same length as readBalance)
	readBalance : numpy array
		Read balance corresponding with the positions in the positions array
	sampleID : str
		Sample name or id to include in the plot title
	output_prefix : str
		Desired prefix (including full path) of the output files
	MarkerSize : float
		Size of markers (matplotlib sizes) to use in the figure
	MarkerAlpha : float
		Transparency (matplotlib values) of markers for the figure
	bamfile_obj : BamFile() object
		Used to get chromosome lengths only

	Returns
	-------
	int
		0

	"""
	chrom_len = bamfile_obj.get_chrom_length(chrom)
	if "x" in chrom.lower():
		Color = "green"
	elif "y" in chrom.lower():
		Color = "blue"
	else:
		Color = "red"
	fig = plt.figure(figsize=(15, 5))
	axes = fig.add_subplot(111)
	axes.scatter(
		positions, readBalance, c=Color, alpha=MarkerAlpha, s=MarkerSize, lw=0)
	axes.set_xlim(0, chrom_len)
	axes.set_title(sampleID)
	axes.set_xlabel("Chromosomal Coordinate")
	axes.set_ylabel("Read Balance")
	# print(len(positions))
	plt.savefig("{}_{}_ReadBalance_GenomicScatter.svg".format(
		output_prefix, chrom))
	plt.savefig("{}_{}_ReadBalance_GenomicScatter.png".format(
		output_prefix, chrom))
	plt.close(fig)
	variants_logger.info("Genomic read balance plot of {} complete.".format(
		chrom))
	return 0


def hist_read_balance(chrom, readBalance, sampleID, output_prefix):
	"""
	Plots a histogram of read balance values between 0.05 and 0.95

	Parameters
	----------

	chrom : str
		Name of the chromosome
	readBalance : numpy array
		Read balance values
	sampleID : str
		Sample name or id to include in the plot title
	output_prefix : str
		Desired prefix (including full path) of the output files

	Returns
	-------

	int
		0 if plotting successful, 1 otherwise.

	"""
	try:
		read_balance = readBalance[0.05 < readBalance < 0.95]
	except IndexError:
		variants_logger.info(
			"No sites on {} to plot histogram. Skipping.".format(chrom))
		return 1
	if "x" in chrom.lower():
		Color = "green"
	elif "y" in chrom.lower():
		Color = "blue"
	else:
		Color = "red"
	fig = plt.figure(figsize=(8, 8))
	axes = fig.add_subplot(111)
	axes.set_title(sampleID)
	axes.set_xlabel("Read Balance")
	axes.set_ylabel("Frequency")
	axes.hist(readBalance, bins=50, color=Color)
	plt.savefig("{}_{}_ReadBalance_Hist.svg".format(output_prefix, chrom))
	plt.savefig("{}_{}_ReadBalance_Hist.png".format(output_prefix, chrom))
	plt.close(fig)
	variants_logger.info(
		"Genomic read balance histogram of {} complete.".format(
			chrom))
	return 0
