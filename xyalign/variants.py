# Part of XYalign
# Functions for calling and processing variants

from __future__ import division
import csv
import gzip
import logging
import subprocess
import time
import xyalign  # from xyalign import utils
import numpy as np
import pandas as pd
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
		rc = subprocess.call([self.bgzip, "-f", self.filepath])
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
		rc = subprocess.call([self.tabix, "-f", "-p", "vcf", self.filepath])
		if rc == 0:
			self.logger.info("Indexing complete. Elapsed time: {} seconds.".format(
				time.time() - index_start))
			return True
		else:
			self.logger.info("Unable to index vcf file: {}. Exiting".format(
				self.filepath))
			logging.shutdown()
			raise RuntimeError("unable to index vcf file. Exiting.")

	def parse_platypus_VCF(self, site_qual, genotype_qual, depth, chrom):
		"""
		Parse vcf generated by Platypus to grab read balance. Note that this
		is hard-coded to Platypus (version 0.8.1) and will not generalize to vcfs
		generated with other programs (and, potentially, other versions of Platypus)

		Parameters
		----------

		site_qual : int
			Minimum (PHRED) site quality at which sites should be included
		genotype_qual : int
			Minimum (PHRED) genotype quality at which sites should be included
		depth : int
			Minimum depth at which sites should be included
		chrom : str
			Name of the chromosome to include

		Returns
		-------
		tuple
			four corresponding arrays of the same length:
				(position across the chromosome, site quality, read balance,
				genotype quality)

		"""
		parse_start = time.time()
		self.logger.info("Parsing {} for read balance.".format(self.filepath))

		positions = []
		quality = []
		readBalance = []
		gqs = []
		dps = []

		# Right now, cyvcf2 (and the htslib have trouble reading platypus vcfs),
		# so currently hard coding parsing.

		# for variant in cyvcf2.VCF(self.filepath):
		# 	pos = variant.start
		# 	qual = variant.QUAL
		# 	if qual < site_qual:
		# 		continue
		# 	if variant.format("GQ")[0] < genotype_qual:
		# 		continue
		# 	TC = variant.INFO.get("TC")
		# 	TR = variant.INFO.get("TR")
		# 	if TR == 0 or TC == 0:
		# 		continue
		# 	if TR + TC < depth:
		# 		continue

		if self.is_bgzipped() is True:
			infile = gzip.open("{}".format(self.filepath), 'rb')
		else:
			infile = open("{}".format(self.filepath), 'r')

		for line in infile:
			cols = line.strip('\n').split('\t')
			if cols[0] != chrom:
				continue
			pos = int(cols[1])
			qual = float(cols[5])

			# Check site quality
			if qual < site_qual:
				continue

			# Check genotype quality
			try:
				GQ = int(cols[9].split(':')[3])
			except IndexError:
				self.logger.error("Error parsing line: {}".format(line))
				continue
			if GQ < genotype_qual:
				continue

			# Grab allele depths and check for depth
			try:
				TR = cols[7].split(';')[17].split('=')[1]
				TC = cols[7].split(';')[14].split('=')[1]
			except IndexError:
				self.logger.error("Error parsing line: {}".format(line))
				continue
			if ',' in TR or ',' in TC:
				continue
			TR = float(TR)
			TC = float(TC)
			if TR == 0 or TC == 0:
				continue
			if TC < depth:
				continue

			# Site passes filters - grab the read ratio
			ReadRatio = TR / TC

			# Add to arrays
			readBalance.append(ReadRatio)
			positions.append(pos)
			quality.append(qual)
			gqs.append(GQ)
			dps.append(TC)

		infile.close()
		self.logger.info("Parsing complete. Elapsed time: {} seconds".format(
			time.time() - parse_start))
		return (positions, quality, readBalance, gqs, dps)

	def plot_variants_per_chrom(
		self, chrom_list, sampleID, output_prefix, site_qual, genotype_qual,
		depth, MarkerSize, MarkerAlpha, bamfile_obj, variant_caller, homogenize,
		dataframe_out, min_count, window_size, x_scale=1000000, target_file=None):
		"""
		Parses a vcf file and plots read balance in separate plots
		for each chromosome in the input list

		Parameters
		----------

		chrom_list : list
			Chromosomes to include
		sampleID : str
			Sample ID (for plot titles)
		output_prefix : str
			Full path to and prefix of desired output plots
		site_qual : int
			Minimum (PHRED) site quality at which sites should be included
		genotype_qual : int
			Minimum (PHRED) genotype quality at which sites should be included
		depth : int
			Minimum depth at which sites should be included
		MarkerSize : float
			Size of markers (matplotlib sizes) to use in the figure
		MarkerAlpha : float
			Transparency (matplotlib values, 0 to 1) of markers
		bamfile_obj : BamFile() object
			Used to get chromosome lengths only
		variant_caller : str
			Variant caller used to generate vcf - currently only "platypus" supported
		homogenize: bool
			If True, all read balance values less than 0.5 will be transformed
			by subtracting the value from 1. For example, the values 0.25 and
			0.75 would be treated as equivalent.
		dataframe_out : str
			Full path of file to write pandas dataframe to. Will overwire if exists
		min_count : int
			Minimum number of variants to include a window for plotting.
		window_size
			If int, the window size to use for sliding window analyses, if None
			intervals from target_file
		x_scale : int
			Divide all x values (including Xlim) by this value. Default is 1000000 (1MB)
		target_file : str
			Path to bed_file containing regions to analyze instead of
			windows of a fixed size. Will only be engaged if window_size is None

		Returns
		-------
		int
			0 if variants to analyze; 1 if no variants to analyze on any chromosome

		"""
		plot_start = time.time()
		self.logger.info("Plotting read balance from {} for chroms: {}".format(
			self.filepath, " ".join(chrom_list)))
		if variant_caller.lower() != "platypus":
			self.logger.error(
				"Error. Only 'platypus' currently supported as variant_caller.")
			logging.shutdown()
			raise RuntimeError(
				"Error. Only 'platypus' currently supported as variant_caller "
				"for plot_variants_per_chrom.")
		no_sites = []
		all_df = []
		for i in chrom_list:
			parse_results = self.parse_platypus_VCF(
				site_qual, genotype_qual, depth, i)
			if len(parse_results[0]) < 1:
				no_sites.append(i)
			else:
				chrom_len = bamfile_obj.get_chrom_length(i)
				plot_read_balance(
					i, parse_results[0], parse_results[2],
					sampleID, output_prefix, MarkerSize,
					MarkerAlpha, homogenize, chrom_len, x_scale)
				hist_read_balance(
					i, parse_results[2], sampleID, homogenize, output_prefix)
				rb_df = read_balance_per_window(
					i, parse_results[0], parse_results[2], sampleID, homogenize,
					chrom_len, window_size, target_file)
				all_df.append(rb_df)
				xyalign.utils.chromosome_wide_plot(
					i, rb_df["start"].values, rb_df["count"], "Window_Variant_Count",
					sampleID, output_prefix, MarkerSize, MarkerAlpha,
					chrom_len, rb_df["count"].max(), x_scale)
				rb_df = rb_df[rb_df["count"] >= min_count]
				xyalign.utils.chromosome_wide_plot(
					i, rb_df["start"].values, rb_df["balance"], "Window_Read_Balance",
					sampleID, output_prefix, MarkerSize, MarkerAlpha,
					chrom_len, 1.0, x_scale)
		if len(all_df) < 1:
			self.logger.error(
				"No chromosomes with any variants present")
			return 1
		else:
			all_concat = pd.concat(all_df)
			all_concat.to_csv(
				dataframe_out, index=False, sep="\t", quoting=csv.QUOTE_NONE)
			if len(no_sites) >= 1:
				self.logger.info(
					"No variants passing filters on the following chromosomes: {}".format(
						" ".join(no_sites)))
			else:
				self.logger.info(
					"All chromosomes had variant sites passing filters.")
			self.logger.info(
				"Read balance plotting complete. Elapsed time: {} seconds".format(
					time.time() - plot_start))
			return 0


def read_balance_per_window(
	chrom, positions, readBalance, sampleID, homogenize, chr_len,
	window_size, target_file=None):
	"""
	Calculates mean read balance per genomic window (defined by size or an
	external target bed file) for a given chromosome. Takes as input an array
	of positions and an array of read balances - the order of which must
	correspond exactly. In addition, the positions are expected to ALL BE ON
	THE SAME CHROMOSOME and be in numerically sorted order (i.e., the output
	of parse_platypus_VCF())

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
	homogenize: bool
		If True, all read balance values less than 0.5 will be transformed
		by subtracting the value from 1. For example, the values 0.25 and
		0.75 would be treated as equivalent.
	chr_len : int
		Length of chromosome. Ignored if target_file is provided.
	window_size
		If int, the window size to use for sliding window analyses, if None
		intervals from target_file
	target_file : str
		Path to bed file containing regions to analyze instead of
		windows of a fixed size. Will only be engaged if window_size is None

	Returns
	-------

	pandas dataframe
		With columns: "chrom", "start", "stop", "balance", and "count"
	"""
	readbalance_start = time.time()
	variants_logger.info(
		"Traversing {} to calculate mean read balance".format(
			chrom))
	positions = np.asarray(positions)
	readBalance = np.asarray(readBalance)

	if window_size is not None:
		variants_logger.info(
			"Using windows size: {}".format(window_size))

		num_windows = chr_len // window_size + 1
		if chr_len % num_windows == 0:
			last_window_len = window_size
		else:
			last_window_len = chr_len % num_windows

		chr_list = []
		start_list = []
		stop_list = []
		balance_list = []
		count_list = []

		window_id = 0
		start = 0
		end = window_size
		window_count = 0
		window_balances = []

		for idx, i in enumerate(positions):
			if i < start:
				variants_logger.info(
					"Position {} is less than window start {}. Check that "
					"positions are sorted numerically. Exiting.")
				raise RuntimeError(
					"Position {} is less than window start {}. Check that "
					"positions are sorted numerically. Exiting.")
			elif i < end:
				window_count += 1
				if homogenize is False:
					window_balances.append(readBalance[idx])
				else:
					if readBalance[idx] > 0.5:
						window_balances.append(readBalance[idx])
					else:
						window_balances.append(1.0 - readBalance[idx])
			else:
				# i is >= end, so window needs to be reset
				while i >= end:
					chr_list.append(chrom)
					start_list.append(start)
					stop_list.append(end)
					if len(window_balances) == 0:
						balance_list.append(0)
					else:
						balance_list.append(np.mean(window_balances))
					count_list.append(window_count)

					window_id += 1
					if window_id > num_windows:
						variants_logger.info(
							"Exhausted windows, but positions still remaining. "
							"Ensure correct chromosome provided and that positions "
							"are numerically sorted.")
						raise RuntimeError(
							"Exhausted windows, but positions still remaining. "
							"Ensure correct chromosome provided and that positions "
							"are numerically sorted.")
					if window_id == num_windows - 1:
						start += window_size
						end += last_window_len
					else:
						start += window_size
						end += window_size
					window_balances = []
					window_count = 0

				window_count += 1
				if homogenize is False:
					window_balances.append(readBalance[idx])
				else:
					if readBalance[idx] > 0.5:
						window_balances.append(readBalance[idx])
					else:
						window_balances.append(1.0 - readBalance[idx])

		# process last window
		chr_list.append(chrom)
		start_list.append(start)
		stop_list.append(end)
		if len(window_balances) == 0:
			balance_list.append(0)
		else:
			balance_list.append(np.mean(window_balances))
		count_list.append(window_count)

	elif target_file is not None:
		variants_logger.info(
			"Using targets from: {}".format(target_file))
		with open(target_file) as f:
			targets = [x.strip() for x in f]
			targets = [x.split() for x in targets]
			targets = [x for x in targets if x[0] == chrom]
			while [""] in targets:
				targets.remove([""])

		num_windows = len(targets)

		chr_list = [chrom] * num_windows
		start_list = []
		stop_list = []
		balance_list = []
		count_list = []

		window_id = 0
		window_count = 0
		window_balances = []

		start = int(targets[0][1])
		end = int(targets[0][2])

		num_pos = len(positions)
		pos_idx = 0

		while window_id < num_windows:
			if pos_idx < num_pos:
				if positions[pos_idx] < start:
					pos_idx += 1
				elif positions[pos_idx] < end:
					window_count += 1
					if homogenize is False:
						window_balances.append(readBalance[pos_idx])
					else:
						if readBalance[pos_idx] > 0.5:
							window_balances.append(readBalance[pos_idx])
						else:
							window_balances.append(1.0 - readBalance[pos_idx])
					pos_idx += 1
				else:
					start_list.append(start)
					stop_list.append(end)
					if len(window_balances) == 0:
						balance_list.append(0)
					else:
						balance_list.append(np.mean(window_balances))
					count_list.append(window_count)

					window_id += 1
					try:
						start = int(targets[window_id][1])
						end = int(targets[window_id][2])
					except IndexError:
						break
					window_count = 0
					window_balances = []
			else:
					start_list.append(start)
					stop_list.append(end)
					if len(window_balances) == 0:
						balance_list.append(0)
					else:
						balance_list.append(np.mean(window_balances))
					count_list.append(window_count)

					window_id += 1
					try:
						start = int(targets[window_id][1])
						end = int(targets[window_id][2])
					except IndexError:
						break
					window_count = 0
					window_balances = []
	else:
		variants_logger.error(
			"Both window_size and target_file set to None. "
			"Cannot proceed with read balance traversal. Exiting.")
		logging.shutdown()
		raise RuntimeError(
			"Both window_size and target_file set to None. "
			"Cannot proceed with read balance traversal. Exiting.")

	# Convert data into pandas data frames
	windows_df = pd.DataFrame({
		"chrom": np.asarray(chr_list),
		"start": np.asarray(start_list),
		"stop": np.asarray(stop_list),
		"balance": np.asarray(balance_list),
		"count": np.asarray(count_list)
	})[["chrom", "start", "stop", "balance", "count"]]

	variants_logger.info("Analysis complete. Elapsed time: {} seconds".format(
		time.time() - readbalance_start))
	return windows_df


def plot_read_balance(
	chrom, positions, readBalance, sampleID, output_prefix, MarkerSize,
	MarkerAlpha, homogenize, chrom_len, x_scale=1000000):
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
	homogenize: bool
		If True, all read balance values less than 0.5 will be transformed
		by subtracting the value from 1. For example, the values 0.25 and
		0.75 would be treated as equivalent.
	chrom_len : int
		Length of chromosome
	x_scale : int
		Divide all x values (including Xlim) by this value. Default is 1000000 (1MB)

	Returns
	-------
	int
		0

	"""
	if homogenize is True:
		readBalance = np.asarray([1.0 - x if x < 0.5 else x for x in readBalance])
	if "x" in chrom.lower():
		Color = "green"
	elif "y" in chrom.lower():
		Color = "blue"
	else:
		Color = "red"
	fig = plt.figure(figsize=(15, 5))
	axes = fig.add_subplot(111)
	positions = np.divide(positions, float(x_scale))
	axes.scatter(
		positions, readBalance, c=Color, alpha=MarkerAlpha, s=MarkerSize, lw=0)
	axes.set_xlim(0, (chrom_len / float(x_scale)))
	axes.set_title(sampleID)
	if x_scale == 1000000:
		scale_label = "(MB)"
	elif x_scale == 1000:
		scale_label = "(KB)"
	elif x_scale == 1:
		scale_label = "(BP)"
	else:
		scale_label = "(divided by {})".formatt(x_scale)
	axes.set_xlabel("Chromosomal Position {}".format(scale_label))
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


def hist_read_balance(
	chrom, readBalance, sampleID, homogenize, output_prefix):
	"""
	Plots a histogram of read balance values between 0.05 and 0.95

	Parameters
	----------

	chrom : str
		Name of the chromosome
	readBalance : list or numpy array
		Read balance values
	sampleID : str
		Sample name or id to include in the plot title
	homogenize: bool
		If True, all read balance values less than 0.5 will be transformed
		by subtracting the value from 1. For example, the values 0.25 and
		0.75 would be treated as equivalent.
	output_prefix : str
		Desired prefix (including full path) of the output files

	Returns
	-------

	int
		0 if plotting successful, 1 otherwise.

	"""
	readBalance = np.asarray(readBalance)
	read_balance = readBalance[
		np.where((readBalance > 0.05) & (readBalance < 1))]
	if len(read_balance) == 0:
		variants_logger.info(
			"No sites on {} to plot histogram. Skipping.".format(chrom))
		return 1
	if homogenize is True:
		read_balance = np.asarray([1.0 - x if x < 0.5 else x for x in read_balance])
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
	axes.hist(read_balance, bins=50, color=Color)
	plt.savefig("{}_{}_ReadBalance_Hist.svg".format(output_prefix, chrom))
	plt.savefig("{}_{}_ReadBalance_Hist.png".format(output_prefix, chrom))
	plt.close(fig)
	variants_logger.info(
		"Genomic read balance histogram of {} complete.".format(
			chrom))
	return 0
