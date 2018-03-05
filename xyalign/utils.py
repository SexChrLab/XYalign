# Part of XYalign
# Collection of miscellaneous functions

from __future__ import division
from __future__ import print_function
import logging
import os
import subprocess
import numpy as np
import pandas as pd
import pybedtools
import time
# Matplotlib needs to be called in this way to set the display variable
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


# Create logger for utils submodule
utils_logger = logging.getLogger("xyalign.utils")


def validate_external_prog(prog_path, prog_name):
	"""
	Checks to see if external program can be called using provided path

	Parameters
	----------

	prog_path: path to call program
	prog_name: name of program

	Returns
	-------

	int
		0

	"""
	try:
		a = subprocess.Popen(
			[prog_path], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
		utils_logger.info(
			"{} successfully called using path: {}".format(prog_name, prog_path))
	except OSError as e:
		utils_logger.error(
			"ERROR: {} not available from path: {}".format(prog_name, prog_path))
		logging.shutdown()
		raise OSError("ERROR: {} not available from path: {}".format(
			prog_name, prog_path))
	return 0


def validate_dir(parent_dir, dir_name):
	"""
	Checks if directory exists and if not, creates it.

	Parameters
	----------

	parent_dir : Parent directory name
	dir_name : Name of the new directory

	Returns
	-------

	bool
		whether the directory existed

	"""
	full_path = os.path.join(parent_dir, dir_name)
	exists = os.path.exists(full_path)
	if not exists:
		os.makedirs(full_path)
	return exists


def check_bam_fasta_compatibility(bam_object, fasta_object):
	"""
	Checks to see if bam and fasta sequence names and lengths are
	equivalent (i.e., if it is likely that the bam file was generated
	using the fasta in question).

	Parameters
	----------

	bam_object : BamFile() object
	fasta_object: RefFasta() object

	Returns
	-------

	bool
		True if sequence names and lengths match. False otherwise.

	"""
	utils_logger.info(
		"Checking compatibility of {} and {}".format(
			bam_object.filepath, fasta_object.filepath))

	b_names = bam_object.chromosome_names()
	f_names = fasta_object.chromosome_names()
	b_lengths = bam_object.chromosome_lengths()
	f_lengths = fasta_object.chromosome_lengths()

	if b_lengths == f_lengths:
		if b_names == f_names:
			utils_logger.info("Bam and Fasta are compatible")
			return True
		else:
			utils_logger.error(
				"Bam and Fasta are incompatible.  Sequence lengths are identical "
				"but names are not.  Check chromosome ids to ensure compatibility. "
				"{} ids are: {}\n"
				"{} ids are: {}".format(
					bam_object.filepath, b_names,
					fasta_object.filepath, f_names))
			return False
	else:
		utils_logger.error(
			"Bam and Fasta are incompatible. Check sequence names and lengths "
			"to ensure correct fasta was used as input. Chromosome ids and "
			"lengths for {} are:\n"
			"{}\n"
			"Chromosome ids and lengths for {} are:\n"
			"{}".format(
				bam_object.filepath, zip(b_names, b_lengths),
				fasta_object.filepath, zip(f_names, f_lengths)))
		return False


def check_compatibility_bam_list(bam_obj_list):
	"""
	Checks to see if bam sequence names and lengths are
	equivalent (i.e., if it is likely that the bam files were generated
	using the same reference genome).

	Parameters
	----------

	bam_obj_list : list
		List of bam.BamFile() objects

	Returns
	-------

	bool
		True if sequence names and lengths match. False otherwise.

	"""
	utils_logger.info(
		"Checking compatibility of: ".format(
			", ".join([x.filepath for x in bam_obj_list])))

	seq_names = [x.chromosome_names() for x in bam_obj_list]
	seq_lengths = [x.chromosome_lengths() for x in bam_obj_list]

	# Checking for ANY incompatibility
	uniq_seq_names = [x for x in seq_names if x != seq_names[0]]
	uniq_seq_lens = [x for x in seq_lengths if x != seq_lengths[0]]

	if len(uniq_seq_names) == 0 and len(uniq_seq_lens) == 0:
		return True
	else:
		utils_logger.info("Bam files contain different headers.")
		return False


def merge_bed_files(output_file, *bed_files):
	"""
	This function simply takes an output_file (full path to desired output file)
	and an arbitrary number of external bed files (including full path),
	and merges the bed files into the output_file

	Parameters
	----------

	output_file : str
		Full path to and name of desired output bed file
	*bed_files
		Variable length argument list of external bed files (include full path)

	Returns
	-------

	str
		path to output_file

	"""
	merged_bed_start = time.time()
	utils_logger.info("Merging bed files: {}".format(
		" ".join(bed_files)))
	a = pybedtools.BedTool(bed_files[0])
	for i in bed_files[1:]:
		b = a.cat(i)
	c = b.sort().merge()
	c.saveas(output_file)
	utils_logger.info(
		"Merging bed files complete. Elapsed time: {} seconds".format(
			time.time() - merged_bed_start))
	return output_file


def make_region_lists_genome_filters(
	depthAndMapqDf, mapqCutoff, min_depth, max_depth):
	"""
	Filters a pandas dataframe for mapq and depth based on using all values
	from across the entire genome

	Parameters
	----------

	depthAndMapqDf : pandas dataframe
		Must have 'depth' and 'mapq' columns
	mapqCutoff : int
		The minimum mapq for a window to be considered high quality
	min_depth : float
		Fraction of mean to set as minimum depth
	max_depth : float
		Multiple of mean to set as maximum depth

	Returns
	-------

	tuple
		(passing dataframe, failing dataframe)
	"""
	make_region_lists_start = time.time()
	depth_mean = depthAndMapqDf["depth"].mean()
	depth_sd = depthAndMapqDf["depth"].std()

	depthMin = depth_mean * min_depth
	depthMax = depth_mean * max_depth

	utils_logger.info(
		"Filtering dataframe for mapq (MAPQ >= {}) "
		"and depth (between {} and {})".format(
			mapqCutoff, depthMin, depthMax))

	good = (
		(depthAndMapqDf["mapq"] >= mapqCutoff) &
		(depthAndMapqDf["depth"] > depthMin) &
		(depthAndMapqDf["depth"] < depthMax))
	dfGood = depthAndMapqDf[good]
	dfBad = depthAndMapqDf[~good]

	utils_logger.info("Filtering complete. Elapsed time: {} seconds".format(
		time.time() - make_region_lists_start))
	return (dfGood, dfBad)


def make_region_lists_chromosome_filters(
	depthAndMapqDf, mapqCutoff, min_depth, max_depth):
	"""
	Filters a pandas dataframe for mapq and depth based on thresholds calculated
	per chromosome

	Parameters
	----------

	depthAndMapqDf : pandas dataframe
		Must have 'depth' and 'mapq' columns
	mapqCutoff : int
		The minimum mapq for a window to be considered high quality
	min_depth : float
		Fraction of mean to set as minimum depth
	max_depth : float
		Multiple of mean to set as maximum depth

	Returns
	-------

	tuple
		(passing dataframe, failing dataframe)
	"""
	make_region_lists_start = time.time()
	utils_logger.info(
		"Calculating filtering thresholds for mapq and depth on a per "
		"chromosome basis")

	indices = np.unique(depthAndMapqDf.chrom, return_index=True)[1]
	ordered_chrom_list = [
		depthAndMapqDf.chrom[index] for index in sorted(indices)]

	good_list = []
	bad_list = []

	for i in ordered_chrom_list:
		df = depthAndMapqDf.loc[depthAndMapqDf["chrom"] == i]
		depth_mean = df["depth"].mean()
		depth_sd = df["depth"].std()

		depthMin = depth_mean * min_depth
		depthMax = depth_mean * max_depth

		utils_logger.info(
			"Filtering chromosome {} for mapq (MAPQ >= {}) "
			"and depth (between {} and {})".format(
				i, mapqCutoff, depthMin, depthMax))

		good = (
			(df["chrom"] == i) &
			(df["mapq"] >= mapqCutoff) &
			(df["depth"] > depthMin) &
			(df["depth"] < depthMax))
		good_list.append(df[good])
		bad_list.append(df[~good])
	dfGood = pd.concat(good_list)
	dfBad = pd.concat(bad_list)

	utils_logger.info("Filtering complete. Elapsed time: {} seconds".format(
		time.time() - make_region_lists_start))
	return (dfGood, dfBad)


def output_bed(outBed, *regionDfs):
	"""
	Concatenate and merges dataframes into an output bed file

	Parameters
	----------

	outBed : str
		The full path to and name of the output bed file
	*regionDfs
		Variable length list of dataframes to be included

	Returns
	-------

	int
		0

	"""
	dfComb = pd.concat(regionDfs)
	regionList = dfComb.ix[:, "chrom":"stop"].values.tolist()
	merge = pybedtools.BedTool(regionList).sort().merge()
	with open(outBed, 'w') as output:
		output.write(str(merge))
	return 0


def output_bed_no_merge(outBed, *regionDfs):
	"""
	Concatenate dataframes into an output bed file. This will preserve all
	columns after the first three as well.

	Parameters
	----------

	outBed : str
		The full path to and name of the output bed file
	*regionDfs
		Variable length list of dataframes to be included

	Returns
	-------

	int
		0

	"""
	dfComb = pd.concat(regionDfs)
	regionList = dfComb.ix[:,:].values.tolist()
	regionList_str = [[str(y) for y in x] for x in regionList]
	sorted = pybedtools.BedTool(regionList_str).sort()
	with open(outBed, 'w') as output:
		output.write(str(sorted))
	return 0


def chromosome_wide_plot(
	chrom, positions, y_value, measure_name, sampleID, output_prefix,
	MarkerSize, MarkerAlpha, Xlim, Ylim, x_scale=1000000):
	"""
	Plots values across a chromosome, where the x axis is the position along the
	chromosome and the Y axis is the value of the measure of interest.

	Parameters
	----------

	chrom : str
		Name of the chromosome
	positions : numpy array
		Genomic coordinates
	y_value : numpy	array
		The values of the measure of interest
	measure_name : str
		The name of the measure of interest (y axis title)
	sampleID : str
		The name of the sample
	output_prefix : str
		Full path to and prefix of desired output plot
	MarkerSize : float
		Size in points^2
	MarkerAlpha : float
		Transparency (0 to 1)
	Xlim : float
		Maximum X value
	Ylim : float
		Maximum Y value
	x_scale : int
		Divide all x values (including Xlim) by this value. Default is 1000000 (1MB)

	Returns
	-------

	int
		0

	"""
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
		positions, y_value, c=Color, alpha=MarkerAlpha, s=MarkerSize, lw=0)
	axes.set_xlim(0, (Xlim / float(x_scale)))
	axes.set_ylim(0, Ylim)
	axes.set_title("%s - %s" % (sampleID, chrom))
	if x_scale == 1000000:
		scale_label = "(MB)"
	elif x_scale == 1000:
		scale_label = "(KB)"
	elif x_scale == 1:
		scale_label = "(BP)"
	else:
		scale_label = "(divided by {})".formatt(x_scale)
	axes.set_xlabel("Chromosomal Position {}".format(scale_label))
	axes.set_ylabel(measure_name)
	plt.savefig("{}_{}_{}_GenomicScatter.pdf".format(
		output_prefix, chrom, measure_name), transparent=True)
	# plt.savefig("{}_{}_{}_GenomicScatter.svg".format(
	# 	output_prefix, chrom, measure_name))
	# plt.savefig("{}_{}_{}_GenomicScatter.png".format(
	# 	output_prefix, chrom, measure_name))
	plt.close(fig)
	return 0


def hist_array(chrom, value_array, measure_name, sampleID, output_prefix):
	"""
	Plots a histogram of an array of values of interest. Intended for mapq and
	depth, but generalizeable.  Separate function from variants.hist_read_balance
	because that function eliminates fixed variants, while this function will
	plot all values.

	Parameters
	----------

	chrom : str
		Name of the chromosome
	value_array : numpy array
		Read balance values
	measure_name : str
		The name of the measure of interest (y axis title)
	sampleID : str
		Sample name or id to include in the plot title
	output_prefix : str
		Desired prefix (including full path) of the output files

	Returns
	-------

	int
		0 if plotting successful, 1 otherwise.

	"""
	if len(value_array) < 1:
		utils_logger.info(
			"No {} values on {} to plot histogram. Skipping.".format(
				measure_name, chrom))
		return 1
	else:
		value_array = value_array[~np.isnan(value_array)]
		if "x" in chrom.lower():
			Color = "green"
		elif "y" in chrom.lower():
			Color = "blue"
		else:
			Color = "red"
		fig = plt.figure(figsize=(8, 8))
		axes = fig.add_subplot(111)
		axes.set_title("{} - {}".format(sampleID, chrom))
		axes.set_xlabel("{}".format(measure_name))
		axes.set_ylabel("Frequency")
		axes.hist(value_array, bins=50, color=Color)
		plt.savefig("{}_{}_{}_Hist.pdf".format(
			output_prefix, chrom, measure_name), transparent=True)
		# plt.savefig("{}_{}_{}_Hist.svg".format(output_prefix, chrom, measure_name))
		# plt.savefig("{}_{}_{}_Hist.png".format(output_prefix, chrom, measure_name))
		plt.close(fig)
		utils_logger.info(
			"{} histogram of {} complete.".format(measure_name, chrom))
		return 0


def plot_depth_mapq(
	window_df, output_prefix, sampleID, chrom_length, MarkerSize,
	MarkerAlpha, x_scale=1000000):
	"""
	Creates histograms and genome-wide plots of various metrics.

	Parameters
	----------

	window_df : pandas dataframe
		Columns must include chrom, start, depth, and mapq (at least)
	output_prefix : str
		Path and prefix of output files to create
	sampleID : str
		Sample ID
	chrom_length: int
		Length of chromosome
	x_scale : int
		Divide all x values (including Xlim) by this value for chromosome_wide_plot.
		Default is 1000000 (1MB)

	Returns
	-------

	int
		0

	"""
	chromosome = window_df["chrom"][1]

	# Create genome-wide plots based on window means
	if window_df is not None:
		# depth plot
		chromosome_wide_plot(
			chromosome, window_df["start"].values, window_df["depth"].values,
			"Depth", sampleID, output_prefix,
			MarkerSize, MarkerAlpha,
			chrom_length, 100, x_scale)
		hist_array(
			chromosome, window_df["depth"], "Depth", sampleID, output_prefix)

		# mapping quality plot
		chromosome_wide_plot(
			chromosome, window_df["start"].values, window_df["mapq"].values,
			"Mapq", sampleID, output_prefix,
			MarkerSize, MarkerAlpha, chrom_length, 80, x_scale)
		hist_array(
			chromosome, window_df["mapq"], "Mapq", sampleID, output_prefix)

	return 0


def before_after_plot(
	chrom, positions, values_before, values_after, measure_name, sampleID,
	output_prefix, MarkerSize, MarkerAlpha, Xlim,
	YMin="auto", YMax="auto", x_scale=1000000, Color="black"):
	"""
	Plots difference between before/after values (after minus before) across
	a chromosome.

	Parameters
	----------

	chrom : str
		Name of the chromosome
	positions : numpy array
		Genomic coordinates
	values_before : numpy array
		The values of the measure of interest in the "before" condidtion
	values_after : numpy array
		The values of the measure of interest in the "after" condidtion
	measure_name : str
		The name of the measure of interest (for y-axis title)
	sampleID : str
		The name of the sample
	output_prefix : str
		Full path to and prefix of desired output plot
	MarkerSize : float
		Size in points^2
	MarkerAlpha : float
		Transparency (0 to 1)
	Xlim : float
		Maximum X value
	YMin : str, int, or float
		If "auto", will allow matplotlib to automatically determine limit. Otherwise,
		will set the y axis minimum to the value provided (int or float)
	YMax : str, int, or float
		If "auto", will allow matplotlib to automatically determine limit. Otherwise,
		will set the y axis maximum to the value provided (int or float)
	x_scale : int
		Divide all x values (including Xlim) by this value. Default is 1000000 (1MB)
	Color : str
		Color to use for points. See matplotlib documentation for acceptable options

	Returns
	-------

	int
		0 if plotting successful, 1 otherwise
	"""
	# Check that lengths are identical
	if len(values_before) != len(values_after):
		utils_logger.error(
			"Error. values_before and values_after must have the same length")
		return 1
	else:
		value_array = np.nan_to_num(values_after) - np.nan_to_num(values_before)
		fig = plt.figure(figsize=(15, 5))
		axes = fig.add_subplot(111)
		positions = np.divide(positions, float(x_scale))
		axes.scatter(
			positions, value_array, c=Color, alpha=MarkerAlpha, s=MarkerSize, lw=0)
		axes.set_xlim(0, (Xlim / float(x_scale)))
		if YMin != "auto":
			if "." in YMin:
				YMin = float(YMin)
			else:
				YMin = int(YMin)
			if YMax != "auto":
				if "." in YMax:
					YMax = float(YMax)
				else:
					YMax = int(YMax)
				axes.ylim((YMin, YMax))
			else:
				axes.ylim(ymin=YMin)
		elif YMax != "auto":
				if "." in YMax:
					YMax = float(YMax)
				else:
					YMax = int(YMax)
				axes.ylim(ymax=YMax)
		axes.set_title("%s - %s" % (sampleID, chrom))
		if x_scale == 1000000:
			scale_label = "(MB)"
		elif x_scale == 1000:
			scale_label = "(KB)"
		elif x_scale == 1:
			scale_label = "(BP)"
		else:
			scale_label = "(divided by {})".formatt(x_scale)
		axes.set_xlabel("Chromosomal Position {}".format(scale_label))
		axes.set_ylabel("Difference in {}".format(measure_name))
		plt.savefig("{}_{}_{}_BeforeAfterScatter.pdf".format(
			output_prefix, chrom, measure_name), transparent=True)
		# plt.savefig("{}_{}_{}_BeforeAfterScatter.svg".format(
		# 	output_prefix, chrom, measure_name))
		# plt.savefig("{}_{}_{}_BeforeAfterScatter.png".format(
		# 	output_prefix, chrom, measure_name))
		plt.close(fig)
		return 0
