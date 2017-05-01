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
import pysam
import sys
import time
# Matplotlib needs to be called in this way to set the display variable
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt


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
		sys.exit(1)
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


def chromosome_bed(bamfile_obj, output_file, chromosome_list):
	"""
	Takes list of chromosomes, uses a BamFile() object to find chromosome length,
	and outputs a bed file with the length of each chromosome on each line
	(e.g., chr1    0   247249719).

	Parameters
	----------

	bamfile_obj : BamFile() object
	output_file : str
		Name of (including full path to) desired output file
	chromosome_list : list
		Chromosome/scaffolds to include

	Returns
	-------

	str
		output_file

	Raises
	------

	RuntimeError
		If chromocomse name is not in bam header.

	"""
	c_bed_start = time.time()
	utils_logger.info("Creating bed file with chromosome lengths for {}".format(
		" ".join(chromosome_list)))
	with open(output_file, "w") as f:
		for i in chromosome_list:
			try:
				lengths = bamfile_obj.get_chrom_length(i)
				f.write("{}\t{}\t{}\n".format(i, "0", lengths))
			except:
				utils_logger.error(
					"Error finding chromosome length in bam file {} "
					"(for bed file)".format(bamfile_obj.filepath))
				logging.shutdown()
				raise RuntimeError(
					"Error finding chromosome length in bam file {}.  Check "
					"chromosome names and bam header.".format(
						bamfile_obj.filepath))
	utils_logger.info(
		"Bed file ({}) created. Elapsed time: {} seconds".format(
			output_file, time.time() - c_bed_start))
	return output_file


def check_chrom_in_bam(bam_object, chromosome_list):
	"""
	Checks to see if all chromosomes in chromosome_list are in bam file

	Parameters
	----------

	bam_object : BamFile() object
	chromosome_list : list
			Chromosomes/scaffolds to check

	Returns
	-------

	list
		List of chromosomes not in bam file
	"""
	utils_logger.info(
		"Checking to ensure all chromosomes are found in {}".format(bam_object.filepath))
	bam_chroms = bam_object.chromosome_names()

	missing_list = []
	for c in chromosome_list:
		if c not in bam_chroms:
			missing_list.append(c)

	if len(missing_list) > 0:
		utils_logger.info(
			"The following chromosomes were not found in {}: {}".format(
				bam_object.filepath, ",".join(missing_list)))
	else:
		utils_logger.info(
			"All chromosomes present in {}".format(bam_object.filepath))
	return missing_list


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


def make_region_lists_genome_filters(depthAndMapqDf, mapqCutoff, depth_thresh):
	"""
	Filters a pandas dataframe for mapq and depth based on using all values
	from across the entire genome

	Parameters
	----------

	depthAndMapqDf : pandas dataframe
		Must have 'depth' and 'mapq' columns
	mapqCutoff : int
		The minimum mapq for a window to be considered high quality
	depth_thresh : float
		Factor to use in filtering regions based on depth. Li (2014) recommends:
		mean_depth +- (depth_thresh * (depth_mean ** 0.5)), where depth_thresh
		is 3 or 4.

	Returns
	-------

	tuple
		(passing dataframe, failing dataframe)
	"""
	make_region_lists_start = time.time()
	depth_mean = depthAndMapqDf["depth"].mean()
	depth_sd = depthAndMapqDf["depth"].std()

	depthMin = depth_mean - (depth_thresh * (depth_mean ** 0.5))
	depthMax = depth_mean + (depth_thresh * (depth_mean ** 0.5))

	utils_logger.info(
		"Filtering dataframe for mapq (MAPQ >= mapqCutoff) "
		"and depth (between depthMin and depthMax)")

	good = (
		(depthAndMapqDf.mapq >= mapqCutoff) &
		(depthAndMapqDf.depth > depthMin) &
		(depthAndMapqDf.depth < depthMax))
	dfGood = depthAndMapqDf[good]
	dfBad = depthAndMapqDf[~good]

	utils_logger.info("Filtering complete. Elapsed time: {} seconds".format(
		time.time() - make_region_lists_start))
	return (dfGood, dfBad)


def make_region_lists_chromosome_filters(
	depthAndMapqDf, mapqCutoff, depth_thresh):
	"""
	Filters a pandas dataframe for mapq and depth based on thresholds calculated
	per chromosome

	Parameters
	----------

	depthAndMapqDf : pandas dataframe
		Must have 'depth' and 'mapq' columns
	mapqCutoff : int
		The minimum mapq for a window to be considered high quality
	depth_thresh : float
		Factor to use in filtering regions based on depth. Li (2014) recommends:
		mean_depth +- (depth_thresh * (depth_mean ** 0.5)), where depth_thresh
		is 3 or 4.

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
		depth_mean = depthAndMapqDf["depth"].mean()
		depth_sd = depthAndMapqDf["depth"].std()

		depthMin = depth_mean - (depth_thresh * (depth_mean ** 0.5))
		depthMax = depth_mean + (depth_thresh * (depth_mean ** 0.5))

		utils_logger.info(
			"Filtering chromosome {} for mapq (MAPQ >= mapqCutoff) "
			"and depth (between depthMin and depthMax)".format(i))

		good = (
			(depthAndMapqDf.chrom == i) &
			(depthAndMapqDf.mapq >= mapqCutoff) &
			(depthAndMapqDf.depth > depthMin) &
			(depthAndMapqDf.depth < depthMax))
		good_list.append(depthAndMapqDf[good])
		bad_list.append(depthAndMapqDf[~good])
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


def chromosome_wide_plot(
	chrom, positions, y_value, measure_name, sampleID, output_prefix,
	MarkerSize, MarkerAlpha, Xlim, Ylim):
	"""
	Plots values across a chromosome, where the x axis is the position along the
	chromosome and the Y axis is the value of the measure of interest.

	Parameters
	----------

	positions : numpy array
		Genomic coordinates
	y_value : numpy	array
		The values of the measure of interest
	measure_name : str
		The name of the measure of interest (y axis title)
	chromosome : str
		The name of the chromosome being plotted
	sampleID : str
		The name of the sample
	MarkerSize : float
		Size in points^2
	MarkerAlpha : float
		Transparency (0 to 1)
	Xlim : float
		Maximum X value
	Ylim : float
		Maximum Y value

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
	axes.scatter(
		positions, y_value, c=Color, alpha=MarkerAlpha, s=MarkerSize, lw=0)
	axes.set_xlim(0, Xlim)
	axes.set_ylim(0, Ylim)
	axes.set_title("%s - %s" % (sampleID, chrom))
	axes.set_xlabel("Chromosomal Position")
	axes.set_ylabel(measure_name)
	plt.savefig("{}_{}_{}_GenomicScatter.svg".format(
		output_prefix, chrom, measure_name))
	plt.savefig("{}_{}_{}_GenomicScatter.png".format(
		output_prefix, chrom, measure_name))
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
		plt.savefig("{}_{}_{}_Hist.svg".format(output_prefix, chrom, measure_name))
		plt.savefig("{}_{}_{}_Hist.png".format(output_prefix, chrom, measure_name))
		plt.close(fig)
		variants_logger.info(
			"{} histogram of {} complete.".format(measure_name, chrom))
		return 0


def plot_depth_mapq(
	data_dict, output_prefix, sampleID, chrom_length, MarkerSize, MarkerAlpha):
	"""
	Creates histograms and genome-wide plots of various metrics.

	Note that the odd import format (a dictionary, in which "windows" is the
	key whose value is the pandas dataframe of interest) is a carryover from
	previous versions of this function that required a variety of pandas
	dataframes.

	Parameters
	----------

	data_dict : dict
		Key must be 'windows' value is a pandas data frame
	output_prefix : str
		Path and prefix of output files to create
	sampleID : str
		Sample ID
	chrom_length: int
		Length of chromosome

	Returns
	-------

	int
		0

	"""

	window_df = None if "windows" not in data_dict else data_dict[
		"windows"]

	chromosome = window_df["chrom"][1]

	# Create genome-wide plots based on window means
	if window_df is not None:
		# depth plot
		chromosome_wide_plot(
			chromosome, window_df["start"].values, window_df["depth"].values,
			"Depth", sampleID, output_prefix,
			MarkerSize, MarkerAlpha,
			chrom_length, 100)
		hist_array(
			chromosome, window_df["depth"], "Depth", sampleID, output_prefix)

		# mapping quality plot
		chromosome_wide_plot(
			chromosome, window_df["start"].values, window_df["mapq"].values,
			"Mapq", sampleID, output_prefix,
			MarkerSize, MarkerAlpha, chrom_length, 80)
		hist_array(
			chromosome, window_df["mapq"], "Mapq", sampleID, output_prefix)

	return 0
