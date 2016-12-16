# Part of XYalign
# Collection of miscellaneous functions

from __future__ import division
from __future__ import print_function
import logging
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


def chromosome_bed(bamfile_obj, output_file, chromosome_list):
	"""
	Takes list of chromosomes and outputs a bed file with the length of each
	chromosome on each line (e.g., chr1    0   247249719).

	args:
		bamfile_obj: a BamFile object for calculating length
		output_file: name of (including full path to) desired output file
		chromosome_list: list of chromosome/scaffolds to include

	returns:
		output_file
	"""
	c_bed_start = time.time()
	utils_logger.info("Creating bed file with chromosome lengths for {}".format(
		" ".join(chromosome_list)))
	with open(output_file, "w") as f:
		for i in chromosome_list:
			try:
				l = bamfile_obj.get_chrom_length(i)
				f.write("{}\t{}\t{}\n".format(i, "0", l))
			except:
				utils_logger.error(
					"Error finding chromosome length in bam file {} "
					"(for bed file)".format(bamfile_obj.filepath))
				sys.exit(1)
	utils_logger.info("Bed file ({}) created. Elapsed time: {} seconds".format(
		output_file, time.time() - c_bed_start))
	return output_file


def check_bam_fasta_compatibility(bam_object, fasta_object):
	"""
	Takes as input BamFile() and RefFasta() objects and checks to see if
	sequence names and lengths are equivalent (i.e., if it is likely that the
	bam file was generated using the fasta in question).

	Returns:
		True if sequence names and lengths match.
		False if names or lengths do not match.
	"""
	utils_logger.info(
		"Checking compatibility of {} and {}".format(bam_object.filepath,
			fasta_object.filepath))

	b_names = bam_object.chromosome_names()
	f_names = fasta_object.chromosome_names()
	b_lengths = bam_object.chromosome_lengths()
	f_lengths = fasta_object.chromosome_lengths()

	if b_lengths == f_lengths:
		if b_names == f_names():
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
				bam_object.filepath, " ".join(zip(b_names, b_lenghts)),
				fasta_object.filepath, " ".join(zip(f_names, f_lenghts))))
		return False


def merge_bed_files(output_file, *bed_files):
	"""
	This function simply takes an output_file (full path to desired output file)
	and an arbitrary number of external bed files (including full path),
	and merges the bed files into the output_file

	args:
		output_file: full path to and name of desired output bed file
		*bed_files: arbitrary number of external bed files (include full path)
	returns:
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


def make_region_lists(depthAndMapqDf, mapqCutoff, depth_thresh):
	"""
	Filters a pandas dataframe for mapq and depth

	depthAndMapqDf is a dataframe with 'depth' and 'mapq' columns
	mapqCutoff is the minimum mapq for a window to be considered high quality
	depth_thresh is the factor to use in filtering regions based on depth:
		Li (2014) recommends:
			mean_depth +- (depth_thresh * (depth_mean ** 0.5)),
				where depth_thresh is 3 or 4.

	Returns:
		A tuple containing two dataframes (passing, failing)
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


def output_bed(outBed, *regionDfs):
	"""
	Takes a list of dataframes to concatenate and merge into an output bed file

	outBed is the full path to and name of the output bed file
	*regionDfs is an arbitrary number of dataframes to be included

	Returns:
		Nothing
	"""
	dfComb = pd.concat(regionDfs)
	regionList = dfComb.ix[:, "chrom":"stop"].values.tolist()
	merge = pybedtools.BedTool(regionList).sort().merge()
	with open(outBed, 'w') as output:
		output.write(str(merge))
	pass


def chromosome_wide_plot(
	chrom, positions, y_value, measure_name, sampleID, output_prefix,
	MarkerSize, MarkerAlpha, Xlim, Ylim):
	"""
	Plots values across a chromosome, where the x axis is the position along the
	chromosome and the Y axis is the value of the measure of interest.

	positions is an array of coordinates
	y_value is an array of the values of the measure of interest
	measure_name is the name of the measure of interest (y axis title)
	chromosome is the name of the chromosome being plotted
	sampleID is the name of the sample
	MarkerSize is the size in points^2
	MarkerAlpha is the transparency (0 to 1)
	Xlim is the maximum X value
	Ylim is the maximum Y value

	Returns:
		Nothing
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
	pass


def plot_depth_mapq(
	data_dict, output_prefix, sampleID, chrom_length, MarkerSize, MarkerAlpha):
	"""
	Takes a dictionary (output from traverseBam) and outputs histograms and
	genome-wide plots of various metrics.
	Args:
		data_dict: Dictionary of pandas data frames
		output_prefix: Path and prefix of output files to create
		sampleID: name/ID of sample
		chrom_length: length of chromosome
	Returns:
		Nothing
	"""

	window_df = None if "windows" not in data_dict else data_dict[
		"windows"]
	depth_hist = None if "depth_freq" not in data_dict else data_dict[
		"depth_freq"]
	readbal_hist = None if "readbal_freq" not in data_dict else data_dict[
		"readbal_freq"]
	mapq_hist = None if "mapq_freq" not in data_dict else data_dict[
		"mapq_freq"]

	chromosome = window_df["chrom"][1]

	# Create genome-wide plots based on window means
	if window_df is not None:
		# depth plot
		# depth_genome_plot_path = os.path.join(
		# 	output_dir, "depth_windows." + chromosome + ".png")
		# depth_genome_plot = sns.lmplot(
		# 	'start', 'depth', data=window_df, fit_reg=False,
		# 	scatter_kws={'alpha': 0.3})
		# depth_genome_plot.savefig(depth_genome_plot_path)
		chromosome_wide_plot(
			chromosome, window_df["start"].values, window_df["depth"].values,
			"Depth", sampleID, output_prefix,
			MarkerSize, MarkerAlpha,
			chrom_length, 100)

		# mapping quality plot
		# mapq_genome_plot_path = os.path.join(
		# 	output_dir, "mapq_windows." + chrom + ".png")
		# mapq_genome_plot = sns.lmplot(
		# 	'start', 'mapq', data=window_df, fit_reg=False)
		# mapq_genome_plot.savefig(mapq_genome_plot_path)
		chromosome_wide_plot(
			chromosome, window_df["start"].values, window_df["mapq"].values,
			"Mapq", sampleID, output_prefix,
			MarkerSize, MarkerAlpha, chrom_length, 80)

	# Create histograms
	# TODO: update filenames dynamically like window_df above
	# TODO: update Count column name
	if depth_hist is not None:
		depth_bar_plot = sns.countplot(
			x='depth', y='Count', data=depth_hist)
		depth_bar_plot.savefig("depth_hist.png")
	if readbal_hist is not None:
		balance_bar_plot = sns.countplot(
			x='ReadBalance', y='Count', data=readbal_hist)
		balance_bar_plot.savefig("readbalance_hist.png")
	if mapq_hist is not None:
		mapq_bar_plot = sns.countplot(
			x='Mapq', y='Count', data=mapq_hist)
		mapq_bar_plot.savefig("mapq_hist.png")
	pass
