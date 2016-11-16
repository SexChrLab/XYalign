from __future__ import division
from __future__ import print_function
import argparse
import os
import subprocess
import sys
import numpy as np
import pandas as pd
import pybedtools
import pysam
import time
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns


def main():

	args = parse_args()

	# Setup output paths
	fastq_path = os.path.join(args.output_dir, "fastq")
	bam_path = os.path.join(args.output_dir, "bam")
	reference_path = os.path.join(args.output_dir, "reference")
	bed_path = os.path.join(args.output_dir, "bed")
	vcf_path = os.path.join(args.output_dir, "vcf")
	plots_path = os.path.join(args.output_dir, "plots")
	results_path = os.path.join(args.output_dir, "results")

	depth_mapq_prefix_noprocessing = os.path.join(
		plots_path, "{}_noprocessing".format(args.sample_id))
	depth_mapq_prefix_postprocessing = os.path.join(
		plots_path, "{}_postprocessing".format(args.sample_id))
	if args.high_quality_bed_out is not None:
		# high_prefix = args.high_quality_bed_out
		print(
			"--high_quality_bed_out is currently unsupported.  Please remove "
			"this flag")
		sys.exit(1)
	else:
		high_prefix = "{}_highquality_preprocessing".format(args.sample_id)
	output_bed_high = os.path.join(
		bed_path, "{}.bed".format(high_prefix))
	if args.low_quality_bed_out is not None:
		# low_prefix = args.low_quality_bed_out
		print(
			"--low_quality_bed_out is currently unsupported.  Please remove "
			"this flag")
	else:
		low_prefix = "{}_lowquality_preprocessing".format(args.sample_id)
	output_bed_low = os.path.join(
		bed_path, "{}.bed".format(low_prefix))

	bam_analysis_start = time.time()
	if args.bam is not None:
		print("Beginning bam analyses on {}\n".format(args.bam))
		samfile = pysam.AlignmentFile(args.bam, "rb")
	else:
		print("Beginning cram analyses on {}\n".format(args.cram))
		samfile = pysam.AlignmentFile(args.cram, "rc")
	pass_df = []
	fail_df = []
	for chromosome in args.chromosomes:
		data = traverse_bam_fetch(samfile, chromosome, args.window_size)
		tup = make_region_lists(
			data["windows"], args.mapq_cutoff, args.depth_filter)
		pass_df.append(tup[0])
		fail_df.append(tup[1])
		plot_depth_mapq(
			data, depth_mapq_prefix_noprocessing, args.sample_id,
			get_length(samfile, chromosome), args.marker_size,
			args.marker_transparency)
	output_bed(output_bed_high, *pass_df)
	output_bed(output_bed_low, *fail_df)
	bam_analysis_end = time.time()
	print("Bam-cram analyses complete. Elapsed time: {} seconds\n".format(
		bam_analysis_end - bam_analysis_start))


def parse_args():

	"""Parse command-line arguments"""
	parser = argparse.ArgumentParser(description="XYalign")

	parser.add_argument(
		"--ref", required=True,
		help="REQUIRED. Path to reference sequence (including file name).")

	parser.add_argument(
		"--output_dir", "-o",
		help="REQUIRED. Output directory. XYalign will create a directory "
		"structure within this directory")

	parser.add_argument(
		"--chromosomes", "-c", nargs="+", default=["chrX", "chrY", "chr19"],
		help="Chromosomes to analyze (names must match reference exactly). "
		"Defaults to chr19, chrX, chrY.")

	parser.add_argument(
		"--sample_id", "-id", default="sample",
		help="Name/ID of sample - for use in plot titles and file naming. "
		"Default is sample")

	# Bam Analysis Flags
	parser.add_argument(
		"--window_size", "-w", type=int, default=50000,
		help="Window size (integer) for sliding window calculations. Default "
		"is 50000.")

	parser.add_argument(
		"--mapq_cutoff", "-mq", type=int, default=20,
		help="Minimum mean mapq threshold for a window to be "
		"considered high quality. Default is 20.")

	parser.add_argument(
		"--depth_filter", "-df", type=float, default=4.0,
		help="Filter for depth (f), where the threshold used is mean_depth +- "
		"(f * square_root(mean_depth)).  See Li 2014 (Bioinformatics 30: "
		"2843-2851) for more information.  Default is 4.")

	parser.add_argument(
		"--high_quality_bed_out", "-hq", default=None,
		help="Prefix of output file for high quality regions. Defaults to "
		"sample_id_highquality")

	parser.add_argument(
		"--low_quality_bed_out", "-lq", default=None,
		help="Prefix of output file for low quality regions. Defaults to "
		"sample_id_lowquality")

	# Plotting flags
	parser.add_argument(
		"--marker_size", type=float, default=10.0,
		help="Marker size for genome-wide plots in matplotlib. Default is 10.")

	parser.add_argument(
		"--marker_transparency", "-mt", type=float, default=0.5,
		help="Transparency of markers in genome-wide plots.  "
		"Alpha in matplotlib.  Default is 0.5")

	# Mutually exclusive group 1 - bam or cram file
	group = parser.add_mutually_exclusive_group(required=True)

	group.add_argument(
		"--bam", help="Input bam file.")

	group.add_argument(
		"--cram", help="Input cram file.")

	args = parser.parse_args()

	# Create directory structure if not already in place
	if not os.path.exists(os.path.join(args.output_dir, "bed")):
		os.makedirs(os.path.join(args.output_dir, "bed"))
	if not os.path.exists(os.path.join(args.output_dir, "plots")):
		os.makedirs(os.path.join(args.output_dir, "plots"))
	if not os.path.exists(os.path.join(args.output_dir, "logfiles")):
		os.makedirs(os.path.join(args.output_dir, "logfiles"))

	# Return arguments namespace
	return args


def get_length(bamfile, chrom):
	"""
	Extract chromosome length from BAM header.

	args:
		bamfile: pysam AlignmentFile object
			- can be bam, cram, or sam, needs to be declared
				in pysam.AlignmentFile call before passing to function
		chrom: chromosome name (string)

	returns:
		Length (int)

	"""
	lengths = dict(zip(bamfile.references, bamfile.lengths))
	return lengths[chrom]


def plot_variants_per_chrom(
	chrom_list, vcf_file, sampleID, output_prefix, qual_cutoff,
	MarkerSize, MarkerAlpha, bamfile):
	"""
	Parses a vcf file and plots read balance in separate plots
	for each chromosome in the input list

	chrom_list is the list of chromosomes to run parse_platypus_VCF and plotting
		functions on
	vcf_file is the file (including path) of platypus vcf to analyze
	sampleID is the sample name (for plot titles)
	output_prefix is the full path to and prefix of desired output plots
	qual_cutoff is the minimum (Phred) quality to consider a site in the vcf
	MarkerSize is the size of markers (matplotlib sizes) to use in the figure
	MarkerAlpha is the transparency (matplotlib values) of markers for the figure
	bamfile is the name of the corresponding bam file (used to get chromosome
		lengths only)

	Returns:
		Nothing
	"""
	for i in chrom_list:
		parse_results = parse_platypus_VCF(vcf_file, qual_cutoff, i)
		plot_read_balance(
			i, parse_results[0], parse_results[2],
			sampleID, output_prefix, MarkerSize, MarkerAlpha, bamfile)
		hist_read_balance(
			i, parse_results[2], sampleID, output_prefix)
	pass


def traverse_bam_fetch(samfile, chrom, window_size):
	"""Analyze the `samfile` BAM (or CRAM) file for various metrics.
	Currently, this function looks at the following metrics across genomic
	windows:
	- Read depth
	- Mapping quality
	The average of each metric will be calculated for each window of
	size `window_size` and stored altogether in a pandas data frame.


	samfile is a pysam AlignmentFile object
	chrom is the chromosome to analyze
	window size is the integer window size to use for sliding window analyses

	Returns:
		A dictionary of pandas data frames with the following key:
			- windows: The averages for each metric for each window
	"""
	chr_len = get_length(samfile, chrom)
	num_windows = chr_len // window_size + 1
	if chr_len % num_windows == 0:
		last_window_len = window_size
	else:
		last_window_len = chr_len % num_windows

	window_id = 0

	chr_list = [chrom] * num_windows
	start_list = []
	stop_list = []
	depth_list = []
	mapq_list = []

	start = 0
	end = window_size
	for window in range(0, num_windows):
		mapq = []
		total_read_length = 0
		for read in samfile.fetch(chrom, start, end):
			if read.is_secondary is False:
				if read.is_supplementary is False:
					total_read_length += read.infer_query_length()
					mapq.append(read.mapping_quality)
		start_list.append(start)
		stop_list.append(end)
		depth_list.append(total_read_length / window_size)
		mapq_list.append(np.mean(np.asarray(mapq)))

		window_id += 1
		if window_id == num_windows - 1:
			start += window_size
			end += last_window_len
		else:
			start += window_size
			end += window_size

		# Print progress
		print("{} out of {} windows processed on {}".format(
			window_id, num_windows, chrom))

	# Convert data into pandas data frames
	windows_df = pd.DataFrame({
		"chrom": np.asarray(chr_list),
		"start": np.asarray(start_list),
		"stop": np.asarray(stop_list),
		"depth": np.asarray(depth_list),
		"mapq": np.asarray(mapq_list)
	})[["chrom", "start", "stop", "depth", "mapq"]]

	results = {"windows": windows_df}
	return results


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
	depth_mean = depthAndMapqDf["depth"].mean()
	depth_sd = depthAndMapqDf["depth"].std()

	depthMin = depth_mean - (depth_thresh * (depth_mean ** 0.5))
	depthMax = depth_mean + (depth_thresh * (depth_mean ** 0.5))

	good = (
		(depthAndMapqDf.mapq >= mapqCutoff) &
		(depthAndMapqDf.depth > depthMin) &
		(depthAndMapqDf.depth < depthMax))
	dfGood = depthAndMapqDf[good]
	dfBad = depthAndMapqDf[~good]

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
	# plt.show()


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

if __name__ == "__main__":
	main()
