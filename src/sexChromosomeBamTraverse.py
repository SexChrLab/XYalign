from __future__ import division
import os
import argparse
from collections import Counter, defaultdict
from itertools import chain
import numpy as np
import pandas as pd
import pybedtools
import pysam
import seaborn as sns


def main():
	"""Main program"""
	args = parse_args()
	samfile = pysam.AlignmentFile(args.bam, "rb")
	pass_df = []
	fail_df = []
	for chromosome in args.chromosomes:
		data = traverse_bam_fetch(samfile, chromosome, args.window_size)
		tup = makeRegionLists(data["windows"], args.mapq_cutoff, args.depth_filter)
		pass_df.append(tup[0])
		fail_df.append(tup[1])
		plot_data(data, args.output_dir)
	outputBed(os.path.join(args.output_dir, args.high_quality_bed), *pass_df)
	outputBed(os.path.join(args.output_dir, args.low_quality_bed), *fail_df)


def parse_args():
	"""Parse and validation command-line arguments."""
	# Parse arguments
	parser = argparse.ArgumentParser(description="Add description")
	parser.add_argument("--bam", required=True, help="Input bam file")
	parser.add_argument("--chromosomes", "-c", nargs="+", default=["chrX", "chrY", "chr19"],
						help="Chromosomes to analyze.")
	parser.add_argument("--window_size", "-w", type=int, default=50000,
						help="Window size (integer) for sliding window calculations.")
	parser.add_argument("--depth_filter", "-df", type=float, default=4.0,
						help="Filter for depth (f), where the threshold used is mean_depth +- (f * square_root(mean_depth)).")
	parser.add_argument("--mapq_cutoff", "-mq", type=int, default=20,
						help="Minimum mean mapq threshold for a window to be considered high quality.")
	parser.add_argument("--high_quality_bed", "-hq", default="highquality.bed",
						help="Name of output file for high quality regions.")
	parser.add_argument("--low_quality_bed", "-lq", default="lowquality.bed",
						help="Name of output file for high quality regions.")
	parser.add_argument("--min_depth", "-d", type=int, default=2,
						help="Minimum depth for a site to be considered.")
	parser.add_argument("--min_minor_depth", "-ar", type=int, default=1,
						help="Minimum number of reads supporting the alternate allele to be considered.")
	parser.add_argument("--min_minor_fraction", "-af", type=float, default=0.1,
						help="Minimum fraction of reads supporting the alternate allele to be considered.")
	parser.add_argument("--output_dir", "-o", help="Output directory")
	args = parser.parse_args()
	# Validate arguments
	if not os.path.exists(args.output_dir):
		os.makedirs(args.output_dir)
	# Return arguments namespace
	return args


def get_length(samfile, chrom):
	"""Extract chromosome length from BAM header.

	Args:
		samfile: pysam AlignmentFile object
		chromosome: Chromosome name (string)

	Returns:
		Length (int)
	"""
	lengths = dict(zip(samfile.references, samfile.lengths))
	return lengths[chrom]


def reset_lists(n, size):
	"""Initialize `n` lists of size `size`."""
	return [[0] * size for i in range(n)]


def reset_counters(n):
	"""Initialize `n` lists of size `size`."""
	return [Counter() for i in range(n)]

def traverse_bam_fetch(samfile, chrom, window_size):
	"""Analyze the `samfile` BAM file for various metrics and statistics.

	Currently, this function looks at the following metrics across genomic windows:
	- Read depth
	- Mapping quality

	The average of each metric will be calculated for each window of
	size `window_size` and stored altogether in a pandas data frame.

	Returns a dictionary of pandas data frames with the following key:
	- windows: The averages for each metric for each window
	"""
	chr_len = get_length(samfile, chrom)
	num_windows = chr_len // window_size + 1
	if chr_len % num_windows == 0:
		last_window_len = window_size
	else:
		last_window_len = chr_len % num_windows
	
	window_id = 0
	win_pos = 0
	
	chr_list = [chrom] * num_windows
	start_list = []
	stop_list = []
	depth_list = []
	mapq_list = []
	
	start = 0
	end = window_size
	for window in range(0, num_windows):
		mapq = []
		num_reads = 0
		total_read_length = 0
		for read in samfile.fetch(chrom, start, end):
			num_reads += 1
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
		print "{} out of {} windows processed on {}".format(window_id, num_windows, chrom)

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


def makeRegionLists(depthAndMapqDf, mapqCutoff, sd_thresh):
	"""
	(pandas.core.frame.DataFrame, int, float) -> (list, list)
	return two lists of regions (keepList, excludeList) based on cutoffs for depth and mapq
	"""
	depth_mean = depthAndMapqDf["depth"].mean()
	depth_sd = depthAndMapqDf["depth"].std()

	depthMin = depth_mean - (sd_thresh * (depth_sd ** 0.5))
	depthMax = depth_mean + (sd_thresh * (depth_sd ** 0.5))

	good = (depthAndMapqDf.mapq > mapqCutoff) & (depthAndMapqDf.depth > depthMin) & (depthAndMapqDf.depth < depthMax)
	dfGood = depthAndMapqDf[good]
	dfBad = depthAndMapqDf[~good]

	return (dfGood, dfBad)


def outputBed(outBed, *regionDfs):
    '''
    (list, list, str) -> bedtoolsObject
    Take two sorted lists.  Each list is a list of tuples (chrom[str], start[int], end[int])
    Return a pybedtools object and output a bed file.
    '''
	dfComb = pd.concat(regionDfs)
	regionList = dfComb.ix[:, "chrom":"stop"].values.tolist()
	merge = pybedtools.BedTool(regionList).sort().merge()
	with open(outBed, 'w') as output:
		output.write(str(merge))
	pass


def traverse_bam_pileup(samfile, chrom, window_size, min_depth, min_minor_depth,
					min_minor_fraction, depth_binsize, readbal_binsize):
	"""Analyze the `samfile` BAM file for various metrics and statistics.

	Currently, this function looks at the following metrics:
	- Read depth
	- Read balance (allelic fractions)
	- Mapping quality

	For each metric, a pandas data frame will be created with frequencies
	across a number of ranges (for plotting a histogram).

	The average of each metric will also be calculated for each window of
	size `window_size` and stored altogether in a pandas data frame.

	Returns a dictionary of pandas data frames with the following keys:
	- depth_freq: Frequencies of read depths
	- readbal_freq: Frequencies of read balances
	- mapq_freq: Frequencies of mapping qualities
	- windows: The averages for each metric for each window
	"""
	chr_len = get_length(samfile, chrom)
	num_windows = chr_len // window_size + 1
	depths, readbals, mapqs = reset_lists(3, window_size)
	depths_binned, readbals_binned, mapqs_binned = reset_lists(3, window_size)
	depths_counter, readbals_counter, mapqs_counter = reset_counters(3)
	window_id = 0
	win_pos = 0
	for chr_pos, column in enumerate(samfile.pileup(chrom, 0, chr_len)):
		base_counter = defaultdict(int)
		col_mapqs = []

		# Iterate over reads in the pileup column
		for read in column.pileups:
			if not read.is_del and not read.is_refskip:
				base = read.alignment.query_sequence[read.query_position]
				base_counter[base] += 1
				col_mapqs.append(read.alignment.mapping_quality)

		# Identify major and minor allele counts
		base_counts = sorted(base_counter.values(), reverse=True)
		num_major = 0 if len(base_counts) < 1 else base_counts[0]
		num_minor = 0 if len(base_counts) < 2 else base_counts[1]
		total_depth = num_major + num_minor

		# Apply threshold and track metrics
		if total_depth >= min_depth:
			depths[win_pos] = total_depth
			mapqs[win_pos] = sum(col_mapqs) / len(col_mapqs)  # Mean
			allele_fraction = num_minor / total_depth
			if num_minor >= min_minor_depth and allele_fraction >= min_minor_fraction:
				readbals[win_pos] = allele_fraction

		# Increment counters
		win_pos += 1
		chr_pos += 1

		# Are we at the end of a window?
		if win_pos == window_size - 1:
			# Bin values
			depths_binned[window_id] = np.mean(np.asarray(depths))
			readbals_binned[window_id] = np.mean(np.asarray(readbals))
			mapqs_binned[window_id] = np.mean(np.asarray(mapqs))
			# Increment counters
			depths_counter.update(depths)
			readbals_counter.update(readbals)
			mapqs_counter.update(mapqs)
			# Reset and increment
			num_passed_depth = len(depths)
			num_passed_minor = len(readbals)
			depths, readbals, mapqs = reset_lists(3, window_size)
			win_pos = 0
			window_id += 1
			# Print progress
			print "{} out of {} positions passed depth threshold.".format(num_passed_depth, window_size)
			print "{} out of {} positions passed minor allele thresholds.".format(num_passed_minor, window_size)
			print "{} out of {} windows processed.".format(window_id, num_windows)

	# Convert data into pandas data frames
	windows_df = pd.DataFrame({
		"depth": np.asarray(depths_binned),
		"readbal": np.asarray(readbals_binned),
		"mapq": np.asarray(mapqs_binned)
	})

	mapq_freq_df = pd.DataFrame(dict(zip(np.asarray(mapqs_counter.keys()),
											np.asarray(mapqs_counter.values()))))

	# TODO: Create depth_freq_df and readbal_freq_df data frames

	results = {
		"windows": windows_df,
		"mapq_freq": mapq_freq_df
	}
	return results

def chromsomeWidePlot(positions,y_value,measure_name,chromosome,sampleID,MarkerSize,MarkerAlpha,Xlim,Ylim):
    '''
    Plots values across a chromosome, where the x axis is the position along the
    chromosome and the Y axis is the value of the measure of interest.
    
    positions is an array of coordinates 
    y_value is an array of the values of the measure of interest
    measure_name is the name of the measure of interest (y axis title)
    chromsome is the name of the chromosome being plotted
    sampleID is the name of the sample
    MarkerSize is the size in points^2
    MarkerAlpha is the transparency (0 to 1)
    Xlim is the maximum X value
    Ylim is the maximum Y value
    '''
    
    if "X" in sampleID:
        Color="green"
    elif "Y" in sampleID:
        Color = "blue"
    else:
    	Color = "red"
    fig = plt.figure(figsize=(15,5))
    axes = fig.add_subplot(111)
    axes.scatter(positions,y_value,c=Color,alpha=MarkerAlpha,s=MarkerSize,lw=0)
    axes.set_xlim(0,Xlim)
    axes.set_ylim(0,Ylim)
    axes.set_title("%s - %s" % (sampleID, chromosome))
    axes.set_xlabel("Chromosomal Position")
    axes.set_ylabel(measure_name)
    plt.savefig("%s_%s_%s_GenomicScatter.svg" % (sampleID, chromsome, measure_name))
    plt.savefig("%s_%s_%s_GenomicScatter.png"% (sampleID, chromsome, measure_name))
	#plt.show()

def plot_data(data_dict, output_dir):
	"""
	Takes a dictionary (output from traverseBam) and outputs histograms and
	genome-wide plots of various metrics.

	Args:
		data_dict: Dictionary of pandas data frames
		output_dir: Directory where the PNG file with the plots will be stored

	Returns:
		None
	"""

	window_df = None if "windows" not in data_dict else data_dict["windows"]
	depth_hist = None if "depth_freq" not in data_dict else data_dict["depth_freq"]
	readbal_hist = None if "readbal_freq" not in data_dict else data_dict["readbal_freq"]
	mapq_hist = None if "mapq_freq" not in data_dict else data_dict["mapq_freq"]

	chrom = window_df["chrom"][1]

	# Create genome-wide plots based on window means
	if window_df is not None:
		# depth plot
		depth_genome_plot_path = os.path.join(output_dir, "depth_windows." + chrom + ".png")
		depth_genome_plot = sns.lmplot('start', 'depth', data=window_df, fit_reg=False,
										scatter_kws={'alpha': 0.3})
		depth_genome_plot.savefig(depth_genome_plot_path)
		# mapping quality plot
		mapq_genome_plot_path = os.path.join(output_dir, "mapq_windows." + chrom + ".png")
		mapq_genome_plot = sns.lmplot('start', 'mapq', data=window_df, fit_reg=False)
		mapq_genome_plot.savefig(mapq_genome_plot_path)

	# Create histograms
	# TODO: update filenames dynamically like window_df above
	# TODO: update Count column name
	if depth_hist is not None:
		depth_bar_plot = sns.countplot(x='depth', y='Count', data=depth_hist)
		depth_bar_plot.savefig("depth_hist.png")
	if readbal_hist is not None:
		balance_bar_plot = sns.countplot(x='ReadBalance', y='Count', data=readbal_hist)
		balance_bar_plot.savefig("readbalance_hist.png")
	if mapq_hist is not None:
		mapq_bar_plot = sns.countplot(x='Mapq', y='Count', data=mapq_hist)
		mapq_bar_plot.savefig("mapq_hist.png")


if __name__ == "__main__":
	main()
