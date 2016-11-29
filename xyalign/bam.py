# Part of XYalign
# Functions for calling and processing variants
from __future__ import division
from __future__ import print_function
import numpy as np
import pandas as pd


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
