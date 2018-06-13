from __future__ import absolute_import
from __future__ import print_function
import argparse
import numpy as np
import pandas as pd
import xyalign.utils as utils
# from xyalign import utils


def parse_args():
	"""
	Parse command-line arguments
	Returns
	-------
	Parser argument namespace
	"""

	parser = argparse.ArgumentParser(
		description="This script takes two dataframes output by XYalign's "
		"BAM_ANALYSIS module and plots the difference in depth and mapq for "
		"a given chromosome. This was designed to plot these metrics before "
		"and after XYalign processing on the same sample, but could conceivably "
		"be used for other purposes.")

	parser.add_argument(
		"--before", type=str, help="Full path to dataframe 1 (will be treated) "
		"as 'before' condition.")

	parser.add_argument(
		"--after", type=str, help="Full path to dataframe 2 (will be treated) "
		"as 'after' condition.")

	parser.add_argument(
		"--color", help="Color of points to use. Consult matplotlib documentation "
		"for available options.")

	parser.add_argument(
		"--chrom", help="Chromosome to analyze. Must match name in dataframes.")

	parser.add_argument(
		"--sample_id", help="Sample ID or other identifier to be used in naming")

	parser.add_argument(
		"--output_prefix", help="Full path to and prefix of desired output plot")

	parser.add_argument(
		"--marker_size", type=float, default=10.0,
		help="Marker size in matplotlib. Default is 10.")

	parser.add_argument(
		"--marker_transparency", "-mt", type=float, default=0.5,
		help="Transparency of markers.  "
		"Alpha in matplotlib.  Default is 0.5")

	parser.add_argument(
		"--coordinate_scale", type=int, default=1000000,
		help="Divide all coordinates by this value."
		"Default is 1000000, which will plot everything in megabases.")

	parser.add_argument(
		"--y_min", default="auto", help="If 'auto', will allow matplotlib "
		"to automatically determine limit. Otherwise, will set the y axis "
		"minimum to the value provided (int or float)")

	parser.add_argument(
		"--y_max", default="auto", help="If 'auto', will allow matplotlib "
		"to automatically determine limit. Otherwise, will set the y axis "
		"maximum to the value provided (int or float)")

	parser.add_argument(
		"--x_limit", help="Max value on x axis. We recommend you use the "
		"chromosome length.")

	parser.add_argument(
		"--log_transform_depth", action="store_true", default=False,
		help="Include flag to plot the absolute value of the log of the depth "
		"difference *IN THE DIRECTION OF THE DIFFERENCE*. For exampe, "
		"if the difference is 2, this would plot that value as abs(log10(2)), "
		"while if the difference is -2, this would plot that value as "
		"-abs(log10(2)). This allows the sign of the difference to remain "
		"intact, while controlling for taking the log of negative numbers "
		"or values between 0 and 1.")

	args = parser.parse_args()

	return args


def transform_depth(numpy_array):
	"""
	Performs custom version of log transformation on a numpy array. Where each
	value is processed to be equal to:
	initial_sign * abs(log10(abs(value)))
	Parameters
	----------
	numpy_array : numpy array
		Array of values without NaNs
	Returns
	-------
	numpy array
	"""
	signs = np.sign(numpy_array)
	step1 = np.absolute(numpy_array)
	id_zeros = step1 != 0
	step2 = np.absolute(np.log10(step1, where=id_zeros))
	return signs * step2


def main():
	args = parse_args()

	df_before = pd.read_csv(
		args.before, sep="\t", header=0)
	df_after = pd.read_csv(
		args.after, sep="\t", header=0)

	df_before = df_before.loc[df_before["chrom"] == args.chrom]
	df_after = df_after.loc[df_after["chrom"] == args.chrom]

	df_before = df_before.fillna(0)
	df_after = df_after.fillna(0)

	print(
		"Mean depth difference = {}".format(
			np.mean(df_after["depth"].values - df_before["depth"].values)))
	print(
		"Mean MAPQ difference = {}".format(
			np.mean(df_after["mapq"].values - df_before["mapq"].values)))

	if args.log_transform_depth is True:
		v_before = transform_depth(df_before["depth"].values)
		v_after = transform_depth(df_after["depth"].values)
	else:
		v_before = df_before["depth"].values
		v_after = df_after["depth"].values
	return_val = utils.before_after_plot(
		chrom=args.chrom,
		positions=df_before["start"].values,
		values_before=v_before,
		values_after=v_after,
		measure_name="depth",
		sampleID=args.sample_id,
		output_prefix=args.output_prefix,
		MarkerSize=args.marker_size,
		MarkerAlpha=args.marker_transparency,
		Xlim=float(args.x_limit),
		YMin=args.y_min,
		YMax=args.y_max,
		x_scale=args.coordinate_scale,
		Color=args.color)

	return_val = utils.before_after_plot(
		chrom=args.chrom,
		positions=df_before["start"].values,
		values_before=df_before["mapq"].values,
		values_after=df_after["mapq"].values,
		measure_name="MAPQ",
		sampleID=args.sample_id,
		output_prefix=args.output_prefix,
		MarkerSize=args.marker_size,
		MarkerAlpha=args.marker_transparency,
		Xlim=float(args.x_limit),
		YMin=args.y_min,
		YMax=args.y_max,
		x_scale=args.coordinate_scale,
		Color=args.color)


if __name__ == "__main__":
	main()
