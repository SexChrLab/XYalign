from __future__ import print_function
import argparse
import numpy as np
import pandas as pd
from xyalign import utils


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

	args = parser.parse_args()

	return args


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

	return_val = utils.before_after_plot(
		chrom=args.chrom,
		positions=df_before["start"].values,
		values_before=df_before["depth"].values,
		values_after=df_after["depth"].values,
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
