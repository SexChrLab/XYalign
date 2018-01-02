from __future__ import division
from __future__ import print_function
import argparse
import csv
import numpy as np
import pandas as pd
import sys
# Matplotlib needs to be called in this way to set the display variable
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt


def parse_args():
	"""
	Parse command-line arguments

	Returns
	-------

	Parser argument namespace
	"""
	parser = argparse.ArgumentParser(
		description="This script will plot scatter plots of 'X/A' and 'Y/A' "
		"ratios from CHROM_STATS output. However, its features are not meant "
		"to be exhaustive. Rather, it is likely best used as a template for "
		"users to customize and adjust as needed.")

	parser.add_argument(
		"--input", type="str", required=True,
		help="Full path to file containing table output by CHROM_STATS")

	parser.add_argument(
		"--meta", type="str", required=True,
		help="Full path to file containing metadata table. This file should "
		"have the following columns separated by tabs: "
		"Sample NameOfVariable1 NameOfVariable2. NameOfVariable1 and NameOfVariable2 "
		"should be the names of whatever you're interested in plotting (e.g., Sex). "
		"NameOfVariable2 is optional. This script handles a max of two variables.")

	parser.add_argument(
		"--output_prefix", type="str", required=True,
		help="'Prefix' of output files. This includes full path to desired file "
		"and desired file name before suffix (suffixes will be .png and .svg).")

	parser.add_argument(
		"--exclude_suffix", type="str", default="",
		help="Text to remove from end of sample names in input file. Default is "
		"to remove nothing. Note that the sample names in the input file have "
		"to match the names in the meta file AFTER they undergo this step.")

	parser.add_argument(
		"--first_chr", type="str", required=True,
		help="Chromosome to use a numerator on X-axis. For example, if one was "
		"comparing chrX and chrY, and using chr19 to normalize, recommended "
		"values would be: --first_chr chrX --second_chr chrY --const_chr chr19.")

	parser.add_argument(
		"--second_chr", type="str", required=True,
		help="Chromosome to use a numerator on Y-axis. For example, if one was "
		"comparing chrX and chrY, and using chr19 to normalize, recommended "
		"values would be: --first_chr chrX --second_chr chrY --const_chr chr19.")

	parser.add_argument(
		"--const_chr", type="str", required=True,
		help="Chromosome to use denominator on both the X- and Y-axis. For "
		"example, if one was comparing chrX and chrY, and using chr19 to "
		"normalize, recommended values would be: --first_chr chrX "
		"--second_chr chrY --const_chr chr19.")

	parser.add_argument(
		"--var1_marker", choices=["color", "shape", "size"], default="color"
		help="Way of designating variable 1 values in plot. Choices are 'color', "
		"'shape', or 'size'. Must be used in conjunction with --var1_marker_vals. "
		"Default is 'color'.")

	parser.add_argument(
		"--var1_marker_vals", nargs="+", default=["red", "blue"],
		help="Marker values to use for variable 1 values. If --var1_marker is "
		"'color', then --var1_marker_vals should be a space-separated list "
		"of Matplotlib colors (e.g., 'red blue green'). If --var1_marker is "
		"'shape' then --var1_marker_vals should be a space-separated list "
		"of Matplotlib scatter markers (e.g., 'x o D' for x, cicle, and Diamond). "
		"Finally, if --var1_marker is 'size', then --var1_marker_vals should "
		"be a space-separated list of Matplotlib marker sizes in units of "
		"points^2 (e.g., '5 10 15'). Default is 'red blue'.")

	parser.add_argument(
		"--var2_marker", choices=["color", "shape", "size", "none"], default="none"
		help="Way of designating variable 2 values in plot. Choices are 'color', "
		"'shape', or 'size'. Must be used in conjunction with --var2_marker_vals. "
		"Default is 'none', which will only process --var1_marker.")

	parser.add_argument(
		"--var2_marker_vals", nargs="*", default=["8", "15"],
		help="Marker values to use for variable 2 values. If --var2_marker is "
		"'color', then --var2_marker_vals should be a space-separated list "
		"of Matplotlib colors (e.g., 'red blue green'). If --var2_marker is "
		"'shape' then --var2_marker_vals should be a space-separated list "
		"of Matplotlib scatter markers (e.g., 'x o D' for x, cicle, and Diamond). "
		"Finally, if --var2_marker is 'size', then --var2_marker_vals should "
		"be a space-separated list of Matplotlib marker sizes in units of "
		"points^2 (e.g., '5 10 15'). Default is 'red blue'.")

	args = parser.parse_args()

	if args.var1_marker == args.var2_marker:
		print(
			"--var1_marker (={}) and --var2_marker (={}) must be different.".format(
				args.var1_marker, args.var2_marker))
		sys.exit(1)

	return args


def main():
	args = parse_args()

	df = pd.read_csv(
		args.input, sep="\t", header=0)

	if len(args.exclude_suffix) > 0:
		col_names = df.values[0]
		col_names = [x[:-(len(args.exclude_suffix))] if args.exclude_suffix in x else x for x in col_names]
		df.loc[0] = col_names
	df2 = df.transpose()
	df2.columns = df2.iloc[0]
	df2.reindex(df2.index.drop(0))
	df2.columns = ["Sample" if x == "chrom" else x for x in df.columns]

	df2["ratiox"] = df2.apply(
		lambda row: row[args.first_chr] / row[args.const_chr], axis=1)
	df2["ratioy"] = df2.apply(
		lambda row: row[args.second_chr] / row[args.const_chr], axis=1)

	meta_df = pd.read_csv(args.meta, sep="\t", header=0)
	merged = pd.concat[df2, meta_df]

	fig = plt.figure(figsize=(10,10))
	ax = plt.subplot(111)

	# One variable
	if len(meta_df.columns) == 2:
		for idx, i in enumerate(merged[meta_df.columns[1]].unique()):
			filtered = merged.loc[merged[args.var1_marker] == i]
			if args.var1_marker == "color":
				ax.scatter(
					filtered["ratiox"], filtered["ratioy"],
					color=args.var1_marker_vals[idx],
					alpha=0.5, label=i)
			elif args.var1_marker = "shape":
				ax.scatter(
					filtered["ratiox"], filtered["ratioy"],
					marker=args.var1_marker_vals[idx], label=i)
			else:
				ax.scatter(
					filtered["ratiox"], filtered["ratioy"],
					s=args.var1_marker_vals[idx], label=i)
	# Two or more variables - only grab first two
	elif len(meta_df.columns) > 2:
		if args.var2_marker == "none":
			print(
				"More than one variables in --meta. Please set --var2_marker "
				"or only keep a single variable (other than 'Sample') in "
				"the --meta file.")
			sys.exit(1)

		for idx, i in enumerate(merged[meta_df.columns[1]].unique()):
			for jdx, j in enumerate(merged[meta_df.columns[1]].unique()):
				filtered = merged.loc[merged[args.var1_marker] == i]
				filtered = filtered.loc[filtered[args.var2_marker] == j]
				if args.var1_marker == "color":
					if args.var2_marker == "shape":
						ax.scatter(
							filtered["ratiox"], filtered["ratioy"],
							color=args.var1_marker_vals[idx],
							marker=args.var2_marker_vals[jdx],
							label="{}_{}".format(i, j))
					elif args.var2_marker == "size":
						ax.scatter(
							filtered["ratiox"], filtered["ratioy"],
							color=args.var1_marker_vals[idx],
							s=args.var2_marker_vals[jdx],
							label="{}_{}".format(i, j))
				elif args.var1_marker == "shape":
					if args.var2_marker == "color":
						ax.scatter(
							filtered["ratiox"], filtered["ratioy"],
							color=args.var2_marker_vals[idx],
							marker=args.var1_marker_vals[jdx],
							label="{}_{}".format(i, j))
					elif args.var2_marker == "size":
						ax.scatter(
							filtered["ratiox"], filtered["ratioy"],
							s=args.var2_marker_vals[idx],
							marker=args.var1_marker_vals[jdx],
							label="{}_{}".format(i, j))
				else:
					if args.var2_marker == "color":
						ax.scatter(
							filtered["ratiox"], filtered["ratioy"],
							color=args.var2_marker_vals[idx],
							s=args.var1_marker_vals[jdx],
							label="{}_{}".format(i, j))
					elif args.var2_marker == "shape":
						ax.scatter(
							filtered["ratiox"], filtered["ratioy"],
							marker=args.var2_marker_vals[idx],
							s=args.var1_marker_vals[jdx],
							label="{}_{}".format(i, j))

	# Need at least one variable
	else:
		print("Need at least one variable other than 'Sample' in --meta file.")
		sys.exit(1)

	ax.set_xlabel("{} / {} Ratio".format(args.first_chr, args.const_chr))
	ax.set_ylabel("{} / {} Ratio".format(args.second_chr, args.const_chr))
	fig.savefig("{}.svg".format(args.output_prefix))
	fig.savefig("{}.png".format(args.output_prefix))

if __name__ == "__main__":
	main()
