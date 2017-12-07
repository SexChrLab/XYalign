import argparse
import pandas as pd
import pybedtools


def parse_args():
	"""
	Parse command-line arguments

	Returns
	-------

	Parser argument namespace
	"""

	parser = argparse.ArgumentParser(
		description="This script takes as input a CSV file of a pandas dataframe "
		"and Platypus VCF file output by the BAM_ANALYSIS module and plots a "
		"histogram of read balances given specified MAPQ and depth thresholds.")

	parser.add_argument(
		"--dataframe", type=str, required=True, help="Full path to csv output of "
		"pandas dataframe from BAM_ANALYSIS module")

	parser.add_argument(
		"--chrom", type=str, default="ALL", help="Name of chromosome to analyze. "
		"Default is 'ALL', which will analyze all chromosomes. Otherwise, will only "
		"plot for chromosome listed.")

	parser.add_argument(
		"--whole_genome_threshold", action="store_true", default=False,
		help="If flag provided, use full dataset to calculate mean for filters. "
		"Otherwise, will calculate mean per chromosome.")

	args = parser.parse_args()

	return args


def main():
	args = parse_args()

	pd_df = pd.read_csv(
		args.dataframe, sep="\t", header=0)
