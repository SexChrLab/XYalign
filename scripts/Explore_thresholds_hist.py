import argparse
import numpy as np
import pandas as pd
import pybedtools
from xyalign import utils, variants


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
		"--vcf", type=str, required=True, help="Full path to Platypus VCF output "
		"from BAM_ANALYSIS module")

	parser.add_argument(
		"--chrom", type=str, default="ALL", help="Name of chromosome to analyze. "
		"Default is 'ALL', which will analyze all chromosomes. Otherwise, will only "
		"plot for chromosome listed.")

	parser.add_argument(
		"--whole_genome_threshold", action="store_true", default=False,
		help="If flag provided, use full dataset to calculate mean for filters. "
		"Otherwise, will calculate mean per chromosome.")

	parser.add_argument(
		"--min_depth_filter", nargs='*', default=[0.0],
		help="Minimum depth threshold for a window to be considered high "
		"quality. Calculated as mean depth * min_depth_filter. So, a "
		"min_depth_filter of 0.2 would require at least a minimum depth "
		"of 2 if the mean depth was 10. Default is 0.0 to consider all windows.")

	parser.add_argument(
		"--max_depth_filter", nargs='*', default=[10000.0],
		help="Maximum depth threshold for a window to be considered high "
		"quality. Calculated as mean depth * max_depth_filter. So, a "
		"max_depth_filter of 4 would require depths to be less than or "
		"equal to 40 if the mean depth was 10. "
		"Default is 10000.0 to consider all windows.")

	parser.add_argument(
		"--mapq_cutoff", nargs='*', default=[20],
		help="Minimum mean mapq threshold for a window to be "
		"considered high quality. Default is 20.")

	parser.add_argument(
		"--variant_site_quality", type=int, default=30,
		help="Consider all SNPs with a site quality (QUAL) greater than or "
		"equal to this value. Default is 30.")

	parser.add_argument(
		"--variant_genotype_quality", type=int, default=30,
		help="Consider all SNPs with a sample genotype quality greater than or "
		"equal to this value. Default is 30.")

	parser.add_argument(
		"--variant_depth", type=int, default=4,
		help="Consider all SNPs with a sample depth greater than or "
		"equal to this value. Default is 4.")

	args = parser.parse_args()

	return args


def main():
	args = parse_args()

	pd_df = pd.read_csv(
		args.dataframe, sep="\t", header=0)

	vcf = variants.VCFFile(args.vcf)

	if args.chrom != "ALL":
		chrom_df = pd_df.loc[pd_df['chrom'] == args.chrom]
	else:
		chrom_df = pd_df

	if args.whole_genome_threshold is True:
		fail_df, pass_df = utils.make_region_lists_genome_filters(
			depthAndMapqDf=chrom_df,
			mapqCutoff=args.mapq_cutoff,
			min_depth=args.min_depth_filter,
			max_depth=args.max_depth_filter)
	else:
		fail_df, pass_df = utils.make_region_lists_chromosome_filters(
			depthAndMapqDf=chrom_df,
			mapqCutoff=args.mapq_cutoff,
			min_depth=args.min_depth_filter,
			max_depth=args.max_depth_filter)

	if args.CHROM == "ALL":
		chrom_dict = {}
		for i in pass_df.chrom.unique():
			chrom_dict[i] = vcf.parse_platypus_VCF(
				site_qual=args.variant_site_quality,
				genotype_qual=args.variant_genotype_quality,
				depth=args.variant_depth,
				chrom=i)
	else:
		parsed = vcf.parse_platypus_VCF(
			site_qual=args.variant_site_quality,
			genotype_qual=args.variant_genotype_quality,
			depth=args.variant_depth,
			chrom=args.chrom)

		start_sites = chrom_df["start"].values
		tmp_df = chrom_df[["start", "stop"]]
		coord_dict = tmp_df.set_index("start").to_dict()

		pos = []
		qual = []
		rb = []

		for idx, i in enumerate(parsed[0]):
			tmp_start = start_sites[np.searchsorted(start_sites, i, side="right")]
			if i < coord_dict[tmp_start]:
				pos.append(i)
				qual.append(parsed[1][idx])
				rb.append(parsed[2][idx])

		rc = variants.hist_read_balance(
			chrom=args.chrom,
			readBalance=np.asarray(rb),
			homogenize=False,
			output_prefix="{}_mapq{}_mindepth{}_max_depth".format(args.output_prefix)
		)
