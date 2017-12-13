from __future__ import print_function
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
		"histogram of read balances given specified MAPQ and depth thresholds. "
		"It either plots the chromosome specified by --chrom or all data "
		"in the input dataframe when --chrom ALL is specified.")

	parser.add_argument(
		"--dataframe", type=str, required=True, help="Full path to csv output of "
		"pandas dataframe from BAM_ANALYSIS module")

	parser.add_argument(
		"--vcf", type=str, required=True, help="Full path to Platypus VCF output "
		"from BAM_ANALYSIS module")

	parser.add_argument(
		"--output_prefix", type=str, help="Full path to and prefix of output files.")

	parser.add_argument(
		"--chrom", type=str, default="ALL", help="Name of chromosome to analyze. "
		"Default is 'ALL', which will analyze all chromosomes in dataframe "
		"together. Otherwise, will only plot for chromosome listed.")

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
		"--max_depth_filter", nargs='*', default=[100000.0],
		help="Maximum depth threshold for a window to be considered high "
		"quality. Calculated as mean depth * max_depth_filter. So, a "
		"max_depth_filter of 4 would require depths to be less than or "
		"equal to 40 if the mean depth was 10. "
		"Default is 100000.0 to consider all windows.")

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

	parser.add_argument(
		"--sample_id", help="Sample ID or other identifier to be used in naming")

	parser.add_argument(
		"--plot_snp_distance", action="store_true", help="If True, will also plot "
		"a histogram of distances between SNPs. Will only run on a single "
		"chromosome. Default is False.")

	args = parser.parse_args()

	return args


def main():
	args = parse_args()

	pd_df = pd.read_csv(
		args.dataframe, sep="\t", header=0)

	vcf = variants.VCFFile(args.vcf)

	if args.chrom != "ALL":
		chrom_df = pd_df.loc[pd_df['chrom'] == args.chrom]
		chrom_df.index = pd.RangeIndex(len(chrom_df.index))
		#chrom_df.reset_index(drop=True)
	else:
		chrom_df = pd_df

	for min_dp in args.min_depth_filter:
		for max_dp in args.max_depth_filter:
			for mapq_thresh in args.mapq_cutoff:
				print(min_dp, max_dp, mapq_thresh)

				if args.whole_genome_threshold is True:
					pass_df, fail_df = utils.make_region_lists_genome_filters(
						depthAndMapqDf=chrom_df,
						mapqCutoff=int(mapq_thresh),
						min_depth=float(min_dp),
						max_depth=float(max_dp))
				else:
					pass_df, fail_df = utils.make_region_lists_chromosome_filters(
						depthAndMapqDf=chrom_df,
						mapqCutoff=int(mapq_thresh),
						min_depth=float(min_dp),
						max_depth=float(max_dp))

				print(np.amin(pass_df["mapq"].values), np.amin(pass_df["depth"].values), np.amax(pass_df["depth"].values))

				if args.chrom == "ALL":
					chrom_dict = {}
					for i in pass_df.chrom.unique():
						chrom_dict[i] = vcf.parse_platypus_VCF(
							site_qual=args.variant_site_quality,
							genotype_qual=args.variant_genotype_quality,
							depth=args.variant_depth,
							chrom=i)

					rb = []
					for i in chrom_dict:
						tmp_chrom_df = chrom_df.loc[chrom_df['chrom'] == i]
						start_sites = chrom_df["start"].values
						tmp_df = chrom_df[["start", "stop"]]
						coord_dict = tmp_df.set_index("start").to_dict()

						for idx, i in enumerate(chrom_dict[i][0]):
							tmp_start = start_sites[np.searchsorted(start_sites, i, side="right")]
							if i < coord_dict[tmp_start]:
								rb.append(parsed[2][idx])

				else:
					parsed = vcf.parse_platypus_VCF(
						site_qual=args.variant_site_quality,
						genotype_qual=args.variant_genotype_quality,
						depth=args.variant_depth,
						chrom=args.chrom)

					print(len(parsed[0]))

					start_sites = pass_df["start"].values
					tmp_df = pass_df[["start", "stop"]]
					coord_dict = tmp_df.set_index("start").to_dict()['stop']

					print(len(start_sites), len(coord_dict))

					rb = []
					dist = []

					last_snp = None
					for idx, i in enumerate(parsed[0]):
						if last_snp is not None:
							dist.append(i - last_snp)
							last_snp = i
						else:
							last_snp = i
						try:
							tmp_start = start_sites[np.searchsorted(start_sites, i, side="right")]
						except IndexError:
							tmp_start = start_sites[-1]
						if i < coord_dict[tmp_start]:
							rb.append(parsed[2][idx])

					if args.plot_snp_distance is True:
						rc = utils.hist_array(
							chrom=args.chrom,
							value_array=np.asarray(dist),
							measure_name="Distance between SNPs",
							sampleID=args.sample_id,
							output_prefix="{}_mapq{}_mindepth{}_maxdepth{}".format(
								args.output_prefix, mapq_thresh, min_dp, max_dp))

					print("Mean distance between SNPs: {} ({} std)".format(
						np.mean(dist), np.std(dist)))

				print(len(rb))

				rc = variants.hist_read_balance(
					chrom=args.chrom,
					readBalance=np.asarray(rb),
					sampleID=args.sample_id,
					homogenize=False,
					output_prefix="{}_mapq{}_mindepth{}_maxdepth{}".format(
						args.output_prefix, mapq_thresh, min_dp, max_dp))

if __name__ == "__main__":
	main()
