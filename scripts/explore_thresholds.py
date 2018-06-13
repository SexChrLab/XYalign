from __future__ import print_function
import argparse
import numpy as np
import os
import pandas as pd
import pybedtools
import random
import string
import xyalign.utils as utils
import xyalign.variants as variants
# from xyalign import utils, variants


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
		"--callable_bed", type=str, default="None", help="Full path to OPTIONAL "
		"external bed file with callable regions. This script will plot based "
		"on filters only, callable sites only, and filters and callable "
		"sites combined. Default is 'None', which will plot based on filters only.")

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


def filter_and_plot(
	whole_genome, dataframe, mapq_thresh, min_dp, max_dp, chromosome, vcf_file,
	variant_site_quality, variant_genotype_quality, variant_depth, sample_id,
	out_pre):
	"""
	"""
	if whole_genome is True:
		pass_df, fail_df = utils.make_region_lists_genome_filters(
			depthAndMapqDf=dataframe,
			mapqCutoff=int(mapq_thresh),
			min_depth=float(min_dp),
			max_depth=float(max_dp))
	else:
		pass_df, fail_df = utils.make_region_lists_chromosome_filters(
			depthAndMapqDf=dataframe,
			mapqCutoff=int(mapq_thresh),
			min_depth=float(min_dp),
			max_depth=float(max_dp))

	print(
		"Checking empirical values. Min_depth={}, Max_depth={}, Min_Mapq={}".format(
			np.amin(pass_df["depth"].values),
			np.amax(pass_df["depth"].values),
			np.amin(pass_df["mapq"].values)))

	if chromosome == "ALL":
		chrom_dict = {}
		for i in pass_df.chrom.unique():
			chrom_dict[i] = vcf_file.parse_platypus_VCF(
				site_qual=int(variant_site_quality),
				genotype_qual=int(variant_genotype_quality),
				depth=int(variant_depth),
				chrom=i)

		rb = []
		for i in chrom_dict:
			tmp_chrom_df = pass_df.loc[pass_df['chrom'] == i]
			start_sites = tmp_chrom_df["start"].values
			tmp_df = tmp_chrom_df[["start", "stop"]]
			coord_dict = tmp_df.set_index("start").to_dict()['stop']

			for idx, i in enumerate(chrom_dict[i][0]):
				try:
					tmp_start = start_sites[np.searchsorted(start_sites, i, side="right")]
				except IndexError:
					tmp_start = start_sites[-1]
				if i < coord_dict[tmp_start]:
					rb.append(chrom_dict[i][2][idx])

	else:
		parsed = vcf_file.parse_platypus_VCF(
			site_qual=variant_site_quality,
			genotype_qual=variant_genotype_quality,
			depth=variant_depth,
			chrom=chromosome)

		print(
			"{} sites with only variant-specific filters".format(len(parsed[0])))

		start_sites = pass_df["start"].values
		tmp_df = pass_df[["start", "stop"]]
		coord_dict = tmp_df.set_index("start").to_dict()['stop']

		print(
			"{} windows passing filters".format(len(start_sites)))

		rb = []

		for idx, i in enumerate(parsed[0]):
			try:
				tmp_start = start_sites[np.searchsorted(start_sites, i, side="right")]
			except IndexError:
				tmp_start = start_sites[-1]
			if i < coord_dict[tmp_start]:
				rb.append(parsed[2][idx])

		# if args.plot_snp_distance is True:
		# 	rc = utils.hist_array(
		# 		chrom=args.chrom,
		# 		value_array=np.asarray(dist),
		# 		measure_name="Distance between SNPs",
		# 		sampleID=args.sample_id,
		# 		output_prefix="{}_mapq{}_mindepth{}_maxdepth{}".format(
		# 			args.output_prefix, mapq_thresh, min_dp, max_dp))

		# print("Mean distance between SNPs: {} ({} std)".format(
		# 	np.mean(dist), np.std(dist)))

	print(
		"{} variant sites after window filtering\n".format(len(rb)))

	rc = variants.hist_read_balance(
		chrom=chromosome,
		readBalance=np.asarray(rb),
		sampleID=sample_id,
		homogenize=False,
		output_prefix="{}_mapq{}_mindepth{}_maxdepth{}".format(
			out_pre, mapq_thresh, min_dp, max_dp))

	return rc


def main():
	pybedtools.set_tempdir(".")
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

	tmp_bed = "tmp_{}_{}.bed".format(
		''.join(
			random.choice(string.ascii_uppercase + string.digits) for _ in range(6)),
		''.join(
			random.choice(string.ascii_uppercase + string.digits) for _ in range(6)))
	if args.callable_bed != "None":
		callable_bed = pybedtools.BedTool(args.callable_bed)
		callable_df = callable_bed.to_dataframe(
			names=["chrom", "start", "stop"])
		# tmp_bed = "tmp_{}_{}.bed".format(
		# 	''.join(
		# 		random.choice(string.ascii_uppercase + string.digits) for _ in range(6)),
		# 	''.join(
		# 		random.choice(string.ascii_uppercase + string.digits) for _ in range(6)))
		utils.output_bed_no_merge(tmp_bed, chrom_df)
		# chrom_df_bed = pybedtools.BedTool().from_dataframe(chrom_df).saveas("test_tmp.bed")
		# print(chrom_df_bed)
		import_tmp = pybedtools.BedTool(tmp_bed)
		intersection = import_tmp.intersect(callable_bed)
		# intersection = chrom_df_bed.intersect(callable_bed)
		intersection_df = intersection.to_dataframe(
			names=["chrom", "start", "stop", "depth", "mapq"])

	for min_dp in args.min_depth_filter:
		for max_dp in args.max_depth_filter:
			for mapq_thresh in args.mapq_cutoff:
				print(
					"Filtering with min_depth={} times mean, "
					"max_depth={} times mean, and min_mapq={}".format(
						min_dp, max_dp, mapq_thresh))

				print("Using filters only")
				filters_only = filter_and_plot(
					whole_genome=args.whole_genome_threshold,
					dataframe=chrom_df,
					mapq_thresh=int(mapq_thresh),
					min_dp=float(min_dp),
					max_dp=float(max_dp),
					chromosome=args.chrom,
					vcf_file=vcf,
					variant_site_quality=int(args.variant_site_quality),
					variant_genotype_quality=int(args.variant_genotype_quality),
					variant_depth=int(args.variant_depth),
					sample_id=args.sample_id,
					out_pre="{}_filtersonly".format(args.output_prefix))

				if args.callable_bed != "None":
					print("Using filters and callable bed file")
					filters_only = filter_and_plot(
						whole_genome=args.whole_genome_threshold,
						dataframe=intersection_df,
						mapq_thresh=int(mapq_thresh),
						min_dp=float(min_dp),
						max_dp=float(max_dp),
						chromosome=args.chrom,
						vcf_file=vcf,
						variant_site_quality=int(args.variant_site_quality),
						variant_genotype_quality=int(args.variant_genotype_quality),
						variant_depth=int(args.variant_depth),
						sample_id=args.sample_id,
						out_pre="{}_filters_and_callablebed".format(args.output_prefix))

	if os.path.exists(tmp_bed):
		os.remove(tmp_bed)


if __name__ == "__main__":
	main()
