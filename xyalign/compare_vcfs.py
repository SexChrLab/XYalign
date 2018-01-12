from __future__ import print_function
import argparse
import numpy as np
import xyalign  # from xyalign import utils, variants


def parse_args():
	"""
	Parse command-line arguments

	Returns
	-------

	Parser argument namespace
	"""

	parser = argparse.ArgumentParser(
		description="This script takes as input two Platypus VCF files and "
		"compares them to find differences in variant presence/absence "
		"and genotype quality.")

	parser.add_argument(
		"--vcf_before", type=str, required=True, help="Full path to first Platypus "
		"VCF output from BAM_ANALYSIS module. This will be treated as the 'before' "
		"condition for comparisons.")

	parser.add_argument(
		"--vcf_after", type=str, required=True, help="Full path to second Platypus "
		"VCF output from BAM_ANALYSIS module. This will be treated as the 'after' "
		"condition for comparisons.")

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
		"--chrom", type=str, required=True, help="Name of chromosome to analyze.")

	args = parser.parse_args()

	return args


def main():
	args = parse_args()

	vcf_before = xyalign.variants.VCFFile(args.vcf_before)
	vcf_after = xyalign.variants.VCFFile(args.vcf_after)

	parsed_before = vcf_before.parse_platypus_VCF(
		site_qual=args.variant_site_quality,
		genotype_qual=args.variant_genotype_quality,
		depth=args.variant_depth,
		chrom=args.chrom)

	parsed_after = vcf_after.parse_platypus_VCF(
		site_qual=args.variant_site_quality,
		genotype_qual=args.variant_genotype_quality,
		depth=args.variant_depth,
		chrom=args.chrom)

	missing_in_before = [[], []]
	before_also_in_after = [[], [], [], []]
	for idx, i in enumerate(parsed_after[0]):
		if i not in parsed_before[0]:
			missing_in_before[0].append(i)
			missing_in_before[1].append(parsed_after[2][idx])
		else:
			before_also_in_after[0].append(parsed_before[0][idx])
			before_also_in_after[1].append(parsed_before[1][idx])
			before_also_in_after[2].append(parsed_before[2][idx])
			before_also_in_after[3].append(parsed_before[3][idx])

	missing_in_after = [[], []]
	after_also_in_before = [[], [], [], []]
	for idx, i in enumerate(parsed_before[0]):
		if i not in parsed_after[0]:
			missing_in_after[0].append(i)
			missing_in_after[1].append(parsed_before[2][idx])
		else:
			after_also_in_before[0].append(parsed_after[0][idx])
			after_also_in_before[1].append(parsed_after[1][idx])
			after_also_in_before[2].append(parsed_after[2][idx])
			after_also_in_before[3].append(parsed_after[3][idx])

	# missing_in_before = [x for x in parsed_after[0] if x not in parsed_before[0]]
	# missing_in_after = [x for x in parsed_before[0] if x not in parsed_after[0]]


if __name__ == "__main__":
	main()
