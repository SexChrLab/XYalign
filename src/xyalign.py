# To-do list
# 1) Add ploidy estimation
# 		- added permutation tests
# 		- need to add likelihood analyses (model fitting)
# 2) Compartmentalize all steps of analysis
# 		- Add flags to make each part of the pipeline optional
# 		- Allow users to call specific parts of the pipeline
# 					(e.g. only vcf plotting)
# 		- Add checkpointing
# 3) Write to a better designed output directory structure
# 4) Better (and unified) naming scheme for plots and output
# 5) Generalize mapping and calling (perhaps by allowing users to
# 		add command lines as  strings)
# 6) Add plotting of high-quality windows (depth, mapq), also after remapping

# 7) Check for behavior when files already exist (e.g., overwrite, quit, etc.?)
# 8) Incorporate mask integration on the fly
# 9) Check with Python 3 and see if any incompatibilities (e.g., printing) can
# 		be solved with from __future__

from __future__ import division
import argparse
import csv
import os
import subprocess
import sys
import numpy as np
import pandas as pd
import pybedtools
import pysam
# Setting the matplotlib display variable requires importing
# 	in exactly the following order:
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns


def main():
	""" Main program"""

	args = parse_args()

	# First round of Platypus calling and plotting
	if args.platypus_calling == "both" or "before":
		if args.bam is not None:
			a = platypus_caller(
				args.bam, args.ref, args.chromosomes, args.cpus,
				args.output_dir + "/{}.noprocessing.vcf".format(
					args.sample_id), None)
			if a != 0:
				print "Error in initial Platypus calling."
				sys.exit(1)
			if args.no_variant_plots is True:
				plot_variants_per_chrom(
					args.chromosomes,
					args.output_dir + "/{}.noprocessing.vcf".format(
						args.sample_id),
					args.sample_id, args.output_dir, "noprocessing",
					args.variant_quality_cutoff, args.marker_size,
					args.marker_transparency, args.bam)
		else:
			a = platypus_caller(
				args.cram, args.ref, args.chromosomes, args.cpus,
				args.output_dir + "/{}.noprocessing.vcf".format(
					args.sample_id), None)
			if a != 0:
				print "Error in initial Platypus calling."
				sys.exit(1)
			if args.no_variant_plots is True:
				plot_variants_per_chrom(
					args.chromosomes,
					args.output_dir + "/{}.noprocessing.vcf".format(
						args.sample_id),
					args.sample_id, args.output_dir, "noprocessing",
					args.variant_quality_cutoff, args.marker_size,
					args.marker_transparency, args.cram)

	# Analyze bam for depth and mapq
	if args.bam is not None:
		samfile = pysam.AlignmentFile(args.bam, "rb")
	else:
		samfile = pysam.AlignmentFile(args.cram, "rc")
	pass_df = []
	fail_df = []
	for chromosome in args.chromosomes:
		data = traverse_bam_fetch(samfile, chromosome, args.window_size)
		tup = make_region_lists(
			data["windows"], args.mapq_cutoff, args.depth_filter)
		pass_df.append(tup[0])
		fail_df.append(tup[1])
		plot_depth_mapq(
			data, args.output_dir, args.sample_id,
			get_length(samfile, chromosome), args.marker_size,
			args.marker_transparency)
	output_bed(os.path.join(args.output_dir, args.high_quality_bed), *pass_df)
	output_bed(os.path.join(args.output_dir, args.low_quality_bed), *fail_df)

	# Infer ploidy (needs to be finished)
	if args.y_present is True:
		y_present = True
	elif args.y_absent is True:
		y_present = False
	else:
		# Replace this with code to infer ploidy, etc.
		# Permutation tests
		if args.no_perm_test is not True:
			sex_chromosomes = args.x_chromosome + args.y_chromosome
			autosomes = [x for x in args.chromosomes not in sex_chromosomes]
			perm_res = []
			for c in autosomes:
				perm_res.append(permutation_test_chromosomes(
					pd.concat(pass_df), c, args.x_chromosome, "chrom",
					"depth", args.num_permutations,
					"{}_{}_permutation_results.txt".format(
						c, args.x_chromosome)))
			sex_perm_res = permutation_test_chromosomes(
				pd.concat(pass_df), args.x_chromosome, args.y_chromosome,
				"chrom", "depth", args.num_permutations,
				"{}_{}_permutation_results.txt".format(
					args.x_chromosome, args.y_chromosome))
			if 0.025 < sex_perm_res[2] < 0.95:
				y_present = True

		# Likelihood analyses

	# Remapping
	if args.no_remapping is True:
		if y_present is True:
			# Isolate sex chromosomes from reference and index new reference
			new_reference = isolate_chromosomes_reference(
				args.ref, "{}/{}.sex_chroms".format(
					args.output_dir, args.sample_id),
				args.x_chromosome + args.y_chromosome)
			# Strip reads from sex chromosomes
			if args.bam is not None:
				new_fastqs = bam_to_fastq(
					args.bam, args.single_end, args.output_dir, args.sample_id,
					args.x_chromosome + args.y_chromosome)
			else:
				new_fastqs = bam_to_fastq(
					args.cram, args.single_end, args.output_dir, args.sample_id,
					args.x_chromosome + args.y_chromosome)
			# Remap
			new_bam = bwa_mem_mapping(
				new_reference, "{}/{}.sex_chroms".format(
					args.output_dir, args.sample_id),
				new_fastqs)
			# Merge bam files
			if args.bam is not None:
				merged_bam = switch_sex_chromosomes_bam(
					args.bam, new_bam, args.x_chromosome + args.y_chromosome,
					args.output_dir, args.sample_id)
			else:
				merged_bam = switch_sex_chromosomes_bam(
					args.cram, new_bam, args.x_chromosome + args.y_chromosome,
					args.output_dir, args.sample_id)

		else:
			# Isolate sex chromosomes from reference and index new reference
			new_reference = isolate_chromosomes_reference(
				args.ref, "{}/{}.sex_chroms".format(
					args.output_dir, args.sample_id),
				args.x_chromosome)
			# Strip reads from sex chromosomes
			if args.bam is not None:
				new_fastqs = bam_to_fastq(
					args.bam, args.single_end, args.output_dir, args.sample_id,
					args.x_chromosome)
			else:
				new_fastqs = bam_to_fastq(
					args.cram, args.single_end, args.output_dir, args.sample_id,
					args.x_chromosome)
			# Remap
			new_bam = bwa_mem_mapping(
				new_reference, "{}/{}.sex_chroms".format(
					args.output_dir, args.sample_id),
				new_fastqs)
			# Merge bam files
			if args.bam is not None:
				merged_bam = switch_sex_chromosomes_bam(
					args.bam, new_bam, args.x_chromosome + args.y_chromosome,
					args.output_dir, args.sample_id)
			else:
				merged_bam = switch_sex_chromosomes_bam(
					args.cram, new_bam, args.x_chromosome + args.y_chromosome,
					args.output_dir, args.sample_id)

	# Final round of calling and plotting
	variant_mask = os.path.join(args.output_dir, args.high_quality_bed)

	if args.platypus_calling == "both" or "after":
		a = platypus_caller(
			merged_bam, args.ref, args.chromosomes, args.cpus,
			args.output_dir + "/{}.postprocessing.vcf".format(args.sample_id),
			variant_mask)
		if a != 0:
			print "Error in initial Platypus calling."
			sys.exit(1)
		if args.no_variant_plots is True:
			plot_variants_per_chrom(
				args.chromosomes,
				args.output_dir + "/{}.postprocessing.vcf".format(
					args.sample_id),
				args.sample_id, args.output_dir, "postprocessing",
				args.variant_quality_cutoff, args.marker_size,
				args.marker_transparency, merged_bam)


def parse_args():
	"""Parse command-line arguments"""
	parser = argparse.ArgumentParser(description="XYalign")

	parser.add_argument(
		"--ref", required=True,
		help="Path to reference sequence (including file name).")

	parser.add_argument(
		"--chromosomes", "-c", nargs="+", default=["chrX", "chrY", "chr19"],
		help="Chromosomes to analyze.")

	parser.add_argument(
		"--x_chromosome", "-x", nargs="+", default=["chrX"],
		help="Names of y-linked scaffolds in reference fasta.")

	parser.add_argument(
		"--y_chromosome", "-y", nargs="+", default=["chrY"],
		help="Names of y-linked scaffolds in reference fasta.")

	parser.add_argument(
		"--sample_id", "-id", default="sample",
		help="Name/ID of sample - for use in plot titles and file naming.")

	parser.add_argument(
		"--single_end", action="store_true", default=False,
		help="Include flag if reads are single-end and NOT paired-end.")

	parser.add_argument(
		"--platypus_calling", default="both",
		help="Platypus calling withing the pipeline "
		"(before processing, after processing, both, "
		"or neither). Options: both, none, before, after.")

	parser.add_argument(
		"--no_variant_plots", action="store_false", default=True,
		help="Include flag to prevent plotting read balance from VCF files.")

	parser.add_argument(
		"--no_remapping", action="store_false", default=True,
		help="Include this flag to prevent remapping sex chromosome reads.")

	parser.add_argument(
		"--variant_quality_cutoff", "-vqc", type=int, default=20,
		help="Consider all SNPs with a quality greater than or "
		"equal to this value. Default is 20.")

	parser.add_argument(
		"--marker_size", type=float, default=10.0,
		help="Marker size for genome-wide plots in matplotlib.")

	parser.add_argument(
		"--marker_transparency", "-mt", type=float, default=0.5,
		help="Transparency of markers in genome-wide plots.  "
		"Alpha in matplotlib.")

	parser.add_argument(
		"--cpus", type=int, default=1,
		help="Number of cores/threads to use.")

	parser.add_argument(
		"--window_size", "-w", type=int, default=50000,
		help="Window size (integer) for sliding window calculations.")

	parser.add_argument(
		"--mapq_cutoff", "-mq", type=int, default=20,
		help="Minimum mean mapq threshold for a window to be "
		"considered high quality.")

	parser.add_argument(
		"--depth_filter", "-df", type=float, default=4.0,
		help="Filter for depth (f), where the threshold used is mean_depth +- "
		"(f * square_root(mean_depth)).")

	parser.add_argument(
		"--high_quality_bed", "-hq", default="highquality.bed",
		help="Name of output file for high quality regions.")

	parser.add_argument(
		"--low_quality_bed", "-lq", default="lowquality.bed",
		help="Name of output file for high quality regions.")

	parser.add_argument(
		"--num_permutations", type=int, default=10000,
		help="Number of permutations to use for permutation analyses")

	parser.add_argument(
		"--no_perm_test", action="store_true", default=False,
		help="Include flag to turn off permutation tests. Requires either "
		"--y_present or --y_absent to also be called")

	parser.add_argument(
		"--output_dir", "-o",
		help="Output directory")

	group = parser.add_mutually_exclusive_group(required=True)

	group.add_argument(
		"--bam", help="Input bam file.")

	group.add_argument(
		"--cram", help="Input cram file.")

	group2 = parser.add_mutually_exclusive_group(required=False)

	group2.add_argument(
		"--y_present", action="store_true", default=False,
		help="Overrides sex chr estimation by XYalign and remaps with Y present.")

	group2.add_argument(
		"--y_absent", action="store_true", default=False,
		help="Overrides sex chr estimation by XY align and remaps with Y absent.")

	args = parser.parse_args()

	# Validate arguments
	if not os.path.exists(args.output_dir):
		os.makedirs(args.output_dir)
	if args.platypus_calling not in ["both", "none", "before", "after"]:
		print "Error. Platypus calling must be both, none, before, or after. ",\
			"Default is both."
		sys.exit(1)
	if args.no_perm_test is True:
		if args.y_present is False and args.y_absent is False:
			print "Error. Either --y_present or --y_absent needs to be ",\
				"included with --no_perm_test"
		sys.exit(1)

	# Return arguments namespace
	return args


def get_length(bamfile, chrom):
	""" Extract chromosome length from BAM header.

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


def permutation_test_chromosomes(
	data_frame, first_chrom, second_chrom, chrom_column,
	value_column, num_perms, output_file=None):
	""" Takes a dataframe and runs a permutation test comparing mean values
	of two chromosomes.
	"""
	first_vals = data_frame[chrom_column == first_chrom].value_column
	second_vals = data_frame[chrom_column == second_chrom].value_column
	combined = np.append(first_vals, second_vals)

	first_mean = np.mean(first_vals)
	second_mean = np.mean(second_vals)

	observed = first_mean - second_mean
	perms = []
	for i in range(0, num_perms):
		np.random.shuffle(combined)
		first = np.mean(combined[:len(first_vals)])
		second = np.mean(combined[-len(second_vals):])
		perms.append(first - second)
	perms = np.asarray(perms)
	sig = len(np.where(perms > observed)) / num_perms
	if output_file is not None:
		a = [
			"{}_mean".format(first_chrom),
			"{}_mean".format(second_chrom),
			"{}_{}_diff".format(first_chrom, second_chrom),
			"p_val_({}_>_{})".format(first_chrom, second_chrom),
			"perm_2.5",
			"perm_50",
			"perm_97.5"]
		b = [
			"{}".format(first),
			"{}".format(second),
			"{}".format(observed),
			"{}".format(sig),
			"{}".format(np.percentile(perms, 2.5)),
			"{}".format(np.percentile(perms, 50)),
			"{}".format(np.percentile(perms, 97.5))]
		with open(output_file, "w") as f:
			w = csv.writer(f, dialect="excel-tab")
			w.writerows([a, b])
	return (first_mean, second_mean, sig)


def bwa_mem_mapping(reference, output_prefix, fastqs, cram=False):
	""" Maps reads to a reference genome using bwa mem.
	"""
	fastqs = ' '.join(fastqs)
	subprocess.call("bwa index {}".format(reference), shell=True)
	if cram is False:
		command_line = "bwa mem {} {} | samtools fixmate -O bam - - | samtools sort -O bam -o {}_sorted.bam -".format(reference, fastqs, output_prefix)
		subprocess.call(command_line, shell=True)
		subprocess.call(
			"samtools index {}_sorted.bam".format(output_prefix), shell=True)
		return "{}_sorted.bam".format(output_prefix)
	else:
		command_line = "bwa mem {} {} | samtools fixmate -O cram - - | samtools sort -O cram -o {}_sorted.cram -".format(reference, fastqs, output_prefix)
		subprocess.call(command_line, shell=True)
		subprocess.call(
			"samtools index {}_sorted.cram".format(output_prefix), shell=True)
		return "{}_sorted.cram".format(output_prefix)


def switch_sex_chromosomes_bam(
	bam_orig, bam_new, sex_chroms, output_directory, output_prefix, cram=False):
	""" Removes sex chromosomes from original bam and merges in remmapped
	sex chromosomes, while retaining the original bam header
	"""
	# Grab original header
	subprocess.call(
		"samtools view -H {} > {}/header.sam".format(
			bam_orig, output_directory), shell=True)
	if cram is False:
		# Remove sex chromosomes from original bam
		samfile = pysam.AlignmentFile(bam_orig, "rb")
		non_sex_scaffolds = filter(
			lambda x: x not in sex_chroms, list(samfile.references))
		subprocess.call(
			"samtools view -h -b {} {} > {}/no_sex.bam".format(
				bam_orig, " ".join(non_sex_scaffolds), output_directory),
			shell=True)
		subprocess.call(
			"samtools index {}/no_sex.bam".format(output_directory), shell=True)

		# Merge bam files
		subprocess.call(
			"samtools merge -h {}/header.sam {}/{}.bam {}/no_sex.bam {}".format(
				output_directory, output_directory, output_prefix,
				output_directory, bam_new), shell=True)
		subprocess.call("samtools index {}/{}.bam".format(
			output_directory, output_prefix), shell=True)

		return "{}/{}.bam".format(output_directory, output_prefix)

	else:
		# Remove sex chromosomes from original bam
		samfile = pysam.AlignmentFile(bam_orig, "rc")
		non_sex_scaffolds = filter(
			lambda x: x not in sex_chroms, list(samfile.references))
		subprocess.call(
			"samtools view -h -b {} {} > {}/no_sex.cram".format(
				bam_orig, " ".join(non_sex_scaffolds), output_directory),
			shell=True)
		subprocess.call(
			"samtools index {}/no_sex.cram".format(output_directory), shell=True)

		# Merge bam files
		subprocess.call(
			"samtools merge -h {}/header.sam {}/{}.cram {}/no_sex.cram {}".format(
				output_directory, output_directory, output_prefix,
				output_directory, bam_new), shell=True)
		subprocess.call("samtools index {}/{}.cram".format(
			output_directory, output_prefix), shell=True)

		return "{}/{}.cram".format(output_directory, output_prefix)


def platypus_caller(bam, ref, chroms, cpus, output_file, regions_file=None):
	""" Uses platypus to make variant calls on provided bam file

	bam is input bam (or cram) file
	ref is path to reference sequence
	chroms is a list of chromosomes to call on, e.g., ["chrX", "chrY", "chr19"]
	cpus is the number of threads/cores to use
	output_file is the name of the output vcf
	"""
	if regions_file is None:
		regions = ','.join(map(str, chroms))
	else:
		regions = regions_file
	command_line = "platypus callVariants --bamFiles {} -o {} --refFile {} --nCPU {} --regions {} --assemble 1".format(bam, output_file, ref, cpus, regions)
	return_code = subprocess.call(command_line, shell=True)
	return return_code


def isolate_chromosomes_reference(reference_fasta, new_ref_prefix, chroms):
	""" Takes a reference fasta file and a list of chromosomes to include
	and outputs a new, indexed reference fasta.
	"""
	outpath = "{}.fa".format(new_ref_prefix)
	if type(chroms) != list:
		chroms = list(chroms)
	command_line = "samtools faidx {} {} > {}".format(reference_fasta, " ".join(chroms), outpath)
	subprocess.call(command_line, shell=True)
	subprocess.call("samtools faidx {}".format(outpath), shell=True)
	return outpath


def bam_to_fastq(bamfile, single, output_directory, output_prefix, regions):
	""" Strips reads from a bam or cram file in provided regions and outputs
	sorted fastqs containing reads.
	"""
	if single is False:
		command_line = "samtools view -b {} {} | samtools bam2fq -1 {}/temp_1.fastq -2 {}/temp_2.fastq -t -n - ".format(bamfile, ' '.join(map(str, regions)), output_directory, output_directory)
		subprocess.call(command_line, shell=True)
		command_line = "repair.sh in1={} in2={} out1={} out2={} overwrite=true".format(
			output_directory + "/temp_1.fastq",
			output_directory + "/temp_2.fastq",
			output_directory + "/" + output_prefix + "_1.fastq",
			output_directory + "/" + output_prefix + "_2.fastq")
		subprocess.call(command_line, shell=True)
		return [output_directory + "/" + output_prefix + "_1.fastq", output_directory + "/" + output_prefix + "_2.fastq"]
	else:
		command_line = "samtools view -b {} {} | samtools bam2fq -t -n - > {}/temp.fastq".format(bamfile, ' '.join(map(str, regions)), output_directory)
		subprocess.call(command_line, shell=True)
		command_line = "repair.sh in={} out={} overwrite=true".format(output_directory + "/temp.fastq", output_directory + "/" + output_prefix + ".fastq")
		return [output_directory + "/temp.fastq", output_directory + "/" + output_prefix + ".fastq"]


def parse_platypus_VCF(filename, qualCutoff, chrom):
	""" Parse vcf generated by Platypus """
	infile = open("{}".format(filename), 'r')
	positions = []
	quality = []
	readBalance = []
	for line in infile:
		cols = line.strip('\n').split('\t')
		if cols[0] != chrom:
			continue
		pos = int(cols[1])
		qual = float(cols[5])
		if qual < qualCutoff:
			continue
		TR = cols[7].split(';')[17].split('=')[1]
		TC = cols[7].split(';')[14].split('=')[1]
		if ',' in TR or ',' in TC:
			continue
		if (float(TR) == 0) or (float(TC) == 0):
			continue
		ReadRatio = float(TR) / float(TC)

		# Add to arrays
		readBalance.append(ReadRatio)
		positions.append(pos)
		quality.append(qual)

	return (positions, quality, readBalance)


def plot_read_balance(
	chrom, positions, readBalance, sampleID, output_prefix, MarkerSize,
	MarkerAlpha, bamfile):
	""" Plots read balance at each SNP along a chromosome """
	if bamfile[-3] == "bam" or bamfile[-3] == "BAM":
		chrom_len = get_length(pysam.AlignmentFile(bamfile, "rb"), chrom)
	else:
		chrom_len = get_length(pysam.AlignmentFile(bamfile, "rc"), chrom)
	if "x" in chrom.lower():
		Color = "green"
	elif "y" in chrom.lower():
		Color = "blue"
	else:
		Color = "red"
	fig = plt.figure(figsize=(15, 5))
	axes = fig.add_subplot(111)
	axes.scatter(
		positions, readBalance, c=Color, alpha=MarkerAlpha, s=MarkerSize, lw=0)
	axes.set_xlim(0, chrom_len)
	axes.set_title(sampleID)
	axes.set_xlabel("Chromosomal Coordinate")
	axes.set_ylabel("Read Balance")
	# print(len(positions))
	plt.savefig("{}_{}_ReadBalance_GenomicScatter.svg".format(
		output_prefix, chrom))
	plt.savefig("{}_{}_ReadBalance_GenomicScatter.png".format(
		output_prefix, chrom))
	# plt.show()


def hist_read_balance(chrom, readBalance, sampleID, output_prefix):
	""" Plot a histogram of read balance """
	if "x" in chrom.lower():
		Color = "green"
	elif "y" in chrom.lower():
		Color = "blue"
	else:
		Color = "red"
	fig = plt.figure(figsize=(8, 8))
	axes = fig.add_subplot(111)
	axes.set_title(sampleID)
	axes.set_xlabel("Read Balance")
	axes.set_ylabel("Frequency")
	axes.hist(readBalance, bins=50, color=Color)
	plt.savefig("{}_{}_ReadBalance_Hist.svg".format(output_prefix, chrom))
	plt.savefig("{}_{}_ReadBalance_Hist.png".format(output_prefix, chrom))
	# plt.show()


def plot_variants_per_chrom(
	chrom_list, vcf_file, sampleID, output_directory, output_prefix, qualCutoff,
	MarkerSize, MarkerAlpha, bamfile):
	""" Parses a vcf file and plots read balance in separate plots
	for each chromosome in the input list
	"""
	for i in chrom_list:
		parse_results = parse_platypus_VCF(vcf_file, qualCutoff, i)
		plot_read_balance(
			i, parse_results[0], parse_results[2],
			sampleID, output_directory + "/{}.{}".format(
				sampleID, output_prefix), MarkerSize, MarkerAlpha, bamfile)
		hist_read_balance(
			i, parse_results[2], sampleID,
			output_directory + "/{}.{}".format(sampleID, output_prefix))
	pass


def traverse_bam_fetch(samfile, chrom, window_size):
	"""Analyze the `samfile` BAM (or CRAM) file for various metrics.
	Currently, this function looks at the following metrics across genomic
	windows:
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
		print "{} out of {} windows processed on {}".format(
			window_id, num_windows, chrom)

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


def make_region_lists(depthAndMapqDf, mapqCutoff, sd_thresh):
	"""
	(pandas.core.frame.DataFrame, int, float) -> (list, list)
	return two lists of regions (keepList, excludeList) based on
	cutoffs for depth and mapq
	"""
	depth_mean = depthAndMapqDf["depth"].mean()
	depth_sd = depthAndMapqDf["depth"].std()

	depthMin = depth_mean - (sd_thresh * (depth_sd ** 0.5))
	depthMax = depth_mean + (sd_thresh * (depth_sd ** 0.5))

	good = (depthAndMapqDf.mapq > mapqCutoff) & (depthAndMapqDf.depth > depthMin) & (depthAndMapqDf.depth < depthMax)
	dfGood = depthAndMapqDf[good]
	dfBad = depthAndMapqDf[~good]

	return (dfGood, dfBad)


def output_bed(outBed, *regionDfs):
	'''
	(list, list, str) -> bedtoolsObject
	Take two sorted lists.  Each list is a list of tuples
	(chrom[str], start[int], end[int])
	Return a pybedtools object and output a bed file.
	'''
	dfComb = pd.concat(regionDfs)
	regionList = dfComb.ix[:, "chrom":"stop"].values.tolist()
	merge = pybedtools.BedTool(regionList).sort().merge()
	with open(outBed, 'w') as output:
		output.write(str(merge))
	pass


def chromosome_wide_plot(
	chrom, positions, y_value, measure_name, sampleID, output_prefix,
	MarkerSize, MarkerAlpha, Xlim, Ylim):
	'''
	Plots values across a chromosome, where the x axis is the position along the
	chromosome and the Y axis is the value of the measure of interest.
	positions is an array of coordinates
	y_value is an array of the values of the measure of interest
	measure_name is the name of the measure of interest (y axis title)
	chromosome is the name of the chromosome being plotted
	sampleID is the name of the sample
	MarkerSize is the size in points^2
	MarkerAlpha is the transparency (0 to 1)
	Xlim is the maximum X value
	Ylim is the maximum Y value
	'''
	if "x" in chrom.lower():
		Color = "green"
	elif "y" in chrom.lower():
		Color = "blue"
	else:
		Color = "red"
	fig = plt.figure(figsize=(15, 5))
	axes = fig.add_subplot(111)
	axes.scatter(
		positions, y_value, c=Color, alpha=MarkerAlpha, s=MarkerSize, lw=0)
	axes.set_xlim(0, Xlim)
	axes.set_ylim(0, Ylim)
	axes.set_title("%s - %s" % (sampleID, chrom))
	axes.set_xlabel("Chromosomal Position")
	axes.set_ylabel(measure_name)
	plt.savefig("{}_{}_{}_GenomicScatter.svg".format(
		output_prefix, chrom, measure_name))
	plt.savefig("{}_{}_{}_GenomicScatter.png".format(
		output_prefix, chrom, measure_name))
	# plt.show()


def plot_depth_mapq(
	data_dict, output_dir, sampleID, chrom_length, MarkerSize, MarkerAlpha):
	"""
	Takes a dictionary (output from traverseBam) and outputs histograms and
	genome-wide plots of various metrics.
	Args:
		data_dict: Dictionary of pandas data frames
		output_dir: Directory where the PNG file with the plots will be stored
		sampleID: name/ID of sample
		chrom_length: length of chromosome
	Returns:
		None
	"""

	window_df = None if "windows" not in data_dict else data_dict["windows"]
	depth_hist = None if "depth_freq" not in data_dict else data_dict["depth_freq"]
	readbal_hist = None if "readbal_freq" not in data_dict else data_dict["readbal_freq"]
	mapq_hist = None if "mapq_freq" not in data_dict else data_dict["mapq_freq"]

	chromosome = window_df["chrom"][1]

	# Create genome-wide plots based on window means
	if window_df is not None:
		# depth plot
		# depth_genome_plot_path = os.path.join(
		# 	output_dir, "depth_windows." + chromosome + ".png")
		# depth_genome_plot = sns.lmplot(
		# 	'start', 'depth', data=window_df, fit_reg=False,
		# 	scatter_kws={'alpha': 0.3})
		# depth_genome_plot.savefig(depth_genome_plot_path)
		chromosome_wide_plot(
			chromosome, window_df["start"].values, window_df["depth"].values,
			"Depth", sampleID, "{}/{}".format(output_dir, sampleID),
			MarkerSize, MarkerAlpha,
			chrom_length, 100)

		# mapping quality plot
		# mapq_genome_plot_path = os.path.join(
		# 	output_dir, "mapq_windows." + chrom + ".png")
		# mapq_genome_plot = sns.lmplot(
		# 	'start', 'mapq', data=window_df, fit_reg=False)
		# mapq_genome_plot.savefig(mapq_genome_plot_path)
		chromosome_wide_plot(
			chromosome, window_df["start"].values, window_df["mapq"].values,
			"Mapq", sampleID, "{}/{}".format(output_dir, sampleID),
			MarkerSize, MarkerAlpha, chrom_length, 80)

	# Create histograms
	# TODO: update filenames dynamically like window_df above
	# TODO: update Count column name
	if depth_hist is not None:
		depth_bar_plot = sns.countplot(
			x='depth', y='Count', data=depth_hist)
		depth_bar_plot.savefig("depth_hist.png")
	if readbal_hist is not None:
		balance_bar_plot = sns.countplot(
			x='ReadBalance', y='Count', data=readbal_hist)
		balance_bar_plot.savefig("readbalance_hist.png")
	if mapq_hist is not None:
		mapq_bar_plot = sns.countplot(
			x='Mapq', y='Count', data=mapq_hist)
		mapq_bar_plot.savefig("mapq_hist.png")

if __name__ == "__main__":
	main()
