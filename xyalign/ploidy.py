# Part of XYalign
# Collection of functions related to ploidy estimation and chromosome depth
# comparison

from __future__ import division
from __future__ import print_function
import csv
import logging
import numpy as np
from scipy.stats import ks_2samp
import time

# Create logger for ploidy submodule
ploidy_logger = logging.getLogger("xyalign.ploidy")


def permutation_test_chromosomes(
	data_frame, first_chrom, second_chrom, chrom_column,
	value_column, num_perms, output_file=None):
	"""
	Takes a dataframe and runs a permutation test comparing mean values
	of two chromosomes.

	data_frame is a pandas dataframe
	first_chrom is the name of the first chromosome in comparison
	second_chrom is the name of the second chromosome in comparison
	chrom_column is the name of the column containing chromosome names
	value_column is the name of the column containing the value of interest
	num_perms is the number of permutations to use
	output_file: if not none, will print results to this file

	Returns:
		A tuple containing (mean of first chrom, mean of second chrom, p-value)
	"""
	perm_start = time.time()
	ploidy_logger.info(
		"Running permutation test ({} reps) comparing ratio of {} over {}".format(
			num_perms, first_chrom, second_chrom))
	first_vals = data_frame[
		data_frame[chrom_column] == first_chrom][value_column]
	second_vals = data_frame[
		data_frame[chrom_column] == second_chrom][value_column]
	combined = np.append(first_vals, second_vals)

	first_mean = np.mean(first_vals)
	second_mean = np.mean(second_vals)

	observed = first_mean / second_mean
	perms = []
	for i in range(0, num_perms):
		np.random.shuffle(combined)
		first = np.mean(combined[:len(first_vals)])
		second = np.mean(combined[-len(second_vals):])
		perms.append(first / second)
	perms = np.asarray(perms)
	sig = len(np.where(perms > observed)) / num_perms
	if output_file is not None:
		a = [
			"{}_mean".format(first_chrom),
			"{}_mean".format(second_chrom),
			"{}_{}_ratio".format(first_chrom, second_chrom),
			"p_val_({}_/_{})".format(first_chrom, second_chrom),
			"perm_dist_2.5",
			"perm_dist_50",
			"perm_dist_97.5"]
		b = [
			"{}".format(first_mean),
			"{}".format(second_mean),
			"{}".format(observed),
			"{}".format(sig),
			"{}".format(np.percentile(perms, 2.5)),
			"{}".format(np.percentile(perms, 50)),
			"{}".format(np.percentile(perms, 97.5))]
		with open(output_file, "w") as f:
			w = csv.writer(f, dialect="excel-tab")
			w.writerows([a, b])
	ploidy_logger.info(
		"Permutations on {} and {} complete. Elapsed time: {} seconds".format(
			first_chrom, second_chrom, time.time() - perm_start))
	return (first_mean, second_mean, sig)


def ks_two_sample(
	data_frame, first_chrom, second_chrom, chrom_column,
	value_column, output_file=None):
	"""
	Takes a dataframe and runs a Two-sample Kolmogorov-Smirnov test on desired
	value column

	data_frame is a pandas dataframe
	first_chrom is the name of the first chromosome in comparison
	second_chrom is the name of the second chromosome in comparison
	chrom_column is the name of the column containing chromosome names
	value_column is the name of the column containing the value of interest
	num_perms is the number of permutations to use
	output_file: if not none, will print results to this file

	Returns:
		A tuple containing (ks_statistic, ks_pvalue)
	"""
	ks_start = time.time()
	ploidy_logger.info(
		"Running KS two sample test on {} and {}".format(
			first_chrom, second_chrom))
	first_vals = data_frame[
		data_frame[chrom_column] == first_chrom][value_column]
	second_vals = data_frame[
		data_frame[chrom_column] == second_chrom][value_column]

	first_mean = np.mean(first_vals)
	second_mean = np.mean(second_vals)
	mean_ratio = first_mean / second_mean

	result = ks_2samp(first_vals, second_vals)

	if output_file is not None:
		a = [
			"{}_mean".format(first_chrom),
			"{}_mean".format(second_chrom),
			"{}_{}_ratio".format(first_chrom, second_chrom),
			"ks_statistic",
			"p_val"]
		b = [
			"{}".format(first_mean),
			"{}".format(second_mean),
			"{}".format(mean_ratio),
			"{}".format(result[0]),
			"{}".format(result[1])]

		with open(output_file, "w") as f:
			w = csv.writer(f, dialect="excel-tab")
			w.writerows([a, b])
	ploidy_logger.info(
		"KS two sample test on {} and {} complete. Elapsed time: {} seconds".format(
			first_chrom, second_chrom, time.time() - ks_start))
	return result


def bootstrap(
	data_frame, first_chrom, second_chrom, chrom_column,
	value_column, num_reps, output_file=None):
	"""
	Takes a dataframe and bootstraps the 95 percent confidence interval
	of the mean ratio of measure for two chromosomes (chrom1 / chrom2).

	data_frame is a pandas dataframe
	first_chrom is the name of the first chromosome in comparison
	second_chrom is the name of the second chromosome in comparison
	chrom_column is the name of the column containing chromosome names
	value_column is the name of the column containing the value of interest
	num_perms is the number of permutations to use
	output_file: if not none, will print results to this file

	Returns:
		A tuple containing (0.025 percentile, 0.5 percentile, 0.975 percentile)
	"""
	boot_start = time.time()

	ks_start = time.time()
	ploidy_logger.info(
		"Bootstrapping mean depth ratio of {} over {}".format(
			first_chrom, second_chrom))
	first_vals = data_frame[
		data_frame[chrom_column] == first_chrom][value_column]
	second_vals = data_frame[
		data_frame[chrom_column] == second_chrom][value_column]

	first_vals = np.asarray(first_vals)
	second_vals = np.asarray(second_vals)

	first_mean = np.mean(first_vals)
	second_mean = np.mean(second_vals)
	mean_ratio = first_vals / second_vals

	samples = []
	for i in range(0, num_reps):
		first_boot = np.random.shuffle(first_vals)
		second_boot = np.random.shuffle(second_vals)
		samples.append(np.mean(first_boot / second_boot))

	samples = np.asarray(samples)

	if output_file is not None:
		a = [
			"{}_mean".format(first_chrom),
			"{}_mean".format(second_chrom),
			"{}_{}_ratio".format(first_chrom, second_chrom),
			"boot_2.5_percentile",
			"boot_50_percentile",
			"boot_97.5_percentile"]
		b = [
			"{}".format(first_mean),
			"{}".format(second_mean),
			"{}".format(mean_ratio),
			"{}".format(np.percentile(samples, 2.5)),
			"{}".format(np.percentile(samples, 50)),
			"{}".format(np.percentile(samples, 97.5))]
		with open(output_file, "w") as f:
			w = csv.writer(f, dialect="excel-tab")
			w.writerows([a, b])
	ploidy_logger.info(
		"Bootstrapping of {} and {} (ratio) complete. "
		"Elapsed time: {} seconds".format(
			first_chrom, second_chrom, time.time() - boot_start))
	return (mean_ratio, np.percentile(samples, 2.5), np.percentile(samples, 97.5))
