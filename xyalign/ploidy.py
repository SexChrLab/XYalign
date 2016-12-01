# Part of XYalign
# Collection of functions related to ploidy estimation and chromosome depth
# comparison

from __future__ import division
from __future__ import print_function
import csv
import numpy as np


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
			"{}_{}_diff".format(first_chrom, second_chrom),
			"p_val_({}_/_{})".format(first_chrom, second_chrom),
			"perm_2.5",
			"perm_50",
			"perm_97.5"]
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
	return (
		first_mean, second_mean, sig, np.percentile(perms, 2.5),
		np.percentile(perms, 97.5))
