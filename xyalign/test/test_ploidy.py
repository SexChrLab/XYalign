import os
import pandas as pd
from xyalign import ploidy

# Get current directory
dir = os.path.dirname(__file__)


def teardown_module(function):
	teardown_files = [
		"perm.txt",
		"boot.txt",
		"ks.txt"]
	for file_name in teardown_files:
		if os.path.exists(os.path.join(dir, file_name)):
			os.remove(os.path.join(dir, file_name))


def test_permutation_test_chromosomes():
	df = pd.read_csv(os.path.join(dir, "dataframe.csv"), sep="\t")
	results = ploidy.permutation_test_chromosomes(
		df, "chr19", "chrX", "chrom", "depth", 100, os.path.join(dir, "perm.txt"))
	assert 9.81 < results[0] < 9.82
	assert 2.83 < results[1] < 2.84
	assert os.path.exists(os.path.join(dir, "perm.txt"))


def test_ks_two_sample():
	df = pd.read_csv(os.path.join(dir, "dataframe.csv"), sep="\t")
	results = ploidy.ks_two_sample(
		df, "chr19", "chrX", "chrom", "depth", os.path.join(dir, "ks.txt"))
	print(results)
	assert 0.23 < results[0] < 0.25
	assert os.path.exists(os.path.join(dir, "ks.txt"))


def test_bootstrap():
	df = pd.read_csv(os.path.join(dir, "dataframe.csv"), sep="\t")
	results = ploidy.bootstrap(
		df, "chr19", "chrX", "chrom", "depth", 100, os.path.join(dir, "boot.txt"))
	assert 3.45 < results[0] < 3.47
	assert os.path.exists(os.path.join(dir, "boot.txt"))
