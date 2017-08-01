import os
import pybedtools
import pytest
import pandas as pd
from xyalign import bam
from xyalign import reftools
from xyalign import utils

# Get current directory
dir = os.path.dirname(__file__)


def read_bed(file):
	with open(file, "r") as f:
		return_list = []
		for line in f:
			return_list.append(line.strip().split())
	return return_list


def teardown_module(function):
	teardown_files = ["merged.bed", "out_test.bed"]
	teardown_dirs = ["temp_test"]
	for file_name in teardown_files:
		if os.path.exists(os.path.join(dir, file_name)):
			os.remove(os.path.join(dir, file_name))
	for direc in teardown_dirs:
		if os.path.exists(os.path.join(dir, direc)):
			os.rmdir(os.path.join(dir, direc))


def test_validate_external_prog():
	assert utils.validate_external_prog("echo", "echo") == 0
	with pytest.raises(OSError):
		utils.validate_external_prog(
			"this_program_should_not_exist", "this_program_should_not_exist")


def test_validate_dir_exist():
	assert utils.validate_dir(dir, "temp_test") is False
	assert utils.validate_dir(dir, "temp_test") is True


def test_check_bam_fasta_compatibility():
	test_bam = bam.BamFile(
		os.path.join(dir, "toy.bam"), "samtools", no_initial_index=True)
	test2_bam = bam.BamFile(
		os.path.join(dir, "tinyheader.bam"), "samtools", no_initial_index=True)
	test_fasta = reftools.RefFasta(
		os.path.join(dir, "toy.fasta"), "samtools", "bwa", no_initial_index=True)
	test2_fasta = reftools.RefFasta(
		os.path.join(dir, "toy3.fasta"), "samtools", "bwa", no_initial_index=True)
	assert utils.check_bam_fasta_compatibility(test_bam, test_fasta) is True
	assert utils.check_bam_fasta_compatibility(test2_bam, test_fasta) is False
	assert utils.check_bam_fasta_compatibility(test_bam, test2_fasta) is False


def test_merge_bed_files():
	path = utils.merge_bed_files(
		os.path.join(dir, "merged.bed"),
		os.path.join(dir, "toy_mask.bed"),
		os.path.join(dir, "toy_mask2.bed"))
	assert path == os.path.join(dir, "merged.bed")
	assert os.path.exists(os.path.join(dir, "merged.bed"))
	assert read_bed(os.path.join(dir, "merged.bed")) == [["seq1", "0", "25"]]


def test_make_region_lists_genome_filters():
	df = pd.read_csv(os.path.join(dir, "dataframe.csv"), sep="\t")
	genome_filters = utils.make_region_lists_genome_filters(
		df, 0, 10, 20)
	df2 = df.loc[df["depth"] >= 38.6496]
	df3 = df2.loc[df2["depth"] <= 77.2992]
	assert df3["depth"].mean() == genome_filters[0]["depth"].mean()


def test_make_region_lists_chromosome_filters():
	df = pd.read_csv(os.path.join(dir, "dataframe.csv"), sep="\t")
	chrom_filters = utils.make_region_lists_chromosome_filters(
		df, 0, 10, 100000)
	chr19 = df.loc[df["chrom"] == "chr19"]
	chrX = df.loc[df["chrom"] == "chrX"]
	chrY = df.loc[df["chrom"] == "chrY"]
	chr19_2 = chr19.loc[df["depth"] >= 98.178]
	chrX_2 = chrX.loc[df["depth"] >= 28.32]
	chrY_2 = chrY.loc[df["depth"] >= 0.194]
	df2 = pd.concat([chr19_2, chrX_2, chrY_2])
	assert df2["depth"].mean() == chrom_filters[0]["depth"].mean()


def test_output_bed():
	df = pd.read_csv(os.path.join(dir, "dataframe.csv"), sep="\t")
	return_val = utils.output_bed(os.path.join(dir, "out_test.bed"), *[df])
	assert return_val == 0
	old_bed = read_bed(os.path.join(dir, "out.bed"))
	new_bed = read_bed(os.path.join(dir, "out_test.bed"))
	assert new_bed == old_bed
