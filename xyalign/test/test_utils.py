import os
import pytest
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
	if os.path.exists(os.path.join(dir, "merged.bed")):
		os.remove(os.path.join(dir, "merged.bed"))
	if os.path.exists(os.path.join(dir, "temp_test")):
		os.rmdir(os.path.join(dir, "temp_test"))


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
