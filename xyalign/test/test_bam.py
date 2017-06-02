import os
import pytest
import subprocess
from xyalign import bam

# Get directory
dir = os.path.dirname(__file__)

# Test if "samtools" available
try:
	a = subprocess.Popen(
		["samtools"], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	samtools_present = True
except OSError:
	samtools_precesnt = False


def teardown_module(function):
	if os.path.exists(os.path.join(dir, "header2.bam.bai")):
		os.remove(os.path.join(dir, "header2.bam.bai"))


def test_BamFile():
	test_header = bam.BamFile(
		os.path.join(dir, "header.bam"), "samtools", no_initial_index=True)
	test_header2 = bam.BamFile(
		os.path.join(dir, "header2.bam"), "samtools", no_initial_index=True)

	assert test_header.filepath == os.path.join(dir, "header.bam")
	assert test_header2.filepath == os.path.join(dir, "header2.bam")
	assert test_header.is_indexed() is True
	assert test_header2.is_indexed() is False


@pytest.mark.skipif(
	samtools_present is False,
	reason="samtools needs too be callable with 'samtools'")
def test_BamFile_setup():
	test_header3 = bam.BamFile(
		os.path.join(dir, "header2.bam"), "samtools", no_initial_index=False)

	assert test_header3.filepath == os.path.join(dir, "header2.bam")
	assert test_header3.is_indexed() is True


@pytest.mark.parametrize(
	"input,expected", [
		("chrM", 16571),
		("chr10", 135534747),
		("chrX", 155270560)])
def test_get_chrom_length(input, expected):
	test_header = bam.BamFile(
		os.path.join(dir, "header.bam"), "samtools", no_initial_index=True)
	assert test_header.get_chrom_length(input) == expected
	with pytest.raises(RuntimeError):
		test_header.get_chrom_length("foo")


def test_chromosome_lengths():
	tiny_header = bam.BamFile(
		os.path.join(dir, "tinyheader.bam"), "samtools", no_initial_index=False)
	assert tiny_header.chromosome_lengths() == (16571, 249250621)


def test_chromosome_names():
	tiny_header = bam.BamFile(
		os.path.join(dir, "tinyheader.bam"), "samtools", no_initial_index=False)
	assert tiny_header.chromosome_names() == ('chrM', 'chr1')
