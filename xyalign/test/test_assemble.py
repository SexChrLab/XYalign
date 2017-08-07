import os
import pysam
import pytest
import subprocess

from xyalign import assemble
from xyalign import bam
from xyalign import reftools

# Get current directory
dir = os.path.dirname(__file__)

# Test if "samtools" available
try:
	a = subprocess.Popen(
		["samtools"], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	samtools_present = True
except OSError:
	samtools_present = False

# Test if "bwa" available
try:
	a = subprocess.Popen(
		["bwa"], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	bwa_present = True
except OSError:
	bwa_present = False

# Test if "sambamba" available
try:
	a = subprocess.Popen(
		["sambamba"], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	sambamba_present = True
except OSError:
	sambamba_present = False


def teardown_module(function):
	teardown_files = [
		"assemble_test.sorted.bam",
		"assemble_test.sorted.bam.bai",
		"test_assemble.sorted.bam",
		"test_assemble.sorted.bam.bai",
		"test_assemble.rg.sorted.bam",
		"test_assemble.rg.sorted.bam.bai"]
	for file_name in teardown_files:
		if os.path.exists(os.path.join(dir, file_name)):
			os.remove(os.path.join(dir, file_name))


@pytest.mark.skipif(
	all([samtools_present, bwa_present, sambamba_present]) is False,
	reason="samtools, bwa, and sambamba need too be callable with 'samtools', "
	"'bwa', and 'sambamba' respectively")
def test_bwa_mem_mapping_sambamba():
	# Test single end, no rg
	toy_fasta = reftools.RefFasta(os.path.join(dir, "toy.fasta"))
	test_bam = assemble.bwa_mem_mapping_sambamba(
		"bwa", "samtools", "sambamba", toy_fasta,
		os.path.join(dir, "assemble_test"), [os.path.join(dir, "toy_1.fastq")], 1,
		"None", [""], cram=False)
	assert os.path.exists(os.path.join(dir, "assemble_test.sorted.bam"))
	assert pysam.idxstats(
		os.path.join(
			dir, "assemble_test.sorted.bam")) == 'seq1\t20\t0\t0\nseq2\t40\t0\t0\n*\t0\t0\t3\n'
	# Test single end, with rg
	toy_fasta = reftools.RefFasta(os.path.join(dir, "toy.fasta"))
	test_bam = assemble.bwa_mem_mapping_sambamba(
		"bwa", "samtools", "sambamba", toy_fasta,
		os.path.join(dir, "test_assemble.rg"), [os.path.join(dir, "toy_1.fastq")], 1,
		"@RG\tID:{}".format("test"), [""], cram=False)
	assert os.path.exists(os.path.join(dir, "test_assemble.rg.sorted.bam"))
	assert pysam.idxstats(
		os.path.join(
			dir, "test_assemble.rg.sorted.bam")) == 'seq1\t20\t0\t0\nseq2\t40\t0\t0\n*\t0\t0\t3\n'
	header = pysam.view("-H", os.path.join(dir, "test_assemble.rg.sorted.bam"))
	assert header.find("@RG\\tID:test") != -1
	# Test raising error for inaccessible reference
	with pytest.raises(RuntimeError):
		fake_fasta = reftools.RefFasta(os.path.join(
			dir, "DOES_NOT_EXIST.fasta"), no_initial_index=True)
		test_bam = assemble.bwa_mem_mapping_sambamba(
			"bwa", "samtools", "sambamba", fake_fasta,
			os.path.join(dir, "test_assemble"),
			[os.path.join(dir, "toy_1.fastq")], 1,
			"None", [""], cram=False)
	# Test raising error for no fastqs
	toy_fasta = reftools.RefFasta(os.path.join(dir, "toy.fasta"))
	with pytest.raises(RuntimeError):
		test_bam = assemble.bwa_mem_mapping_sambamba(
			"bwa", "samtools", "sambamba", toy_fasta,
			os.path.join(dir, "test_assemble"), [], 1,
			"None", [""], cram=False)
	# Test raising error for inaccessible fastq
	with pytest.raises(RuntimeError):
		toy_fasta = reftools.RefFasta(os.path.join(dir, "toy.fasta"))
		test_bam = assemble.bwa_mem_mapping_sambamba(
			"bwa", "samtools", "sambamba", toy_fasta,
			os.path.join(dir, "test_assemble"),
			[os.path.join(dir, "toy_FOO_fake.fastq")], 1,
			"None", [""], cram=False)
