import gzip
import os
import pytest
import subprocess
from xyalign import bam

# Get current directory
dir = os.path.dirname(__file__)

# Test if "samtools" available
try:
	a = subprocess.Popen(
		["samtools"], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	samtools_present = True
except OSError:
	samtools_present = False

# Test if "repair.sh" available
try:
	a = subprocess.Popen(
		["repair.sh"], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	repairsh_present = True
except OSError:
	repairsh_present = False

# Test if "shuffle.sh" available
try:
	a = subprocess.Popen(
		["shuffle.sh"], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	shufflesh_present = True
except OSError:
	shufflesh_present = False


def read_bed(file):
	with open(file, "r") as f:
		return_list = []
		for line in f:
			return_list.append(line.strip().split())
	return return_list


def read_file_to_string(file, zipped):
	if zipped is False:
		with open(file, "r") as f:
			output = ""
			for line in f:
				output += line.strip()
	else:
		with gzip.open(file, "rb") as f:
			output = ""
			for line in f:
				output += line.strip()
	return output


def teardown_module(function):
	teardown_files = [
		"header2.bam.bai",
		"test1.bed",
		"test2.bed",
		"toy_test_rg.fastq",
		"toy_test_rg.fastq.gz",
		"toy.fastq",
		"toy.full_rg.list",
		"toy.header_lines_rg.list",
		"toy.rg_fastq_key.list",
		"toy.temp.fastq",
		"toy_test_rg.fastq",
		"toy_test_rg.fastq.gz",
		"toy.test_rg.temp.fastq",
		"toy_none.fastq",
		"toy_none.full_rg.list",
		"toy_none.header_lines_rg.list",
		"toy_none.rg_fastq_key.list",
		"toy_none.temp.fastq",
		"xmx_none_compress.fastq.gz",
		"xmx_none_compress.full_rg.list",
		"xmx_none_compress.header_lines_rg.list",
		"xmx_none_compress.rg_fastq_key.list",
		"xmx_none.fastq",
		"xmx_none.full_rg.list",
		"xmx_none.header_lines_rg.list",
		"xmx_none.rg_fastq_key.list"]
	for file_name in teardown_files:
		if os.path.exists(os.path.join(dir, file_name)):
			os.remove(os.path.join(dir, file_name))


def test_BamFile():
	test_header = bam.BamFile(
		os.path.join(dir, "header.bam"), "samtools", no_initial_index=True)
	test_header2 = bam.BamFile(
		os.path.join(dir, "header2.bam"), "samtools", no_initial_index=True)
	test_header3 = bam.BamFile(
		os.path.join(dir, "tinyheader2.bam"), "samtools", no_initial_index=True)
	assert test_header.filepath == os.path.join(dir, "header.bam")
	assert test_header2.filepath == os.path.join(dir, "header2.bam")
	assert test_header3.filepath == os.path.join(dir, "tinyheader2.bam")
	assert test_header.is_indexed() is True
	assert test_header2.is_indexed() is False
	assert test_header3.is_indexed() is True


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
	lengths = tiny_header.chromosome_lengths()
	assert type(lengths) == tuple
	assert lengths == (16571, 249250621)


def test_chromosome_names():
	tiny_header = bam.BamFile(
		os.path.join(dir, "tinyheader.bam"), "samtools", no_initial_index=False)
	names = tiny_header.chromosome_names()
	assert type(names) == tuple
	assert names == ('chrM', 'chr1')


def test_chromosome_bed():
	test_header = bam.BamFile(
		os.path.join(dir, "header.bam"), "samtools", no_initial_index=True)
	path1 = test_header.chromosome_bed(
		os.path.join(dir, "test1.bed"), ["chrX"])
	path2 = test_header.chromosome_bed(
		os.path.join(dir, "test2.bed"), ["chr1", "chrX"])
	assert path1 == os.path.join(dir, "test1.bed")
	assert path2 == os.path.join(dir, "test2.bed")
	assert read_bed(
		os.path.join(dir, "test1.bed")) == [["chrX", "0", "155270560"]]
	assert read_bed(
		os.path.join(dir, "test2.bed")) == [
			["chr1", "0", "249250621"], ["chrX", "0", "155270560"]]
	with pytest.raises(RuntimeError):
		path3 = test_header.chromosome_bed(
			os.path.join(dir, "test1.bed"), ["foo"])


@pytest.mark.parametrize(
	"input,expected", [
		(["chrM"], []),
		(["chr10", "chr11", "chrY"], []),
		(["chrX", "foo"], ["foo"]),
		(["foo"], ["foo"])])
def test_check_chrom_in_bam(input, expected):
	test_header = bam.BamFile(
		os.path.join(dir, "header.bam"), "samtools", no_initial_index=True)
	assert test_header.check_chrom_in_bam(input) == expected


@pytest.mark.skipif(
	all([samtools_present, repairsh_present, shufflesh_present]) is False,
	reason="samtools and repair.sh need too be callable with 'samtools'")
def test_strip_reads():
	no_rg_bam = bam.BamFile(
		os.path.join(dir, "toy.bam"), "samtools", no_initial_index=True)
	assert no_rg_bam.strip_reads(
		"repair.sh", "shuffle.sh", True, dir, "toy", [], "None", 0, False, "None") == [
			os.path.join(dir, "toy.rg_fastq_key.list"), None]
	assert read_file_to_string(
		os.path.join(dir, "toy.fastq"), False) == '@read1GGCCCC+FFFFFF@read2TTTTTA+FGGGGG@read3TTTGGG+BBBBBB'

	addrg_test = no_rg_bam.strip_reads(
		"repair.sh", "shuffle.sh", True, dir, "toy", [], "None", 0, False, "test_rg")
	assert addrg_test == [
		os.path.join(dir, "toy.rg_fastq_key.list"), os.path.join(
			dir, "toy.header_lines_rg.list")]
	assert os.path.exists(os.path.join(dir, "toy_test_rg.fastq")) is True
	assert read_file_to_string(
		os.path.join(dir, "toy_test_rg.fastq"), False) == '@read1GGCCCC+FFFFFF@read2TTTTTA+FGGGGG@read3TTTGGG+BBBBBB'

	compression_test = no_rg_bam.strip_reads(
		"repair.sh", "shuffle.sh", True, dir, "toy", [], "None", 2, False, "test_rg")
	assert os.path.exists(os.path.join(dir, "toy_test_rg.fastq.gz")) is True
	assert read_file_to_string(
		os.path.join(dir, "toy_test_rg.fastq.gz"), True) == '@read1GGCCCC+FFFFFF@read2TTTTTA+FGGGGG@read3TTTGGG+BBBBBB'

	none_rg_test = no_rg_bam.strip_reads(
		"repair.sh", "shuffle.sh", True, dir, "toy_none", [], "None", 0, False, "None")
	assert none_rg_test == [
		os.path.join(dir, "toy_none.rg_fastq_key.list"), None]
	assert os.path.exists(os.path.join(dir, "toy_none.fastq")) is True
	assert read_file_to_string(
		os.path.join(dir, "toy_none.fastq"), False) == '@read1GGCCCC+FFFFFF@read2TTTTTA+FGGGGG@read3TTTGGG+BBBBBB'

	xmx_single_test = no_rg_bam.strip_reads(
		"repair.sh", "shuffle.sh", True, dir, "xmx_none", [], "100m", 0, True, "None")
	assert xmx_single_test == [
		os.path.join(dir, "xmx_none.rg_fastq_key.list"), None]
	assert os.path.exists(os.path.join(dir, "xmx_none.fastq")) is True
	assert read_file_to_string(
		os.path.join(dir, "xmx_none.fastq"), False) == '@read1GGCCCC+FFFFFF@read2TTTTTA+FGGGGG@read3TTTGGG+BBBBBB'

	xmx_single_test_compression = no_rg_bam.strip_reads(
		"repair.sh", "shuffle.sh", True, dir, "xmx_none_compress", [], "100m", 2, True, "None")
	assert xmx_single_test_compression == [
		os.path.join(dir, "xmx_none_compress.rg_fastq_key.list"), None]
	assert os.path.exists(os.path.join(dir, "xmx_none_compress.fastq.gz")) is True
	assert read_file_to_string(
		os.path.join(dir, "xmx_none_compress.fastq.gz"), True) == '@read1GGCCCC+FFFFFF@read2TTTTTA+FGGGGG@read3TTTGGG+BBBBBB'
