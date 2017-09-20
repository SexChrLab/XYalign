import os
import pytest
import subprocess
from xyalign import reftools

# Get current directory
dir = os.path.dirname(__file__)

# Test if "samtools" and "bwa" are available
try:
	a = subprocess.Popen(
		["samtools"], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	b = subprocess.Popen(
		["bwa"], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	tools_present = True
except OSError:
	tools_present = False


def read_file_to_list(input_file):
	file_contents = []
	with open(input_file, "r") as f:
		for line in f:
			file_contents.append(line.strip())
	return file_contents


def read_bed(file):
	with open(file, "r") as f:
		return_list = []
		for line in f:
			return_list.append(line.strip().split())
	return return_list


def teardown_module(function):
	teardown_files = [
		"toy2.fasta.amb",
		"toy2.fasta.ann",
		"toy2.fasta.bwt",
		"toy2.fasta.fai",
		"toy2.fasta.pac",
		"toy2.fasta.sa",
		"toy2.fasta.dict",
		"toy.masked.fasta",
		"toy.masked.fasta.fai",
		"newfasta.fa",
		"newfasta.fa.fai",
		"newfasta.masked.fa",
		"newfasta.masked.fa.fai",
		"reftest1.bed",
		"reftest2.bed"]
	for file_name in teardown_files:
		if os.path.exists(os.path.join(dir, file_name)):
			os.remove(os.path.join(dir, file_name))


def test_RefFasta():
	test_fasta = reftools.RefFasta(
		os.path.join(dir, "toy.fasta"), "samtools", "bwa", no_initial_index=True)
	test_fasta2 = reftools.RefFasta(
		os.path.join(dir, "toy2.fasta"), "samtools", "bwa", no_initial_index=True)
	os.utime(os.path.join(dir, "toy.dict"), None)
	assert test_fasta.filepath == os.path.join(dir, "toy.fasta")
	assert test_fasta.is_faidxed() is True
	assert test_fasta.check_bwa_index() is True
	assert test_fasta.check_seq_dict() is True
	assert test_fasta2.filepath == os.path.join(dir, "toy2.fasta")
	assert test_fasta2.is_faidxed() is False
	assert test_fasta2.check_bwa_index() is False
	assert test_fasta2.check_seq_dict() is False


@pytest.mark.skipif(
	tools_present is False,
	reason="samtools needs too be callable with 'samtools' "
	"and bwa needs to be callable with 'bwa'")
def test_RefFasta_setup():
	test_fasta2 = reftools.RefFasta(
		os.path.join(dir, "toy2.fasta"), "samtools", "bwa", no_initial_index=False)
	test_fasta2.conditional_index_bwa()
	test_fasta2.conditional_seq_dict()
	assert os.path.exists(os.path.join(dir, "toy2.fasta.amb"))
	assert os.path.exists(os.path.join(dir, "toy2.fasta.ann"))
	assert os.path.exists(os.path.join(dir, "toy2.fasta.bwt"))
	assert os.path.exists(os.path.join(dir, "toy2.fasta.fai"))
	assert os.path.exists(os.path.join(dir, "toy2.fasta.pac"))
	assert os.path.exists(os.path.join(dir, "toy2.fasta.sa"))
	assert os.path.exists(os.path.join(dir, "toy2.fasta.dict"))


def test_mask_reference():
	test_fasta = reftools.RefFasta(
		os.path.join(dir, "toy.fasta"), "samtools", "bwa", no_initial_index=True)
	new_path = test_fasta.mask_reference(
		os.path.join(dir, "toy_mask.bed"), os.path.join(dir, "toy.masked.fasta"))
	assert new_path == os.path.join(dir, "toy.masked.fasta")
	assert os.path.exists(os.path.join(dir, "toy.masked.fasta.fai"))
	contents = read_file_to_list(os.path.join(dir, "toy.masked.fasta"))
	assert contents[1] == "NNNNNNNNNNNNNNNNNNNN"


def test_isolate_chroms():
	test_fasta = reftools.RefFasta(
		os.path.join(dir, "toy.fasta"), "samtools", "bwa", no_initial_index=True)
	seq1 = test_fasta.isolate_chroms(
		os.path.join(dir, "newfasta"), ["seq1"], bed_mask=None)
	seq1_masked = test_fasta.isolate_chroms(
		os.path.join(
			dir, "newfasta"), ["seq1"], bed_mask=os.path.join(
				dir, "toy_mask.bed"))
	assert seq1 == os.path.join(dir, "newfasta.fa")
	assert seq1_masked == os.path.join(dir, "newfasta.masked.fa")
	assert os.path.exists(os.path.join(dir, "newfasta.fa.fai"))
	assert os.path.exists(os.path.join(dir, "newfasta.masked.fa.fai"))
	contents1 = read_file_to_list(os.path.join(dir, "newfasta.fa"))
	contents2 = read_file_to_list(os.path.join(dir, "newfasta.masked.fa"))
	assert contents1[1] == "AAAATTTTAAAATTTTGGGG"
	assert contents2[1] == "NNNNNNNNNNNNNNNNNNNN"


def test_chromosome_lengths():
	test_fasta = reftools.RefFasta(
		os.path.join(dir, "toy.fasta"), "samtools", "bwa", no_initial_index=True)
	lengths = test_fasta.chromosome_lengths()
	assert type(lengths) == tuple
	assert lengths == (20, 40)


def test_chromosome_names():
	test_fasta = reftools.RefFasta(
		os.path.join(dir, "toy.fasta"), "samtools", "bwa", no_initial_index=True)
	names = test_fasta.chromosome_names()
	assert type(names) == tuple
	assert names == ("seq1", "seq2")


def test_get_chrom_length():
	test_fasta = reftools.RefFasta(
		os.path.join(dir, "toy.fasta"), "samtools", "bwa", no_initial_index=True)
	lengths = test_fasta.get_chrom_length("seq1")
	assert type(lengths) == int
	assert lengths == 20
	with pytest.raises(RuntimeError):
		lengths = test_fasta.get_chrom_length("foo")


def test_chromosome_bed():
	test_fasta = reftools.RefFasta(
		os.path.join(dir, "toy.fasta"), "samtools", "bwa", no_initial_index=True)
	path1 = test_fasta.chromosome_bed(
		os.path.join(dir, "reftest1.bed"), ["seq1"])
	path2 = test_fasta.chromosome_bed(
		os.path.join(dir, "reftest2.bed"), ["seq1", "seq2"])
	assert path1 == os.path.join(dir, "reftest1.bed")
	assert path2 == os.path.join(dir, "reftest2.bed")
	assert read_bed(
		os.path.join(dir, "reftest1.bed")) == [["seq1", "0", "20"]]
	assert read_bed(
		os.path.join(dir, "reftest2.bed")) == [
			["seq1", "0", "20"], ["seq2", "0", "40"]]
	with pytest.raises(RuntimeError):
		path3 = test_fasta.chromosome_bed(
			os.path.join(dir, "reftest1.bed"), ["foo"])
