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


def teardown_module(function):
	for file in [
		os.path.join(dir, "toy2.fasta.amb"),
		os.path.join(dir, "toy2.fasta.ann"),
		os.path.join(dir, "toy2.fasta.bwt"),
		os.path.join(dir, "toy2.fasta.fai"),
		os.path.join(dir, "toy2.fasta.pac"),
		os.path.join(dir, "toy2.fasta.sa"),
		os.path.join(dir, "toy2.fasta.dict"),
		os.path.join(dir, "toy.masked.fasta"),
		os.path.join(dir, "toy.masked.fasta.fai"),
		os.path.join(dir, "newfasta.fa"),
		os.path.join(dir, "newfasta.fa.fai"),
		os.path.join(dir, "newfasta.masked.fa"),
		os.path.join(dir, "newfasta.masked.fa.fai")]:
		if os.path.exists(file):
			os.remove(file)


def test_RefFasta():
	test_fasta = reftools.RefFasta(
		os.path.join(dir, "toy.fasta"), "samtools", "bwa", no_initial_index=True)
	test_fasta2 = reftools.RefFasta(
		os.path.join(dir, "toy2.fasta"), "samtools", "bwa", no_initial_index=True)
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
