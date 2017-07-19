import gzip
import os
import pytest
import subprocess
from xyalign import variants

# Get current directory
dir = os.path.dirname(__file__)

# Test if "tabix" available
try:
	a = subprocess.Popen(
		["tabix"], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	tabix_present = True
except OSError:
	tabix_present = False

# Test if "bgzip" available
try:
	a = subprocess.Popen(
		["bgzip"], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	bgzip_present = True
except OSError:
	bgzip_present = False

# Test if "platypus" available
try:
	a = subprocess.Popen(
		["platypus"], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	platypus_present = True
except OSError:
	platypus_present = False


def gunzip_file(input_file):
	with open(gzip.GzipFile(input_file, "rb")) as f:
		contents = f.read()
	with open(input_file[-3], "w") as o:
		o.write(contents)


def teardown_module(function):
	teardown_files = [
		os.path.join(dir, "platypus.vcf.gz.tbi")]
	for file_name in teardown_files:
		if os.path.exists(os.path.join(dir, file_name)):
			os.remove(os.path.join(dir, file_name))
	os.system(
		"gunzip {}".format(os.path.join(dir, "platypus.vcf.gz")))
	# gunzip_file(os.path.join(dir, "platypus.vcf.gz"))


def test_isbgzipped():
	test_vcf = variants.VCFFile(
		os.path.join(dir, "platypus.vcf"), no_initial_compress=True)
	assert test_vcf.is_bgzipped() is False


@pytest.mark.skipif(
	bgzip_present is False or tabix_present is False,
	reason="samtools needs too be callable with 'samtools'")
def test_VCFFile():
	test_vcf = variants.VCFFile(
		os.path.join(dir, "platypus.vcf"))
	assert os.path.exists(os.path.join(dir, "platypus.vcf.gz"))
	assert os.path.exists(os.path.join(dir, "platypus.vcf.gz.tbi"))
	assert test_vcf.is_bgzipped() is True


@pytest.mark.skipif(
	bgzip_present is False or tabix_present is False,
	reason="samtools needs too be callable with 'samtools'")
def test_read_balance_per_window():
	test_vcf2 = variants.VCFFile(
		os.path.join(dir, "platypus.vcf.gz"))
	pass
