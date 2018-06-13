import gzip
import os
import pytest
import subprocess
import numpy as np
from xyalign import bam, variants

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
		"platypus.vcf.gz.tbi",
		"plot_test_chr19_ReadBalance_GenomicScatter.pdf",
		"plot_test_chr19_ReadBalance_Hist.pdf",
		"scatter_test_chr19_ReadBalance_GenomicScatter.pdf",
		"hist_test_chr19_ReadBalance_Hist.pdf",
		"hg19_header_rb.csv",
		"hg19_header_rb_target_test.csv",
		"plot_test_chr19_Window_Read_balance_GenomicScatter.pdf",
		"plot_test_chr19_Window_Variant_Count_GenomicScatter.pdf"]
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
	a = test_vcf.plot_variants_per_chrom(
		["chr19"], "test", os.path.join(dir, "plot_test"), 30, 30, 8, 4, 0.5,
		bam.BamFile(os.path.join(dir, "hg19_header.bam")), "platypus", True,
		os.path.join(dir, "hg19_header_rb.csv"), 1, 10000, 1)
	assert os.path.exists(
		os.path.join(dir, "plot_test_chr19_ReadBalance_GenomicScatter.pdf"))
	assert os.path.exists(
		os.path.join(dir, "plot_test_chr19_ReadBalance_Hist.pdf"))
	assert os.path.exists(
		os.path.join(dir, "hg19_header_rb.csv"))
	parsed = test_vcf.parse_platypus_VCF(30, 30, 8, "chr19")
	positions = [
		156497, 308662, 311825, 312020, 312143, 325481, 326481,
		327683, 327998, 334170, 335946, 362283, 362411, 362432,
		366804, 366846, 367313, 371396, 372551, 372867, 373618,
		374453, 374492, 474588, 474607, 501580, 501701, 501738,
		535634, 536038, 536068, 536088, 536090, 536107, 536878,
		536900, 541685, 547224, 547491, 547555, 549492, 549509,
		549595, 549678, 580268, 580333, 580353, 580665, 580793,
		580810, 582253, 582468, 582505]
	qualities = [
		1802.0, 138.0, 1467.0, 648.0, 538.0, 699.0, 850.0, 528.0,
		2965.0, 1370.0, 1350.0, 2019.0, 1071.0, 383.0, 386.0, 1668.0,
		2330.0, 512.0, 1735.0, 610.0, 493.0, 1412.0, 2505.0, 196.0,
		153.0, 359.0, 394.0, 36.0, 714.0, 608.0, 600.0, 531.0, 588.0,
		201.0, 1238.0, 1479.0, 546.0, 874.0, 199.0, 486.0, 517.0,
		235.0, 986.0, 2965.0, 196.0, 420.0, 338.0, 176.0, 303.0,
		372.0, 320.0, 375.0, 257.0]
	read_balances = [
		0.35323383084577115, 0.6363636363636364, 0.48672566371681414,
		0.95, 0.9473684210526315, 0.45454545454545453, 0.64, 1.0,
		1.0, 0.84375, 1.0, 0.5070422535211268, 0.6166666666666667,
		0.36, 0.4838709677419355, 1.0, 1.0, 0.6363636363636364, 0.53125,
		0.6470588235294118, 0.4166666666666667, 0.45112781954887216,
		0.9764705882352941, 0.55, 0.4782608695652174, 0.46875,
		0.16049382716049382, 0.1411764705882353, 1.0, 0.17391304347826086,
		0.3235294117647059, 0.5333333333333333, 0.5, 0.18181818181818182,
		0.9743589743589743, 0.9791666666666666, 0.9444444444444444, 1.0,
		0.6, 0.7142857142857143, 1.0, 0.4090909090909091, 0.9722222222222222,
		1.0, 0.5714285714285714, 0.6666666666666666, 0.5483870967741935,
		0.391304347826087, 0.43333333333333335, 0.27586206896551724,
		0.45714285714285713, 0.7, 0.56]
	gqs = [
		99, 75, 99, 51, 48, 99, 99, 45, 99, 87, 99, 99, 99, 99,
		99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 32,
		63, 99, 99, 99, 99, 99, 99, 99, 45, 75, 99, 99, 45, 99,
		90, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99]
	dps = [
		201.0, 11.0, 113.0, 20.0, 19.0, 66.0, 50.0, 17.0, 98.0,
		32.0, 42.0, 142.0, 60.0, 50.0, 31.0, 49.0, 68.0, 33.0,
		128.0, 34.0, 36.0, 133.0, 85.0, 20.0, 23.0, 32.0, 81.0,
		85.0, 22.0, 23.0, 34.0, 30.0, 30.0, 22.0, 39.0, 48.0,
		18.0, 28.0, 15.0, 14.0, 16.0, 22.0, 36.0, 106.0, 14.0,
		21.0, 31.0, 23.0, 30.0, 29.0, 35.0, 20.0, 25.0]
	assert parsed[0] == positions
	assert parsed[1] == qualities
	assert parsed[2] == read_balances
	assert parsed[3] == gqs
	assert parsed[4] == dps
	val = variants.plot_read_balance(
		"chr19", positions, read_balances, "scatter_test",
		os.path.join(dir, "scatter_test"), 4, 0.5, True, bam.BamFile(
			os.path.join(dir, "hg19_header.bam")).get_chrom_length("chr19"))
	assert val == 0
	assert os.path.exists(
		os.path.join(dir, "scatter_test_chr19_ReadBalance_GenomicScatter.pdf"))
	variants.hist_read_balance(
		"chr19", read_balances, "hist_test", True, os.path.join(dir, "hist_test"))
	assert os.path.exists(
		os.path.join(dir, "hist_test_chr19_ReadBalance_Hist.pdf"))
	val = variants.hist_read_balance(
		"chr19", [0, 0, 1], "hist_test_broken", True,
		os.path.join(dir, "hist_test_broken"))
	assert val == 1
	assert os.path.exists(
		os.path.join(dir, "hist_test_broken_chr19_ReadBalance_Hist.pdf")) is False
	target_test = variants.read_balance_per_window(
		"chr19", positions, read_balances, "target_test", True, 1,
		None, os.path.join(dir, "platypus.bed"))
	assert len(target_test) == 27


def check_hist_error():
	val = variants.hist_read_balance(
		"chr19", np.asarray([0, 0, 1]), "hist_test_broken",
		os.path.join(dir, "hist_test_broken"))
	assert val == 1
	assert os.path.exists(
		os.path.join(dir, "hist_test_broken_chr19_ReadBalance_Hist.pdf")) is False


@pytest.mark.skipif(
	bgzip_present is False or tabix_present is False,
	reason="samtools needs too be callable with 'samtools'")
def test_read_balance_per_window():
	test_vcf2 = variants.VCFFile(
		os.path.join(dir, "platypus.vcf.gz"))
	pass
