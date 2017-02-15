# Part of XYalign
# Collection of functions for mapping reads, processing bams, etc.
from __future__ import print_function
import logging
import os
import subprocess
import time


# Create logger for assemble submodule
assemble_logger = logging.getLogger("xyalign.assemble")


def bwa_mem_mapping_sambamba(
	bwa_path, samtools_path, sambamba_path, reference, output_prefix, fastqs,
	threads, read_group_line, bwa_params, cram=False):
	"""
	Maps reads to a reference genome using bwa mem.  If output is in bam format,
	will sort using sambamba, else will sort with samtools

	Parameters
	----------

	bwa_path : str
		The path to bwa
	samtools_path : str
		The path to samtools
	sambamba_path : str
		The path to sambamba
	reference : str
		The path to the reference genome (in fasta format)
	output_prefix : str
		The prefix (including path) to the desired output files
	fastqs : list
		Fastqs, e.g. ['sample_1.fastq', 'sample_2.fastq']
	threads : int
		The number of threads/cpus to use
	read_group_line : str
		Read group info for bwa to add
	bwa_params : list
		Bwa parameters to be joined into a string and inserted into the command
	cram : bool
		If True, will output a sorted cram, else a sorted bam. Default is False.

	Returns
	-------

	str
		Path to output bam file (indexed)

	"""
	map_start = time.time()

	fastqs = ' '.join(fastqs)
	ref_time = os.path.getmtime("{}".format(reference))
	assemble_logger.info(
		"Beginning steps mapping fastqs ({}) to reference ({}) "
		"using bwa_mem_mapping_sambamba".format(
			fastqs, reference))
	# Check that bwa index is not newer than reference (and re-index if it is)
	try:
		amb = os.path.getmtime("{}.amb".format(reference))
		ann = os.path.getmtime("{}.ann".format(reference))
		bwt = os.path.getmtime("{}.bwt".format(reference))
		pac = os.path.getmtime("{}.pac".format(reference))
		sa = os.path.getmtime("{}.sa".format(reference))
		if not all(x > ref_time for x in (amb, ann, bwt, pac, sa)):
			assemble_logger.info(
				"BWA indices older than reference file. Indexing")
			subprocess.call([bwa_path, "index", reference])
			assemble_logger.info(
				"BWA indexing of {} successful".format(reference))
	except:
		assemble_logger.info(
			"Could not find all BWA indices. Indexing")
		subprocess.call([bwa_path, "index", reference])
		assemble_logger.info(
			"BWA indexing of {} successful".format(reference))

	# Check that .fai is not newer than reference (and re-index if it is)
	try:
		faidx = os.path.getmtime("{}.fai".format(reference))
		if ref_time >= faidx:
			assemble_logger.info(
				"Reference index (fai) is older than reference. Indexing.")
			subprocess.call([samtools_path, "faidx", reference])
			assemble_logger.info(
				"Faidx indexing complete for {}".format(reference))
	except:
		assemble_logger.info("Could not find .fai index. Indexing")
		subprocess.call([samtools_path, "faidx", reference])
		assemble_logger.info("Faidx indexing complete for {}".format(reference))

	# BAM mapping
	if cram is False:
		command_line = "{} mem -t {} -R {} {} {} {} | {} fixmate -O bam - - | "\
			"{} sort -t {} -o {}_sorted.bam /dev/stdin".format(
				bwa_path, threads, repr(read_group_line), " ".join(bwa_params),
				reference, fastqs, samtools_path,
				sambamba_path, threads, output_prefix)
		assemble_logger.info(
			"Mapping reads with the command: {}".format(command_line))
		subprocess.call(command_line, shell=True)
		assemble_logger.info(
			"Indexing {}_sorted.bam".format(output_prefix))
		subprocess.call(
			[sambamba_path, "index", "-t", str(threads),
				"{}_sorted.bam".format(output_prefix)])
		assemble_logger.info(
			"Completed mapping for fastqs ({}) to reference ({}). "
			"Elapsed time: {} seconds".format(
				fastqs, reference, time.time() - map_start))
		return "{}_sorted.bam".format(output_prefix)

	# CRAM mapping
	# 	- note, doesn't currently support sambamba sorting (samtools instead)
	else:
		command_line = "{} mem -t {} -R {} {} {} {} | {} fixmate -O cram - - | "\
			"{} sort -O cram -o {}_sorted.cram -".format(
				bwa_path, threads, repr(read_group_line), " ".join(bwa_params),
				reference, fastqs, samtools_path,
				samtools_path, output_prefix)
		assemble_logger.info(
			"Mapping reads with the command: {}".format(command_line))
		subprocess.call(command_line, shell=True)
		assemble_logger.info(
			"Indexing {}_sorted.cram".format(output_prefix))
		subprocess.call(
			[samtools_path, "index", "{}_sorted.cram".format(output_prefix)])
		assemble_logger.info(
			"Completed mapping for fastqs ({}) to reference ({}). "
			"Elapsed time: {} seconds".format(
				fastqs, reference, time.time() - map_start))
		return "{}_sorted.cram".format(output_prefix)
