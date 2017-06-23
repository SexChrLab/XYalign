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
		Read group info for bwa to add. If 'None', will not add read group.
	bwa_params : list
		Bwa parameters to be joined into a string and inserted into the command
	cram : bool
		If True, will output a sorted cram, else a sorted bam. Default is False.

	Returns
	-------

	str
		Path to output bam file (indexed)

	Raises
	------

	RuntimeError
		If fastq or reference files cannot be found

	"""
	map_start = time.time()
	fastq_status = [os.path.exists(x) for x in fastqs]
	if all(fastq_status) is False:
		assemble_logger.error(
			"One or more fastq files cannot be found. Check paths.")
		raise RuntimeError("One or more fastq files cannot be found. Check paths.")
	if os.path.exists(reference) is False:
		assemble_logger.error("Reference file cannot be found. Check path.")
		raise RuntimeError("Reference file cannot be found. Check path.")
	if len(fastqs) == 2:
		single = False
	elif len(fastqs) == 1:
		single = True
	else:
		assemble_logger.error(
			"{} fastq files provided. Please provide only one or two".format(
				len(fastqs)))
		raise RuntimeError(
			"{} fastq files provided. Please provide only one or two".format(
				len(fastqs)))
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
		output_file = "{}.sorted.bam".format(output_prefix)
		if single is False:
			if read_group_line != "None":
				command_line = "{} mem -t {} -R {} {} {} {} | {} fixmate -O bam - - | "\
					"{} sort -t {} -o {} /dev/stdin".format(
						bwa_path, threads, repr(read_group_line), " ".join(bwa_params),
						reference, fastqs, samtools_path,
						sambamba_path, threads, output_file)
			else:
				command_line = "{} mem -t {} {} {} {} | {} fixmate -O bam - - | "\
					"{} sort -t {} -o {} /dev/stdin".format(
						bwa_path, threads, " ".join(bwa_params),
						reference, fastqs, samtools_path,
						sambamba_path, threads, output_file)
		else:
			if read_group_line != "None":
				command_line = "{} mem -t {} -R {} {} {} {} | {} view -hb - | "\
					"{} sort -t {} -o {} /dev/stdin".format(
						bwa_path, threads, repr(read_group_line), " ".join(bwa_params),
						reference, fastqs, samtools_path,
						sambamba_path, threads, output_file)
			else:
				command_line = "{} mem -t {} {} {} {} | {} view -hb - | "\
					"{} sort -t {} -o {} /dev/stdin".format(
						bwa_path, threads, " ".join(bwa_params),
						reference, fastqs, samtools_path,
						sambamba_path, threads, output_file)
		assemble_logger.info(
			"Mapping reads with the command: {}".format(command_line))
		subprocess.call(command_line, shell=True)
		assemble_logger.info(
			"Indexing {}".format(output_file))
		subprocess.call(
			[sambamba_path, "index", "-t", str(threads),
				"{}".format(output_file)])
		assemble_logger.info(
			"Completed mapping for fastqs ({}) to reference ({}). "
			"Elapsed time: {} seconds".format(
				fastqs, reference, time.time() - map_start))
		return output_file

	# CRAM mapping
	# 	- note, doesn't currently support sambamba sorting (samtools instead)
	else:
		output_file = "{}.sorted.cram".format(output_prefix)
		if read_group_line != "None":
			command_line = "{} mem -t {} -R {} {} {} {} | {} fixmate -O cram - - | "\
				"{} sort -O cram -o {} -".format(
					bwa_path, threads, repr(read_group_line), " ".join(bwa_params),
					reference, fastqs, samtools_path,
					samtools_path, output_file)
		else:
			command_line = "{} mem -t {} {} {} {} | {} fixmate -O cram - - | "\
				"{} sort -O cram -o {} -".format(
					bwa_path, threads, " ".join(bwa_params),
					reference, fastqs, samtools_path,
					samtools_path, output_file)
		assemble_logger.info(
			"Mapping reads with the command: {}".format(command_line))
		subprocess.call(command_line, shell=True)
		assemble_logger.info(
			"Indexing {}".format(output_file))
		subprocess.call(
			[samtools_path, "index", "{}".format(output_file)])
		assemble_logger.info(
			"Completed mapping for fastqs ({}) to reference ({}). "
			"Elapsed time: {} seconds".format(
				fastqs, reference, time.time() - map_start))
		return output_file
