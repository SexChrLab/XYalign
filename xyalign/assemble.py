# Part of XYalign
# Collection of functions for mapping reads, processing bams, etc.
from __future__ import print_function
import os
import subprocess


def bwa_mem_mapping_sambamba(
	bwa_path, samtools_path, sambamba_path, reference, output_prefix, fastqs,
	threads, read_group_line, cram=False):
	"""
	Maps reads to a reference genome using bwa mem.  If output is in bam format,
	will sort using sambamba, else will sort with samtools

	bwa_path is the path to bwa
	samtools_path is the path to samtools
	sambamba_path is the path to sambamba
	reference is the path to the reference genome (in fasta format)
	output_prefix is the path and prefix to the desired output files
	fastqs is a list of fastqs, e.g. ['sample_1.fastq', 'sample_2.fastq']
	threads is the number of threads/cpus to use
	read_group_line is a string containing read group info for bwa to add
	cram (default is False) - if True, will output a sorted cram file
	"""
	fastqs = ' '.join(fastqs)
	ref_time = os.path.getmtime("{}".format(reference))

	# Check that bwa index is not newer than reference (and re-index if it is)
	try:
		amb = os.path.getmtime("{}.amb".format(reference))
		ann = os.path.getmtime("{}.ann".format(reference))
		bwt = os.path.getmtime("{}.bwt".format(reference))
		pac = os.path.getmtime("{}.pac".format(reference))
		sa = os.path.getmtime("{}.sa".format(reference))
		if not all(x > ref_time for x in (amb, ann, bwt, pac, sa)):
			subprocess.call([bwa_path, "index", reference])
	except:
		subprocess.call([bwa_path, "index", reference])

	# Check that .fai is not newer than reference (and re-index if it is)
	try:
		faidx = os.path.getmtime("{}.fai".format(reference))
		if ref_time >= faidx:
			subprocess.call([samtools_path, "faidx", reference])
	else:
		subprocess.call([samtools_path, "faidx", reference])

	# BAM mapping
	if cram is False:
		command_line = "{} mem -t {} -R {} {} {} | {} fixmate -O bam - - | "\
			"{} sort -t {} -o {}_sorted.bam /dev/stdin".format(
				bwa_path, threads, repr(read_group_line), reference, fastqs, samtools_path,
				sambamba_path, threads, output_prefix)
		print("\Mapping reads with the command: {}\n".format(command_line))
		subprocess.call(command_line, shell=True)
		subprocess.call(
			[sambamba_path, "index", "-t", str(threads),
				"{}_sorted.bam".format(output_prefix)])
		return "{}_sorted.bam".format(output_prefix)

	# CRAM mapping
	else:
		command_line = "{} mem -t {} -R {} {} {} | {} fixmate -O cram - - | "\
			"{} sort -O cram -o {}_sorted.cram -".format(
				bwa_path, threads, repr(read_group_line), reference, fastqs, samtools_path,
				samtools_path, output_prefix)
		print("\Mapping reads with the command: {}\n".format(command_line))
		subprocess.call(command_line, shell=True)
		subprocess.call(
			[samtools_path, "index", "{}_sorted.cram".format(output_prefix)])
		return "{}_sorted.cram".format(output_prefix)


def bwa_mem_mapping(
	bwa_path, samtools_path, reference, output_prefix, fastqs,
	threads, cram=False):
	"""
	Maps reads to a reference genome using bwa mem.  Sorting done with samtools,
	and does not allow for adding read groups.

	bwa_path is the path to bwa
	samtools_path is the path to samtools
	reference is the path to the reference genome (in fasta format)
	output_prefix is the path and prefix to the desired output files
	fastqs is a list of fastqs, e.g. ['sample_1.fastq', 'sample_2.fastq']
	threads is the number of threads/cpus to use
	read_group_line is a string containing read group info for bwa to add
	cram (default is False) - if True, will output a sorted cram file
	"""
	fastqs = ' '.join(fastqs)
	subprocess.call(
		[bwa_path, "index", reference])
	if cram is False:
		command_line = "{} mem -t {} {} {} | {} fixmate -O bam - - | "\
			"{} sort -O bam -o {}_sorted.bam -".format(
				bwa_path, threads, reference, fastqs, samtools_path,
				samtools_path, output_prefix)
		print("\Mapping reads with the command: {}\n".format(command_line))
		subprocess.call(command_line, shell=True)
		subprocess.call(
			[samtools_path, "index", "{}_sorted.bam".format(output_prefix)])
		return "{}_sorted.bam".format(output_prefix)
	else:
		command_line = "{} mem -t {} {} {} | {} fixmate -O cram - - | "\
			"{} sort -O cram -o {}_sorted.cram -".format(
				bwa_path, threads, reference, fastqs, samtools_path,
				samtools_path, output_prefix)
		print("\Mapping reads with the command: {}\n".format(command_line))
		subprocess.call(command_line, shell=True)
		subprocess.call(
			[samtools_path, "index", "{}_sorted.cram".format(output_prefix)])
		return "{}_sorted.cram".format(output_prefix)
