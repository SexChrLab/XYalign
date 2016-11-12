# Part of XYalign
# Collection of functions for mapping reads, processing bams, etc.
import subprocess

def bwa_mem_mapping_sambamba(
	bwa_path, samtools_path, sambamba_path, reference, output_prefix, fastqs,
	threads, read_group_line, cram=False):
	""" Maps reads to a reference genome using bwa mem.
	"""
	fastqs = ' '.join(fastqs)
	subprocess.call("{} index {}".format(bwa_path, reference), shell=True)
	if cram is False:
		command_line = "{} mem -t {} -R {} {} {} | {} fixmate -O bam - - | "\
			"{} sort -t {} -o {}_sorted.bam /dev/stdin".format(
				bwa_path, threads, repr(read_group_line), reference, fastqs, samtools_path,
				sambamba_path, threads, output_prefix)
		subprocess.call(command_line, shell=True)
		subprocess.call(
			"{} index -t {} {}_sorted.bam".format(
				sambamba_path, threads, output_prefix), shell=True)
		return "{}_sorted.bam".format(output_prefix)
	else:
		command_line = "{} mem -t {} -R {} {} {} | {} fixmate -O cram - - | "\
			"{} sort -O cram -o {}_sorted.cram -".format(
				bwa_path, threads, repr(read_group_line), reference, fastqs, samtools_path,
				samtools_path, output_prefix)
		subprocess.call(command_line, shell=True)
		subprocess.call(
			"{} index {}_sorted.cram".format(
				samtools_path, output_prefix), shell=True)
		return "{}_sorted.cram".format(output_prefix)
