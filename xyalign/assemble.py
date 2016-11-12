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
		p1 = subprocess.Popen([
            bwa_path, "mem", "-t", threads, "-R", repr(read_group_line),
            reference, fastqs], stdout=subprocess.PIPE)
        p2 = subprocess.Popen([
            samtools_path, "fixmate", "-O", "bam", "-", "-"],
            stdin=p1.stdout, stdout=subprocess.PIPE)
        p3 = subprocess.Popen([
            sambamba_path, "sort", "-t", threads, "-o", "{}_sorted.bam".format(
                output_prefix), "/dev/stdin"], stdin=p2.stdout)
        p3.communicate()
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
