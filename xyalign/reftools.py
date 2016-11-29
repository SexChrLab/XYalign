# Part of XYalign
# Functions, etc. related to reference fasta processing

import subprocess
import pybedtools


def isolate_chromosomes_reference(
	samtools_path, reference_fasta, new_ref_prefix, chroms, bed_mask):
	"""
	Takes a reference fasta file and a list of chromosomes to include
	and outputs a new, indexed (and optionally masked) reference fasta.

	samtools_path is the path to samtools
	reference_fasta is the path to the reference genome (in fasta format)
	new_ref_prefix is the desired path to and prefix of the output files
	chroms should be a list of chromosomes to include in the output fasta
	bed_mask is a bed file of regions to mask (as N) in the new reference

	Returns:
		Path to new, indexed (optionally masked) fasta
	"""
	outpath = "{}.fa".format(new_ref_prefix)
	if type(chroms) != list:
		chroms = list(chroms)
	if bed_mask is not None:
		maskedpath = "{}.masked.fa".format(new_ref_prefix)
		with open(outpath, "w") as f:
			subprocess.call(
				[samtools_path, "faidx", reference_fasta, "{}".format(
					" ".join(chroms))], stdout=f)
		subprocess.call(
			[samtools_path, "faidx", outpath])
		b_fasta = pybedtools.BedTool(outpath)
		b_tool = pybedtools.BedTool(bed_mask)
		b = b_tool.mask_fasta(fi=b_fasta, fo=maskedpath)
		subprocess.call(
			[samtools_path, "faidx", "{}".format(maskedpath)])
	else:
		with open(outpath, "w") as f:
			subprocess.call(
				[samtools_path, "faidx", reference_fasta, "{}".format(
					" ".join(chroms))], stdout=f)
		subprocess.call([samtools_path, "faidx", outpath])
		return outpath


def create_masked_reference(
	samtools_path, reference_fasta, new_ref_prefix, bed_mask):
	"""
	Creates a new masked references by hardmasking regions included
	in the bed_mask

	samtools_path is the path to samtools
	reference_fasta is the path to the reference genome (in fasta format)
	new_ref_prefix is the desired path to and prefix of the output files
	bed_mask is a bed file of regions to mask (as N) in the new reference

	Returns:
		Path to new, indexed, masked) fasta
	"""
	maskedpath = "{}.masked.fa".format(new_ref_prefix)
	b_fasta = pybedtools.BedTool(reference_fasta)
	b_tool = pybedtools.BedTool(bed_mask)
	b = b_tool.mask_fasta(fi=b_fasta, fo=maskedpath)
	subprocess.call(
		[samtools_path, "faidx", "{}".format(maskedpath)])
	return maskedpath
