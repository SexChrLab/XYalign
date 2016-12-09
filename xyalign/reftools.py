# Part of XYalign
# Functions, etc. related to reference fasta processing
import logging
import os
import subprocess
import time
import pybedtools

# Create logger for reftools submodule
reftools_logger = logging.getLogger("xyalign.reftools")


class RefFasta():
	"""
	A class for working with external reference fasta files
	"""
	def __init__(self, filepath, samtools="samtools"):
		""" Initiate object with associated filepath """
		self.filepath = filepath
		self.samtools = samtools
		self.logger = logging.getLogger("xyalign.reftools.RefFasta")
		self.logger.info("Creating a RefFasta instance for {}".format(
			self.filepath))
		if self.is_faidxed() is False:
			self.index_fai()

	def is_faidxed(self):
		"""
		Checks that fai index exists, is not empty and is newer than reference.
		If any case is False, return False.  Other wise, return True.
		"""
		self.logger.info("Checking fai indexing of {}".format(self.filepath))
		if os.path.exists("{}.fai".format(self.filepath)):
			if os.stat("{}.fai".format(self.filepath)).st_size != 0:
				idx_stamp = os.path.getmtime("{}.fai".format(self.filepath))
			else:
				self.logger.info("Fai index empty")
				return False
		else:
			self.logger.info("No fai index detected")
			return False

		ref_stamp = os.path.getmtime(self.filepath)
		if ref_stamp < idx_stamp:
			self.logger.info("Fai index present and newer than ref file")
			return True
		else:
			self.logger.info("Fai index is older than ref file")
			return False

	def index_fai(self):
		"""
		Create fai index for reference using samtools
		"""
		self.logger.info("Creating fai index for: {}".format(self.filepath))
		idx_start = time.time()
		rc = subprocess.call([self.samtools, "faidx", self.filepath])
		if rc == 0:
			self.logger.info(
				"Fai indexing complete. Elapsed time: {} seconds".format(
					time.time() - idx_start))
			return True
		else:
			self.logger.error("Unable to create fai index for {}. Exiting".format(
				self.filepath))
			raise RuntimeError("Unable to create faidx. Exiting")

	def index_bwa(self, bwa="bwa"):
		"""
		Index reference using bwa

		bwa is path to bwa program (default is 'bwa')
		"""
		self.logger.info("Creating bwa indices for: {}".format(
			self.filepath))
		bwa_idx_start = time.time()
		rc = subprocess.call([bwa, "index", self.filepath])
		if rc == 0:
			self.logger.info(
				"BWA indexing complete. Elapsed time: {} seconds".format(
					time.time() - bwa_idx_start))
			return True
		else:
			self.logger.error(
				"Unable to create bwa indices for {}. Exiting".format(
					self.filepath))
			raise RuntimeError("Unable to create bwa indicies. Exiting")

	def seq_dict(self, out_dict=None):
		"""
		Create sequence dictionary .dict file using samtools

		out_dict is the desired file name for the sequence dictionary -
			defaults to adding '.dict' to the fasta name
		"""
		self.logger.info("Creating sequence dictionary for: {}".format(
			self.filepath))
		if out_dict is None:
			out_dict = "{}.dict".format(self.filepath)
		dict_start = time.time()
		rc = subprocess.call(
			[self.samtools, "dict", "-o", out_dict, self.filepath])
		if rc == 0:
			self.logger.info(
				"Sequence dictionary complete. "
				"Elapsed time: {} seconds".format(
					time.time() - dict_start))
			return True
		else:
			self.logger.error(
				"Unable to create sequence dictionary for {}. "
				"Exiting".format(
					self.filepath))
			raise RuntimeError(
				"Unable to create sequence dictionary. Exiting")

	def mask_reference(new_ref_prefix, bed_mask):
		"""
		Creates a new masked references by hardmasking regions included
		in the bed_mask

		new_ref_prefix is the desired path to and prefix of the output files
		bed_mask is a bed file of regions to mask (as N) in the new reference

		Returns:
			Path to new, indexed, masked) fasta
		"""
		mask_start = time.time()
		reftools_logger.info("Masking {} using regions in {}".format(
			self.filepath, bed_mask))
		maskedpath = "{}.masked.fa".format(new_ref_prefix)
		b_fasta = pybedtools.BedTool(self.filepath)
		b_tool = pybedtools.BedTool(bed_mask)
		b = b_tool.mask_fasta(fi=b_fasta, fo=maskedpath)
		reftools_logger.info("Creating fai index for {}".format(maskedpath))
		subprocess.call(
			[self.samtools, "faidx", "{}".format(maskedpath)])
		reftools_logger.info(
			"Masking complete. Masked fasta output to {}. "
			"Elapsed time: {} seconds".format(
				maskedpath, time.time() - mask_start))
		return maskedpath

	def isolate_chroms(new_ref_prefix, chroms, bed_mask):
		"""
		Takes a reference fasta file and a list of chromosomes to include
		and outputs a new, indexed (and optionally masked) reference fasta.

		new_ref_prefix is the desired path to and prefix of the output files
		chroms should be a list of chromosomes to include in the output fasta
		bed_mask is a bed file of regions to mask (as N) in the new reference

		Returns:
			Path to new, indexed (optionally masked) fasta
		"""
		iso_start = time.time()
		reftools_logger.info("Isolating chromosomes ({}) from {}.".format(
			" ".join(chroms), self.filepath))
		outpath = "{}.fa".format(new_ref_prefix)
		if type(chroms) != list:
			chroms = list(chroms)
		if bed_mask is not None:
			maskedpath = "{}.masked.fa".format(new_ref_prefix)
			with open(outpath, "w") as f:
				subprocess.call(
					[self.samtools, "faidx", self.filepath, "{}".format(
						" ".join(chroms))], stdout=f)
			subprocess.call(
				[self.samtools, "faidx", outpath])
			b_fasta = pybedtools.BedTool(outpath)
			b_tool = pybedtools.BedTool(bed_mask)
			b = b_tool.mask_fasta(fi=b_fasta, fo=maskedpath)
			subprocess.call(
				[self.samtools, "faidx", "{}".format(maskedpath)])
			reftools_logger.info(
				"Isolating (masked) chromosomes complete. "
				"Elapsed time: {} seconds".format(time.time() - iso_start))
			return maskedpath
		else:
			with open(outpath, "w") as f:
				subprocess.call(
					[self.samtools, "faidx", self.filepath, "{}".format(
						" ".join(chroms))], stdout=f)
			subprocess.call([self.samtools, "faidx", outpath])
			reftools_logger.info(
				"Isolating (un-masked) chromosomes complete. "
				"Elapsed time: {} seconds".format(time.time() - iso_start))
			return outpath


# Legacy functions - saved in case useful in future
# def isolate_chromosomes_reference(
# 	samtools_path, reference_fasta, new_ref_prefix, chroms, bed_mask):
# 	"""
# 	Takes a reference fasta file and a list of chromosomes to include
# 	and outputs a new, indexed (and optionally masked) reference fasta.
#
# 	samtools_path is the path to samtools
# 	reference_fasta is the path to the reference genome (in fasta format)
# 	new_ref_prefix is the desired path to and prefix of the output files
# 	chroms should be a list of chromosomes to include in the output fasta
# 	bed_mask is a bed file of regions to mask (as N) in the new reference
#
# 	Returns:
# 		Path to new, indexed (optionally masked) fasta
# 	"""
# 	outpath = "{}.fa".format(new_ref_prefix)
# 	if type(chroms) != list:
# 		chroms = list(chroms)
# 	if bed_mask is not None:
# 		maskedpath = "{}.masked.fa".format(new_ref_prefix)
# 		with open(outpath, "w") as f:
# 			subprocess.call(
# 				[samtools_path, "faidx", reference_fasta, "{}".format(
# 					" ".join(chroms))], stdout=f)
# 		subprocess.call(
# 			[samtools_path, "faidx", outpath])
# 		b_fasta = pybedtools.BedTool(outpath)
# 		b_tool = pybedtools.BedTool(bed_mask)
# 		b = b_tool.mask_fasta(fi=b_fasta, fo=maskedpath)
# 		subprocess.call(
# 			[samtools_path, "faidx", "{}".format(maskedpath)])
# 	else:
# 		with open(outpath, "w") as f:
# 			subprocess.call(
# 				[samtools_path, "faidx", reference_fasta, "{}".format(
# 					" ".join(chroms))], stdout=f)
# 		subprocess.call([samtools_path, "faidx", outpath])
# 		return outpath
#
#
# def create_masked_reference(
# 	samtools_path, reference_fasta, new_ref_prefix, bed_mask):
# 	"""
# 	Creates a new masked references by hardmasking regions included
# 	in the bed_mask
#
# 	samtools_path is the path to samtools
# 	reference_fasta is the path to the reference genome (in fasta format)
# 	new_ref_prefix is the desired path to and prefix of the output files
# 	bed_mask is a bed file of regions to mask (as N) in the new reference
#
# 	Returns:
# 		Path to new, indexed, masked) fasta
# 	"""
# 	maskedpath = "{}.masked.fa".format(new_ref_prefix)
# 	b_fasta = pybedtools.BedTool(reference_fasta)
# 	b_tool = pybedtools.BedTool(bed_mask)
# 	b = b_tool.mask_fasta(fi=b_fasta, fo=maskedpath)
# 	subprocess.call(
# 		[samtools_path, "faidx", "{}".format(maskedpath)])
# 	return maskedpath
