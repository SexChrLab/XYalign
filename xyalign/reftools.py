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
	def __init__(self, filepath, samtools="samtools", bwa="bwa"):
		""" Initiate object with associated filepath """
		self.filepath = filepath
		self.samtools = samtools
		self.bwa = bwa
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

	def index_bwa(self):
		"""
		Index reference using bwa

		bwa is path to bwa program (default is 'bwa')
		"""
		self.logger.info("Creating bwa indices for: {}".format(
			self.filepath))
		bwa_idx_start = time.time()
		rc = subprocess.call([self.bwa, "index", self.filepath])
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

	def check_bwa_index(self):
		"""
		Checks to see if bwa indices are newer than fasta.  Returns True if so,
		and False if any of the indices are the same age as the fasta or
		older.
		"""
		self.logger.info("Checking bwa indexing of {}".format(self.filepath))
		if (
			os.path.exists("{}.amb".format(self.filepath)) and
			os.path.exists("{}.ann".format(self.filepath)) and
			os.path.exists("{}.bwt".format(self.filepath)) and
			os.path.exists("{}.pac".format(self.filepath)) and
			os.path.exists("{}.sa".format(self.filepath))):
			if (
				os.stat("{}.amb".format(self.filepath)).st_size != 0 and
				os.stat("{}.ann".format(self.filepath)).st_size != 0 and
				os.stat("{}.bwt".format(self.filepath)).st_size != 0 and
				os.stat("{}.pac".format(self.filepath)).st_size != 0 and
				os.stat("{}.sa".format(self.filepath)).st_size != 0):
					bwa_idx_stamp = [
						os.path.getmtime("{}.amb".format(self.filepath)),
						os.path.getmtime("{}.ann".format(self.filepath)),
						os.path.getmtime("{}.bwt".format(self.filepath)),
						os.path.getmtime("{}.pac".format(self.filepath)),
						os.path.getmtime("{}.sa".format(self.filepath))]
			else:
				self.logger.info("One or more BWA indices empty")
				return False
		else:
			self.logger.info("Bwa index incomplete or absent")
			return False

		ref_stamp = os.path.getmtime(self.filepath)
		if all(x > ref_stamp for x in bwa_idx_stamp) is True:
			self.logger.info("BWA index present and newer than ref file")
			return True
		else:
			self.logger.info("BWA index is older than ref file")
			return False

	def conditional_index_bwa(self, bwa="bwa"):
		"""
		Like index_bwa, but only indexes if indices are the same age or older
		than the fasta.  Use index_bwa to force indexing.

		bwa is path to bwa program (default is 'bwa')
		"""
		if self.check_bwa_index is False:
			self.index_bwa()

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

	def mask_reference(bed_mask, output_fasta=None):
		"""
		Creates a new masked references by hardmasking regions included
		in the bed_mask

		bed_mask is a bed file of regions to mask (as N) in the new reference
		output_fasta is the full path to and filename of the output fasta

		Returns:
			Path to new, indexed, masked) fasta
		"""
		mask_start = time.time()
		reftools_logger.info("Masking {} using regions in {}".format(
			self.filepath, bed_mask))
		maskedpath = output_fasta
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
