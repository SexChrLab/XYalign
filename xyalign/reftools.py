# Part of XYalign
# Functions, etc. related to reference fasta processing
from __future__ import division
import logging
import os
import subprocess
import time
import pybedtools
import pysam


# Create logger for reftools submodule
reftools_logger = logging.getLogger("xyalign.reftools")


class RefFasta():
	"""
	A class for working with external reference fasta files

	Attributes
	----------

	filepath : str
		Full path to external bam file.
	samtools : str
		Full path to samtools. Default = 'samtools'
	bwa : str
		Full path to bwa. Default = 'bwa'

	"""
	def __init__(
		self, filepath, samtools="samtools", bwa="bwa", no_initial_index=False):
		self.filepath = filepath
		self.samtools = samtools
		self.bwa = bwa
		self.logger = logging.getLogger("xyalign.reftools.RefFasta")
		self.logger.info("Creating a RefFasta instance for {}".format(
			self.filepath))
		if no_initial_index is False:
			if self.is_faidxed() is False:
				self.index_fai()

	def is_faidxed(self):
		"""
		Checks that fai index exists, is not empty and is newer than reference.

		Returns
		-------

		bool
			True if fai index exists and is newer than fasta, False otherwise.

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
		Create fai index for reference using samtools ('samtools faidx ref.fa')

		Returns
		-------

		bool
			True if successful

		Raises
		------

		RuntimeError
			If return code from external call is not 0

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
			logging.shutdown()
			raise RuntimeError("Unable to create faidx. Exiting")

	def index_bwa(self):
		"""
		Index reference using bwa

		Returns
		-------

		bool
			True if successful

		Raises
		------

		RuntimeError
			If return code from external call is not 0

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
			logging.shutdown()
			raise RuntimeError("Unable to create bwa indicies. Exiting")

	def check_bwa_index(self):
		"""
		Checks to see if bwa indices are newer than fasta.

		Returns
		-------

		bool
			True if indices exist and are newer than fasta. False otherwise.

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
		Indexes if indices are the same age or older than the fasta.
		Use RefFasta.index_bwa() to force indexing.

		Parameters
		----------

		bwa : str
			Path to bwa program (default is 'bwa')

		"""
		if self.check_bwa_index() is False:
			self.index_bwa()

	def check_seq_dict(self):
		"""
		Checks that sequence dictionary exists, is not empty and
		is newer than reference.

		Returns
		-------

		bool
			True if seq dict exists and is newer than fasta, False otherwise.

		"""
		self.logger.info(
			"Checking for sequence dictionary of {}".format(self.filepath))

		filename, file_extension = os.path.splitext(self.filepath)

		if os.path.exists("{}.dict".format(self.filepath)):
			if os.stat("{}.dict".format(self.filepath)).st_size != 0:
				idx_stamp = os.path.getmtime("{}.dict".format(self.filepath))
			else:
				self.logger.info("Seq dict empty")
				return False
		elif os.path.exists("{}.dict".format(filename)):
			if os.stat("{}.dict".format(filename)).st_size != 0:
				idx_stamp = os.path.getmtime("{}.dict".format(filename))
			else:
				self.logger.info("Seq dict empty")
				return False
		else:
			self.logger.info("No seq dict detected")
			return False

		ref_stamp = os.path.getmtime(self.filepath)
		if ref_stamp < idx_stamp:
			self.logger.info("Sequence dictionary present and newer than ref file")
			return True
		else:
			self.logger.info("Sequence dictionary is older than ref file")
			return False

	def seq_dict(self, out_dict=None):
		"""
		Create sequence dictionary .dict file using samtools

		Parameters
		----------

		out_dict : str
			The desired file name for the sequence dictionary -
			defaults to adding '.dict' to the fasta name

		Returns
		-------

		bool
			True if exit code of external call is 0.

		Raises
		------

		RuntimeError
			If external call exit code is not 0.

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
			logging.shutdown()
			raise RuntimeError(
				"Unable to create sequence dictionary. Exiting")

	def conditional_seq_dict(self):
		"""
		Creates sequence dictionary if .dict the same age or older than the fasta,
		or doesn't exist.

		Use RefFasta.seq_dict() to force creation.

		"""
		if self.check_seq_dict() is False:
			self.seq_dict()

	def mask_reference(self, bed_mask, output_fasta):
		"""
		Creates a new masked references by hardmasking regions included
		in the bed_mask

		Parameters
		----------

		bed_mask : str
			Bed file of regions to mask (as N) in the new reference
		output_fasta : str
			The full path to and filename of the output fasta

		Returns
		-------

		str
			Path to new (indexed and masked) fasta

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

	def isolate_chroms(self, new_ref_prefix, chroms, bed_mask=None):
		"""
		Takes a reference fasta file and a list of chromosomes to include
		and outputs a new, indexed (and optionally masked) reference fasta.

		Parameters
		----------

		new_ref_prefix : str
			The desired path to and prefix of the output files
		chroms : list
			Chromosomes to include in the output fasta
		bed_mask : str
			Bed file of regions to mask (as N) in the new reference

		Returns
		-------

		str
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

	def get_chrom_length(self, chrom):
		"""
		Extract chromosome length from fasta.

		Parameters
		----------

		chrom : str
			The name of the chromosome or scaffold.

		Returns
		-------

		length : int
			The length (integer) of the chromsome/scaffold

		Raises
		------

		RuntimeError
			If chromosome name not present in bam header

		"""
		fastafile = pysam.FastaFile(self.filepath)
		lengths = dict(zip(fastafile.references, fastafile.lengths))
		try:
			lens = lengths[chrom]
			fastafile.close()
			return lens
		except:
			self.logger.error(
				"{} not present in {}. Exiting.".format(
					chrom, self.filepath))
			logging.shutdown()
			raise RuntimeError(
				"Chromosome name not present in fasta. Exiting")

	def chromosome_bed(self, output_file, chromosome_list):
		"""
		Takes list of chromosomes and outputs a bed file with the
		length of each chromosome on each line
		(e.g., chr1    0   247249719).

		Parameters
		----------

		output_file : str
			Name of (including full path to) desired output file
		chromosome_list : list
			Chromosome/scaffolds to include

		Returns
		-------

		str
			output_file

		Raises
		------

		RuntimeError
			If chromosome name is not in fasta.

		"""
		c_bed_start = time.time()
		self.logger.info("Creating bed file with chromosome lengths for {}".format(
			" ".join(chromosome_list)))
		with open(output_file, "w") as f:
			for i in chromosome_list:
				try:
					lengths = self.get_chrom_length(i)
					f.write("{}\t{}\t{}\n".format(i, "0", lengths))
				except:
					self.logger.error(
						"Error finding chromosome length in {} "
						"(for bed file)".format(self.filepath))
					logging.shutdown()
					raise RuntimeError(
						"Error finding chromosome length in {}.  Check "
						"chromosome names and fasta IDs.".format(
							self.filepath))
		self.logger.info(
			"Bed file ({}) created. Elapsed time: {} seconds".format(
				output_file, time.time() - c_bed_start))
		return output_file

	def chromosome_lengths(self):
		"""
		Returns
		-------

		tuple
			Chromosome lengths ordered by sequence order in fasta

		"""
		fastafile = pysam.FastaFile(self.filepath)
		lengths = fastafile.lengths
		fastafile.close()
		# pysam's .lengths does not return a tuple (despite what is in the docs),
		# so, convert to tuple before returning.
		return tuple(lengths)

	def chromosome_names(self):
		"""
		Returns
		-------

		tuple
			Chromosome names ordered by sequence order in fasta

		"""
		fastafile = pysam.FastaFile(self.filepath)
		names = fastafile.references
		fastafile.close()
		# pysam's .lengths does not return a tuple (despite what is in the docs),
		# so, convert to tuple before returning.
		return tuple(names)
