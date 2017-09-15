# Part of XYalign
# Functions for calling and processing variants
from __future__ import division
from __future__ import print_function
import logging
import numpy as np
import os
import pandas as pd
import pysam
import subprocess
import time


# Create logger for bam submodule
bam_logger = logging.getLogger("xyalign.bam")


class BamFile():
	"""
	A class for working with external bam files

	Attributes
	----------

	filepath : str
		Full path to external bam file.
	samtools : str
		Full path to samtools. Default = 'samtools'

	"""
	def __init__(self, filepath, samtools="samtools", no_initial_index=False):
		self.filepath = filepath
		self.samtools = samtools
		self.logger = logging.getLogger("xyalign.bam.BamFile")
		self.logger.info("Creating a BamFile instance for {}".format(
			self.filepath))
		if no_initial_index is False:
			if self.is_indexed() is False:
				self.index_bam()

	def is_indexed(self):
		"""
		Checks that bam index exists, is not empty, and is newer than bam.

		Returns
		-------

		bool
			True if bam index exists and is newer than bam, False otherwise.

		"""
		self.logger.info("Checking indexing of {}".format(self.filepath))
		if os.path.exists("{}.bai".format(self.filepath)):
			if os.stat("{}.bai".format(self.filepath)).st_size != 0:
				idx_stamp = os.path.getmtime("{}.bai".format(self.filepath))
			else:
				self.logger.info("Bam index empty")
				return False
		elif os.path.exists("{}.bai".format(self.filepath[:-4])):
			if os.stat("{}.bai".format(self.filepath[:-4])).st_size != 0:
				idx_stamp = os.path.getmtime("{}.bai".format(self.filepath[:-4]))
			else:
				self.logger.info("Bam index empty")
				return False
		else:
			self.logger.info("No bam index detected")
			return False

		bam_stamp = os.path.getmtime(self.filepath)
		if bam_stamp < idx_stamp:
			self.logger.info("Bam index present and newer than bam file")
			return True
		else:
			self.logger.info("Bam index is older than bam file")
			return False

	def index_bam(self):
		"""
		Indexes a bam using samtools ('samtools index file.bam').

		Returns
		-------

		bool
			True if successful.

		Raises
		------

		RuntimeError
			If return code from external call is not 0.

		"""
		self.logger.info("Indexing bam file: {}".format(self.filepath))
		idx_start = time.time()
		rc = subprocess.call([self.samtools, "index", self.filepath])
		if rc == 0:
			self.logger.info("Indexing complete. Elapsed time: {} seconds".format(
				time.time() - idx_start))
			return True
		else:
			self.logger.error("Unable to index bamfile {}. Exiting".format(
				self.filepath))
			logging.shutdown()
			raise RuntimeError("Unable to index bamfile. Exiting")

	def get_chrom_length(self, chrom):
		"""
		Extract chromosome length from BAM header.

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
		bamfile = pysam.AlignmentFile(self.filepath, "rb")
		lengths = dict(zip(bamfile.references, bamfile.lengths))
		try:
			lens = lengths[chrom]
			bamfile.close()
			return lens
		except:
			self.logger.error(
				"{} not present in bam header for {}. Exiting.".format(
					chrom, self.filepath))
			logging.shutdown()
			raise RuntimeError(
				"Chromosome name not present in bam header. Exiting")

	def chromosome_lengths(self):
		"""
		Returns
		-------

		tuple
			chromosome lengths ordered by sequence order in bam header

		"""
		bamfile = pysam.AlignmentFile(self.filepath, "rb")
		lengths = bamfile.lengths
		bamfile.close()
		return lengths

	def chromosome_names(self):
		"""
		Returns
		-------

		tuple
			chromosome names ordered by sequence order in bam header

		"""
		bamfile = pysam.AlignmentFile(self.filepath, "rb")
		names = bamfile.references
		bamfile.close()
		return names

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
			If chromosome name is not in bam header.

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
						"Error finding chromosome length in bam file {} "
						"(for bed file)".format(self.filepath))
					logging.shutdown()
					raise RuntimeError(
						"Error finding chromosome length in bam file {}.  Check "
						"chromosome names and bam header.".format(
							self.filepath))
		self.logger.info(
			"Bed file ({}) created. Elapsed time: {} seconds".format(
				output_file, time.time() - c_bed_start))
		return output_file

	def check_chrom_in_bam(self, chromosome_list):
		"""
		Checks to see if all chromosomes in chromosome_list are in bam file

		Parameters
		----------

		chromosome_list : list
				Chromosomes/scaffolds to check

		Returns
		-------

		list
			List of chromosomes not in bam file
		"""
		self.logger.info(
			"Checking to ensure all chromosomes are found in {}".format(self.filepath))
		bam_chroms = self.chromosome_names()

		missing_list = []
		for c in chromosome_list:
			if c not in bam_chroms:
				missing_list.append(c)

		if len(missing_list) > 0:
			self.logger.info(
				"The following chromosomes were not found in {}: {}".format(
					self.filepath, ",".join(missing_list)))
		else:
			self.logger.info(
				"All chromosomes present in {}".format(self.filepath))
		return missing_list

	def sort_bam(
		self, sorted_bam, query_name=False):
		"""
		Sorts bam file by coordinate (query_name=False) or
		query name (query_name=True)

		Parameters
		----------

		sorted_bam : str
			Full path to (including desired name of) output bam file
		query_name : bool
			If True, sort by query name (read ID), else sort by coordinate

		Returns
		-------

		BamFile() object
			BamFile() object of new, sorted bam file
		"""
		if query_name is True:
			command_line = "{} sort -n -o {} {}".format(
				self.samtools, sorted_bam, self.filepath)
			self.logger.info(
				"Sorting {} by query name with command: {}".format(
					self.filepath, command_line))
		else:
			command_line = "{} sort -o {} {}".format(
				self.samtools, sorted_bam, self.filepath)
			self.logger.info(
				"Sorting {} by coordinate with command: {}".format(
					self.filepath, command_line))
		subprocess.call(command_line, shell=True)
		sorted_bam_obj = BamFile(sorted_bam, self.samtools, no_initial_index=True)
		return sorted_bam_obj

	def extract_regions(
		self, regions, new_bam, rg_id=None):
		"""
		Extracts regions from a bam file into new bam file.

		Parameters
		----------

		regions : list
			regions from which reads will be stripped
		new_bam : str
			Full path to and desired name of output bam file
		rg_id : str or None
			Path to text file containing read group ids to use when isolating regions.
			If None, all reads from regions will be extracted.

		Returns
		-------

		BamFile() object
			BamFile() object of new bam file (containing extracted regions)
		"""
		if rg_id is None:
			command_line = "{} view -bh -o {} {} {}".format(
				self.samtools, new_bam, self.filepath, ' '.join(map(str, regions)))
			self.logger.info(
				"Extracting regions from {} into {} using the command: {}".format(
					self.filepath, new_bam, command_line))
		else:
			command_line = "{} view -bh -R {} -o {} {} {}".format(
				self.samtools, rg_id, new_bam, self.filepath, ' '.join(map(str, regions)))
			self.logger.info(
				"Extracting reads with RG ID {} from regions ({}) into {} "
				"using the command: {}".format(
					self.filepath, rg_id, new_bam, command_line))
		subprocess.call(command_line, shell=True)
		new_bam_obj = BamFile(new_bam, self.samtools, no_initial_index=True)
		return new_bam_obj

	def extract_read_group(
		self, new_bam, rg_id):
		"""
		Extracts all reads belonging to a given RG ID from a
		bam file into new bam file.

		Parameters
		----------

		new_bam : str
			Full path to and desired name of output bam file
		rg_id : str
			Path to text file containing read group ids to use when isolating regions.

		Returns
		-------

		BamFile() object
			BamFile() object of new bam file (containing extracted regions)
		"""
		command_line = "{} view -bh -R {} -o {} {}".format(
			self.samtools, rg_id, new_bam, self.filepath)
		self.logger.info(
			"Extracting reads with RG ID {} from regions ({}) into {} "
			"using the command: {}".format(
				self.filepath, rg_id, new_bam, command_line))
		subprocess.call(command_line, shell=True)
		new_bam_obj = BamFile(new_bam, self.samtools)
		return new_bam_obj

	def strip_reads(
		self, repairsh, shufflesh, single, output_directory,
		output_prefix, regions, repair_xmx, compression,
		cleanup=True, default_rg="None"):
		"""
		Strips reads from a bam or cram file in provided regions and outputs
		sorted fastqs containing reads, one set of fastq files per read group.

		Parameters
		----------

		repairsh : str
			Path to repair.sh (from BBmap)
		shufflesh : str
			Path to shuffle.sh (from BBmap)
		single : bool
			If true output single-end fastq, otherwise output paired-end fastqs
		output_directory : str
			The directory for ALL output (including temporary files)
		output_prefix : str
			The name (without path) to use for prefix to output fastqs
		regions : list
			regions from which reads will be stripped
		repair_xmx : str
			If "None", repair.sh will allocate its own memory. Otherwise value
			will be provided in the form of -Xmx4g, where 4g is the value provided
			as repair_xmx
		compression : int
			Desired compression level (0-9) for output fastqs. If 0, fastqs
			will be uncompressed.
		cleanup : bool
			If true, will clean up temporary files.
		default_rg : str
			If "None", no default read group will be created. Otherwise, default
			read group will be string provided.  This read group will consist
			exclusively of an ID.

		Returns
		-------

		list
			A two-item list containing the path to a text file pairing read group
			names with associated output fastqs, and a text file containing a
			list of @RG lines associated with each read group

		"""
		# Collect RGs
		rg_start = time.time()
		rg_list = os.path.join(
			output_directory, output_prefix + ".full_rg.list")
		command_line = """{} view -H {} | awk '$1=="\x40RG"' | """\
			"""awk {} """\
			"""| cut -d':' -f 2 > {}""".format(
				self.samtools, self.filepath,
				repr('{for(i=1;i<=NF;i++){if (substr($i,1,2) ~ /ID/){print $i}}}'),
				rg_list)
		self.logger.info("Grabbing read groups from {} with the command: {}".format(
			self.filepath, command_line))
		subprocess.call(command_line, shell=True)
		rg_header_lines = os.path.join(
			output_directory, output_prefix + ".header_lines_rg.list")
		command_line = """{} view -H {} | awk '$1=="\x40RG"' > {}""".format(
			self.samtools, self.filepath, rg_header_lines)
		self.logger.info(
			"Grabbing RG header lines from {} with the command: {}".format(
				self.filepath, command_line))
		subprocess.call(command_line, shell=True)
		temporary_fastqs = []
		temporary_rg_files = []
		temporary_bams = []
		if os.stat(rg_list).st_size == 0:
			if default_rg != "None":
				rg = default_rg
				with open(rg_header_lines, "w") as rg_head:
					rg_head.write("@RG\tID:{}".format(default_rg))
			else:
				rg = "None"
				rg_header_lines = None
			out_rg_table = os.path.join(
				output_directory, output_prefix + ".rg_fastq_key.list")
			with open(out_rg_table, "w") as ortab:
				self.logger.info(
					"No read group information found in {}.  Will therefore treat "
					"all reads as coming from the same read group and write all "
					"reads to same fastqs.".format(self.filepath))
				if single is False:
					if rg == "None":
						temp1 = os.path.join(
							output_directory, "{}.temp_1.fastq".format(output_prefix))
						temp2 = os.path.join(
							output_directory, "{}.temp_2.fastq".format(output_prefix))
						out1 = os.path.join(
							output_directory, "{}_1.fastq".format(output_prefix))
						out2 = os.path.join(
							output_directory, "{}_2.fastq".format(output_prefix))
						iso_bam = os.path.join(
							output_directory, "{}_extracted.bam".format(output_prefix))
						iso_bam_sort = os.path.join(
							output_directory, "{}_extracted.sorted.bam".format(output_prefix))
					else:
						temp1 = os.path.join(
							output_directory, "{}.{}.temp_1.fastq".format(output_prefix, rg))
						temp2 = os.path.join(
							output_directory, "{}.{}.temp_2.fastq".format(output_prefix, rg))
						out1 = os.path.join(
							output_directory, "{}_{}_1.fastq".format(output_prefix, rg))
						out2 = os.path.join(
							output_directory, "{}_{}_2.fastq".format(output_prefix, rg))
						iso_bam = os.path.join(
							output_directory, "{}_{}_extracted.bam".format(output_prefix, rg))
						iso_bam_sort = os.path.join(
							output_directory, "{}_{}_extracted.sorted.bam".format(
								output_prefix, rg))
					if compression != 0:
						temp1 += ".gz"
						temp2 += ".gz"
						out1 += ".gz"
						out2 += ".gz"
					if len(regions) < len(self.chromosome_names()):
						# Isolate chromosomes
						iso_bam_obj = self.extract_regions(regions, iso_bam)
						# Sort isolated bam
						iso_bam_sort_obj = iso_bam_obj.sort_bam(iso_bam_sort, query_name=True)
						temporary_bams.extend([
							iso_bam, iso_bam + ".bai",
							iso_bam_sort, iso_bam_sort + ".bai"])
						# Strip reads
						if compression == 0:
							command_line = "{} fastq -n -t -1 {} -2 {} {}".format(
								self.samtools, temp1, temp2, iso_bam_sort_obj.filepath)
						else:
							command_line = "{} fastq -n -t -c {} -1 {} -2 {} {}".format(
								self.samtools, compression, temp1, temp2, iso_bam_sort_obj.filepath)
						self.logger.info(
							"Stripping reads with command: {}".format(
								command_line))
						subprocess.call(command_line, shell=True)
						temporary_fastqs.extend([temp1, temp2])
					else:
						# Strip all reads
						self.logger.info(
							"Length of region list same or greater than number of "
							"chroms in bam, so stripping all reads")
						if compression == 0:
							command_line = "{} fastq -n -t -1 {} -2 {} {}".format(
								self.samtools, temp1, temp2, self.filepath)
						else:
							command_line = "{} fastq -n -t -c {} -1 {} -2 {} {}".format(
								self.samtools, compression, temp1, temp2, self.filepath)
						self.logger.info(
							"Stripping reads with command: {}".format(
								command_line))
						subprocess.call(command_line, shell=True)
						temporary_fastqs.extend([temp1, temp2])
					if repair_xmx == "None":
						if compression == 0:
							command_line = "{} in1={} in2={} out1={} out2={} " \
								"overwrite=true".format(
									repairsh, temp1, temp2, out1, out2)
						else:
							command_line = "{} in1={} in2={} out1={} out2={} ziplevel={} " \
								"overwrite=true".format(
									repairsh, temp1, temp2, out1, out2, compression)
					else:
						if compression == 0:
							command_line = "{} -Xmx{} in1={} in2={} out1={} out2={} " \
								"overwrite=true".format(
									repairsh, repair_xmx, temp1, temp2, out1, out2)
						else:
							command_line = "{} -Xmx{} in1={} in2={} out1={} out2={} ziplevel={} " \
								"overwrite=true".format(
									repairsh, repair_xmx, temp1, temp2, out1, out2, compression)
					self.logger.info(
						"Sorting reads with command: {}".format(
							command_line))
					subprocess.call(command_line, shell=True)
					ortab.write("{}\t{}\t{}\n".format(
						rg, out1, out2))
				else:
					if rg == "None":
						temp1 = os.path.join(
							output_directory, "{}.temp.fastq".format(output_prefix))
						out1 = os.path.join(
							output_directory, "{}.fastq".format(output_prefix))
						iso_bam = os.path.join(
							output_directory, "{}_extracted.bam".format(output_prefix))
						iso_bam_sort = os.path.join(
							output_directory, "{}_extracted.sorted.bam".format(output_prefix))
					else:
						temp1 = os.path.join(
							output_directory, "{}.{}.temp.fastq".format(output_prefix, rg))
						out1 = os.path.join(
							output_directory, "{}_{}.fastq".format(output_prefix, rg))
						iso_bam = os.path.join(
							output_directory, "{}_{}_extracted.bam".format(output_prefix, rg))
						iso_bam_sort = os.path.join(
							output_directory, "{}_{}_extracted.sorted.bam".format(
								output_prefix, rg))
					if compression != 0:
						temp1 += ".gz"
						out1 += ".gz"
					if len(regions) < len(self.chromosome_names()):
						# Isolate chromosomes
						iso_bam_obj = self.extract_regions(regions, iso_bam)
						# Sort isolated bam
						iso_bam_sort_obj = iso_bam_obj.sort_bam(iso_bam_sort, query_name=True)
						temporary_bams.extend([
							iso_bam, iso_bam + ".bai",
							iso_bam_sort, iso_bam_sort + ".bai"])
						# Strip reads
						if compression == 0:
							command_line = "{} fastq -n -t -0 {} {}".format(
								self.samtools, temp1, iso_bam_sort_obj.filepath)
						else:
							command_line = "{} fastq -n -t -c {} -0 {} {}".format(
								self.samtools, compression, temp1, iso_bam_sort_obj.filepath)
						self.logger.info(
							"Stripping reads with command: {}".format(
								command_line))
						subprocess.call(command_line, shell=True)
						temporary_fastqs.append(temp1)
					else:
						# Strip all reads
						self.logger.info(
							"Length of region list same or greater than number of "
							"chroms in bam, so stripping all reads")
						if compression == 0:
							command_line = "{} fastq -n -t -0 {} {} > {}".format(
								self.samtools, temp1, self.filepath)
						else:
							command_line = "{} fastq -n -t -c {} -0 {} {} > {}".format(
								self.samtools, compression, temp1, self.filepath)
						self.logger.info(
							"Stripping reads with command: {}".format(
								command_line))
						subprocess.call(command_line, shell=True)
						temporary_fastqs.append(temp1)
					if repair_xmx == "None":
						if compression == 0:
							command_line = "{} in={} out={} name overwrite=true".format(
								shufflesh, temp1, out1)
						else:
							command_line = "{} in={} out={} name ziplevel={} overwrite=true".format(
								shufflesh, temp1, out1, compression)
					else:
						if compression == 0:
							command_line = "{} -Xmx{} in={} out={} name overwrite=true".format(
								shufflesh, repair_xmx, temp1, out1)
						else:
							command_line = "{} -Xmx{} in={} out={} name ziplevel={} " \
								"overwrite=true".format(
									shufflesh, repair_xmx, temp1, out1, compression)
					self.logger.info(
						"Sorting reads with command: {}".format(
							command_line))
					subprocess.call(command_line, shell=True)
					ortab.write("{}\t{}\n".format(
						rg, out1))
		else:
			with open(rg_list, "r") as f:
				out_rg_table = os.path.join(
					output_directory, output_prefix + ".rg_fastq_key.list")
				self.logger.info(
					"Iteratively removing reads by read group. Writing table of RGs and "
					"fastqs to {}".format(out_rg_table))
				with open(out_rg_table, "w") as ortab:
					for line in f:
						rg = line.strip()
						if rg != "":
							self.logger.info("Removing reads from group: {}".format(rg))
							tmp_out = "{}/{}.txt".format(output_directory, rg)
							temporary_rg_files.append(tmp_out)
							with open(tmp_out, "w") as o:
								o.write(rg)
							if single is False:
								temp1 = os.path.join(
									output_directory, "{}.{}.temp_1.fastq".format(output_prefix, rg))
								temp2 = os.path.join(
									output_directory, "{}.{}.temp_2.fastq".format(output_prefix, rg))
								out1 = os.path.join(
									output_directory, "{}_{}_1.fastq".format(output_prefix, rg))
								out2 = os.path.join(
									output_directory, "{}_{}_2.fastq".format(output_prefix, rg))
								iso_bam = os.path.join(
									output_directory, "{}_{}_extracted.bam".format(output_prefix, rg))
								iso_bam_sort = os.path.join(
									output_directory, "{}_{}_extracted.sorted.bam".format(
										output_prefix, rg))
								if compression != 0:
									temp1 += ".gz"
									temp2 += ".gz"
									out1 += ".gz"
									out2 += ".gz"
								if len(regions) < len(self.chromosome_names()):
									# Isolate chromosomes
									iso_bam_obj = self.extract_regions(regions, iso_bam, tmp_out)
									# Sort isolated bam
									iso_bam_sort_obj = iso_bam_obj.sort_bam(iso_bam_sort, query_name=True)
									temporary_bams.extend([
										iso_bam, iso_bam + ".bai",
										iso_bam_sort, iso_bam_sort + ".bai"])
									# Strip reads
									if compression == 0:
										command_line = "{} fastq -n -t -1 {} -2 {} {}".format(
											self.samtools, temp1, temp2, iso_bam_sort_obj.filepath)
									else:
										command_line = "{} fastq -n -t -c {} -1 {} -2 {} {}".format(
											self.samtools, compression, temp1, temp2, iso_bam_sort_obj.filepath)
									self.logger.info(
										"Stripping reads with command: {}".format(
											command_line))
									subprocess.call(command_line, shell=True)
									temporary_fastqs.extend([temp1, temp2])
								else:
									# Isolate RG ID
									iso_bam_obj = self.extract_read_group(iso_bam, tmp_out)
									# Sort isolated bam
									iso_bam_sort_obj = iso_bam_obj.sort_bam(iso_bam_sort, query_name=True)
									temporary_bams.extend([
										iso_bam, iso_bam + ".bai",
										iso_bam_sort, iso_bam_sort + ".bai"])
									# Strip all reads
									self.logger.info(
										"Length of region list same or greater than number of "
										"chroms in bam, so stripping all reads")
									if compression == 0:
										command_line = "{} fastq -n -t -1 {} -2 {} {}".format(
											self.samtools, temp1, temp2, iso_bam_sort_obj.filepath)
									else:
										command_line = "{} fastq -n -t -c {} -1 {} -2 {} {}".format(
											self.samtools, compression, temp1, temp2, iso_bam_sort_obj.filepath)
									self.logger.info(
										"Stripping reads with command: {}".format(
											command_line))
									subprocess.call(command_line, shell=True)
									temporary_fastqs.extend([temp1, temp2])
								if repair_xmx == "None":
									if compression == 0:
										command_line = "{} in1={} in2={} out1={} out2={} " \
											"overwrite=true".format(
												repairsh, temp1, temp2, out1, out2)
									else:
										command_line = "{} in1={} in2={} out1={} out2={} ziplevel={} " \
											"overwrite=true".format(
												repairsh, temp1, temp2, out1, out2, compression)
								else:
									if compression == 0:
										command_line = "{} -Xmx{} in1={} in2={} out1={} out2={} " \
											"overwrite=true".format(
												repairsh, repair_xmx, temp1, temp2, out1, out2)
									else:
										command_line = "{} -Xmx{} in1={} in2={} out1={} out2={} " \
											"ziplevel={} overwrite=true".format(
												repairsh, repair_xmx, temp1, temp2, out1, out2, compression)
								self.logger.info(
									"Sorting reads with command: {}".format(
										command_line))
								subprocess.call(command_line, shell=True)
								ortab.write("{}\t{}\t{}\n".format(
									rg, out1, out2))
							else:
								temp1 = os.path.join(
									output_directory, "{}.{}.temp.fastq".format(output_prefix, rg))
								out1 = os.path.join(
									output_directory, "{}_{}.fastq".format(output_prefix, rg))
								iso_bam = os.path.join(
									output_directory, "{}_{}_extracted.bam".format(output_prefix, rg))
								iso_bam_sort = os.path.join(
									output_directory, "{}_{}_extracted.sorted.bam".format(
										output_prefix, rg))
								if compression != 0:
									temp1 += ".gz"
									out1 += ".gz"
								if len(regions) < len(self.chromosome_names()):
									# Isolate chromosomes
									iso_bam_obj = self.extract_regions(regions, iso_bam, tmp_out)
									# Sort isolated bam
									iso_bam_sort_obj = iso_bam_obj.sort_bam(iso_bam_sort, query_name=True)
									temporary_bams.extend([
										iso_bam, iso_bam + ".bai",
										iso_bam_sort, iso_bam_sort + ".bai"])
									# Strip reads
									if compression == 0:
										command_line = "{} fastq -n -t -0 {}".format(
											self.samtools, temp1, temp1, iso_bam_sort_obj.filepath)
									else:
										command_line = "{} fastq -n -t -c {} -0 {} {}".format(
											self.samtools, compression, temp1, iso_bam_sort_obj.filepath)
									self.logger.info(
										"Stripping reads with command: {}".format(
											command_line))
									subprocess.call(command_line, shell=True)
									temporary_fastqs.append(temp1)
								else:
									# Isolate RG ID
									iso_bam_obj = self.extract_read_group(iso_bam, tmp_out)
									# Sort isolated bam
									iso_bam_sort_obj = iso_bam_obj.sort_bam(iso_bam_sort, query_name=True)
									temporary_bams.extend([
										iso_bam, iso_bam + ".bai",
										iso_bam_sort, iso_bam_sort + ".bai"])
									# Strip all reads
									self.logger.info(
										"Length of region list same or greater than number of "
										"chroms in bam, so stripping all reads")
									if compression == 0:
										command_line = "{} fastq -n -t -0 {} {}".format(
											self.samtools, temp1, self.filepath)
									else:
										command_line = "{} fastq -n -t -c {} -0 {} {} > {}".format(
											self.samtools, compression, temp1, self.filepath)
									self.logger.info(
										"Stripping reads with command: {}".format(
											command_line))
									subprocess.call(command_line, shell=True)
									temporary_fastqs.append(temp1)
								if repair_xmx == "None":
									if compression == 0:
										command_line = "{} in={} out={} name overwrite=true".format(
											shufflesh, temp1, out1)
									else:
										command_line = "{} in={} out={} name ziplevel={} " \
											"overwrite=true".format(
												shufflesh, temp1, out1, compression)
								else:
									if compression == 0:
										command_line = "{} -Xmx{} in={} out={} name overwrite=true".format(
											shufflesh, repair_xmx, temp1, out1)
									else:
										command_line = "{} -Xmx{} in={} out={} name ziplevel={} " \
											"overwrite=true".format(
												shufflesh, repair_xmx, temp1, out1, compression)
								self.logger.info(
									"Sorting reads with command: {}".format(
										command_line))
								subprocess.call(command_line, shell=True)
								# write line
								ortab.write("{}\t{}\n".format(
									rg, out1))
		files_to_remove = temporary_fastqs + temporary_rg_files + temporary_bams
		if cleanup is True:
			self.logger.info("Cleanup is True. Removing temporary files: {}".format(
				",".join(files_to_remove)))
			for file_path in files_to_remove:
				if os.path.exists(file_path):
					os.remove(file_path)
				else:
					self.logger.error("Unable to remove file: {}".format(file_path))
		else:
			self.logger.info(
				"Cleanup is False. The following temporary files will be left "
				"intact: {}".format(", ".join(files_to_remove)))
		self.logger.info(
			"Stripping and sorting reads complete. Elapsed time: {} seconds".format(
				time.time() - rg_start))
		return [out_rg_table, rg_header_lines]

	def analyze_bam(
		self, chrom, duplicates, exact, window_size, target_file=None):
		"""
		Analyze BAM (or CRAM) file for depth and mapping quality across genomic
		windows.

		Parameters
		----------

		chrom : str
			The name of the chromosome to analyze
		duplicates : bool
			If True, duplicates included in analyses.
		exact : bool
			If True, mean depth is calculated exactly within each window.
			If False, an accurate (and much faster) approximation is used
		window_size
			If int, the window size to use for sliding window analyses, if None
			intervals from target_file
		target_file : str
			Path to bed_file containing regions to analyze instead of
			windows of a fixed size. Will only be engaged if window_size is None

		Returns
		-------

		pandas dataframe
			pandas data frame with the columns: "chrom", "start", "stop",
			"depth", "mapq"

		"""
		self.logger.info(
			"Traversing {} in {} to analyze depth and mapping quality".format(
				chrom, self.filepath))

		# exact not yet implemented
		if exact is True:
			self.logger.error(
				"Exact is not yet implemented. Exiting")
			logging.shutdown()
			raise RuntimeError(
				"Exact is not yet implemented. Exiting")

		analyze_start = time.time()
		samfile = pysam.AlignmentFile(self.filepath, "rb")
		chr_len = self.get_chrom_length(chrom)
		no_query_len = 0

		if window_size is not None:
			self.logger.info(
				"Using windows size: {}".format(window_size))

			num_windows = chr_len // window_size + 1
			if chr_len % num_windows == 0:
				last_window_len = window_size
			else:
				last_window_len = chr_len % num_windows

			window_id = 0

			chr_list = [chrom] * num_windows
			start_list = []
			stop_list = []
			depth_list = []
			mapq_list = []

			start = 0
			end = window_size
			if duplicates is True:
				for window in range(0, num_windows):
					mapq = []
					total_read_length = 0
					for read in samfile.fetch(chrom, start, end):
						if read.is_secondary is False:
							if read.is_supplementary is False:
								try:
									total_read_length += read.infer_query_length()
								except TypeError:
									no_query_len += 1
									continue
								mapq.append(read.mapping_quality)
					start_list.append(start)
					stop_list.append(end)
					depth_list.append(total_read_length / window_size)
					mapq_list.append(np.mean(np.asarray(mapq)))

					window_id += 1
					if window_id == num_windows - 1:
						start += window_size
						end += last_window_len
					else:
						start += window_size
						end += window_size

					# Print progress
					if window_id % 1000 == 0:
						self.logger.info("{} out of {} windows processed on {}".format(
							window_id, num_windows, chrom))

			else:
				for window in range(0, num_windows):
					mapq = []
					total_read_length = 0
					for read in samfile.fetch(chrom, start, end):
						if read.is_secondary is False:
							if read.is_supplementary is False:
								if read.is_duplicate is False:
									try:
										total_read_length += read.infer_query_length()
									except TypeError:
										no_query_len += 1
										continue
									mapq.append(read.mapping_quality)
					start_list.append(start)
					stop_list.append(end)
					depth_list.append(total_read_length / window_size)
					mapq_list.append(np.mean(np.asarray(mapq)))

					window_id += 1
					if window_id == num_windows - 1:
						start += window_size
						end += last_window_len
					else:
						start += window_size
						end += window_size

					# Print progress
					if window_id % 1000 == 0:
						self.logger.info("{} out of {} windows processed on {}".format(
							window_id, num_windows, chrom))

		elif target_file is not None:
			self.logger.info(
				"Using targets from: {}".format(target_file))
			with open(target_file) as f:
				targets = [x.strip() for x in f]
				targets = [x.split() for x in targets]
				targets = [x for x in targets if x[0] == chrom]
				while [""] in targets:
					targets.remove([""])

			num_windows = len(targets)

			window_id = 0

			chr_list = []
			start_list = []
			stop_list = []
			depth_list = []
			mapq_list = []

			if duplicates is True:
				for window in targets:
					mapq = []
					total_read_length = 0
					start = int(window[1])
					end = int(window[2])
					window_size = end - start
					for read in samfile.fetch(window[0], start, end):
						if read.is_secondary is False:
							if read.is_supplementary is False:
								try:
									total_read_length += read.infer_query_length()
								except TypeError:
									no_query_len += 1
									continue
								mapq.append(read.mapping_quality)
					chr_list.append(window[0])
					start_list.append(start)
					stop_list.append(end)
					depth_list.append(total_read_length / window_size)
					mapq_list.append(np.mean(np.asarray(mapq)))

					window_id += 1

					# Print progress
					if window_id % 1000 == 0:
						self.logger.info("{} out of {} targets processed on {}".format(
							window_id, num_windows, chrom))
			else:
				for window in targets:
					mapq = []
					total_read_length = 0
					start = int(window[1])
					end = int(window[2])
					window_size = end - start
					for read in samfile.fetch(window[0], start, end):
						if read.is_secondary is False:
							if read.is_supplementary is False:
								if read.is_duplicate is False:
									try:
										total_read_length += read.infer_query_length()
									except TypeError:
										no_query_len += 1
										continue
									mapq.append(read.mapping_quality)
					chr_list.append(window[0])
					start_list.append(start)
					stop_list.append(end)
					depth_list.append(total_read_length / window_size)
					mapq_list.append(np.mean(np.asarray(mapq)))

					window_id += 1

					# Print progress
					if window_id % 1000 == 0:
						self.logger.info("{} out of {} targets processed on {}".format(
							window_id, num_windows, chrom))

		else:
			self.logger.error(
				"Both window_size and target_file set to None. "
				"Cannot proceed with bam traversal. Exiting.")
			logging.shutdown()
			raise RuntimeError(
				"Both window_size and target_file set to None. "
				"Cannot proceed with bam traversal. Exiting.")

		# Convert data into pandas data frames
		windows_df = pd.DataFrame({
			"chrom": np.asarray(chr_list),
			"start": np.asarray(start_list),
			"stop": np.asarray(stop_list),
			"depth": np.asarray(depth_list),
			"mapq": np.asarray(mapq_list)
		})[["chrom", "start", "stop", "depth", "mapq"]]

		samfile.close()
		self.logger.info(
			"Analysis complete. {} reads ignored because infer_query_length() "
			"returned None. Elapsed time: {} seconds".format(
				no_query_len, time.time() - analyze_start))
		return windows_df

	def chrom_stats(self, chrom, duplicates):
		"""
		Analyze BAM (or CRAM) file for depth and mapping quality across a
		single chromosome.

		Parameters
		----------

		chrom : str
			The name of the chromosome to analyze
		duplicates : bool
			If True, duplicates included in analyses.

		Returns
		-------

		tuple
			(mean_depth, mean_mapq)

		"""
		self.logger.info(
			"Traversing {} (whole chromosome) in {} to analyze depth and "
			"mapping quality".format(
				chrom, self.filepath))

		analyze_start = time.time()
		samfile = pysam.AlignmentFile(self.filepath, "rb")
		chr_len = self.get_chrom_length(chrom)

		no_query_len = 0

		start = 0
		end = chr_len
		if duplicates is True:
			mapq = []
			total_read_length = 0
			for read in samfile.fetch(chrom, start, end):
				if read.is_secondary is False:
					if read.is_supplementary is False:
						try:
							total_read_length += read.infer_query_length()
						except TypeError:
							no_query_len += 1
							continue
						mapq.append(read.mapping_quality)

		else:
			mapq = []
			total_read_length = 0
			for read in samfile.fetch(chrom, start, end):
				if read.is_secondary is False:
					if read.is_supplementary is False:
						if read.is_duplicate is False:
							try:
								total_read_length += read.infer_query_length()
							except TypeError:
								no_query_len += 1
								continue
							mapq.append(read.mapping_quality)

		samfile.close()
		self.logger.info(
			"Analysis complete. {} reads ignored because infer_query_length() "
			"returned None. Elapsed time: {} seconds".format(
				no_query_len, time.time() - analyze_start))
		return ((float(total_read_length) / chr_len), np.mean(np.asarray(mapq)))

	def chrom_counts(self):
		"""
		Get read counts per chrom from a bamfile
		"""
		self.logger.info(
			"Using index of {} to get read counts per chromosome".format(
				self.filepath))

		analyze_start = time.time()
		samfile = pysam.AlignmentFile(self.filepath, "rb")
		idx_out = samfile.get_index_statistics()
		samfile.close()
		self.logger.info("Index analysis complete. Elapsed time: {} seconds".format(
			time.time() - analyze_start))
		return idx_out

	def platypus_caller(
		self, platypus_path, log_path, ref, chroms, cpus, output_file,
		regions_file=None):
		"""
		Uses platypus to make variant calls on provided bam file

		Parameters
		----------

		platypus_path : str
			Path to platypus
		log_path : str
			Path to and name of desired log file for platypus
		ref : str
			Path to reference sequence
		chroms : list
			Chromosomes to call variants on, e.g., ["chrX", "chrY", "chr19"]
		cpus : int
			Number of threads/cores to use
		output_file : path
			Path to and name of the output vcf
		regions_file : {str, None}
			If not None, must be path to bed file containing regions to call variants
			in.  If None, calls in call regions of provided chromosomes. Default =
			None.

		Returns
		-------

		int
			Exit code of the platypus call

		"""
		platy_start = time.time()
		if regions_file is None:
			regions = ','.join(map(str, chroms))
		else:
			regions = regions_file
		command_line = [
			platypus_path, "callVariants", "--bamFiles", self.filepath, "-o",
			output_file, "--refFile", ref, "--nCPU", str(cpus), "--regions", regions,
			"--assemble", "1", "--logFileName", log_path]
		self.logger.info("Calling variants with command line: {}".format(
			" ".join(command_line)))
		return_code = subprocess.call(command_line)
		self.logger.info(
			"Variant calling complete. Elapsed time: {} seconds".format(
				time.time() - platy_start))
		return return_code

###############################################################################
# Functions associated with bam files


def switch_sex_chromosomes_sambamba(
	samtools_path, sambamba_path, bam_orig, bam_new, sex_chroms,
	output_directory, output_prefix, threads, pg_header_dict, cram=False):
	"""
	Removes sex chromosomes from original bam and merges in remmapped
	sex chromosomes, while retaining the original bam header (and adding new
	@PG line)

	Parameters
	----------

	samtools_path : str
		The path to samtools
	sambamba_path :
		The path to sambamba
	bam_orig : str
		The path to the original full bam file
	bam_new : str
		The path to the bam file containing the remapped sex chromosomes
	sex_chroms : list
		Sex chromosomes (to be removed from bam_orig)
	output_directory : str
		The path to directory where all files (inc. temp) will be output
	output_prefix : str
		The name (without path) to use as prefix for all files
	threads : int
		The number of threads/cpus to use
	pg_header_dict : dict
		dictionary with information to be included in the new @PG line
			- must contain:
				Key = 'CL', value = list of command line values
				Key = 'ID', value = string of program ID
			- optional:
				Key = 'VN', value = string of program version
	cram : bool
		If True, will treat input as cram files and output cram files.
		Otherwise, will treat input as bam.  Defaule is False. True is currently
		unsupported.

	Returns
	-------

	str
		Bam or cram file path with new sex chromosomes, but all others intact.

	Raises
	------

	RuntimeError
		If cram is not False.

	"""
	switch_start = time.time()
	bam_logger.info(
		"Beginning process of swapping sex chromosomes into full assembly")
	# Grab original header
	command_line = [samtools_path, "view", "-H", bam_orig]
	bam_logger.info(
		"Isolating header with command: {} > {}/{}.header.sam".format(
			" ".join(command_line), output_directory, output_prefix))
	with open(
		"{}/{}.header.sam".format(output_directory, output_prefix), "w") as f:
		subprocess.call(command_line, stdout=f)
	# Reheader new bam (for merge)
	command_line = [samtools_path, "reheader", "-P", "{}/{}.header.sam".format(
		output_directory, output_prefix), bam_new]
	bam_logger.info(
		"Reheading {} with command: {} > {}/{}.reheadered.temp.new.bam".format(
			bam_new, " ".join(command_line), output_directory, output_prefix))
	with open(
		"{}/{}.reheadered.temp.new.bam".format(
			output_directory, output_prefix), "w") as f:
			subprocess.call(command_line, stdout=f)
	bam_logger.info("Indexing {}".format("{}/{}.reheadered.temp.new.bam".format(
		output_directory, output_prefix)))
	subprocess.call(
		[samtools_path, "index", "{}/{}.reheadered.temp.new.bam".format(
			output_directory, output_prefix)])
	# Add XYalign @PG line to header
	bam_logger.info("Adding new PG line to header")
	cl_string = " ".join(pg_header_dict["CL"])
	if "VN" in pg_header_dict:
		pg_line = [
			"@PG", "ID:{}".format(pg_header_dict["ID"]), "VN:{}".format(
				pg_header_dict["VN"]), "CL:{}".format(cl_string)]
	else:
		pg_line = [
			"@PG", "ID:{}".format(pg_header_dict["ID"]), "CL:{}".format(cl_string)]
	subprocess.call("echo '{}' >> {}/{}.header.sam".format(
		"\t".join(pg_line), output_directory, output_prefix), shell=True)
	if cram is False:
		# Remove sex chromosomes from original bam and merge
		samfile = pysam.AlignmentFile(bam_orig, "rb")
		non_sex_scaffolds = filter(
			lambda x: x not in sex_chroms, list(samfile.references))
		command_line = "{} view -h -t {} -f bam -o " \
			"{}/{}.temp.nosexchr.bam {} {}".format(
				sambamba_path, threads, output_directory, output_prefix, bam_orig,
				" ".join(non_sex_scaffolds))
		bam_logger.info(
			"Removing sex chromosomes with command: {}".format(command_line))
		subprocess.call(command_line, shell=True)
		# subprocess.call(
		# 	[sambamba_path, "view", "-h", "-t", "{}".format(threads), "-f",
		# 		"bam", "-o", "{}/temp.nosexchr.bam".format(output_directory),
		# 		bam_orig, "{}".format(" ".join(non_sex_scaffolds))])
		bam_logger.info("Indexing {}".format("{}/{}.temp.nosexchr.bam".format(
			output_directory, output_prefix)))
		subprocess.call(
			[samtools_path, "index", "{}/{}.temp.nosexchr.bam".format(
				output_directory, output_prefix)])
		command_line = [
			samtools_path, "merge", "-cp", "-@", "{}".format(threads), "-h",
			"{}/{}.header.sam".format(output_directory, output_prefix), "-f",
			"{}/{}.merged.bam".format(output_directory, output_prefix),
			"{}/{}.temp.nosexchr.bam".format(
				output_directory, output_prefix), bam_new]
		bam_logger.info("Merging bam files with command: {}".format(
			" ".join(command_line)))
		subprocess.call(command_line)
		bam_logger.info("Indexing {}".format(
			"{}/{}.merged.bam".format(
				output_directory, output_prefix)))
		subprocess.call(
			[samtools_path, "index", "{}/{}.merged.bam".format(
				output_directory, output_prefix)])
		samfile.close()
		bam_logger.info("Swapping complete.  Elapsed time: {} seconds".format(
			time.time() - switch_start))
		return "{}/{}.merged.bam".format(output_directory, output_prefix)

	else:
		# Remove sex chromosomes from original bam
		# samfile = pysam.AlignmentFile(bam_orig, "rc")
		# non_sex_scaffolds = filter(
		# 	lambda x: x not in sex_chroms, list(samfile.references))
		# with open("{}/no_sex.cram".format(output_directory), "w") as f:
		# 	subprocess.call(
		# 		[samtools_path, "view", "-h", "-b",
		# 			bam_orig, "{}".format(" ".join(non_sex_scaffolds))],
		# 		stdout=f)
		# subprocess.call(
		# 	[samtools_path, "index", "{}/no_sex.cram".format(output_directory)])
		#
		# # Merge bam files
		# subprocess.call(
		# 	[samtools_path, "merge", "-h",
		# 		"{}/header.sam".format(output_directory), "{}/{}.cram".format(
		# 			output_directory, output_prefix), "{}/no_sex.cram".format(
		# 				output_directory), bam_new])
		# subprocess.call(
		# 	[samtools_path, "index", "{}/{}.cram".format(
		# 		output_directory, output_prefix)])
		# return "{}/{}.cram".format(output_directory, output_prefix)
		bam_logger.error("This function does not currently handle cram files")
		logging.shutdown()
		raise RuntimeError(
			"switch_sex_chromosomes_sambamba does not currently support cram files")


def samtools_merge(samtools_path, bam_list, output_prefix, threads):
	"""
	Merges bam files using samtools.

	Parameters
	----------

	samtools_path : str
		The path to samtools
	bam_list : list
		Bam files to be merged. Merging order will match order of this list.
	output_prefix : str

	Returns
	-------
	str
		path to merged bam

	"""
	merge_start = time.time()
	command_line = [
		samtools_path, "merge", "-cp", "-@", "{}".format(threads), "-f",
		"{}.merged.bam".format(output_prefix)] + bam_list
	bam_logger.info(
		"Merging {} with command {}".format(
			" ".join(bam_list), " ".join(command_line)))
	subprocess.call(command_line)
	subprocess.call([
		samtools_path, "index", "{}.merged.bam".format(output_prefix)])
	bam_logger.info("Merging complete. Elapsed time: {} seconds".format(
		time.time() - merge_start))
	return "{}.merged.bam".format(output_prefix)
