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

	Attributes:
		filepath is full path to external bam file
		samtools path is path to samtools (defaults to "samtools")
	"""
	def __init__(self, filepath, samtools="samtools"):
		""" Initiate object with an associated filepath """
		self.filepath = filepath
		self.samtools = samtools
		self.logger = logging.getLogger("xyalign.bam.BamFile")
		self.logger.info("Creating a BamFile instance for {}".format(
			self.filepath))
		if self.is_indexed() is False:
			self.index_bam()

	def is_indexed(self):
		"""
		Checks that bam index exists, is not empty, and is newer than bam.
		If either cases is false, returns False.  Otherwise, returns True.
		"""
		self.logger.info("Checking indexing of {}".format(self.filepath))
		if os.path.exists("{}.bai".format(self.filepath)):
			if os.stat("{}.bai".format(self.filepath)).st_size != 0:
				idx_stamp = os.path.getmtime("{}.bai".format(self.filepath))
			else:
				self.logger.info("Bam index empty")
				return False
		elif os.path.exists("{}.bai".format(self.filepath[:-3])):
			if os.stat("{}.bai".format(self.filepath[:-3])).st_size != 0:
				idx_stamp = os.path.getmtime("{}.bai".format(self.filepath[:-3]))
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
		Indexes a bam using samtools
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
			raise RuntimeError("Unable to index bamfile. Exiting")

	def get_chrom_length(self, chrom):
		"""
		Extract chromosome length from BAM header.

		args:
			bamfile: pysam AlignmentFile object
				- can be bam, cram, or sam, needs to be declared
					in pysam.AlignmentFile call before passing to function
			chrom: chromosome name (string)

		returns:
			Length (int)

		"""
		bamfile = pysam.AlignmentFile(self.filepath, "rb")
		lengths = dict(zip(bamfile.references, bamfile.lengths))
		try:
			return lengths[chrom]
		except:
			self.logger.error(
				"{} not present in bam header for {}. Exiting.".format(
					chrom, self.filepath))
			raise RuntimeError(
				"Chromosome name not present in bam header. Exiting")

	def strip_reads(
		self, repairsh, single, output_directory,
		output_prefix, regions):
		"""
		Strips reads from a bam or cram file in provided regions and outputs
		sorted fastqs containing reads, one set of fastq files per read group.

		repairsh is the path to repair.sh (from BBmap)
		single is either True or False; if true will output single-end fastq file,
			if False, will output paired-end fastq files
		output_directory is the directory for ALL output (including temporary files)
		output_prefix is the name (without path) to use for prefix to output fastqs
		regions is a list of regions from which reads will be stripped

		Returns:
			A two-item list containing the path to a text file pairing read group
				names with associated output fastqs, and a text file containing a
				list of @RG lines associated with each read group
		"""
		# Collect RGs
		rg_start = time.time()
		rg_list = output_directory + "/" + output_prefix + ".full_rg.list"
		command_line = """{} view -H {} | awk '$1=="\x40RG"' | """\
			"""awk {} """\
			"""| cut -d':' -f 2 > {}""".format(
				self.samtools, self.filepath,
				repr('{for(i=1;i<=NF;i++){if (substr($i,1,2) ~ /ID/){print $i}}}'),
				rg_list)
		self.logger.info("Grabbing read groups from {} with the command: {}".format(
			self.filepath, command_line))
		subprocess.call(command_line, shell=True)
		rg_header_lines = output_directory + "/" + output_prefix + ".header_lines_rg.list"
		command_line = """{} view -H {} | awk '$1=="\x40RG"' > {}""".format(
			self.samtools, self.filepath, rg_header_lines)
		self.logger.info(
			"Grabbing RG header lines from {} with the command: {}".format(
				self.filepath, command_line))
		subprocess.call(command_line, shell=True)
		with open(rg_list, "r") as f:
			out_rg_table = output_directory + "/" + output_prefix + ".rg_fastq_key.list"
			self.logger.info(
				"Iteratively removing reads by read group. Writing table of RGs and "
				"fastqs to {}".format(out_rg_table))
			with open(out_rg_table, "w") as ortab:
				for line in f:
					rg = line.strip()
					if rg != "":
						self.logger.info("Removing reads from group: {}".format(rg))
						tmp_out = "{}/{}.txt".format(output_directory, rg)
						with open(tmp_out, "w") as o:
							o.write(rg)
						if single is False:
							command_line = "{} view -R {} -b {} {} | {} bam2fq -1 {}/{}.temp_1.fastq "\
								"-2 {}/{}.temp_2.fastq -t -n - ".format(
									self.samtools, tmp_out, self.filepath, ' '.join(map(str, regions)),
									self.samtools, output_directory, output_prefix,
									output_directory, output_prefix)
							self.logger.info(
								"Stripping paired reads with command: {}".format(
									command_line))
							subprocess.call(command_line, shell=True)
							command_line = "{} in1={} in2={} out1={} out2={} overwrite=true".format(
								repairsh,
								output_directory + "/{}.temp_1.fastq".format(output_prefix),
								output_directory + "/{}.temp_2.fastq".format(output_prefix),
								output_directory + "/" + output_prefix + "_" + rg + "_1.fastq",
								output_directory + "/" + output_prefix + "_" + rg + "_2.fastq")
							self.logger.info(
								"Sorting reads with command: {}".format(
									command_line))
							subprocess.call(command_line, shell=True)
							ortab.write("{}\t{}\t{}\n".format(
								rg,
								output_directory + "/" + output_prefix + "_" + rg + "_1.fastq",
								output_directory + "/" + output_prefix + "_" + rg + "_2.fastq"))
						else:
							command_line = "{} view -R {} -b {} {} | {} bam2fq -t -n - > "\
								"{}/{}.temp.fastq".format(
									self.samtools, tmp_out, self.filepath, ' '.join(map(
										str, regions)), self.samtools, output_directory, output_prefix)
							self.logger.info(
								"Stripping single-end reads with command: {}".format(
									command_line))
							subprocess.call(command_line, shell=True)
							command_line = "{} in={} out={} overwrite=true".format(
								repairsh,
								output_directory + "/{}.temp.fastq".format(output_prefix),
								output_directory + "/" + output_prefix + "_" + rg + ".fastq")
							self.logger.info(
								"Sorting reads with command: {}".format(
									command_line))
							subprocess.call(command_line, shell=True)
							# write line
							ortab.write("{}\t{}\n".format(
								rg,
								output_directory + "/" + output_prefix + "_" + rg + ".fastq"))
		self.logger.info(
			"Stripping and sorting reads complete. Elapsed time: {} seconds".format(
				time.time() - rg_start))
		return [out_rg_table, rg_header_lines]

	def analyze_bam_fetch(self, chrom, window_size, target_file=None):
		"""Analyze BAM (or CRAM) file for various metrics.
		Currently, this function looks at the following metrics across genomic
		windows:
		- Read depth
		- Mapping quality
		The average of each metric will be calculated for each window of
		size `window_size` and stored altogether in a pandas data frame.

		chrom is the chromosome to analyze
		window_size is the integer window size to use for sliding window analyses
			- if set to None, will use intervals from target_file
		target_file is a bed_file containing regions to analyze instead of
			windows of a fixed size
				- will only be engaged if window_size is None

		Returns:
			A dictionary of pandas data frames with the following key:
				- windows: The averages for each metric for each window
		"""
		self.logger.info(
			"Traversing {} in {} to analyze depth and mapping quality".format(
				chrom, self.filepath))
		analyze_start = time.time()
		samfile = pysam.AlignmentFile(self.filepath, "rb")
		chr_len = self.get_chrom_length(chrom)

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
			for window in range(0, num_windows):
				mapq = []
				total_read_length = 0
				for read in samfile.fetch(chrom, start, end):
					if read.is_secondary is False:
						if read.is_supplementary is False:
							total_read_length += read.infer_query_length()
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

			for window in targets:
				mapq = []
				total_read_length = 0
				start = int(window[1])
				end = int(window[2])
				window_size = end - start
				for read in samfile.fetch(window[0], start, end):
					if read.is_secondary is False:
						if read.is_supplementary is False:
							total_read_length += read.infer_query_length()
							mapq.append(read.mapping_quality)
				chr_list.append(window[0])
				start_list.append(start)
				stop_list.append(end)
				depth_list.append(total_read_length / window_size)
				mapq_list.append(np.mean(np.asarray(mapq)))

				window_id += 1

				# Print progress
				self.logger.info("{} out of {} targets processed on {}".format(
					window_id, num_windows, chrom))

		else:
			self.logger.error(
				"Both window_size and target_file set to None. "
				"Cannot proceed with bam traversal. Exiting.")
			sys.exit(1)

		# Convert data into pandas data frames
		windows_df = pd.DataFrame({
			"chrom": np.asarray(chr_list),
			"start": np.asarray(start_list),
			"stop": np.asarray(stop_list),
			"depth": np.asarray(depth_list),
			"mapq": np.asarray(mapq_list)
		})[["chrom", "start", "stop", "depth", "mapq"]]

		results = {"windows": windows_df}
		self.logger.info("Analysis complete. Elapsed time: {} seconds".format(
			time.time() - analyze_start))
		return results

###############################################################################
# Functions associated with bam files


def switch_sex_chromosomes_sambamba(
	samtools_path, sambamba_path, bam_orig, bam_new, sex_chroms,
	output_directory, output_prefix, threads, pg_header_dict, cram=False):
	"""

	Removes sex chromosomes from original bam and merges in remmapped
	sex chromosomes, while retaining the original bam header (and adding new
	@PG line)

	samtools_path is the path to samtools
	sambamba_path is the path to sambamba
	bam_orig is the original full bam file
	bam_new is the bam containing the sex chromosomes
	sex_chroms is a list of sex chromosomes (to be removed from bam_orig)
	output_directory is the path to directory where all files (inc. temp) will
			be output
	output_prefix is the name (without path) to use for prefix for all files
	threads is the number of threads/cpus to use
	pg_header_dict is a dictionary with information to be included in the new
		@PG line
			- must contain:
				Key = 'CL', value = list of command line values
				Key = 'ID', value = string of program ID
			- optional:
				Key = 'VN', value = string of program version
	cram (default is False) - if True, will treat input as cram files and
		output cram files.  Right now slower, with more intermediate/temp files

	Returns:
		New bam or cram file path with original header (plus new @PG line), but sex
			chromosomes swapped out
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
	subprocess.call("echo '{}' >> {}/{}.header.sam".format(
		"\t".join(pg_line), output_directory, output_prefix), shell=True)
	if cram is False:
		# Remove sex chromosomes from original bam and merge
		samfile = pysam.AlignmentFile(bam_orig, "rb")
		non_sex_scaffolds = filter(
			lambda x: x not in sex_chroms, list(samfile.references))
		command_line = "{} view -h -t {} -f bam -o {}/{}.temp.nosexchr.bam {} {}".format(
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
			samtools_path, "merge", "-@", "{}".format(threads), "-h",
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
		return None


def sambamba_merge(sambamba_path, bam_list, output_prefix, threads):
	"""
	Takes a list of bam files, e.g., [bam1,bam2,bam3,...], and merges them
	using sambamba

	Returns:
		path to merged bam
	"""
	merge_start = time.time()
	bam_logger.info("Merging {}".format(" ".join(bam_list)))
	subprocess.call(
		[sambamba_path, "merge", "-t", str(threads), output_prefix, "{}".format(
			" ".join(bam_list))])
	subprocess.call([
		sambamba_path, "index", "{}.merged.bam".format(output_prefix)])
	bam_logger.info("Merging complete. Elapsed time: {} seconds".format(
		time.time() - merge_start))
	return "{}.merged.bam".format(output_prefix)
