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
		if self.is_index() is False:
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
			self.logger.info("Indexing complete. Elapsed time: {}".format(
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
		bamfile = pysam.Alignmentfile(self.filepath, "rb")
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
		self, repairsh_path, single, output_directory,
		output_prefix, regions):
		"""
		Strips reads from a bam or cram file in provided regions and outputs
		sorted fastqs containing reads, one set of fastq files per read group.

		repairsh_path is the path to repair.sh (from BBmap)
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
		rg_list = output_directory + "/" + "full_rg.list"
		command_line = """{} view -H {} | awk '$1=="\x40RG"' | """\
			"""awk {} """\
			"""| cut -d':' -f 2 > {}""".format(
				self.samtools, self.filepath,
				repr('{for(i=1;i<=NF;i++){if (substr($i,1,2) ~ /ID/){print $i}}}'),
				rg_list)
		self.logger.info("Grabbing read groups from {} with the command: {}".format(
			self.filepath, command_line))
		subprocess.call(command_line, shell=True)
		rg_header_lines = output_directory + "/" + "header_lines_rg.list"
		command_line = """{} view -H {} | awk '$1=="\x40RG"' > {}""".format(
			self.samtools, self.filepath, rg_header_lines)
		self.logger.info(
			"Grabbing RG header lines from {} with the command: {}".format(
				self.filepath, command_line))
		subprocess.call(command_line, shell=True)
		with open(rg_list, "r") as f:
			out_rg_table = output_directory + "/" + "rg_fastq_key.list"
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
							command_line = "{} view -R {} -b {} {} | {} bam2fq -1 {}/temp_1.fastq "\
								"-2 {}/temp_2.fastq -t -n - ".format(
									self.samtools, tmp_out, self.filepath, ' '.join(map(str, regions)),
									self.samtools, output_directory, output_directory)
							self.logger.info(
								"Stripping paired reads with command: {}".format(
									command_line))
							subprocess.call(command_line, shell=True)
							command_line = "{} in1={} in2={} out1={} out2={} overwrite=true".format(
								repairsh_path,
								output_directory + "/temp_1.fastq",
								output_directory + "/temp_2.fastq",
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
								"{}/temp.fastq".format(
									self.samtools, tmp_out, self.filepath, ' '.join(map(
										str, regions)), self.samtools, output_directory)
							self.logger.info(
								"Stripping single-end reads with command: {}".format(
									command_line))
							subprocess.call(command_line, shell=True)
							command_line = "{} in={} out={} overwrite=true".format(
								repairsh_path,
								output_directory + "/temp.fastq",
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
			"Stripping and sorting reads complete. Elapsed time: {}".format(
				time.time() - rg_start))
		return [out_rg_table, rg_header_lines]

	def analyze_bam_fetch(self, chrom, window_size):
		"""Analyze BAM (or CRAM) file for various metrics.
		Currently, this function looks at the following metrics across genomic
		windows:
		- Read depth
		- Mapping quality
		The average of each metric will be calculated for each window of
		size `window_size` and stored altogether in a pandas data frame.

		chrom is the chromosome to analyze
		window size is the integer window size to use for sliding window analyses

		Returns:
			A dictionary of pandas data frames with the following key:
				- windows: The averages for each metric for each window
		"""
		self.logger.info(
			"Analyzing {} in {} for depth and mapq using windows of size {}".format(
				chrom, self.filepath, window_size))
		analyze_start = time.time()
		samfile = self.filepath
		chr_len = self.get_chrom_length(chrom)
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
			print("{} out of {} windows processed on {}".format(
				window_id, num_windows, chrom))

		# Convert data into pandas data frames
		windows_df = pd.DataFrame({
			"chrom": np.asarray(chr_list),
			"start": np.asarray(start_list),
			"stop": np.asarray(stop_list),
			"depth": np.asarray(depth_list),
			"mapq": np.asarray(mapq_list)
		})[["chrom", "start", "stop", "depth", "mapq"]]

		results = {"windows": windows_df}
		self.logger.info("Analysis complete. Elapsed time: {}".format(
			time.time() - analyze_start))
		return results

###############################################################################
# Functions associated with bam files


def switch_sex_chromosomes_bam_sambamba_output_temps(
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
	bam_logger.info(
		"Isolating header with command: {}".format(
			" ".join(
				[samtools_path, "view", "-H", bam_orig,
					">", "{}/header.sam".format(output_directory)])))
	with open("{}/header.sam".format(output_directory), "w") as f:
		subprocess.call(
			[samtools_path, "view", "-H", bam_orig], stdout=f)
	# Reheader new bam (for merge)
	bam_logger.info(
		"Reheading {} with command: {}".format(
			bam_new, " ".join(
				[samtools_path, "reheader", "-P", "{}/header.sam".format(
					output_directory),
					bam_new, ">", "{}/reheadered.temp.new.bam".format(
						output_directory)])))
	with open("{}/reheadered.temp.new.bam".format(output_directory), "w") as f:
		subprocess.call(
			[samtools_path, "reheader", "-P", "{}/header.sam".format(
				output_directory), bam_new], stdout=f)
	bam_logger.info("Indexing {}".format("{}/reheadered.temp.new.bam".format(
		output_directory)))
	subprocess.call(
		[samtools_path, "index", "{}/reheadered.temp.new.bam".format(
			output_directory)])
	# Add XYalign @PG line to header
	bam_logger.info("Adding new PG line to header")
	cl_string = " ".join(pg_header_dict["CL"])
	if "VN" in pg_header_dict:
		pg_line = [
			"@PG", "ID:{}".format(pg_header_dict["ID"]), "VN:{}".format(
				pg_header_dict["VN"]), "CL:{}".format(cl_string)]
	subprocess.call("echo '{}' >> {}/header.sam".format(
		"\t".join(pg_line), output_directory), shell=True)
	if cram is False:
		# Remove sex chromosomes from original bam and merge
		samfile = pysam.AlignmentFile(bam_orig, "rb")
		non_sex_scaffolds = filter(
			lambda x: x not in sex_chroms, list(samfile.references))
		command_line = "{} view -h -t {} -f bam -o {}/temp.nosexchr.bam {} {}".format(
			sambamba_path, threads, output_directory, bam_orig,
			" ".join(non_sex_scaffolds))
		bam_logger.info(
			"Removing sex chromosomes with command: {}".format(command_line))
		subprocess.call(command_line, shell=True)
		# subprocess.call(
		# 	[sambamba_path, "view", "-h", "-t", "{}".format(threads), "-f",
		# 		"bam", "-o", "{}/temp.nosexchr.bam".format(output_directory),
		# 		bam_orig, "{}".format(" ".join(non_sex_scaffolds))])
		bam_logger.info("Indexing {}".format("{}/temp.nosexchr.bam".format(
			output_directory)))
		subprocess.call(
			[samtools_path, "index", "{}/temp.nosexchr.bam".format(
				output_directory)])
		bam_logger.info("Merging bam files with command: {}".format(
			[samtools_path, "merge", "-@", "{}".format(threads), "-h",
				"{}/header.sam".format(output_directory), "-f",
				"{}/{}.merged.bam".format(output_directory, output_prefix),
				"{}/temp.nosexchr.bam".format(output_directory), bam_new]))
		subprocess.call(
			[samtools_path, "merge", "-@", "{}".format(threads), "-h",
				"{}/header.sam".format(output_directory), "-f",
				"{}/{}.merged.bam".format(output_directory, output_prefix),
				"{}/temp.nosexchr.bam".format(output_directory), bam_new])
		bam_logger.info("Indexing {}".format(
			"{}/{}.merged.bam".format(
				output_directory, output_prefix)))
		subprocess.call(
			[samtools_path, "index", "{}/{}.merged.bam".format(
				output_directory, output_prefix)])
		bam_logger.info("Swapping complete.  Elapsed time: {}".format(
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
	bam_logger.info("Merging complete. Elapsed time: {}".format(
		time.time() - merge_start))
	return "{}.merged.bam".format(output_prefix)

# Legacy functions (keeping them so they remain callable if needed later)
#
#
# def get_length(bamfile, chrom):
# 	"""
# 	Extract chromosome length from BAM header.
#
# 	args:
# 		bamfile: pysam AlignmentFile object
# 			- can be bam, cram, or sam, needs to be declared
# 				in pysam.AlignmentFile call before passing to function
# 		chrom: chromosome name (string)
#
# 	returns:
# 		Length (int)
#
# 	"""
# 	lengths = dict(zip(bamfile.references, bamfile.lengths))
# 	return lengths[chrom]
#
#
# def traverse_bam_fetch(samfile, chrom, window_size):
# 	"""Analyze the `samfile` BAM (or CRAM) file for various metrics.
# 	Currently, this function looks at the following metrics across genomic
# 	windows:
# 	- Read depth
# 	- Mapping quality
# 	The average of each metric will be calculated for each window of
# 	size `window_size` and stored altogether in a pandas data frame.
#
#
# 	samfile is a pysam AlignmentFile object
# 	chrom is the chromosome to analyze
# 	window size is the integer window size to use for sliding window analyses
#
# 	Returns:
# 		A dictionary of pandas data frames with the following key:
# 			- windows: The averages for each metric for each window
# 	"""
# 	chr_len = get_length(samfile, chrom)
# 	num_windows = chr_len // window_size + 1
# 	if chr_len % num_windows == 0:
# 		last_window_len = window_size
# 	else:
# 		last_window_len = chr_len % num_windows
#
# 	window_id = 0
#
# 	chr_list = [chrom] * num_windows
# 	start_list = []
# 	stop_list = []
# 	depth_list = []
# 	mapq_list = []
#
# 	start = 0
# 	end = window_size
# 	for window in range(0, num_windows):
# 		mapq = []
# 		total_read_length = 0
# 		for read in samfile.fetch(chrom, start, end):
# 			if read.is_secondary is False:
# 				if read.is_supplementary is False:
# 					total_read_length += read.infer_query_length()
# 					mapq.append(read.mapping_quality)
# 		start_list.append(start)
# 		stop_list.append(end)
# 		depth_list.append(total_read_length / window_size)
# 		mapq_list.append(np.mean(np.asarray(mapq)))
#
# 		window_id += 1
# 		if window_id == num_windows - 1:
# 			start += window_size
# 			end += last_window_len
# 		else:
# 			start += window_size
# 			end += window_size
#
# 		# Print progress
# 		print("{} out of {} windows processed on {}".format(
# 			window_id, num_windows, chrom))
#
# 	# Convert data into pandas data frames
# 	windows_df = pd.DataFrame({
# 		"chrom": np.asarray(chr_list),
# 		"start": np.asarray(start_list),
# 		"stop": np.asarray(stop_list),
# 		"depth": np.asarray(depth_list),
# 		"mapq": np.asarray(mapq_list)
# 	})[["chrom", "start", "stop", "depth", "mapq"]]
#
# 	results = {"windows": windows_df}
# 	return results
#
#
# def bam_to_fastq(
# 	samtools_path, repairsh_path, bamfile, single, output_directory,
# 	output_prefix, regions):
# 	"""
# 	Strips reads from a bam or cram file in provided regions and outputs
# 	sorted fastqs containing reads, one set of fastq files per read group.
#
# 	samtools_path is the path to samtools
# 	repairsh_path is the path to repair.sh (from BBmap)
# 	bamfile is the input bam (including path)
# 	single is either True or False; if true will output single-end fastq file,
# 		if False, will output paired-end fastq files
# 	output_directory is the directory for ALL output (including temporary files)
# 	output_prefix is the name (without path) to use for prefix to output fastqs
# 	regions is a list of regions from which reads will be stripped
#
# 	Returns:
# 		A two-item list containing the path to a text file pairing read group
# 			names with associated output fastqs, and a text file containing a
# 			list of @RG lines associated with each read group
# 	"""
# 	# Collect RGs
# 	rg_list = output_directory + "/" + "full_rg.list"
# 	command_line = """{} view -H {} | awk '$1=="\x40RG"' | """\
# 		"""awk {} """\
# 		"""| cut -d':' -f 2 > {}""".format(
# 			samtools_path, bamfile,
# 			repr('{for(i=1;i<=NF;i++){if (substr($i,1,2) ~ /ID/){print $i}}}'),
# 			rg_list)
# 	subprocess.call(command_line, shell=True)
# 	rg_header_lines = output_directory + "/" + "header_lines_rg.list"
# 	command_line = """{} view -H {} | awk '$1=="\x40RG"' > {}""".format(
# 		samtools_path, bamfile, rg_header_lines)
# 	subprocess.call(command_line, shell=True)
# 	with open(rg_list, "r") as f:
# 		out_rg_table = output_directory + "/" + "rg_fastq_key.list"
# 		with open(out_rg_table, "w") as ortab:
# 			for line in f:
# 				rg = line.strip()
# 				if rg != "":
# 					tmp_out = "{}/{}.txt".format(output_directory, rg)
# 					with open(tmp_out, "w") as o:
# 						o.write(rg)
# 					if single is False:
# 						command_line = "{} view -R {} -b {} {} | {} bam2fq -1 {}/temp_1.fastq "\
# 							"-2 {}/temp_2.fastq -t -n - ".format(
# 								samtools_path, tmp_out, bamfile, ' '.join(map(str, regions)),
# 								samtools_path, output_directory, output_directory)
# 						subprocess.call(command_line, shell=True)
# 						command_line = "{} in1={} in2={} out1={} out2={} overwrite=true".format(
# 							repairsh_path,
# 							output_directory + "/temp_1.fastq",
# 							output_directory + "/temp_2.fastq",
# 							output_directory + "/" + output_prefix + "_" + rg + "_1.fastq",
# 							output_directory + "/" + output_prefix + "_" + rg + "_2.fastq")
# 						subprocess.call(command_line, shell=True)
# 						ortab.write("{}\t{}\t{}\n".format(
# 							rg,
# 							output_directory + "/" + output_prefix + "_" + rg + "_1.fastq",
# 							output_directory + "/" + output_prefix + "_" + rg + "_2.fastq"))
# 					else:
# 						command_line = "{} view -R {} -b {} {} | {} bam2fq -t -n - > "\
# 							"{}/temp.fastq".format(samtools_path, tmp_out, bamfile, ' '.join(map(
# 								str, regions)), samtools_path, output_directory)
# 						subprocess.call(command_line, shell=True)
# 						command_line = "{} in={} out={} overwrite=true".format(
# 							repairsh_path,
# 							output_directory + "/temp.fastq",
# 							output_directory + "/" + output_prefix + "_" + rg + ".fastq")
# 						# write line
# 						ortab.write("{}\t{}\n".format(
# 							rg,
# 							output_directory + "/" + output_prefix + "_" + rg + ".fastq"))
# 	return [out_rg_table, rg_header_lines]
#
# def switch_sex_chromosomes_bam_sambamba(
# 	samtools_path, sambamba_path, bam_orig, bam_new, sex_chroms,
# 	output_directory, output_prefix, threads, pg_header_dict, cram=False):
# 	"""
# 	Removes sex chromosomes from original bam and merges in remmapped
# 	sex chromosomes, while retaining the original bam header (and adding new
# 	@PG line)
#
# 	samtools_path is the path to samtools
# 	sambamba_path is the path to sambamba
# 	bam_orig is the original full bam file
# 	bam_new is the bam containing the sex chromosomes
# 	sex_chroms is a list of sex chromosomes (to be removed from bam_orig)
# 	output_directory is the path to directory where all files (inc. temp) will
# 			be output
# 	threads is the number of threads/cpus to use
# 	pg_header_dict is a dictionary with information to be included in the new
# 		@PG line
# 			- must contain:
# 				Key = 'CL', value = list of command line values
# 				Key = 'ID', value = string of program ID
# 			- optional:
# 				Key = 'VN', value = string of program version
# 	cram (default is False) - if True, will treat input as cram files and
# 		output cram files.  Right now slower, with more intermediate/temp files
#
# 	Returns:
# 		New bam or cram file path with original header (plus new @PG line), but sex
# 			chromosomes swapped out
# 	"""
# 	# Grab original header
# 	with open("{}/header.sam".format(output_directory), "w") as f:
# 		subprocess.call(
# 			[samtools_path, "view", "-H", bam_orig], stdout=f)
# 	# Reheader new bam (for merge)
# 	with open("{}/reheadered.temp.new.bam".format(output_directory), "w") as f:
# 		subprocess.call(
# 			[samtools_path, "reheader", "-P", "{}/header.sam".format(
# 				output_directory), bam_new], stdout=f)
# 	# Add XYalign @PG line to header
# 	cl_string = " ".join(pg_header_dict["CL"])
# 	if "VN" in pg_header_dict:
# 		pg_line = [
# 			"@PG", "ID:{}".format(pg_header_dict["ID"]), "VN:{}".format(
# 				pg_header_dict["VN"]), "CL:{}".format(cl_string)]
# 	subprocess.call("echo '{}' >> {}/header.sam".format(
# 		"\t".join(pg_line), output_directory), shell=True)
# 	if cram is False:
# 		# Remove sex chromosomes from original bam and merge
# 		samfile = pysam.AlignmentFile(bam_orig, "rb")
# 		non_sex_scaffolds = filter(
# 			lambda x: x not in sex_chroms, list(samfile.references))
# 		subprocess.call(
# 			"{} view -h -t {} -f bam -o /dev/stdout {} {} | "
# 			"{} merge -t {} {}/{}.merged.bam /dev/stdin {}".format(
# 				sambamba_path, threads, bam_orig, " ".join(non_sex_scaffolds),
# 				sambamba_path, threads, output_directory, output_prefix,
# 				bam_new), shell=True)
#
# 		return "{}/{}.merged.bam".format(output_directory, output_prefix)
#
# 	else:
# 		# Remove sex chromosomes from original bam
# 		samfile = pysam.AlignmentFile(bam_orig, "rc")
# 		non_sex_scaffolds = filter(
# 			lambda x: x not in sex_chroms, list(samfile.references))
# 		with open("{}/no_sex.cram".format(output_directory), "w") as f:
# 			subprocess.call(
# 				[samtools_path, "view", "-h", "-b",
# 					bam_orig, "{}".format(" ".join(non_sex_scaffolds))],
# 				stdout=f)
# 		subprocess.call(
# 			[samtools_path, "index", "{}/no_sex.cram".format(output_directory)])
#
# 		# Merge bam files
# 		subprocess.call(
# 			[samtools_path, "merge", "-h",
# 				"{}/header.sam".format(output_directory), "{}/{}.cram".format(
# 					output_directory, output_prefix), "{}/no_sex.cram".format(
# 						output_directory), bam_new])
# 		subprocess.call(
# 			[samtools_path, "index", "{}/{}.cram".format(
# 				output_directory, output_prefix)])
# 		return "{}/{}.cram".format(output_directory, output_prefix)
#
#
# def switch_sex_chromosomes_bam(
# 	samtools_path, bam_orig, bam_new, sex_chroms, output_directory,
# 	output_prefix, cram=False):
# 	"""
# 	Removes sex chromosomes from original bam and merges in remmapped
# 	sex chromosomes, while retaining the original bam header
#
#
# 		samtools_path is the path to samtools
# 		bam_orig is the original full bam file
# 		bam_new is the bam containing the sex chromosomes
# 		sex_chroms is a list of sex chromosomes (to be removed from bam_orig)
# 		output_directory is the path to directory where all files (inc. temp) will
# 				be output
# 		cram (default is False) - if True, will treat input as cram files and
# 			output cram files.  Right now slower, with more intermediate/temp files
#
# 		Returns:
# 			New bam or cram file path with original header (plus new @PG line), but sex
# 				chromosomes swapped out
# 	"""
# 	# Grab original header
# 	subprocess.call(
# 		"{} view -H {} > {}/header.sam".format(
# 			samtools_path, bam_orig, output_directory), shell=True)
# 	if cram is False:
# 		# Remove sex chromosomes from original bam
# 		samfile = pysam.AlignmentFile(bam_orig, "rb")
# 		non_sex_scaffolds = filter(
# 			lambda x: x not in sex_chroms, list(samfile.references))
# 		subprocess.call(
# 			"{} view -h -b {} {} > {}/no_sex.bam".format(
# 				samtools_path, bam_orig, " ".join(non_sex_scaffolds),
# 				output_directory),
# 			shell=True)
# 		subprocess.call(
# 			"{} index {}/no_sex.bam".format(
# 				samtools_path, output_directory), shell=True)
#
# 		# Merge bam files
# 		subprocess.call(
# 			"{} merge -h {}/header.sam {}/{}.bam {}/no_sex.bam {}".format(
# 				samtools_path, output_directory, output_directory,
# 				output_prefix, output_directory, bam_new), shell=True)
# 		subprocess.call("{} index {}/{}.bam".format(
# 			samtools_path, output_directory, output_prefix), shell=True)
#
# 		return "{}/{}.bam".format(output_directory, output_prefix)
#
# 	else:
# 		# Remove sex chromosomes from original bam
# 		samfile = pysam.AlignmentFile(bam_orig, "rc")
# 		non_sex_scaffolds = filter(
# 			lambda x: x not in sex_chroms, list(samfile.references))
# 		subprocess.call(
# 			"{} view -h -b {} {} > {}/no_sex.cram".format(
# 				samtools_path, bam_orig, " ".join(non_sex_scaffolds),
# 				output_directory),
# 			shell=True)
# 		subprocess.call(
# 			"{} index {}/no_sex.cram".format(
# 				samtools_path, output_directory), shell=True)
#
# 		# Merge bam files
# 		subprocess.call(
# 			"{} merge -h {}/header.sam {}/{}.cram {}/no_sex.cram {}".format(
# 				samtools_path, output_directory, output_directory,
# 				output_prefix, output_directory, bam_new), shell=True)
# 		subprocess.call("{} index {}/{}.cram".format(
# 			samtools_path, output_directory, output_prefix), shell=True)
#
# 		return "{}/{}.cram".format(output_directory, output_prefix)
