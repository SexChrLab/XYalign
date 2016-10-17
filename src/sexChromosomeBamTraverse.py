import sys
import argparse
import numpy as np
import pysam
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def traverseBam(samfile, chrom, start, stop, window_size, minDepth, afThresh, minAltReads):
	"""Analyze a BAM file for various metrics and statistics.

	Returns a tuple of pandas data frames:
	- Statistics for each window
	- Frequencies for mapping quality
	- Frequencies for read depth
	- Frequencies for read balance
	"""
	depthCollector = []
	readBalanceCollector = []
	counter = 0
	histDict = {}
	
	coord = []
	depth = []
	readBalance = []
	windowCounter = 0
	for pileupcolumn in samfile.pileup(chrom, start, stop):
		
		baseCountDict = {'A':0, 'C':0, 'G':0, 'T':0, 'N': 0}
		for pileupread in pileupcolumn.pileups:
			if not pileupread.is_del and not pileupread.is_refskip:
				baseCountDict[pileupread.alignment.query_sequence[pileupread.query_position]] += 1
		baseCounts = sorted(baseCountDict.values())
		numRef = baseCounts[-1]
		numAlt = baseCounts[-2]
		totalDepth = numRef + numAlt
		
		if totalDepth >= minDepth:
			depth.append(totalDepth)
			alleleFraction = float(numAlt) / totalDepth
			if numAlt >= minAltReads and alleleFraction >= afThresh:
				readBalance.append(alleleFraction)
				
		windowCounter += 1
		if windowCounter == window_size:
			### process window
			depthCollector.append(np.mean(np.asarray(depth)))
			readBalanceCollector.append(np.mean(np.asarray(readBalance)))
			
			### Reset counters
			windowCounter == 0
			depth = []
			readBalance = []
			coord.append(int(pileupcolumn.pos))

		counter += 1
		if counter % 100000 == 0:
			print "%d sites processed, %d of which passed filters" % (counter, passed)
	return (np.asarray(readBalance), np.asarray(depth))


def plot_data(dataDict):
	"""
	Takes a dictionary (output from traverseBam) and outputs histograms and 
	genome-wide plots of various metrics
	
	"""
	
	window_df = dataDict[windows]
	depth_hist = dataDict[depth_freq]
	readbal_hist = dataDict[readbal_freq]
	mapq_hist = dataDict[mapq_freq]
	
	## Create genome-wide plots based on window means
	depth_genome_plot = sns.lmplot('Position','Depth', data=window_df, fit_reg=False)
	depth_genome_plot.savefig("depth_windows.png")
	
	balance_genome_plot = sns.lmplot('Position','ReadBalance', data=window_df, fit_reg=False)
	balance_genome_plot.savefig("balance_windows.png")
	
	mapq_genome_plot = sns.lmplot('Position','ReadBalance', data=window_df, fit_reg=False)
	mapq_genome_plot.savefig("mapq_windows.png")
	
	## Create histograms
	depth_bar_plot = sns.countplot(x='Depth', y='Count', data=depth_hist)
	depth_bar_plot.savefig("depth_hist.png")
	
	balance_bar_plot = sns.countplot(x='ReadBalance', y='Count', data=balance_hist)
	balance_bar_plot.savefig("readbalance_hist.png")
	
	mapq_bar_plot = sns.countplot(x='Mapq', y='Count', data=mapq_hist)
	mapq_bar_plot.savefig("mapq_hist.png")
	
	pass
	


def main():
	########################################
	######## Parse the command line ########
	########################################
	parser = argparse.ArgumentParser(description="Add description")
	
	# Parse command lines
	parser.add_argument("--bam", required=True, 
						help="REQUIRED. Input bam file")
	# parser.add_argument("--fai", required=True,
	# 					help="REQUIRED. Fasta index file (.fai) from reference genome.")
	parser.add_argument("--chrX_name", default="chrX",
						help="DEFAULT is chrX.  The name of the X chromosome scaffold.")
	parser.add_argument("--chrY_name", default="chrY",
						help="DEFAULT is chrY.  The name of the Y chromosome scaffold.")
	parser.add_argument("--chromosomes", "-c", nargs="+", default=["chrX", "chrY", "chr19"],
	                    help="Chromosomes to analyze.")
	parser.add_argument("--window_size", type=int, default=50000,
						help="DEFAULT is 50000.  Window size for sliding window calculations. Integer.")
	parser.add_argument("--minimum_depth", type=int, default=2,
						help="DEFAULT is 2. Minimum depth for a site to be considered. Integer.")
	parser.add_argument("--minimum_alt_allele_reads", type=int, default=1,
						help="DEFAULT is 1. Minimum number of reads with minor allele for read balance to be considered. Integer.")
	parser.add_argument("--minor_allele_threshold", type=float, default=0.1,
						help="DEFAULT is 0.1. Minimum fraction of reads with minor allele for read balance to be considered. Float.")

	#Print help/usage if no arguments are supplied
	if len(sys.argv)==1:
		parser.print_usage()
		sys.exit(1)	

	args = parser.parse_args()

	########################################
	########### Traverse bam file ##########
	########################################

	samfile = pysam.AlignmentFile(args.bam, "rb")

	all_lengths = zip(samfile.references, samfile.lengths)
	selected_lengths = dict(filter(lambda x: x[0] in args.chromosomes, all_lengths))

	for chromosome, length in selected_lengths.items():
		data = traverseBam(samfile, chromosome, 0, length, args.minimum_depth, args.minimum_alt_allele_reads, minor_allele_fraction)
		plot_data(data)
    
    
if __name__ == "__main__":
	main()
