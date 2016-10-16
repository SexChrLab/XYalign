import argparse
import numpy as np
import pysam
import sys
from ggplot import *


def traverseBam(samfile, chrom, start, stop, minDepth, afThresh, minAltReads):
	depth = []
	readBalance = []
	counter = 0
	passed = 0
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
			coord = int(pileupcolumn.pos)
			depth.append(totalDepth)
			alleleFraction = float(numAlt) / totalDepth
			if numAlt >= minAltReads and alleleFraction >= afThresh:
				readBalance.append(alleleFraction)
				passed += 1
			else:
				readBalance.append(-1)
		counter += 1
		if counter % 100000 == 0:
			print "%d sites processed, %d of which passed filters" % (counter, passed)
	return (np.asarray(readBalance), np.asarray(depth))
	
def main():
	########################################
	######## Parse the command line ########
	########################################
	parser = argparse.ArgumentParser(description="Add description")
	
	#Print help/usage if no arguments are supplied
	if len(sys.argv)==1:
		parser.print_usage()
	sys.exit(1)	
	
	# Parse command lines
	parser.add_argument("--bam", required=True, 
						help="REQUIRED. Input bam file")
	parser.add_argument("--fai", required=True,
						help="REQUIRED. Fasta index file (.fai) from reference genome.")
	parser.add_argument("--chrX_name", default="chrX",
						help="DEFAULT is chrX.  The name of the X chromosome scaffold.")
	parser.add_argument("--chrY_name", default="chrY",
						help="DEFAULT is chrY.  The name of the Y chromosome scaffold.")
	parser.add_argument("--minimum_depth", type=int, default=2,
						help="DEFAULT is 2. Minimum depth for a site to be considered. Integer.")
	parser.add_argument("--minimum_alt_allele_reads", type=int, default=1,
						help="DEFAULT is 1. Minimum number of reads with minor allele for read balance to be considered. Integer.")
	parser.add_argumber("--minor_allele_threshold", type=float, default=0.1,
						help="DEFAULT is 0.1. Minimum fraction of reads with minor allele for read balance to be considered. Float.")
	args = parser.parse_args()

	# Grab X and Y chromosome coordinates
	with open(args.fai,"r") as f:
		for i in csv.reader(f,delimiter="\t"):
			if i[0] == args.chrX_name:
				x_length = int(i[1])
			if i[0] == args.chrY_name:
				y_length = int(i[1])

	########################################
	########### Traverse bam file ##########
	########################################
    samfile = pysam.AlignmentFile(args.bam, "rb")
    results = traverseBam(samfile, args.chrX_name, 0, x_length, args.minimum_depth, args.minimum_alt_allele_reads, minor_allele_fraction))
    
    

if __name__ == "__main__":
	main()
		
		
