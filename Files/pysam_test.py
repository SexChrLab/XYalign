import numpy as np
import pysam
import sys


def functionName(samfile, chrom, start, stop, minDepth = 2, afThresh = 0.02, minAltReads = 1):
	depth = []
	readBalance = []
	counter = 0
	for pileupcolumn in samfile.pileup(chrom, start, stop):
		baseCountDict = {'A':0, 'C':0, 'G':0, 'T':0, 'N': 0}
		for pileupread in pileupcolumn.pileups:
			if not pileupread.is_del and not pileupread.is_refskip:
				baseCountDict[pileupread.alignment.query_sequence[pileupread.query_position]] += 1
		baseCounts = sorted(baseCountDict.values())
		numRef = baseCounts[-1]
		numAlt = baseCounts[-2]
		totalDepth = numRef + numAlt
		#alleleFraction = float(numAlt) / totalDepth
		if totalDepth >= minDepth and numAlt >= minAltReads and (float(numAlt) / totalDepth) >= afThresh:
			depth.append(totalDepth)
			#readBalance.append(alleleFraction)
			readBalance.append(float(numAlt) / totalDepth)
		counter += 1
		if counter % 1000 == 0:
			print "%d sites processed, %d of which passed filters" % (counter, length(depth))
	return np.mean(np.asarray(readBalance))
	
def main():
    ##usage:  module load python/2.7.8;python DownSampleRefRegion.py in.bam out.bam chrom start end
    inBam = sys.argv[-1]
    samfile = pysam.AlignmentFile(inBam, "rb")
    print(functionName(samfile, 'chrX', 1, 10000000, 3, 0.02, 1))

if __name__ == "__main__":
    main()
			
