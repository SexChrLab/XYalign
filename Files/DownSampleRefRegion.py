import pysam
import random
import sys
import gzip


def makeBaseCountDict(samfile, chrom, testPos):
    '''
    (str, int) -> dict
    '''
    baseCountDict = {'A':0, 'C':0, 'G':0, 'T':0, 'N': 0}
    for pileupcolumn in samfile.pileup(chrom, testPos - 5, testPos + 5):
        if int(pileupcolumn.pos) == testPos - 1:
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    base = pileupread.alignment.query_sequence[pileupread.query_position]
                    baseCountDict[base] += 1
    return baseCountDict




def makePosList(vcfFile = '/DCEG/CGF/Infinium/Clustering/METH/NonMethProjects/Eric/RbVariantCalling/mosaicVariants/RB1.ploidy10.combined.annotated.vcf.gz'):
    posList = []
    with gzip.open(vcfFile) as f:
        for line in f:
            if line[0] != '#':
                line_list = line.split()
                (chrom, pos, x, ref, alt) = line_list[:5]
                pos = int(pos)
                if ',' not in alt and len(alt) == 1 and len(ref) == 1:
                    posList.append((chrom, pos, ref, alt))
    return posList



def makeReadRemoveList(samfile, output, chrom, pos, afThresh = 0.02, minAltReads = 3):
    '''
    (pysam.alignment, pysam.alignment, str, int, str, str, float) -> None
    '''
    readsToRemove = []
    baseCountDict = makeBaseCountDict(samfile, chrom, pos)
    baseCounts = sorted(baseCountDict.values())
    numRef = baseCounts[-1]
    if numRef == 0:
        return readsToRemove
    numAlt = baseCounts[-2]
    if numAlt < minAltReads or (float(numAlt)/(numAlt + numRef) < afThresh):
        return readsToRemove
    ref = 'N'
    alt = 'N'
    for base in baseCountDict.keys():
        if baseCountDict[base] == numRef:
            ref = base
        elif baseCountDict[base] == numAlt:
            alt = base
    if ref == 'N' or alt == 'N':
        return readsToRemove
    if numRef == 0 or float(numAlt)/(numAlt + numRef) < afThresh or numAlt < minAltReads:
        refToKeepList = xrange(numRef)
    else:
        if numAlt + 2 <= numRef:
            refToKeepList = random.sample(xrange(numRef), numAlt + 2)
        else:
            refToKeepList = random.sample(xrange(numRef), numAlt)
    refToKeepDict = {}
    for r in refToKeepList:
        refToKeepDict[r] =1 
    refCount = 0
    for pileupcolumn in samfile.pileup(chrom, pos - 5, pos + 5):
        if int(pileupcolumn.pos) == pos - 1:
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    base = pileupread.alignment.query_sequence[pileupread.query_position]
                    if base == ref:
                        if not refToKeepDict.get(refCount):
                            readsToRemove.append(pileupread.alignment.query_name)
                        refCount += 1

    return readsToRemove

def outputNewBam(inBam, outBam, chrom, start, end):
    '''
    (str, str, str, int, int) -> None
    '''
    samfile = pysam.AlignmentFile(inBam, "rb")
    output = pysam.AlignmentFile(outBam, "wb", template=samfile)
    allReadRemoveDict = {}
    for pos in range(start, end+1):
        readsToRemove = makeReadRemoveList(samfile, output, chrom, pos)
        for read in readsToRemove:
            allReadRemoveDict[read] = 1
    for read in samfile.fetch(chrom, start - 1000, end + 1000):
        if not allReadRemoveDict.get(read.query_name):
            output.write(read)
    samfile.close()
    output.close()




def main():
    ##usage:  module load python/2.7.8;python DownSampleRefRegion.py in.bam out.bam chrom start end
    args = sys.argv[1:]
    outputNewBam(args[0], args[1], args[2], int(args[3]), int(args[4]))
    
    
if __name__ == "__main__":
    main()
