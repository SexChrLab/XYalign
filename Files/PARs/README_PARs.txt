This directory contains files with coordinates for the pseudoautosomal regions (PARs)
for GRCh37, GRCh38, hg19, and hg38.  Files ending in Ymask.bed contain the assembled PAR coordinates,
while files ending in Ymask_startEnd.bed contain the assumed full PARs (including unassembled sequence
denoted by Ns in the assembly).  We recommend using the latter files (ending in Ymask_startEnd.bed) as
masks for XYalign.

Location of GRCh38/hg38 and GRCh37/hg19 PARs: https://www.ncbi.nlm.nih.gov/grc/human

Total length of GRCh38/hg38 chromosomes:
https://www.ncbi.nlm.nih.gov/grc/human/data

Total length of GRCh37/hg19 chromosomes:
https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes

Update on Dec 12, 2016: Coordinates updated to reflect 0 based indexing in bed files

Update on Dec 16, 2016: Fixed issues with chromosome names for different assemblies:
	hg19 and hg38 use "chrY"
	GRCh37 and GRCh38 use "Y"
