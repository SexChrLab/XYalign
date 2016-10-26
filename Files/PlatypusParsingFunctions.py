import sys, os
from matplotlib import pyplot as plt



#This function will take in a Platypus VCF, and a quality cutoff for variants, and return three arrays: 
#1. Position of variants
#2. Quality of variants
#3. Read Balance of variants


def ParseVCF(filename,qualCutoff):
    infile = open("%s"%filename,'r')
    positions = []
    quality = []
    readBalance = []
    for line in infile:
        if line[0]=='#':
            continue
        cols=line.strip('\n').split('\t')
        pos = int(cols[1])
        qual = float(cols[5])
        if qual < qualCutoff:
            continue
        TR = cols[7].split(';')[17].split('=')[1]
        TC = cols[7].split(';')[14].split('=')[1]
        if ',' in TR or ',' in TC:
            continue
        if (float(TR)==0) or (float(TC) == 0):
            continue    
        ReadRatio = float(TR)/float(TC)
        
        # Add to arrays
        readBalance.append(ReadRatio)
        positions.append(pos)
        quality.append(qual)
        
    
    return positions,quality,readBalance
  
  
  
#This Function will plot the read balance (variant allele / total reads mapped), which is the array output (readBalance, positions) from the above functin
# This function also takes in the markerSize and markerAlpha, and x-limit.  Play with these, but examples are given below
  
def PlotReadBalance(positions,readBalance,sampleID,MarkerSize,MarkerAlpha,Xlim):
    if "X" in sampleID:
        Color="green"
    elif "Y" in sampleID:
        Color = "blue"
    fig = plt.figure(figsize=(15,5))
    axes = fig.add_subplot(111)
    axes.scatter(positions,readBalance,c=Color,alpha=MarkerAlpha,s=MarkerSize,lw=0)
    axes.set_xlim(0,Xlim)
    axes.set_title(sampleID)
    axes.set_xlabel("Chromosomal Coordinate")
    axes.set_ylabel("Read Balance")
    #print(len(positions))
    plt.savefig("%s_ReadBalance_GenomicScatter.svg"%sampleID)
    plt.savefig("%s_ReadBalance_GenomicScatter.png"%sampleID)
    plt.show()
    
    
#This Function will plot the read balance (variant allele / total reads mapped), which is the array output (readBalance) from the above functin
    
def HistReadBalance(readBalance,sampleID):
    if "X" in sampleID:
        Color="green"
    elif "Y" in sampleID:
        Color = "blue"
    fig = plt.figure(figsize=(8,8))
    axes = fig.add_subplot(111)
    axes.set_title(sampleID)
    axes.set_xlabel("Read Balance")
    axes.set_ylabel("Frequency")
    axes.hist(readBalance,bins=50,color=Color)
    plt.savefig("%s_ReadBalance_Hist.svg"%sampleID)
    plt.savefig("%s_ReadBalance_Hist.png"%sampleID)
    plt.show()
    


# Example code:
# NA20845 
#positions,quality,readbalance = ParseVCF("/Users/philliprichmond/Desktop/HACKATHON/HG00096.chrX_Platypus.vcf",200)
#PlotReadBalance(positions,readbalance,"HG00096-X",1,.2,156000000)
#HistReadBalance(readbalance,"HG00096-X")

#positions,quality,readbalance = ParseVCF("/Users/philliprichmond/Desktop/HACKATHON/HG00096.chrY_Platypus.vcf",200)
#PlotReadBalance(positions,readbalance,"HG00096-Y",5,1,60000000)
#HistReadBalance(readbalance,"HG00096-Y")
