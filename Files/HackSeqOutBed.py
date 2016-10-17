import pybedtools
import pandas as pd
import numpy as np



def makeRegionLists(depthAndMapqDf, mapqCutoff, depthMin, depthMax):
    '''
    (pandas.core.frame.DataFrame, float, float, float) -> (list, list)
    return two lists of regions (keepList, excludeList) based on cutoffs for depth and mapq
    '''
    good = (depthAndMapqDf.mapq > mapqCutoff) & (depthAndMapqDf.depth > depthMin) & (depthAndMapqDf.depth < depthMax)
    dfGood = depthAndMapqDf[good]
    dfBad = depthAndMapqDf[~good]
    goodList = dfGood.ix[:, 'chrom':'end'].values.tolist()
    badList = dfBad.ix[:, 'chrom':'end'].values.tolist()
    return (goodList, badList)
    



def outputBed(regionList, outBed):
    '''
    (list, list, str) -> bedtoolsObject
    Take two sorted lists.  Each list is a list of tuples (chrom[str], start[int], end[int])
    Return a pybedtools object and output a bed file.
    '''
#    mqBed = pybedtools.BedTool(mqList)
#    dpBed = pybedtools.BedTool(depthList)
#    inter = pybedtools.bedtool.BedTool.intersect(mqBed, dpBed)
    regionBed = pybedtools.BedTool(regionList)
    merge = pybedtools.bedtool.BedTool.merge(regionBed)
    with open(outBed, 'w') as output:
        output.write(str(merge))
    return merge


