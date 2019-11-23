import sys, bisect
import numpy as np
import allel
import scipy.stats
import shicstats
from allel.model.ndarray import SortedIndex
from allel.util import asarray_ndim

def getSnpsOverflowingChr(newPositions, totalPhysLen):
    overflowers = []
    for i in reversed(range(len(newPositions))):
        if newPositions[i] > totalPhysLen:
            overflowers.append(newPositions[i])
    return overflowers

def fillInSnpSlotsWithOverflowers(newPositions, totalPhysLen, overflowers):
    posH = {}
    for pos in newPositions:
        posH[pos]=1
    for i in range(len(overflowers)):
        del newPositions[-1]
    for pos in reversed(range(1, totalPhysLen+1)):
        if not pos in posH:
            bisect.insort_left(newPositions, pos)
            overflowers.pop()
            if len(overflowers) == 0:
                break

def msPositionsToIntegerPositions(positions, totalPhysLen):
    snpNum = 1
    prevPos = -1
    prevIntPos = -1
    newPositions = []
    for position in positions:
        assert position >= prevPos
        origPos = position
        if position == prevPos:
            position += 0.000001
        prevPos = origPos

        intPos = int(totalPhysLen*position)
        if intPos == 0:
            intPos = 1
        if intPos <= prevIntPos:
            intPos = prevIntPos + 1
        prevIntPos = intPos
        newPositions.append(intPos)
    overflowers = getSnpsOverflowingChr(newPositions, totalPhysLen)
    if overflowers:
        fillInSnpSlotsWithOverflowers(newPositions, totalPhysLen, overflowers)
    #sys.stderr.write("number of snps: %d (new) vs. %d (old)\n" %(len(newPositions), len(positions)))
    assert len(newPositions) == len(positions)
    assert all(newPositions[i] <= newPositions[i+1] for i in range(len(newPositions)-1))
    assert newPositions[-1] <= totalPhysLen
    return newPositions

def msRepToHaplotypeArrayIn(samples, positions, totalPhysLen):
    for i in range(len(samples)):
        assert len(samples[i]) == len(positions)

    positions = msPositionsToIntegerPositions(positions, totalPhysLen)

    hapArrayIn = []
    for j in range(len(positions)):
        hapArrayIn.append([])
        for i in range(len(samples)):
            hapArrayIn[j].append(samples[i][j])
    return hapArrayIn, positions

def calcAndAppendStatVal(alleleCounts, snpLocs, statName, subWinStart, subWinEnd, statVals, subWinIndex, hapsInSubWin, unmasked):
    if statName == "tajD":
        statVals[statName].append(allel.stats.diversity.tajima_d(alleleCounts, pos=snpLocs, start=subWinStart, stop=subWinEnd))
    elif statName == "pi":
        statVals[statName].append(allel.stats.diversity.sequence_diversity(snpLocs, alleleCounts, start=subWinStart, stop=subWinEnd, is_accessible=unmasked))
    elif statName == "thetaW":
        statVals[statName].append(allel.stats.diversity.watterson_theta(snpLocs, alleleCounts, start=subWinStart, stop=subWinEnd, is_accessible=unmasked))
    elif statName == "thetaH":
        statVals[statName].append(thetah(snpLocs, alleleCounts, start=subWinStart, stop=subWinEnd, is_accessible=unmasked))
    elif statName == "fayWuH":
        statVals[statName].append(statVals["thetaH"][subWinIndex]-statVals["pi"][subWinIndex])
    elif statName == "HapCount":
        statVals[statName].append(len(hapsInSubWin.distinct()))
    elif statName == "maxFDA":
        statVals[statName].append(maxFDA(snpLocs, alleleCounts, start=subWinStart, stop=subWinEnd, is_accessible=unmasked))
    elif statName == "H1":
        h1, h12, h123, h21 = allel.stats.selection.garud_h(hapsInSubWin)
        statVals["H1"].append(h1)
        if "H12" in statVals:
            statVals["H12"].append(h12)
        if "H123" in statVals:
            statVals["H123"].append(h123)
        if "H2/H1" in statVals:
            statVals["H2/H1"].append(h21)
    elif statName == "ZnS":
        r2Matrix = shicstats.computeR2Matrix(hapsInSubWin)
        statVals["ZnS"].append(shicstats.ZnS(r2Matrix)[0])
        statVals["Omega"].append(shicstats.omega(r2Matrix)[0])
    elif statName == "RH":
        rMatrixFlat = allel.stats.ld.rogers_huff_r(hapsInSubWin.to_genotypes(ploidy=2).to_n_alt())
        rhAvg = rMatrixFlat.mean()
        statVals["RH"].append(rhAvg)
        r2Matrix = squareform(rMatrixFlat ** 2)
        statVals["Omega"].append(shicstats.omega(r2Matrix)[0])
    elif statName == "distVar":
        dists = shicstats.pairwiseDiffs(hapsInSubWin)/float(unmasked[subWinStart-1:subWinEnd].count(True))
        statVals["distVar"].append(np.var(dists, ddof=1))
        statVals["distSkew"].append(scipy.stats.skew(dists))
        statVals["distKurt"].append(scipy.stats.kurtosis(dists))
    elif statName in ["H12", "H123", "H2/H1", "Omega", "distVar", "distSkew", "distKurt"]:
        assert len(statVals[statName]) == subWinIndex+1

def appendStatValsForMonomorphic(statName, statVals, subWinIndex):
    if statName == "tajD":
        statVals[statName].append(0.0)
    elif statName == "pi":
        statVals[statName].append(0.0)
    elif statName == "thetaW":
        statVals[statName].append(0.0)
    elif statName == "thetaH":
        statVals[statName].append(0.0)
    elif statName == "fayWuH":
        statVals[statName].append(0.0)
    elif statName == "maxFDA":
        statVals[statName].append(0.0)
    elif statName == "nDiplos":
        statVals[statName].append(1)
    elif statName == "HapCount":
        statVals[statName].append(1)
    elif statName in ["H1"]:
        statVals["H1"].append(1.0)
        if "H12" in statVals:
            statVals["H12"].append(1.0)
        if "H123" in statVals:
            statVals["H123"].append(1.0)
        if "H2/H1" in statVals:
            statVals["H2/H1"].append(0.0)
    elif statName == "ZnS":
        statVals["ZnS"].append(0.0)
        statVals["Omega"].append(0.0)
    elif statName == "RH":
        statVals["RH"].append(0.0)
        statVals["Omega"].append(0.0)
    elif statName == "iHSMean":
        statVals["iHSMean"].append(0.0)
    elif statName == "nSLMean":
        statVals["nSLMean"].append(0.0)
    elif statName == "iHSMax":
        statVals["iHSMax"].append(0.0)
    elif statName == "nSLMax":
        statVals["nSLMax"].append(0.0)
    elif statName in ["H12", "H123", "H2/H1", "Omega"]:
        #print(statName, statVals[statName], subWinIndex+1)
        assert len(statVals[statName]) == subWinIndex+1
    else:
        statVals[statName].append(0.0)

def getSnpIndicesInSubWins(subWinBounds, snpLocs):
    snpIndicesInSubWins = []
    for subWinIndex in range(len(subWinBounds)):
        snpIndicesInSubWins.append([])

    subWinIndex = 0
    for i in range(len(snpLocs)):
        while not (snpLocs[i] >= subWinBounds[subWinIndex][0] and snpLocs[i] <= subWinBounds[subWinIndex][1]):
            subWinIndex += 1
        snpIndicesInSubWins[subWinIndex].append(i)
    return snpIndicesInSubWins

def getSubWinBounds(subWinLen, totalPhysLen): # get inclusive subwin bounds
    subWinStart = 1
    subWinEnd = subWinStart + subWinLen - 1
    subWinBounds = [(subWinStart, subWinEnd)]
    numSubWins = totalPhysLen//subWinLen
    for i in range(1, numSubWins-1):
        subWinStart += subWinLen
        subWinEnd += subWinLen
        subWinBounds.append((subWinStart, subWinEnd))
    subWinStart += subWinLen
    # if our subwindows are 1 bp too short due to rounding error, the last window picks up all of the slack
    subWinEnd = totalPhysLen
    subWinBounds.append((subWinStart, subWinEnd))
    return subWinBounds

def calculateWindowedStats(continuousSnpPositions, gameteStrs, totalPhysLen, statNames, numSubWins=11):
    continuousSnpPositions = [x/float(totalPhysLen) for x in continuousSnpPositions]
    hapArrayIn, discretePositions = msRepToHaplotypeArrayIn(gameteStrs, continuousSnpPositions, totalPhysLen)
    haps = allel.HaplotypeArray(hapArrayIn, dtype='i1')
    if haps.shape[1] % 2 == 1:
        haps = haps[:,:-1]
    genos = haps.to_genotypes(ploidy=2)

    unmasked = [True]*totalPhysLen
    statVals = {}
    for statName in statNames:
        statVals[statName] = []

    subWinBounds = getSubWinBounds(int(totalPhysLen/numSubWins), totalPhysLen)
    ac = genos.count_alleles()
    unmaskedSnpIndices = [i for i in range(len(discretePositions)) if unmasked[discretePositions[i]-1]]
    positionArrayUnmaskedOnly = [discretePositions[i] for i in unmaskedSnpIndices]
    alleleCountsUnmaskedOnly = allel.AlleleCountsArray(np.array([ac[i] for i in unmaskedSnpIndices]))
    snpIndicesInSubWins = getSnpIndicesInSubWins(subWinBounds, discretePositions)

    for subWinIndex in range(numSubWins):
        subWinStart, subWinEnd = subWinBounds[subWinIndex]
        snpIndicesInSubWinUnmasked = [x for x in snpIndicesInSubWins[subWinIndex] if unmasked[discretePositions[x]-1]]
        if len(snpIndicesInSubWinUnmasked) > 0:
            hapsInSubWin = haps.subset(sel0=snpIndicesInSubWinUnmasked)
            genosInSubWin = genos.subset(sel0=snpIndicesInSubWinUnmasked)
            for statName in statNames:
                calcAndAppendStatVal(alleleCountsUnmaskedOnly, positionArrayUnmaskedOnly, statName, subWinStart, subWinEnd, statVals, subWinIndex, hapsInSubWin, unmasked)
        else:
            for statName in statNames:
                appendStatValsForMonomorphic(statName, statVals, subWinIndex)
    return statVals

def getSubWinCountFromWindowedStats(windowedStats):
    numSubWinsH = {}
    for statName in windowedStats:
        numSubWinsH[len(windowedStats[statName])] = 1
    assert len(numSubWinsH) == 1
    return list(numSubWinsH.keys())[0]

def writeStats(windowedStats, statNames, statFileName):
    with open(statFileName, "wt") as statFile:
        outLine = ["subWinIndex"]
        for statName in statNames:
            outLine.append(statName)
        statFile.write("\t".join(outLine) + "\n")
        for subWinIndex in range(getSubWinCountFromWindowedStats(windowedStats)):
            outLine = [str(subWinIndex)]
            for statName in statNames:
                outLine.append(str(windowedStats[statName][subWinIndex]))
            statFile.write("\t".join(outLine) + "\n")

def normalizeFeatureVec(statVec):
    minVal = min(statVec)
    if minVal < 0:
        statVec = [x-minVal for x in statVec]
    normStatVec = []
    statSum = float(sum(statVec))
    if statSum == 0:
        normStatVec = [1.0/len(statVec)]*len(statVec)
    else:
        for k in range(len(statVec)):
            normStatVec.append(statVec[k]/statSum)
    return normStatVec

def writeFvec(windowedStats, statNames, fvecFileName):
    header = []
    for statName in statNames:
        for i in range(getSubWinCountFromWindowedStats(windowedStats)):
            header.append("%s_win%d" %(statName, i))
    header = "\t".join(header)

    with open(fvecFileName, "wt") as fvecFile:
        fvecFile.write(header + "\n")
        outVec = []
        for statName in statNames:
            outVec += normalizeFeatureVec(windowedStats[statName])
        fvecFile.write("\t".join([str(x) for x in outVec]) + "\n")

def maxFDA(pos, ac, start=None, stop=None, is_accessible=None):
    # check inputs
    if not isinstance(pos, SortedIndex):
        pos = SortedIndex(pos, copy=False)
    ac = asarray_ndim(ac, 2)
    is_accessible = asarray_ndim(is_accessible, 1, allow_none=True)

    # deal with subregion
    if start is not None or stop is not None:
        loc = pos.locate_range(start, stop)
        pos = pos[loc]
        ac = ac[loc]
    if start is None:
        start = pos[0]
    if stop is None:
        stop = pos[-1]

    # calculate values of the stat
    dafs = []
    for i in range(len(ac)):
        p1 = ac[i, 1]
        n = p1+ac[i, 0]
        dafs.append(p1/float(n))
    return max(dafs)

'''
WARNING: this code assumes that the second column of ac gives the derived alleles;
please ensure that this is the case (and that you are using polarized data)!!
'''
def thetah(pos, ac, start=None, stop=None, is_accessible=None):
    # check inputs
    if not isinstance(pos, SortedIndex):
        pos = SortedIndex(pos, copy=False)
    ac = asarray_ndim(ac, 2)
    is_accessible = asarray_ndim(is_accessible, 1, allow_none=True)

    # deal with subregion
    if start is not None or stop is not None:
        loc = pos.locate_range(start, stop)
        pos = pos[loc]
        ac = ac[loc]
    if start is None:
        start = pos[0]
    if stop is None:
        stop = pos[-1]

    # calculate values of the stat
    h = 0
    for i in range(len(ac)):
        p1 = ac[i, 1]
        n = p1+ac[i, 0]
        if n > 1:
            h += (p1*p1)/(n*(n-1.0))
    h *= 2

    # calculate value per base
    if is_accessible is None:
        n_bases = stop - start + 1
    else:
        n_bases = np.count_nonzero(is_accessible[start-1:stop])

    h = h / n_bases
    return h
