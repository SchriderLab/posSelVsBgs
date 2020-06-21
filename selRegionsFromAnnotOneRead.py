import sys
import gzip
import numpy as np
from overlapper import overlap

def readRecRegionsFromWig(recRateFileName):
    if recRateFileName.endswith(".gz"):
        openFunc = gzip.open
    else:
        openFunc = open

    recRates = {}
    with openFunc(recRateFileName, "rt") as recRateFile:
        for line in recRateFile:
            if not line.startswith("#"):
                c, s, e, ratePerBp = line.strip().split()
                s, e = int(s)+1, int(e)
                if not c in recRates:
                    recRates[c] = []
                recRates[c].append((s, e, ratePerBp))
    return recRates

def readRecRegionsInWin(recRates, winC, winStart, winEnd, rRescale=1.0):
    totalRecRateInMorgans = 0
    rregionCoords = []
    for s, e, ratePerBp in recRates[winC]:
        overRange = overlap(s, e, winStart, winEnd)
        if overRange:
            overS, overE = overRange
            overS, overE = overS-winStart, overE-winStart
            rateInRegion = float(ratePerBp)*rRescale*(overE-overS+1)
            rregionCoords.append((overS, overE, rateInRegion))
            totalRecRateInMorgans += rateInRegion
    return rregionCoords, totalRecRateInMorgans

def readSelRegionsFromGtf(geneAnnotFileName, phastConsFileName=None):
    if geneAnnotFileName.endswith(".gz"):
        openFunc = gzip.open
    else:
        openFunc = open

    phastConsCoords = {}
    if phastConsFileName:
        if phastConsFileName.endswith(".gz"):
            openFunc = gzip.open
        else:
            openFunc = open
        if phastConsFileName.rstrip(".gz").endswith(".bed"):
            coordIndices = (0, 1, 2)
        else:
            coordIndices = (1, 2, 3)
        with openFunc(phastConsFileName, "rt") as phastConsFile:
            for line in phastConsFile:
                if not line.startswith("#"):
                    line = line.strip().split()
                    c, s, e = [line[x] for x in coordIndices]
                    if not c in phastConsCoords:
                        phastConsCoords[c] = []
                    phastConsCoords[c].append((int(s)+1, int(e)))

    exonCoords = {}
    with openFunc(geneAnnotFileName, "rt") as geneAnnotFile:
        for line in geneAnnotFile:
            if not line.startswith("#"):
                c, source, annotType, s, e = line.strip().split()[:5]
                if annotType == "exon":
                    if not c in exonCoords:
                        exonCoords[c] = []
                    exonCoords[c].append((s, e))
    return exonCoords, phastConsCoords

def readSelRegionsInWin(exonCoords, winC, winStart, winEnd, cncSelRatio, phastConsCoords=None, codingSelFrac=0.75, nonCodingSelFrac=0.75):
    L=winEnd-winStart+1
    isSel=[0]*L

    if phastConsCoords:
        for s, e in phastConsCoords[winC]:
            overRange = overlap(s, e, winStart, winEnd)
            if overRange:
                overS, overE = overRange
                overS, overE = overS-winStart, overE-winStart
                for pos in range(overS, overE+1):
                    isSel[pos-1] = 1

    for s, e in exonCoords[winC]:
        overRange = overlap(int(s), int(e), winStart, winEnd)
        if overRange:
            overS, overE = overRange
            overS, overE = overS-winStart, overE-winStart
            for pos in range(overS, overE+1):
                isSel[pos-1] = 2

    nregionCoords = []
    sregionCoords = []
    prevState = isSel[0]
    runStart = 1
    totSelRegionSize = 0
    for i in range(1,L):
        if isSel[i] == prevState:
            pass
        else:
            runEnd = i
            if prevState == 2:
                sregionCoords.append((runStart, runEnd, codingSelFrac, 1))
                totSelRegionSize += (runEnd-runStart+1)*codingSelFrac
                nregionCoords.append((runStart, runEnd, 1-codingSelFrac))
            elif prevState == 1:
                sregionCoords.append((runStart, runEnd, nonCodingSelFrac, cncSelRatio))
                totSelRegionSize += (runEnd-runStart+1)*nonCodingSelFrac
                nregionCoords.append((runStart, runEnd, 1-nonCodingSelFrac))
            else:
                nregionCoords.append((runStart, runEnd, 1))
            prevState = isSel[i]
            runStart = i+1
    runEnd = L
    if prevState == 2:
        sregionCoords.append((runStart, runEnd, codingSelFrac, 1))
        totSelRegionSize += (runEnd-runStart+1)*codingSelFrac
        nregionCoords.append((runStart, runEnd, 1-codingSelFrac))
    elif prevState == 1:
        sregionCoords.append((runStart, runEnd, nonCodingSelFrac, cncSelRatio))
        totSelRegionSize += (runEnd-runStart+1)*nonCodingSelFrac
        nregionCoords.append((runStart, runEnd, 1-nonCodingSelFrac))
    else:
        nregionCoords.append((runStart, runEnd, 1))
    return sregionCoords, nregionCoords, totSelRegionSize
