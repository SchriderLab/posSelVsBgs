import sys, random, time
import fwdpy11
import fwdpy11.wright_fisher as fwdpy11wf
import fwdpy11.model_params as mparams
import numpy as np
import scipy.stats
import selRegionsFromAnnot
import discoalParseFuncs
import miscFwdpyFuncs
import summStatFuncs

discoalCmdFileName, prespecifiedCoords, winL, stepSize, chrLenFileName, gapFileName, geneAnnotFileName, phastConsFileName, recRateFileName, totalL, dfeMean, dfeShape, dominance, cncSelRatio, uMean, rMean, codingSelFrac, nonCodingSelFrac, popSizeRescaleFactor, statFilePrefix, fvecFilePrefix, currRepNum = sys.argv[1:]
if phastConsFileName.lower() == "none":
    phastConsFileName = None

npSeed = random.randint(0,(2**31)-1)
np.random.seed(npSeed)

cncSelRatio = float(cncSelRatio)

statNames = ["pi", "thetaW", "tajD", "thetaH", "fayWuH", "maxFDA", "HapCount", "H1", "H12", "H2/H1", "ZnS", "Omega", "distVar", "distSkew", "distKurt"]

if gapFileName.lower() in ["none", "false"]:
    gapFileName = None

popSizeRescaleFactor = float(popSizeRescaleFactor)
currRepNum = int(currRepNum)
dfeMean, dfeShape, dominance = float(dfeMean)/2, float(dfeShape), float(dominance)*2
uMean, rMean = float(uMean)/popSizeRescaleFactor, float(rMean)/popSizeRescaleFactor
codingSelFrac, nonCodingSelFrac = float(codingSelFrac), float(nonCodingSelFrac)
totalL=int(totalL)
winL=int(winL)
stepSize=int(stepSize)
burnTime=10
sys.stderr.write("reading chromosome lengths, selected sites, and recombination map\n")
if prespecifiedCoords.lower() in ["none", "false"]:
    chrLens = selRegionsFromAnnot.readChrLens(chrLenFileName)
    winC, winS, winE, recregionCoords, totalRecRateInMorgans = selRegionsFromAnnot.pickRandomWindow(totalL, chrLens, recRateFileName, rMean, gapFileName=gapFileName, state=np.random.get_state())
else:
    winC, winSE = prespecifiedCoords.split(":")
    winS, winE = [int(x) for x in winSE.split("-")]
    recregionCoords, totalRecRateInMorgans = selRegionsFromAnnot.readRecRegionsInWinFromWig(recRateFileName, winC, winS, winE, rRescale=rMean)
sys.stderr.write("modeling simulated region after %s:%d-%d\n" %(winC, winS, winE))
sregionCoords, nregionCoords, totSelRegionSize = selRegionsFromAnnot.readSelRegionsInWinFromGtf(geneAnnotFileName, winC, winS, winE, cncSelRatio, phastConsFileName=phastConsFileName, codingSelFrac=codingSelFrac, nonCodingSelFrac=nonCodingSelFrac)
sys.stderr.write("\n")

sys.stderr.write("will sample the following windows:\n")
windowsToSample = []
firstWinStart = winS
lastWinStart = winE-winL
for sampWinStart in range(firstWinStart, lastWinStart+1, stepSize):
    sampWinStart -= winS
    sampWinEnd = sampWinStart + winL
    windowsToSample.append((sampWinStart, sampWinEnd))
    sys.stderr.write("\tfrom position %s-%s along the simulated chromosome (zero-based half-open)\n" %(sampWinStart, sampWinEnd))
sys.stderr.write("\n")

sys.stderr.write("reading discoal command\n")
sampleSize, nlist, thetaRange, rhoMean, rhoMax = discoalParseFuncs.getParamBoundsAndPopSizeChangesFromDiscoalCmdFile(discoalCmdFileName, uMean, rMean, winL, burnTime)
thetaLow, thetaHigh = thetaRange
theta = np.random.uniform(thetaLow, thetaHigh)
u = theta/(winL*nlist[-1]*4)
#u = uMean

startTime=time.clock()
sys.stderr.write("initializing regions and population\n")
nregions = [fwdpy11.Region(beg=neutRegionStart, end=neutRegionEnd, weight=neutRegionWeight) for neutRegionStart, neutRegionEnd, neutRegionWeight in nregionCoords]
recregions = [fwdpy11.Region(beg=recRegionStart, end=recRegionEnd, weight=recRegionWeight) for recRegionStart, recRegionEnd, recRegionWeight in recregionCoords]
if dfeMean == 0:
    sregions = [fwdpy11.ConstantS(beg=selRegionStart, end=selRegionEnd, weight=selRegionWeight, s=0.0, h=dominance) for selRegionStart, selRegionEnd, selRegionWeight, selRegionSScale in sregionCoords] #used for debugging
else:
    sregions = [fwdpy11.GammaS(beg=selRegionStart, end=selRegionEnd, shape=dfeShape, mean=dfeMean*selRegionSScale/popSizeRescaleFactor, weight=selRegionWeight, h=dominance) for selRegionStart, selRegionEnd, selRegionWeight, selRegionSScale in sregionCoords]

seed = random.randint(0,(2**31)-1)
rng = fwdpy11.GSLrng(seed)
pop = fwdpy11.SlocusPop(nlist[0])

sys.stderr.write("\n")
sys.stderr.write("ready to rock!\n")
sys.stderr.write("currRepNum: %d; u: %g; r (total Morgans): %g; Nanc: %d; N0: %d; numGens: %d; totSelRegionSize: %d; fwdpy seed: %d; np seed: %d\n" %(currRepNum, u, totalRecRateInMorgans, nlist[0], nlist[-1], len(nlist), totSelRegionSize, seed, npSeed))
params=mparams.SlocusParams(nregions=nregions, sregions=sregions, recregions=recregions, demography=np.array(nlist[0:], dtype=np.uint32), rates=(u*(totalL-totSelRegionSize), u*totSelRegionSize, totalRecRateInMorgans))
fwdpy11wf.evolve(rng, pop, params)
sys.stderr.write("simulation completed; took %f seconds\n" %(time.clock()-startTime))
sys.stderr.write("\n")

sys.stderr.write("calculating and writing stats for sliding windows:\n")
for i in range(len(windowsToSample)):
    currWinStart, currWinEnd = windowsToSample[i]
    numReps=1
    print("fwdpy %d %d %d etc\n%d %d" %(sampleSize, numReps, winL, seed, npSeed))
    mutPositions, gameteStrs, allelesFound = miscFwdpyFuncs.sampleMutsFromDiploids(pop, sampleSize, state=np.random.get_state(), reportStart=currWinStart, reportEnd=currWinEnd)
    snpPositions, gameteStrs = miscFwdpyFuncs.filterMonomorphicSites(mutPositions, gameteStrs, allelesFound)
    assert len(gameteStrs) == sampleSize
    print("\n//\nsegsites: %d\npositions: %s" %(len(snpPositions), " ".join([str(x/winL) for x in snpPositions])))
    print("\n".join(gameteStrs))

    statFileName=statFilePrefix + ".win.%d.stats" %(i)
    fvecFileName=fvecFilePrefix + ".win.%d.stats" %(i)

    windowedStats = summStatFuncs.calculateWindowedStats(snpPositions, gameteStrs, winL, statNames)
    summStatFuncs.writeStats(windowedStats, statNames, statFileName)

    if not fvecFileName.lower() in ["none", "false"]:
        summStatFuncs.writeFvec(windowedStats, statNames, fvecFileName)
    sys.stderr.write("Summary stats/feature vectors calculated and written for window %d (%d-%d)\n" %(i, currWinStart, currWinEnd))
