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

discoalCmdFileName, prespecifiedCoords, chrLenFileName, gapFileName, geneAnnotFileName, phastConsFileName, recRateFileName, L, dfeMean, dfeShape, dominance, cncSelRatio, uMean, rMean, codingSelFrac, nonCodingSelFrac, popSizeRescaleFactor, statFileName, smallStatFileName, fvecFileName, currRepNum = sys.argv[1:]

npSeed = random.randint(0,(2**31)-1)
np.random.seed(npSeed)

cncSelRatio = float(cncSelRatio)

statNames = ["pi", "thetaW", "tajD", "thetaH", "fayWuH", "maxFDA", "HapCount", "H1", "H12", "H2/H1", "ZnS", "Omega", "distVar", "distSkew", "distKurt"]
smallStatNames = ["pi", "thetaW", "thetaH", "tajD", "fayWuH", "ZnS", "Omega"]

if gapFileName.lower() in ["none", "false"]:
    gapFileName = None

popSizeRescaleFactor = float(popSizeRescaleFactor)
currRepNum = int(currRepNum)
dfeMean, dfeShape, dominance = float(dfeMean)/2, float(dfeShape), float(dominance)*2
uMean, rMean = float(uMean)/popSizeRescaleFactor, float(rMean)/popSizeRescaleFactor
codingSelFrac, nonCodingSelFrac = float(codingSelFrac), float(nonCodingSelFrac)
L=int(L)
burnTime=10
sys.stderr.write("reading chromosome lengths, selected sites, and recombination map\n")
if prespecifiedCoords.lower() in ["none", "false"]:
    chrLens = selRegionsFromAnnot.readChrLens(chrLenFileName)
    winC, winS, winE = selRegionsFromAnnot.pickRandomWindow(L, chrLens, gapFileName=gapFileName, state=np.random.get_state())
else:
    winC, winSE = prespecifiedCoords.split(":")
    winS, winE = [int(x) for x in winSE.split("-")]
sys.stderr.write("modeling simulated region after %s:%d-%d\n" %(winC, winS, winE))
sregionCoords, nregionCoords, totSelRegionSize = selRegionsFromAnnot.readSelRegionsInWinFromGtf(geneAnnotFileName, winC, winS, winE, cncSelRatio, phastConsFileName=phastConsFileName, codingSelFrac=codingSelFrac, nonCodingSelFrac=nonCodingSelFrac)
recregionCoords, totalRecRateInMorgans = selRegionsFromAnnot.readRecRegionsInWinFromWig(recRateFileName, winC, winS, winE, rRescale=rMean)

sys.stderr.write("reading discoal command\n")
sampleSize, nlist, thetaRange, rhoMean, rhoMax = discoalParseFuncs.getParamBoundsAndPopSizeChangesFromDiscoalCmdFile(discoalCmdFileName, uMean, rMean, L, burnTime)
thetaLow, thetaHigh = thetaRange
theta = np.random.uniform(thetaLow, thetaHigh)
u = theta/(L*nlist[-1]*4)

startTime=time.clock()
sys.stderr.write("initializing regions and population\n")
nregions = [fwdpy11.Region(beg=neutRegionStart, end=neutRegionEnd, weight=neutRegionWeight) for neutRegionStart, neutRegionEnd, neutRegionWeight in nregionCoords]
recregions = [fwdpy11.Region(beg=recRegionStart, end=recRegionEnd, weight=recRegionWeight) for recRegionStart, recRegionEnd, recRegionWeight in recregionCoords]
if dfeMean == 0:
    sregions = [fwdpy11.ConstantS(beg=selRegionStart, end=selRegionEnd, weight=selRegionWeight, s=0.0, h=1.0) for selRegionStart, selRegionEnd, selRegionWeight, selRegionSScale in sregionCoords] #used for debugging
else:
    sregions = [fwdpy11.GammaS(beg=selRegionStart, end=selRegionEnd, shape=dfeShape, mean=dfeMean*selRegionSScale/popSizeRescaleFactor, weight=selRegionWeight, h=1.0) for selRegionStart, selRegionEnd, selRegionWeight, selRegionSScale in sregionCoords]

seed = random.randint(0,(2**31)-1)
rng = fwdpy11.GSLrng(seed)
pop = fwdpy11.SlocusPop(nlist[0])

sys.stderr.write("ready to rock!\n")
sys.stderr.write("currRepNum: %d; u: %g; r (total Morgans): %g; Nanc: %d; N0: %d; numGens: %d; totSelRegionSize: %d; fwdpy seed: %d; np seed: %d\n" %(currRepNum, u, totalRecRateInMorgans, nlist[0], nlist[-1], len(nlist), totSelRegionSize, seed, npSeed))
params=mparams.SlocusParams(nregions=nregions, sregions=sregions, recregions=recregions, demography=np.array(nlist[0:], dtype=np.uint32), rates=(u*(L-totSelRegionSize), u*totSelRegionSize, totalRecRateInMorgans))
fwdpy11wf.evolve(rng, pop, params)
sys.stderr.write("simulation completed; took %f seconds\n" %(time.clock()-startTime))
numReps=1
print("fwdpy %d %d %d etc\n%d %d" %(sampleSize, numReps, L, seed, npSeed))
mutPositions, gameteStrs, allelesFound = miscFwdpyFuncs.sampleMutsFromDiploids(pop, sampleSize, state=np.random.get_state())
sys.stderr.write("sampling completed\n")
snpPositions, gameteStrs = miscFwdpyFuncs.filterMonomorphicSites(mutPositions, gameteStrs, allelesFound)
sys.stderr.write("filtering completed\n")
assert len(gameteStrs) == sampleSize
print("\n//\nsegsites: %d\npositions: %s" %(len(snpPositions), " ".join([str(x/L) for x in snpPositions])))
print("\n".join(gameteStrs))
sys.stderr.write("printing completed\n")
sys.stderr.write("calculating summary stats\n")
windowedStats = summStatFuncs.calculateWindowedStats(snpPositions, gameteStrs, L, statNames)
summStatFuncs.writeStats(windowedStats, statNames, statFileName)
smallWindowedStats = summStatFuncs.calculateWindowedStats(snpPositions, gameteStrs, L, smallStatNames, numSubWins=1100)
summStatFuncs.writeStats(smallWindowedStats, smallStatNames, smallStatFileName)
if not fvecFileName.lower() in ["none", "false"]:
    summStatFuncs.writeFvec(windowedStats, statNames, fvecFileName)
sys.stderr.write("Summary stats written\n")
