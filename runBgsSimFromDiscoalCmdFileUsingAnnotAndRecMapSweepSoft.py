import sys, random, time, math, pickle
import fwdpy11
import fwdpy11.wright_fisher as fwdpy11wf
import fwdpy11.model_params as mparams
import numpy as np
import scipy.stats
import selRegionsFromAnnot
import discoalParseFuncs
import miscFwdpyFuncs
import summStatFuncs

discoalCmdFileName, prespecifiedCoords, chrLenFileName, gapFileName, geneAnnotFileName, phastConsFileName, recRateFileName, L, dfeMean, dfeShape, dominance, cncSelRatio, uMean, rMean, codingSelFrac, nonCodingSelFrac, popSizeRescaleFactor, statFileName, smallStatFileName, fvecFileName, currRepNum, minTimeSinceSweepSel, maxTimeSinceSweepSel, maxTimeSinceFixation, minSelCoeff, maxSelCoeff, minInitialSelFreq, maxInitialSelFreq = sys.argv[1:]
if phastConsFileName == "None":
    phastConsFileName = None

pySeed = random.randint(0, 2**32 - 1)
random.seed(pySeed)
npSeed = random.randint(0, 2**32 - 1)
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
minTimeSinceSweepSel, maxTimeSinceSweepSel = math.floor(float(minTimeSinceSweepSel)*popSizeRescaleFactor), math.ceil(float(maxTimeSinceSweepSel)*popSizeRescaleFactor)
maxTimeSinceFixation = float(maxTimeSinceFixation)*popSizeRescaleFactor
minInitialSelFreq, maxInitialSelFreq = float(minInitialSelFreq), float(maxInitialSelFreq)
minSelCoeff, maxSelCoeff = float(minSelCoeff)/popSizeRescaleFactor, float(maxSelCoeff)/popSizeRescaleFactor
lossCheckTime = math.floor(100*popSizeRescaleFactor)# time period after which we check to see if our sweeping mut was lost and quickly start over
burnTime=40

sys.stderr.write("reading chromosome lengths, selected sites, and recombination map\n")
if prespecifiedCoords.lower() in ["none", "false"]:
    chrLens = selRegionsFromAnnot.readChrLens(chrLenFileName)
    winC, winS, winE, recregionCoords, totalRecRateInMorgans = selRegionsFromAnnot.pickRandomWindow(L, chrLens, recRateFileName, rMean, gapFileName=gapFileName, state=np.random.get_state())
else:
    winC, winSE = prespecifiedCoords.split(":")
    winS, winE = [int(x) for x in winSE.split("-")]
    recregionCoords, totalRecRateInMorgans = selRegionsFromAnnot.readRecRegionsInWinFromWig(recRateFileName, winC, winS, winE, rRescale=rMean)
sys.stderr.write("modeling simulated region after %s:%d-%d\n" %(winC, winS, winE))
sregionCoords, nregionCoords, totSelRegionSize = selRegionsFromAnnot.readSelRegionsInWinFromGtf(geneAnnotFileName, winC, winS, winE, cncSelRatio, phastConsFileName=phastConsFileName, codingSelFrac=codingSelFrac, nonCodingSelFrac=nonCodingSelFrac)

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
    sregions = [fwdpy11.ConstantS(beg=selRegionStart, end=selRegionEnd, weight=selRegionWeight, s=0.0, h=dominance) for selRegionStart, selRegionEnd, selRegionWeight, selRegionSScale in sregionCoords] #used for debugging
else:
    sregions = [fwdpy11.GammaS(beg=selRegionStart, end=selRegionEnd, shape=dfeShape, mean=dfeMean*selRegionSScale/popSizeRescaleFactor, weight=selRegionWeight, h=dominance) for selRegionStart, selRegionEnd, selRegionWeight, selRegionSScale in sregionCoords]

seed = random.randint(0,(2**31)-1)
rng = fwdpy11.GSLrng(seed)
pop = fwdpy11.SlocusPop(nlist[0])
if maxTimeSinceSweepSel >= len(nlist):
    sys.stderr.write("WARNING: maximum time since the sweeping mutation occurred (%d gen ago) is greater than duration of the simulation (%d gen). Truncating to %d\n" %(maxTimeSinceSweepSel, len(nlist), len(nlist)-1))
    maxTimeSinceSweepSel = len(nlist)-1

sys.stderr.write("drawing sweepSelCoeff from [%s, %s]\n" %(minSelCoeff, maxSelCoeff))
sweepSelCoeff = np.random.uniform(minSelCoeff, maxSelCoeff)

sys.stderr.write("ready to rock!\n")
sys.stderr.write("currRepNum: %d; u: %g; r (total Morgans): %g; Nanc: %d; N0: %d; numGens: %d; totSelRegionSize: %d; fwdpy seed: %d; np seed: %d; python random seed: %d\n" %(currRepNum, u, totalRecRateInMorgans, nlist[0], nlist[-1], len(nlist), totSelRegionSize, seed, npSeed, pySeed))
initialPhaseEnd = maxTimeSinceSweepSel+1
sys.stderr.write("running up until maxTimeSinceSweepSel+1 generations ago (%d gen forward, up until %d gen ago)\n" %(len(nlist)-initialPhaseEnd, maxTimeSinceSweepSel))
params=mparams.SlocusParams(nregions=nregions, sregions=sregions, recregions=recregions, demography=np.array(nlist[:-initialPhaseEnd], dtype=np.uint32), rates=(u*(L-totSelRegionSize), u*totSelRegionSize, totalRecRateInMorgans))
fwdpy11wf.evolve(rng, pop, params)
checkpoint = pickle.dumps(pop)

desiredSweepPos = miscFwdpyFuncs.findSweepPosNearTarget(L/2, pop)
acceptRep=False
numAttempts = 0
acceptProb=1.0
while not acceptRep:
    timeSinceSweepSel = np.random.randint(minTimeSinceSweepSel, maxTimeSinceSweepSel+1)
    sys.stderr.write("drawing timeSinceSweepSel from [%s, %s]\n" %(minTimeSinceSweepSel, maxTimeSinceSweepSel))
    assert lossCheckTime < timeSinceSweepSel
    sys.stderr.write("sweep selection phase onset time: %d generations ago; selection coefficient: %f; desired sweep mut position: %g\n" %(timeSinceSweepSel, sweepSelCoeff, desiredSweepPos))
    sys.stderr.write("simulating up to our sweep mut time (evolving from %d to %d)\n" %(len(nlist)-initialPhaseEnd, len(nlist)-timeSinceSweepSel))
    params=mparams.SlocusParams(nregions=nregions, sregions=sregions, recregions=recregions, demography=np.array(nlist[-initialPhaseEnd:-timeSinceSweepSel], dtype=np.uint32), rates=(u*L, u*totSelRegionSize, totalRecRateInMorgans))
    fwdpy11wf.evolve(rng, pop, params)
    sweepPos, initialSweepSelFreq = miscFwdpyFuncs.findSegsiteNearTarget(desiredSweepPos, pop, minInitialSelFreq, maxInitialSelFreq)
    miscFwdpyFuncs.changeSelCoeffOfMutationAtPos(sweepPos, pop, sweepSelCoeff, tol=1e-8)
    sys.stderr.write("chose a mutation at position %d and frequency %g to be selected\n" %(sweepPos, initialSweepSelFreq))

    newSeed = random.randint(0,(2**31)-1)
    rng = fwdpy11.GSLrng(newSeed)
    sys.stderr.write("preparing to attempt sweep phase of simulation with fwdpy seed of %d\n" %(newSeed))
    sys.stderr.write("simulating for %d gen then checking to see if it was lost (evolving from %d to %d)\n" %(lossCheckTime, len(nlist)-timeSinceSweepSel, len(nlist)-(timeSinceSweepSel-lossCheckTime)))
    params=mparams.SlocusParams(nregions=nregions, sregions=sregions, recregions=recregions, demography=np.array(nlist[-timeSinceSweepSel:-(timeSinceSweepSel-lossCheckTime)], dtype=np.uint32), rates=(u*L, u*totSelRegionSize, totalRecRateInMorgans))
    fwdpy11wf.evolve(rng, pop, params)

    if miscFwdpyFuncs.sweepMutLost(pop, sweepPos):
        sys.stderr.write("sweep lost before our check after %d gen on attempt %d\n" %(lossCheckTime, numAttempts))
    else:
        sys.stderr.write("sweep mut still segregating; simulating for the remaining %d gen (evolving from %d to %d)\n" %(timeSinceSweepSel-lossCheckTime, len(nlist)-(timeSinceSweepSel-lossCheckTime), len(nlist)))
        params=mparams.SlocusParams(nregions=nregions, sregions=sregions, recregions=recregions, demography=np.array(nlist[-(timeSinceSweepSel-lossCheckTime):], dtype=np.uint32), rates=(u*L, u*totSelRegionSize, totalRecRateInMorgans))
        fwdpy11wf.evolve(rng, pop, params)

        fixationTime = miscFwdpyFuncs.sweepFixationTimeIfAfterDate(pop, sweepPos, len(nlist)-maxTimeSinceFixation)
        if fixationTime != None:
            if np.random.random() <= acceptProb:
                sys.stderr.write("sweep fixed and accepted on attempt %d\n" %(numAttempts))
                acceptRep = True
            else:
                sys.stderr.write("sweep fixed but rejected due to imporance sampling on attempt %d\n" %(numAttempts))
        else:
            sys.stderr.write("sweep didn't make it to fixation on attempt %d\n" %(numAttempts))
    numAttempts += 1

    if not acceptRep:
        pop = pickle.loads(checkpoint)

sys.stderr.write("simulation completed; took %f seconds and %d attempts\n" %(time.clock()-startTime, numAttempts))
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
sys.stderr.write("Final sweep characteristics (re-scaled back to orig): (sweepPos\\tsweepSelCoeff\\tinitialSweepSelFreq\\ttimeSinceSweepSel\\ttimeSinceFixation)\n")
sys.stderr.write("%s\t%s\t%s\t%s\t%s\n" %(sweepPos, sweepSelCoeff*popSizeRescaleFactor, initialSweepSelFreq, timeSinceSweepSel/popSizeRescaleFactor, (len(nlist)-fixationTime)/popSizeRescaleFactor))
