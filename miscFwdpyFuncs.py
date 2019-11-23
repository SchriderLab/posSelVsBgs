import sys, math
import numpy as np
import fwdpy11
import fwdpy11.sampling

def filterMonomorphicSites(snps, gameteStrs, allelesFound):
    newSnps = []
    newGameteStrs = [""]*len(gameteStrs)
    for i in range(len(snps)):
        if len(allelesFound[i]) > 1:
            for j in range(len(gameteStrs)):
                newGameteStrs[j] += gameteStrs[j][i]
            newSnps.append(snps[i])
    return newSnps, newGameteStrs

def getMutPositions(pop):
    return [mut.pos for mut in pop.mutations]

def selectSweepPos(desiredPos, pop, step=0.0001):
    mutPositions = getMutPositions(pop)

    if not desiredPos in mutPositions:
        return desiredPos
    else:
        leftDesiredPos = desiredPos-step
        rightDesiredPos = desiredPos+step
        while 1:
            if not leftDesiredPos in mutPositions:
                return leftDesiredPos
            elif not rightDesiredPos in mutPositions:
                return rightDesiredPos
            leftDesiredPos -= step
            rightDesiredPos += step

def sampleMutsFromDiploids(pop, sampleSize, state=None):
    if state:
        state = np.random.set_state(state)
    sample = pop.sample(individuals = np.random.choice(pop.N, int(math.ceil(sampleSize/2.0)), replace=False), separate=False)

    mutPositions = []
    allelesFound = {}
    gameteStrs = [""]*sampleSize
    for j in range(len(sample)): #iterating over segregating sites
        pos, genotypeStr = sample[j]
        mutPositions.append(float(pos))
        for i in range(sampleSize): #iterating over haplotypes (we may skip the last one if odd number requested)
            gameteStrs[i] += genotypeStr[i]
            if not j in allelesFound:
                allelesFound[j] = {}
            allelesFound[j][genotypeStr[i]] = 1

    lens = {}
    for gameteStr in gameteStrs:
        lens[len(gameteStr)] = 1

    assert len(lens) == 1 and len(gameteStr) == len(sample)
    return mutPositions, gameteStrs, allelesFound

def isSelectedMut(mutation):
    return mutation.s > 0

def getSelectedFixations(pop):
    assert len(pop.fixations) == len(pop.fixation_times)
    selFixations, selFixationTimes = [], []
    for i in range(len(pop.fixations)):
        if isSelectedMut(pop.fixations[i]):
            selFixations.append(pop.fixations[i])
            selFixationTimes.append(pop.fixation_times[i])
    return selFixations, selFixationTimes

def calcSweepMutRejectionProb(sweepMutTimeGA, sweepMutTimeGALow, sweepMutTimeGAHigh, nlist):
    nMut = nlist[-sweepMutTimeGA]
    nMax = max(nlist[-int(round(sweepMutTimeGAHigh)):-int(round(sweepMutTimeGALow))])
    return nMut/float(nMax)
