import sys, math
import numpy as np
import fwdpy11
import fwdpy11.sampling
import fwdpy11.util

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

def findSegsiteNearTarget(desiredPos, pop, minFreq, maxFreq):
    np.arange(pop.N)
    sampleN, sampleS = pop.sample(individuals = np.arange(pop.N), separate=True)

    dists = []
    for j in range(len(sampleN)): #iterating over neutral polymorphisms
        pos, genotypeStr = sampleN[j]
        freq =  sum([int(x) for x in genotypeStr])/len(genotypeStr)
        if freq >= minFreq and freq <= maxFreq:
            dist = abs(pos-desiredPos)
            dists.append((dist, pos, freq))
    dists.sort()

    if (len(dists)) == 0:
        sys.exit("Error: failed to find any neutral polymorphisms with freq within [%g, %g]. Let me die. ARRGGGHHHH!!!!\n" %(minFreq, maxFreq))

    smallestDist, nearestMutPos, nearestMutFreq = dists[0]
    return nearestMutPos, nearestMutFreq

def changeSelCoeffOfMutationAtPos(desiredPos, pop, selCoeff, tol=1e-8):
    for i in range(len(pop.mutations)):
        if abs(pop.mutations[i].pos-desiredPos) < tol:
            fwdpy11.util.change_effect_size(pop, i, selCoeff, 2.0)

def findSweepPosNearTarget(desiredPos, pop, step=0.0001):
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

def addMutAtPosition(pop, rng, sweepPos, sweepSelCoeff, dominance):
    mutTuple = (sweepPos, sweepSelCoeff, dominance)
    fwdpy11.util.add_mutation(rng, pop, 1, mutTuple, 0)

def sampleMutsFromDiploids(pop, sampleSize, state=None, reportStart=None, reportEnd=None):
    if state:
        np.random.set_state(state)
    sample = pop.sample(individuals = np.random.choice(pop.N, math.ceil(sampleSize/2.0), replace=False), separate=False)

    mutPositions = []
    allelesFound = {}
    gameteStrs = [""]*sampleSize
    for j in range(len(sample)): #iterating over segregating sites
        pos, genotypeStr = sample[j]
        if reportStart == None or reportEnd == None:
            mutPositions.append(float(pos))
        elif pos >= reportStart and pos < reportEnd:
            mutPositions.append(float(pos)-reportStart)
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

def getSweepFreq(pop, mutPos, tol=1e-8):
    """
    Gets the frequency of a segregating beneficial mutation if present in pop.
    Our population should have at most one such mutation.
    """
    np.arange(pop.N)
    sampleN, sampleS = pop.sample(individuals = np.arange(pop.N), separate=True)

    sweepFreqs = []
    for j in range(len(sampleS)): #iterating over selected mutations
        pos, genotypeStr = sampleS[j]
        if abs(pos-mutPos) < tol:
            sweepFreqs.append(sum([int(x) for x in genotypeStr])/len(genotypeStr))

    if len(sweepFreqs) == 0:
        return 0
    assert len(sweepFreqs) == 1
    return sweepFreqs[0]

def sweepMutLost(pop, sweepPos):
    selFixations, selFixationTimes = getSelectedFixations(pop)
    assert len(selFixations) <= 1
    if len(selFixations) == 0:
        return getSweepFreq(pop, sweepPos) == 0
    else:
        return False

def sweepFixationTimeIfAfterDate(pop, sweepPos, minFixationDate, verbose=True):
    selFixations, selFixationTimes = getSelectedFixations(pop)
    assert len(selFixations) <= 1
    if len(selFixations) == 0:
        if verbose:
            freq = getSweepFreq(pop, sweepPos)
            if freq > 0:
                sys.stderr.write("sweeping mutation didn't fix (freq=%g). Bummer!\n" %(freq))
            else:
                sys.stderr.write("sweeping mutation was lost! Bummer!\n")
        return None
    else:
        selFixationTime = selFixationTimes[0]
        if selFixationTime >= minFixationDate:
            if verbose:
                sys.stderr.write("sweep mut fixed within desired range (fixed at gen %g; min cutoff is %d)\n" %(selFixationTime, minFixationDate))
            return selFixationTime
        else:
            if verbose:
                sys.stderr.write("sweep mut fixed outside of range (fixed at gen %g; min cutoff is %d)\n" %(selFixationTime, minFixationDate))
            return None

def calcSweepMutAcceptanceProb(sweepMutTimeGA, sweepMutTimeGALow, sweepMutTimeGAHigh, nlist, verbose=True):
    nMut = nlist[-sweepMutTimeGA]
    nMax = max(nlist[-int(round(sweepMutTimeGAHigh)):-int(round(sweepMutTimeGALow))])
    sys.stderr.write("population size at sweep mut time was {0}; max in range was {1}, so accept prob is {0}/{1}={2}\n".format(nMut, nMax, nMut/float(nMax)))
    return nMut/float(nMax)
