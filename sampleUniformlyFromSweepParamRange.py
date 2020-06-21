import sys, os
import scipy.stats
import numpy as np
import buildHeatmapForClassifierAndSummarizeRegions

np.random.seed(123003807) #randomly selected seed hard-coded for repeatability

figFileName="plots/covfefe/sweepOnBgs.pdf"
scenarios = [("tennessenEuro", "boykoDFE", "realNeut"), ("tennessenEuro", "boykoDFE", "realBgs")]
sweepTypeToNumBins = {"hard":3, "soft":2}
dataDir = "simData"

classToIndex = {'hard':0, 'linkedHard':1, 'soft':2, 'linkedSoft':3, 'neutral':4}
classLabels = ['hard', 'linkedHard', 'soft', 'linkedSoft', 'neutral']
predictedClassNames = ["Hard sweep", "Hard-linked", "Soft sweep", "Soft-linked", "Neutral"]

def readSweepLogInfo(logFileName):
    with open(logFileName) as f:
        lines = f.readlines()
        if 'rror' in lines[-1] or not lines[-2].startswith("Final sweep characteristics"):
            return None
        paramNames = lines[-2].rstrip("\n").split(": (")[-1].rstrip(")").split("\\t")
        paramVals = [float(x) for x in lines[-1].split()]

    dat = {}
    for i in range(len(paramNames)):
        if not paramNames[i] in ["sweepPos","timeSinceSweepMut","timeSinceSweepSel"]:
            dat[paramNames[i]] = paramVals[i]
    return dat

def updateSweepParams(logInfo, sweepParams):
    if sweepParams:
        for name in logInfo:
            sweepParams[name].append(logInfo[name])
    else:
        for name in logInfo:
            sweepParams[name] = [logInfo[name]]

def getBinInfo(paramValLs, numBins):
    minVal, maxVal = min(paramValLs), max(paramValLs)
    totalWidth = maxVal - minVal
    return min(paramValLs), totalWidth/numBins

def getBinIndex(val, binInfo, nBins):
    index = int((val-binInfo[0]) / binInfo[1])
    if index > 0 and val-binInfo[0] == (index * binInfo[1]):
        index -= 1
    elif index > nBins-1:
        assert val - ((index * binInfo[1]) + binInfo[0]) < 1e-6
        index -= 1
    assert index <= nBins-1
    return index

def getBinIndicesForInstance(sweepParams, i, binInfoAll, nBins):
    binIndices = []
    for name in sweepParams:
        currBinIndex = getBinIndex(sweepParams[name][i], binInfoAll[name], nBins)
        binIndices.append(currBinIndex)
    return binIndices

def updateSweepParamsWithInstances(downsampledSweepParams, sweepParams, indicesToKeep):
    for i in indicesToKeep:
        for name in sweepParams:
            if not name in downsampledSweepParams:
                downsampledSweepParams[name] = []
            downsampledSweepParams[name].append(sweepParams[name][i])
    return downsampledSweepParams

def numInstances(sweepParams):
    if sweepParams:
        name = list(sweepParams.keys())[0]
        return len(sweepParams[name])
    else:
        return 0

def binSweepParams(sweepParams, nBins):
    binnedParams = np.empty([nBins]*len(sweepParams), dtype=object)
    binnedParams.fill([])
    binnedParams = np.frompyfunc(list,1,1)(binnedParams)

    binInfoAll = {}
    for name in sweepParams:
        binInfoAll[name] = getBinInfo(sweepParams[name], nBins)
    return binnedParams, binInfoAll

def downsampleUniformly(sweepParams, binnedParams, binInfoAll, nBins):
    nInstances = numInstances(sweepParams)
    for i in range(nInstances):
        binIndices = getBinIndicesForInstance(sweepParams, i, binInfoAll, nBins)
        binnedParams[tuple(binIndices)].append(i)

    flattennedBinParams = binnedParams.flatten()
    nRepsPerBin = min([len(x) for x in flattennedBinParams])

    downsampledSweepParams = {}
    allIndicesToKeep = []
    for instanceIndexLs in flattennedBinParams:
        indicesToKeep = np.random.choice(instanceIndexLs, size=nRepsPerBin, replace=False)
        allIndicesToKeep.extend(indicesToKeep)
        updateSweepParamsWithInstances(downsampledSweepParams, sweepParams, indicesToKeep)

    indexCounts = {}
    for i in allIndicesToKeep:
        indexCounts[i] = indexCounts.get(i, 0) + 1
    assert len(indexCounts) == len(allIndicesToKeep)

    return downsampledSweepParams, allIndicesToKeep

def textHistogram(ls, binInfo, nBins, totStars=100):
    currHistStr = ""
    binnedData = {}
    for datum in ls:
        binIndex = getBinIndex(datum, binInfo, nBins)
        if not binIndex in binnedData:
            binnedData[binIndex] = []
        binnedData[binIndex].append(datum)

    minBinIndex = min(binnedData)
    maxBinIndex = max(binnedData)
    for binIndex in range(minBinIndex, maxBinIndex+1):
        numStars = round(totStars*(len(binnedData.get(binIndex, []))/float(len(ls))))
        maxVal = binInfo[0] + binInfo[1]*(binIndex+1)
        currHistStr += "{:<10.10s}:{}{}\n".format(str(maxVal), "*"*(numStars), len(binnedData.get(binIndex, [])))
    return currHistStr

def textHistograms(sweepParams, binInfoAll, nBins):
    histStr = ""
    for name in sweepParams:
        histStr += "histogram for {}:\n".format(name)
        histStr += textHistogram(sweepParams[name], binInfoAll[name], nBins)
        histStr += "\n"
    return histStr

resultsH = {}
for demog, dfe, bgsType in scenarios:
    if not (demog, dfe) in resultsH:
        resultsH[(demog,dfe)] = {"realBgs": {}, "realNeut": {}}
    for sweepType in sweepTypeToNumBins:
        nBins = sweepTypeToNumBins[sweepType]
        sys.stderr.write("parsing logs for {}.{}.{}\n".format(demog, dfe, sweepType))
        sweepParams = {}
        goodReps=0
        repIndexMapping = {}
        for i in range(1000):
            logFileName = dataDir + "{0}/randomExamplesSweep/{4}_{1}/simLogs/{2}/{2}_{3}.log".format(dfe, sweepType, demog, i, bgsType)
            logInfo = readSweepLogInfo(logFileName)
            if logInfo:
                updateSweepParams(logInfo, sweepParams)
                repIndexMapping[goodReps] = i
                goodReps += 1
        sys.stderr.write("found {} completed reps.\n".format(goodReps))
        sys.stderr.write("Text histogram before downsampling:\n")
        binnedParams, binInfoAll = binSweepParams(sweepParams, nBins)
        sys.stderr.write("Downsampling...\n")
        sweepParamsDownsampled, allIndicesToKeep = downsampleUniformly(sweepParams, binnedParams, binInfoAll, nBins)
        sys.stderr.write("Done downsampling; {} examples remaining.\n".format(numInstances(sweepParamsDownsampled)))
        sys.stderr.write("Text histogram after downsampling:\n")
        sys.stderr.write("\n\n")

        fvecLines = []
        for goodRepI in allIndicesToKeep:
            i = repIndexMapping[goodRepI]
            fvecFileName = dataDir + "{0}/randomExamplesSweep/{4}_{1}/featureVectors/{2}/{2}_{3}.fvec".format(dfe, sweepType, demog, i, bgsType)
            with open(fvecFileName) as f:
                header, fvec = f.readlines()
                fvecLines.append(fvec)

        sys.stderr.write("writing pred file for {}-{}-{}\n".format(dfe, demog, sweepType))
        os.system("mkdir -p {0}{1}/randomExamplesSweep/{4}_{2}/sweepIndicesKept/{3}".format(dataDir, dfe, sweepType, demog, bgsType))
        fvPredFName = dataDir + "{0}/randomExamplesSweep/{3}_{1}/sweepIndicesKept/{2}/{2}.fvec".format(dfe, sweepType, demog, bgsType)
        with open(fvPredFName, 'wt') as fvPredF:
            fvPredF.write(header)
            for fvecLine in fvecLines:
                fvPredF.write(fvecLine)

        sys.stderr.write("writing accepted rep index file for {}-{}-{}\n".format(dfe, demog, sweepType))
        indexFName = dataDir + "{0}/randomExamplesSweep/{3}_{1}/sweepIndicesKept/{2}/{2}.indices".format(dfe, sweepType, demog, bgsType)
        with open(indexFName, 'wt') as indexF:
            for goodRepI in allIndicesToKeep:
                i = repIndexMapping[goodRepI]
                indexF.write(str(i) + "\n")

        modelDir = dataDir + "/sweeps/{}/classifier/".format(demog)
        sys.stderr.write("running preds for {}-{}-{}-{} using clf in {}\n".format(dfe, demog, bgsType, sweepType, modelDir))
        predFracs, preds = buildHeatmapForClassifierAndSummarizeRegions.predictSingleFile(modelDir + "clf.mod.json", modelDir + "clf.mod.weights.hdf5", fvPredFName)
        print("\t".join(classLabels))
        print("\t".join([str(predFracs[classToIndex[x]]) for x in classLabels]))
        result = ([predFracs[classToIndex[x]] for x in classLabels], predFracs[classToIndex['hard']], predFracs[classToIndex['soft']], len(preds))
        resultsH[(demog,dfe)][bgsType][sweepType] = result
        print("\t".join([str(predFracs[classToIndex[x]] * len(preds)) for x in classLabels]))

def compareFracs(sweepCountNeut, totalCountNeut, sweepCountBgs, totalCountBgs, demog, dfe, sweepType, classType):
    sweepFracNeut = sweepCountNeut/totalCountNeut
    sweepFracBgs = sweepCountBgs/totalCountBgs
    print("{}-{} realNeut vs realBgs; Fraction of {} sweeps classified as {}: {} vs {}".format(demog, dfe, sweepType, classType, sweepFracNeut, sweepFracBgs))
    fetPval = scipy.stats.fisher_exact([[sweepCountNeut, totalCountNeut-sweepCountNeut], [sweepCountBgs, totalCountBgs-sweepCountBgs]])[1]
    print("counts and FET p-val: {} vs {}; p={}\n".format(sweepCountNeut, sweepCountBgs, fetPval))

allPredFracs = []
trueClasses = []
for demog, dfe in resultsH:
    for sweepType in ["hard", "soft"]:
        bgsFracLs, hardFracBgs, softFracBgs, totalCountBgs = resultsH[(demog, dfe)]["realBgs"][sweepType]
        hardCountBgs = int(round(hardFracBgs*totalCountBgs))
        softCountBgs = int(round(softFracBgs*totalCountBgs))
        
        neutFracLs, hardFracNeut, softFracNeut, totalCountNeut = resultsH[(demog, dfe)]["realNeut"][sweepType]
        hardCountNeut = int(round(hardFracNeut*totalCountNeut))
        softCountNeut = int(round(softFracNeut*totalCountNeut))

        compareFracs(hardCountNeut, totalCountNeut, hardCountBgs, totalCountBgs, demog, dfe, sweepType, "hard")
        compareFracs(softCountNeut, totalCountNeut, softCountBgs, totalCountBgs, demog, dfe, sweepType, "soft")
        sweepCountNeut = hardCountNeut + softCountNeut
        sweepCountBgs = hardCountBgs + softCountBgs
        compareFracs(sweepCountNeut, totalCountNeut, sweepCountBgs, totalCountBgs, demog, dfe, sweepType, "sweep")

        trueClasses.append("{} sweep on neutral background".format(sweepType.capitalize()))
        trueClasses.append("{} sweep under real BGS model".format(sweepType.capitalize()))
        
        allPredFracs.append(neutFracLs)
        allPredFracs.append(bgsFracLs)
buildHeatmapForClassifierAndSummarizeRegions.makeConfusionMatrixHeatmap(allPredFracs, "", trueClasses, predictedClassNames, figFileName, figSize="wide")
