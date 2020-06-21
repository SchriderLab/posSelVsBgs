import sys
import scipy.stats
import buildHeatmapForClassifier

dfeLs = ["boykoDFE", "boykoDFE", "flyHubDFE"]
demogLabels = ["humanEquilib", "tennessenEuro", "sheehanSong"]
demogLabelToName = {"humanEquilib": "Human (Equilibrium; ", "tennessenEuro": "Human (European Model; ", "sheehanSong": r"$Drosophila$ ("}

def readTestResFile(fname):
    winIndices = {}
    first = True
    data = []
    with open(fname) as f:
        for line in f:
            if first:
                header = line
                first = False
            else:
                line = line.strip().split()
                winIndex = int(line[1])
                winIndices[winIndex]=1
                data.append(line)
    return data, header, len(winIndices)/2

def getCentralCounts(resultLs, centralWinI):
    "repIndex\twinIndex\tpredClass\tprob(neutral)\tprob(likedSoft)\tprob(linkedHard)\tprob(soft)\tprob(hard)\n"
    centralCounts, nonCentralCounts = {}, {}
    for result in resultLs:
        if int(result[1]) == centralWinI:
            if not result[2] in centralCounts:
                centralCounts[result[2]] = 0
            centralCounts[result[2]] += 1
        else:
            if not result[2] in nonCentralCounts:
                nonCentralCounts[result[2]] = 0
            nonCentralCounts[result[2]] += 1
    return centralCounts, nonCentralCounts

def getSweepCounts(centralCounts, nonCentralCounts):
    centralSweepCount = centralCounts.get("hard", 0)+centralCounts.get("soft", 0)
    centralTotal = sum(centralCounts.values())
    nonCentralSweepCount = nonCentralCounts.get("hard", 0)+nonCentralCounts.get("soft", 0)
    nonCentralTotal = sum(nonCentralCounts.values())
    return centralSweepCount, centralTotal-centralSweepCount, nonCentralSweepCount, nonCentralTotal-nonCentralSweepCount

def fisherTest(centralSweepCount, centralNonSweepCount, nonCentralSweepCount, nonCentralNonSweepCount):
    oR, pVal = scipy.stats.fisher_exact(((centralSweepCount, centralNonSweepCount), (nonCentralSweepCount, nonCentralNonSweepCount)))
    return pVal

def plotCovfefe(covfefeData, demogLabels, demogLabelToName, figFileName):
    predictedClassOrder = ["hard", "linkedHard", "soft", "linkedSoft", "neutral"]
    predictedClassNames = ["Hard sweep", "Hard-linked", "Soft sweep", "Soft-linked", "Neutral"]
    trueClasses = []
    allPreds = []
    for demogLabel in demogLabels:
        central, nonCentral = covfefeData[demogLabel]
        trueClasses.append(demogLabelToName[demogLabel] + "central window)")
        trueClasses.append(demogLabelToName[demogLabel] + "non-central)")
        allPreds.append([central.get(x, 0) for x in predictedClassOrder])
        allPreds.append([nonCentral.get(x, 0) for x in predictedClassOrder])
    buildHeatmapForClassifier.makeConfusionMatrixHeatmap(allPreds, "", trueClasses, predictedClassNames, figFileName, figSize="wide")

def headerError(fname, header, headers):
    sys.exit("Error reading {}: header mismatch:\n{}\n{}\n: ".format(fname),str(header),str(list(headers)[0]))

basePath = "simData/{}/biglyTestResults/{}.out"
figFileName = "plots/covfefe/bigly.pdf"
headers = set()
covfefeData = {}
for i in range(len(demogLabels)):
    fname = basePath.format(dfeLs[i], demogLabels[i])
    print(fname)
    resultLs, header, centralWinI = readTestResFile(fname)
    if headers and not header in headers:
        headerError(fname, header, headers)
    headers.add(header)

    centralCounts, nonCentralCounts = getCentralCounts(resultLs, centralWinI)
    centralSweepCount, centralNonSweepCount, nonCentralSweepCount, nonCentralNonSweepCount = getSweepCounts(centralCounts, nonCentralCounts)
    pval = fisherTest(centralSweepCount, centralNonSweepCount, nonCentralSweepCount, nonCentralNonSweepCount)
    nonCentralSweepFrac = nonCentralSweepCount / (nonCentralSweepCount + nonCentralNonSweepCount)
    centralSweepFrac = centralSweepCount / (centralSweepCount + centralNonSweepCount)
    print("{}; central vs non-central sweep fracs: {} vs {}; p={}".format(demogLabels[i], centralSweepFrac, nonCentralSweepFrac, pval))
    covfefeData[demogLabels[i]] = (centralCounts, nonCentralCounts)
plotCovfefe(covfefeData, demogLabels, demogLabelToName, figFileName)
