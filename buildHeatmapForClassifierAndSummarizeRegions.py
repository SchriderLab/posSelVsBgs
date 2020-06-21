import matplotlib as mpl
mpl.use("Agg")
import sys, os, uuid
import numpy as np
import matplotlib.pyplot as plt
import selRegionsFromAnnotOneRead
import overlapper
from sklearn.preprocessing import normalize

def makeConfusionMatrixHeatmap(data, title, trueClassOrderLs, predictedClassOrderLs, figFileName, figSize="normal"):
    data = np.array(list(reversed(data)))
    data = normalize(data, axis=1, norm='l1')
    fig, ax = plt.subplots()
    heatmap = ax.pcolor(data, cmap=plt.cm.Blues, vmin=0.0, vmax=1.0)

    for i in reversed(range(len(trueClassOrderLs))):
        for j in range(len(predictedClassOrderLs)):
            val = 100*data[i, j]
            if val > 50:
                c = '0.9'
            else:
                c = 'black'
            ax.text(j + 0.5, i + 0.5, '%.2f%%' % val, horizontalalignment='center', verticalalignment='center', color=c, fontsize=9)

    cbar = plt.colorbar(heatmap, cmap=plt.cm.Blues, ax=ax)
    cbar.set_label("Fraction of simulations assigned to class", rotation=270, labelpad=20, fontsize=11)

    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(data.shape[1]) + 0.5, minor=False)
    ax.set_yticks(np.arange(data.shape[0]) + 0.5, minor=False)
    ax.tick_params(axis='y', which='both',length=0)
    ax.axis('tight')

    plt.tick_params(axis='x', which='both', top='off')
    plt.tick_params(axis='y', which='both', right='off')
    #lebels
    ax.set_xticklabels(predictedClassOrderLs, minor=False, fontsize=9, rotation=45)
    ax.set_yticklabels(reversed(trueClassOrderLs), minor=False, fontsize=9)
    ax.set_xlabel("Predicted class")
    ax.set_ylabel("True class")

    #plt.suptitle(title, fontsize=14, y=1.03)
    plt.tight_layout()
    if figSize == "bigly":
        fig.set_size_inches(4.5, 3.5)
    elif figSize == "wide":
        fig.set_size_inches(7, 2.5)
    else:
        fig.set_size_inches(3.5, 2.5)
    plt.savefig(figFileName, bbox_inches='tight')

def generatePredFileName(fvecFileName):
    return fvecFileName + "." + uuid.uuid4().hex + ".tmp"

def readPredsIntoVector(tmpPredFileName, classToIndex = {'hard':0, 'linkedHard':1, 'soft':2, 'linkedSoft':3, 'neutral':4}):
    currPreds = np.loadtxt(tmpPredFileName, skiprows=1, dtype="str")
    predFracs = [0]*len(classToIndex)
    for predLine in currPreds:
        predFracs[classToIndex[predLine[4]]] += 1
    denom = sum(predFracs)
    for i in range(len(predFracs)):
        predFracs[i] /= denom
    return predFracs, currPreds

def writeTmpFvecInForPreds(fvecFileName, tmpFvecInFileName):
    with open(tmpFvecInFileName, "wt") as wFile:
        with open(fvecFileName, "rt") as rFile:
            for line in rFile:
                wFile.write("chrom\tclassifiedWinStart\tclassifiedWinEnd\tbigWinRange\t" + line)

def predictSingleFile(modelArch, modelWeights, fvecFileName):
    tmpPredFileName = generatePredFileName(fvecFileName)
    tmpFvecInFileName = tmpPredFileName + ".input"
    writeTmpFvecInForPreds(fvecFileName, tmpFvecInFileName)
    cmd = "python ~/diploSHIC/diploSHIC.py predict %s %s %s %s" %(modelArch, modelWeights, tmpFvecInFileName, tmpPredFileName)
    os.system(cmd)
    predFracs, preds = readPredsIntoVector(tmpPredFileName)
    os.system("rm %s %s" %(tmpPredFileName, tmpFvecInFileName))
    return predFracs, preds

def predictShicTestSets(modelArch, modelWeights, shicTestDir):
    fileNames = ["hard.fvec", "linkedHard.fvec", "soft.fvec", "linkedSoft.fvec", "neut.fvec"]
    predFracs = []
    for fileName in fileNames:
        predFracs.append(predictSingleFile(modelArch, modelWeights, shicTestDir + "/" + fileName)[0])
    return predFracs

def runPreds(modelArch, modelWeights, shicTestSetDir, absurdBgsFvecFileName, realBgsFvecFileName, realBgsWeakCNCFvecFileName, realNeutFvecFileName):
    predictedClassOrderLs = ["Hard sweep", "Hard-linked", "Soft sweep", "Soft-linked", "Neutral"]
    trueClassOrderLs = predictedClassOrderLs + ["Central BGS", "Real BGS", "Real BGS (Weak CNE)"]
    predFracs = predictShicTestSets(modelArch, modelWeights, shicTestSetDir)[:-1] #only need to use the fwdpy neut sims, coal ones redundant (and also slightly diff model)

    currPredFracs, realNeutPreds = predictSingleFile(modelArch, modelWeights, realNeutFvecFileName)
    predFracs.append(currPredFracs)

    currPredFracs, absurdBgsPreds = predictSingleFile(modelArch, modelWeights, absurdBgsFvecFileName)
    predFracs.append(currPredFracs)

    currPredFracs, realBgsPreds = predictSingleFile(modelArch, modelWeights, realBgsFvecFileName)
    predFracs.append(currPredFracs)

    currPredFracs, realBgsWeakCNCPreds = predictSingleFile(modelArch, modelWeights, realBgsWeakCNCFvecFileName)
    predFracs.append(currPredFracs)

    return predFracs, [realBgsPreds, realBgsWeakCNCPreds, realNeutPreds], trueClassOrderLs, predictedClassOrderLs

def logFileToCoordsAndRecRate(logFileName):
    with open(logFileName) as f:
        for line in f:
            if line.startswith("modeling simulated region after"):
                coords = line.strip().split()[-1]
                c,se = coords.split(":")
                s, e = [int(x) for x in se.split("-")]
            elif "r (total Morgans)" in line:
                #currRepNum: 0; u: 2.80089e-07; r (total Morgans): 0.489661; Nanc: 6467; N0: 6359; numGens: 74670; totSelRegionSize: 55659; fwdpy seed: 179791843; np seed: 1337341592
                recRate = float(line.strip().split()[7].rstrip(";"))
    return c, s, e, recRate

def calcRegionSummaries(c, s, e, sregionCoords, totSelRegionSize, rregionCoords, totalRecRateInMorgans, numSubWins=11):
    winLen = e-s+1
    centralWinLen = winLen/numSubWins/2
    winCenter = 1+(winLen/2)
    centralSubWinStart = winCenter - centralWinLen/2
    centralSubWinEnd = centralSubWinStart + centralWinLen - 1
    leftFlankStart = 1
    leftFlankEnd = centralSubWinStart - 1
    rightFlankStart = centralSubWinEnd + 1
    rightFlankEnd = winLen

    centralSelTotal = 0
    centralSelExon = 0
    centralSelCnc = 0

    flankSelTotal = 0
    flankSelExon = 0
    flankSelCnc = 0

    totalSelTotal = 0
    totalSelExon = 0
    totalSelCnc = 0

    for currS, currE, currSelFrac, currSelRatio in sregionCoords:
        selLen = currE - currS + 1
        centralOverRange = overlapper.overlap(currS, currE, centralSubWinStart, centralSubWinEnd)
        if centralOverRange:
            centralOverLen = centralOverRange[1] - centralOverRange[0] + 1
            centralSelTotal += centralOverLen
            if currSelRatio == 1:
                centralSelExon += centralOverLen
            else:
                centralSelCnc += centralOverLen

        leftOverRange = overlapper.overlap(currS, currE, leftFlankStart, leftFlankEnd)
        leftOverLen, rightOverLen = 0, 0
        if leftOverRange:
            leftOverLen = leftOverRange[1] - leftOverRange[0] + 1
        rightOverRange = overlapper.overlap(currS, currE, rightFlankStart, rightFlankEnd)
        if rightOverRange:
            rightOverLen = rightOverRange[1] - rightOverRange[0] + 1
        flankOverLen = leftOverLen + rightOverLen

        flankSelTotal += flankOverLen
        totalSelTotal += selLen
        if currSelRatio == 1:
            flankSelExon += flankOverLen
            totalSelExon += selLen
        else:
            flankSelCnc += flankOverLen
            totalSelCnc += selLen
    return centralSelTotal/centralWinLen, centralSelExon/centralWinLen, centralSelCnc/centralWinLen,\
           flankSelTotal/(winLen-centralWinLen), flankSelExon/(winLen-centralWinLen), flankSelCnc/(winLen-centralWinLen),\
           totalSelTotal/winLen, totalSelExon/winLen, totalSelCnc/winLen,\
           totalRecRateInMorgans

def writePredAndRegionSummaries(preds, logDir, summaryOutDir, scenarioName, geneAnnotFileName, phastConsFileName, recRateFileName, rMean):
    logFileNames = os.listdir(logDir)
    logFileNames.sort(key=lambda x: int(x.split("_")[1].split(".log")[0]))
    summaryOutFile = summaryOutDir + "/" + scenarioName + ".txt"

    cncSelRatio = 0
    with open(summaryOutFile, "wt") as f:        
        header = "sim", "hardProb", "hardLinkedProb", "softProb", "softLinkedProb", "neutProb", "centralSelTotal", "centralSelExon", "centralSelCnc",\
           "flankSelTotal", "flankSelExon", "flankSelCnc",\
           "totalSelTotal", "totalSelExon", "totalSelCnc",\
           "totalRecRateInMorgans"
        header = "\t".join(header)
        f.write(header + "\n")

        for i in range(len(logFileNames)):
            winC, winStart, winEnd, totalRecRateInMorgans = logFileToCoordsAndRecRate(logDir + "/" + logFileNames[i])
            sys.stderr.write("generating summaries for {} rep {} (based on {}:{}-{})\n".format(scenarioName, i, winC, winStart, winEnd))
            exonCoords, phastConsCoords = selRegionsFromAnnotOneRead.readSelRegionsFromGtf(geneAnnotFileName, phastConsFileName=phastConsFileName)
            sregionCoords, nregionCoords, totSelRegionSize = selRegionsFromAnnotOneRead.readSelRegionsInWin(exonCoords, winC, winStart, winEnd, cncSelRatio, phastConsCoords=phastConsCoords)
            recRates = selRegionsFromAnnotOneRead.readRecRegionsFromWig(recRateFileName)
            rregionCoords, totalRecRateInMorgans2 = selRegionsFromAnnotOneRead.readRecRegionsInWin(recRates, winC, winStart, winEnd, rRescale=rMean)
            summaries = calcRegionSummaries(winC, winStart, winEnd, sregionCoords, totSelRegionSize, rregionCoords, totalRecRateInMorgans)
            chrom, cWinStart, cWinEnd, bigWin, predClass, pNeutral, pLSoft, pLHard, pSoft, pHard = preds[i]
            f.write(logFileNames[i] + "\t" + "\t".join([str(x) for x in [pHard, pLHard, pSoft, pLSoft, pNeutral]]) + "\t" + "\t".join([str(x) for x in summaries]) + "\n")

def main():
    modelArch, modelWeights, shicTestSetDir, absurdBgsFvecFileName, realBgsFvecFileName, realBgsWeakCNCFvecFileName, realNeutFvecFileName, realBgsLogDir, realBgsWeakCNCLogDir, realNeutLogDir, geneAnnotFileName, phastConsFileName, recRateFileName, rMean, summaryOutDir, figFileName = sys.argv[1:]
    rMean = float(rMean)
    allPredFracs, (realBgsPreds, realBgsWeakCNCPreds, realNeutPreds), trueClassOrderLs, predictedClassOrderLs = runPreds(modelArch, modelWeights, shicTestSetDir, absurdBgsFvecFileName, realBgsFvecFileName, realBgsWeakCNCFvecFileName, realNeutFvecFileName)

    writePredAndRegionSummaries(realBgsPreds, realBgsLogDir, summaryOutDir, "realBgs", geneAnnotFileName, phastConsFileName, recRateFileName, rMean)
    writePredAndRegionSummaries(realBgsWeakCNCPreds, realBgsWeakCNCLogDir, summaryOutDir, "realBgsWeakCNC", geneAnnotFileName, phastConsFileName, recRateFileName, rMean)
        
    makeConfusionMatrixHeatmap(allPredFracs, "Predshictions", trueClassOrderLs, predictedClassOrderLs, figFileName, figSize="bigly")

if __name__ == "__main__":
    main()
