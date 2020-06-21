import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import buildHeatmapForClassifierAndSummarizeRegions
import scipy.stats

def demogToTitle(pop):
    if pop == "humanEquilib":
        return "Human (constant population size)"
    elif pop == "tennessenEuro":
        return "Human (Tennessen European model)"
    elif pop == "sheehanSong":
        return "Drosophila (Sheehan and Song model)"
    else:
        sys.exit("Bad input file name. ARRRGHHHHHHHHHHHH!!!!!!!!!")

### The following code taken/modified from paulgb/binom_interval.py on github
from scipy.stats import beta

def binom_interval(outcomes, confint=0.95):
    successes = outcomes.count(1)
    total = len(outcomes)
    frac = successes/total
    quantile = (1 - confint) / 2.
    lower = beta.ppf(quantile, successes, total - successes + 1)
    upper = beta.ppf(1 - quantile, successes + 1, total - successes)
    if successes <= 0.0:
        lower = 0.0
    return (frac-lower, upper-frac)

def plotBinnedMisclassRateData(panelData, figFileName):
    fig, ax = plt.subplots(3,1, figsize=(14, 15))
    i = 0
    for binnedMisclassifications, title in panelData:
        yerr = [binom_interval(misclassForBin) for misclassForBin in binnedMisclassifications]
        lowers = [x[0] for x in yerr]
        uppers = [x[1] for x in yerr]
        bar = ax[i].bar(np.arange(len(binnedMisclassifications)), [x.count(1)/len(x) for x in binnedMisclassifications], yerr=(lowers, uppers))

        ax[i].set_title(title, fontsize=16)
        ax[i].set_ylabel("Fraction of simulations misclassified as sweep", fontsize=14)
        ax[i].set_xticks(np.arange(len(binnedMisclassifications)))
        ax[i].set_xticklabels(['BGS; no recomb\n'+r'$n$={}'.format(len(binnedMisclassifications[0]))] +
                              ['BGS; recomb bin {}\n'.format(i+1)+r'$n$={}'.format(len(binnedMisclassifications[i+1])) for i in range(len(binnedMisclassifications)-2)] +
                              ['Neutral\n'+r'$n$={}'.format(len(binnedMisclassifications[-1]))], fontsize=14)
        ax[i].set_ylim((0, 0.22))
        ax[i].tick_params(axis="y", labelsize=14)

        j = 0
        for rect in bar[:-1]:
            table = [(binnedMisclassifications[j].count(1), binnedMisclassifications[j].count(0)),
                    (binnedMisclassifications[-1].count(1), binnedMisclassifications[-1].count(0))]
            oR, pVal = scipy.stats.fisher_exact(table)
            print("FET result for recomb bin {} vs neutrality: OR={}; P={}".format(j, oR, pVal))
            height = rect.get_height()

            pvalStr = r'$P$ = {:.2g}'.format(pVal)
            if float(pVal) < 0.05:
                pvalStr =  pvalStr+"*"
            ax[i].text(rect.get_x() + rect.get_width()/2.0, rect.get_height(), pvalStr, ha='right', va='bottom', fontsize=12)

            j += 1

        i += 1

    fig.tight_layout()
    fig.savefig(figFileName)

classSummBaseStr = "simData/sweeps/{}/classificationSummaries/{}.txt"

def binMisclassifications(preds, recRates):
    nonZeroRecRateIndices = [i for i in range(len(recRates)) if recRates[i] != 0]

    binnedData = [[]]
    binIndex = 0
    totBins = 4
    doneCount = 0
    for i in sorted(nonZeroRecRateIndices, key=lambda x: recRates[x]):
        if doneCount >= len(nonZeroRecRateIndices)*((binIndex+1)/totBins):
            binIndex += 1
            binnedData.append([])
        binnedData[binIndex].append((preds[i], recRates[i]))
        doneCount += 1
    assert len(binnedData) == totBins

    newBinnedData = [[(preds[i], recRates[i]) for i in range(len(recRates)) if recRates[i] == 0]] + binnedData

    binnedMisclassifications = []
    for b in newBinnedData:
        print(min([x[1] for x in b]), np.mean([x[1] for x in b]), max([x[1] for x in b]), len(b))
        print(np.mean([x[0] for x in b]))
        print("")
        binnedMisclassifications.append([x[0] for x in b])
    return binnedMisclassifications

def readClassificationSummaryFile(fname):
    preds, recRates = [], []

    first = True
    with open(fname) as f:
        for line in f:
            if first:
                first = False
            else:
                sim, hardProb, hardLinkedProb, softProb, softLinkedProb, \
                neutProb, centralSelTotal, centralSelExon, centralSelCnc, \
                flankSelTotal, flankSelExon, flankSelCnc, totalSelTotal, \
                totalSelExon, totalSelCnc, totalRecRateInMorgans = line.strip().split()
                nonSweepProbs = [float(hardLinkedProb), float(softLinkedProb), float(neutProb)]
                if all([float(hardProb) > x for x in nonSweepProbs]) or all([float(softProb) > x for x in nonSweepProbs]):
                    preds.append(1)
                else:
                    preds.append(0)
                recRates.append(float(totalRecRateInMorgans))
    return preds, recRates

def classifyNeutFileAndRecordMisclassifications(fname, modelArch, modelWeights):
    predFracs, preds = buildHeatmapForClassifierAndSummarizeRegions.predictSingleFile(modelArch, modelWeights, fname)
    misClassifications = []
    for predLine in preds:
        if predLine[4] in ["hard","soft"]:
            misClassifications.append(1)
        else:
            misClassifications.append(0)
    return misClassifications

demogDfeLs = [("humanEquilib", "boykoDFE"), ("tennessenEuro", "boykoDFE"), ("sheehanSong", "flyHubDFE")]
neutMisclassifications = {}
for demog, dfe in demogDfeLs:
    neutFName = 'simData/{}/randomExamples/realNeut/shicFeatureVecs/{}.fvec'.format(dfe, demog)
    modelArch = 'simData/sweeps/{}/classifier/clf.mod.json'.format(demog)
    modelWeights = 'simData/sweeps/{}/classifier/clf.mod.weights.hdf5'.format(demog)
    misclassifications = classifyNeutFileAndRecordMisclassifications(neutFName, modelArch, modelWeights)
    neutMisclassifications[demog] = misclassifications

for selType in ['realBgs', 'realBgsWeakCNC']:
    panelData = []
    for demog, dfe in demogDfeLs:
        print(demog, dfe, selType)
        preds, recRates = readClassificationSummaryFile(classSummBaseStr.format(demog, selType))
        binnedMisclassifications = binMisclassifications(preds, recRates)
        binnedMisclassifications.append(neutMisclassifications[demog])
        panelData.append((binnedMisclassifications, demogToTitle(demog)))
        print("\n")
    figFileName = "plots/recRate_{}.pdf".format(selType)
    plotBinnedMisclassRateData(panelData, figFileName)
