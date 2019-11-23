import sys, os, gzip
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import miscPlottingFuncs
import selRegionsFromAnnot

statDir, baseDir, geneAnnotFileName, cncAnnotFileName, winSize, cncWinSize, cncMax, plotFileName = sys.argv[1:]

targC, targSE = plotFileName.split("/")[-1].split(".")[0].split(":")
targS, targE = [int(x) for x in targSE.split("-")]
if not targC.startswith("chr"):
    targC = "chr" + targC

cncMax = float(cncMax)

demogScenarioName = statDir.rstrip("/").split("/")[-1]
winSize = int(winSize)
cncWinSize = int(cncWinSize)
additionalWindowingFactor = 1

def getSpecies(statDir):
    demogScenarioName = statDir.rstrip("/").split("/")[-1]
    if demogScenarioName in ["tennessenEuro", "humanEquilib"]:
        return "human"
    elif demogScenarioName in ["sheehanSong"]:
        return "drosophila"
    else:
        return "?"

def readStatMeansFromStatDir(statDir, winSize):
    statMeans = {}
    prevSubWinMidpts = None

    statVecs, subWinMidpts = miscPlottingFuncs.readStatsFromDir(statDir, winSize)
    if prevSubWinMidpts:
        assert prevSubWinMidpts == subWinMidpts
    prevSubWinMidpts = subWinMidpts
    for statName in statVecs:
        statMeans[statName] = []
        for i in range(len(subWinMidpts)):
            statMeans[statName].append(np.mean([x[i] for x in statVecs[statName]]))

    return statMeans, subWinMidpts

def plotBigStatsOneSet(stats, randomStats, subWinMidpts, color, marker, additionalWindowingFactor, statNames, cncWins, exonDensities, cncDensities, plotFileName):
    fig, ax = plt.subplots(4, 3, figsize=(20, 15))

    for statIndex in range(len(statNames)):
        i = int(statIndex / 3)
        j = statIndex % 3

        assert len(stats[statNames[statIndex]]) == len(subWinMidpts)

        if additionalWindowingFactor > 1:
            newStats = running_mean(stats[statNames[statIndex]], 10)
            p = ax[i,j].plot(miscPlottingFuncs.trimMidpoints(subWinMidpts, 10), newStats, color=color, lw=1, marker=marker)
        else:
            p = ax[i,j].plot(subWinMidpts, stats[statNames[statIndex]], color=color, lw=1, marker=marker)

        overallMin, overallMax = miscPlottingFuncs.overallMinAndMax(randomStats, statNames[statIndex])

        if getSpecies(statDir) == "human":
            ylimFunc = miscPlottingFuncs.getStatYLimsHuman
        elif getSpecies(statDir) == "drosophila":
            ylimFunc = miscPlottingFuncs.getStatYLimsDrosophila
        else:
            sys.exit("Couldn't figure out which species we are plotting!! AARRRGHHHHH")
        ylim = ylimFunc(statNames[statIndex], overallMin, overallMax)

        print(statNames[statIndex], overallMin, overallMax, ylim)
        ax[i,j].set_ylim(ylim)
        ax[i,j].set_ylabel(miscPlottingFuncs.formatStatName(statNames[statIndex]), fontsize=20)

        curr = ax[i,j].twinx()
        curr.stackplot(cncWins, [cncDensities, exonDensities], colors=["gray", "black"], alpha=0.5)
        curr.set_ylabel("Conserved Element Density", color="gray", fontsize=16)
        curr.spines["right"].set_color("gray")
        curr.yaxis.set_ticks_position("right")
        [(t.set_color("gray"), t.set_fontsize(14)) for t in curr.yaxis.get_ticklabels()]
        curr.spines["left"].set_visible(False)
        curr.tick_params(axis='y', colors="gray")
        curr.set_ylim(bottom=0, top=cncMax)

        plt.setp(ax[i, j].get_xticklabels(), fontsize=12)
        plt.setp(ax[i, j].get_yticklabels(), fontsize=18)
    fig.tight_layout()
    fig.savefig(plotFileName)

colors = ['black','red','blue','violet','orange','cyan','gray','brown']
markers = ['o', 'v', '^', 'x', 's', '+', 'D', '']

statMeans, subWinMidpts = readStatMeansFromStatDir(statDir, winSize)
randomStatMeans, randomSubWinMidpts = miscPlottingFuncs.readStatMeansFromBaseDir(baseDir, demogScenarioName, winSize)

cncCoords = selRegionsFromAnnot.readCncCoordsInTargRegion(cncAnnotFileName, targC, targS, targE)
cncWins, cncDensities = selRegionsFromAnnot.coordsToWindowedDensities(cncCoords, targS, targE, winSize=cncWinSize)
cncWins=[midpt-targS for midpt in cncWins]
cncDensities = [x/cncWinSize for x in cncDensities]

exonCoords = selRegionsFromAnnot.readGeneCoordsInTargRegion(geneAnnotFileName, targC, targS, targE)
exonWins, exonDensities = selRegionsFromAnnot.coordsToWindowedDensities(exonCoords, targS, targE, winSize=cncWinSize)
exonDensities = [x/cncWinSize for x in exonDensities]

if additionalWindowingFactor <= 1:
    sys.stderr.write("Not windowing our stats any further.\n")
else:
    sys.stderr.write("Plotting moving averages of our input windows (avg across %d windows)\n" %(additionalWindowingFactor))

statsToPlot = 'pi tajD fayWuH maxFDA HapCount H12 H2/H1 ZnS Omega distVar distSkew distKurt'.split()
plotBigStatsOneSet(statMeans, randomStatMeans, subWinMidpts, colors[0], markers[0], additionalWindowingFactor, statsToPlot, cncWins, exonDensities, cncDensities, plotFileName)
