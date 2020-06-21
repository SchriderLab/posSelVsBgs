import sys, os, gzip
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import miscPlottingFuncs

statPrefix, demogScenarioName, winSize, repRange, logPrefix, plotPrefix = sys.argv[1:]
winSize = int(winSize)
repRangeLow, repRangeHigh = [int(x) for x in repRange.split("-")]
if logPrefix == "None":
    logPrefix = False

def readCoordsFromLog(logFileName):
    with open(logFileName, "rt") as f:
        for line in f:
            if line.startswith("modeling simulated region after "):
                coordStr = line.split("modeling simulated region after ")[1]
    return coordStr

def prefixToPlotTitle(statPrefix):
    if "tennessenEuro" in statPrefix:
        title = "Human (Tennessen Euro)"
    elif "humanEquilib" in statPrefix:
        title = "Human (Equilibrium)"
    elif "sheehanSong" in statPrefix:
        title = "Drosophila (Sheehan Song)"
    else:
        raise ValueError

    if "simHard" in statPrefix:
        title += "; Hard sweep"
    elif "simSoft" in statPrefix:
        title += "; Soft sweep"
    elif "realBgsWeak" in statPrefix:
        title += "; Real Bgs (Weak CNE)"
    elif "realBgs" in statPrefix or "dominanceRealBgsRandom" in statPrefix:
        title += "; Real Bgs"
    elif "absurdBgs" in statPrefix:
        title += "; Central Bgs"
    elif "realNeut" in statPrefix:
        title += "; Neutral"
    else:
        raise ValueError
    return title

def getSpecies(demogScenarioName):
    if demogScenarioName in ["tennessenEuro", "humanEquilib"]:
        return "human"
    elif demogScenarioName in ["sheehanSong"]:
        return "drosophila"
    else:
        return "?"

def plotBigStatsOneSim(stats, species, overallLowHigh, subWinMidpts, statNames, repIndex, plotTitle, plotFileName):
    fig, ax = plt.subplots(4, 3, figsize=(20, 15))

    fig.suptitle(plotTitle, fontsize=20)
    for statIndex in range(len(statNames)):
        i = int(statIndex / 3)
        j = statIndex % 3
        assert len(stats[statNames[statIndex]]) == len(subWinMidpts)
        p = ax[i,j].plot(subWinMidpts, stats[statNames[statIndex]], color='black', lw=1, marker='o')

        overallLow, overallHigh = overallLowHigh[statNames[statIndex]]
        currLow, currHigh = ax[i,j].get_ylim()
        ax[i,j].set_ylim((min(overallLow*0.9, currLow), max(overallHigh*1.1, currHigh)))

        ax[i,j].set_ylabel(miscPlottingFuncs.formatStatName(statNames[statIndex]), fontsize=20)
        plt.setp(ax[i, j].get_xticklabels(), fontsize=14)
        plt.setp(ax[i, j].get_yticklabels(), fontsize=18)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(plotFileName)
    plt.close()

statsToPlot = 'pi tajD fayWuH maxFDA HapCount H12 H2/H1 ZnS Omega distVar distSkew distKurt'.split()
titlePrefix = prefixToPlotTitle(statPrefix)

allStats = []
statsDump = {}
for repIndex in range(repRangeLow, repRangeHigh):
    statFileName = statPrefix + "{}.stats".format(repIndex)
    sys.stderr.write("reading {}\n".format(statFileName))
    stats, subWinMidpts = miscPlottingFuncs.readStatsFromFile(statFileName, winSize)
    allStats.append(stats)

    for statIndex in range(len(statsToPlot)):
        if not statsToPlot[statIndex] in statsDump:
            statsDump[statsToPlot[statIndex]] = []
        statsDump[statsToPlot[statIndex]] += stats[statsToPlot[statIndex]]

overallLowHigh = {}
for statIndex in range(len(statsToPlot)):
    statsDump[statsToPlot[statIndex]].sort()
    lowIndex = int(len(statsDump[statsToPlot[statIndex]])*0.025)
    highIndex = int(len(statsDump[statsToPlot[statIndex]])*0.975)
    overallLowHigh[statsToPlot[statIndex]] = statsDump[statsToPlot[statIndex]][lowIndex], statsDump[statsToPlot[statIndex]][highIndex] 

for repIndex in range(repRangeLow, repRangeHigh):
    if logPrefix:
        logFileName = logPrefix + "{}.log".format(repIndex)
        coordStr = readCoordsFromLog(logFileName)
    
    statFileName = statPrefix + "{}.stats".format(repIndex)
    sys.stderr.write("plotting {}\n".format(statFileName))

    if logPrefix:
        plotTitle = titlePrefix + "; {}".format(coordStr)
    else:
        plotTitle = titlePrefix + "; rep {}".format(repIndex)
    plotFileName = plotPrefix + "_{}.pdf".format(repIndex)
    plotBigStatsOneSim(allStats[repIndex], getSpecies(demogScenarioName), overallLowHigh, subWinMidpts, statsToPlot, repIndex, plotTitle, plotFileName)
