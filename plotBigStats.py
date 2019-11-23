import sys, os, gzip
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import miscPlottingFuncs

baseDir, sweepBaseDir, demogScenarioName, winSize, plotFileName = sys.argv[1:]
winSize = int(winSize)
additionalWindowingFactor = 1

def getSpecies(demogScenarioName):
    if demogScenarioName in ["tennessenEuro", "humanEquilib"]:
        return "human"
    elif demogScenarioName in ["sheehanSong"]:
        return "drosophila"
    else:
        return "?"

def plotBigStats(stats, neutExpectationsCoal, neutExpectationsFwd, selTypes, subWinMidpts, colors, markers, additionalWindowingFactor, statNames, plotFileName):
    fig, ax = plt.subplots(4, 3, figsize=(20, 15))

    for statIndex in range(len(statNames)):
        i = int(statIndex / 3)
        j = statIndex % 3
        for selTypeIndex in range(len(selTypes)):
            assert len(stats[selTypes[selTypeIndex]][statNames[statIndex]]) == len(subWinMidpts)

            if additionalWindowingFactor > 1:
                newStats = running_mean(stats[selTypes[selTypeIndex]][statNames[statIndex]], 10)
                p = ax[i,j].plot(miscPlottingFuncs.trimMidpoints(subWinMidpts, 10), newStats, color=colors[selTypeIndex], lw=1, marker=markers[selTypeIndex], label=miscPlottingFuncs.scenarioName(selTypes[selTypeIndex]))
            else:
                p = ax[i,j].plot(subWinMidpts, stats[selTypes[selTypeIndex]][statNames[statIndex]], color=colors[selTypeIndex], lw=1, marker=markers[selTypeIndex], label=miscPlottingFuncs.scenarioName(selTypes[selTypeIndex]))
        ax[i,j].axhline(neutExpectationsCoal[statNames[statIndex]], color="gray", ls="--", lw=1)
        ax[i,j].axhline(neutExpectationsFwd[statNames[statIndex]], color="black", ls="-", lw=1)
        ax[i,j].set_ylabel(miscPlottingFuncs.formatStatName(statNames[statIndex]), fontsize=20)
        overallMin, overallMax = miscPlottingFuncs.overallMinAndMax(stats, statNames[statIndex])
        if getSpecies(demogScenarioName) == "human":
            ylimFunc = miscPlottingFuncs.getStatYLimsHuman
        elif getSpecies(demogScenarioName) == "drosophila":
            ylimFunc = miscPlottingFuncs.getStatYLimsDrosophila
        else:
            sys.exit("Couldn't figure out which species we are plotting!! AARRRGHHHHH")
        ylim = ylimFunc(statNames[statIndex], overallMin, overallMax)
        print(statNames[statIndex], overallMin, overallMax, ylim)
        ax[i,j].set_ylim(ylim)
        plt.setp(ax[i, j].get_xticklabels(), fontsize=14)
        plt.setp(ax[i, j].get_yticklabels(), fontsize=18)
        if (i,j) == (0,0):
            ax[i,j].legend(loc='lower left', fontsize=12)
    fig.tight_layout()
    fig.savefig(plotFileName)

colors = ['black','red','blue','violet','orange','cyan','gray','brown']
markers = ['o', 'v', '^', 'x', 's', '+', 'D', '']

statMeans, subWinMidpts = miscPlottingFuncs.readStatMeansFromBaseDir(baseDir, demogScenarioName, winSize)
hardStatMeans, hardSubWinMidpts = miscPlottingFuncs.readStatMeansFromDir("%s/%s/testDataFwdpyStats/simHard_5/" %(sweepBaseDir, demogScenarioName), winSize)
statMeans["hardSweep"] = hardStatMeans
softStatMeans, softSubWinMidpts = miscPlottingFuncs.readStatMeansFromDir("%s/%s/testDataFwdpyStats/simSoft_5/" %(sweepBaseDir, demogScenarioName), winSize)
statMeans["softSweep"] = softStatMeans
neutExpectationsCoal = miscPlottingFuncs.readOverallStatMeansFromDir("%s/%s/testDataFwdpyStats/simNeut/" %(sweepBaseDir, demogScenarioName), winSize)
neutExpectationsFwd = miscPlottingFuncs.readOverallStatMeansFromDir(baseDir + "/realNeut/stats/" + demogScenarioName, winSize)

if additionalWindowingFactor <= 1:
    sys.stderr.write("Not windowing our stats any further.\n")
else:
    sys.stderr.write("Plotting moving averages of our input windows (avg across %d windows)\n" %(additionalWindowingFactor))

statsToPlot = 'pi tajD fayWuH maxFDA HapCount H12 H2/H1 ZnS Omega distVar distSkew distKurt'.split()
scenariosToPlot = ["hardSweep", "softSweep", "absurdBgs", "realBgs"]
plotBigStats(statMeans, neutExpectationsCoal, neutExpectationsFwd, scenariosToPlot, subWinMidpts, colors, markers, additionalWindowingFactor, statsToPlot, plotFileName)
