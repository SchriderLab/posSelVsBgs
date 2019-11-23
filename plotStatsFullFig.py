import sys, os, gzip
import matplotlib
matplotlib.use('Agg')
import miscPlottingFuncs
import matplotlib.pyplot as plt
import numpy as np

baseDir, sweepBaseDir, demogScenarioName, winSize, additionalWindowingFactor, plotFileName = sys.argv[1:]
winSize = int(winSize)
additionalWindowingFactor = int(additionalWindowingFactor)

# copied from matplotlib.org example
def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

def directionOpp(d):
    if d == "left":
        return "right"
    elif d == "right":
        return "left"
    else:
        raise Exception

def plotStats(host, stats, statNames, subWinMidpts, neutExpectations, colors, additionalWindowingFactor, axisTitle):
    host.set_xlabel("Position along simulated window")
    directions = ["left", "right", "right"]
    curr = host
    plots = []
    offset = 0
    for statIndex in range(len(statNames)):
        statName = statNames[statIndex]
        assert len(stats[statName]) == len(subWinMidpts)

        direction = directions[statIndex]

        if statIndex > 0:
            curr = host.twinx()

        if additionalWindowingFactor > 1:
            newStats = miscPlottingFuncs.running_mean(stats[statName], additionalWindowingFactor)
            p, = curr.plot(miscPlottingFuncs.trimMidpoints(subWinMidpts, additionalWindowingFactor), newStats, color=colors[statIndex], label=miscPlottingFuncs.formatStatName(statName), lw=1)
        else:
            p, = curr.plot(subWinMidpts, stats[statName], color=colors[statIndex], label=miscPlottingFuncs.formatStatName(statName), lw=1)
        curr.axhline(neutExpectations[statName], color=colors[statIndex], ls="--", lw=1)
        plots.append(p)

        if statIndex > 1:
            offset += 65
            curr.spines[direction].set_position(('outward', offset))
        curr.spines[direction].set_color(colors[statIndex])
        curr.set_ylabel(miscPlottingFuncs.formatStatName(statName), color=colors[statIndex])
        curr.yaxis.set_ticks_position(direction)
        curr.set_title(axisTitle)
        [t.set_color(colors[statIndex]) for t in curr.yaxis.get_ticklabels()]
        curr.spines[directionOpp(direction)].set_visible(False)
        curr.tick_params(axis='y', colors=colors[statIndex])
        curr.set_ylim(bottom=0, top=1.2*max([neutExpectations[statName] for statName in statNames]))
    host.legend(handles=plots, loc='best')

colors = ['black','red','blue']
statsToPlot = ["pi", "thetaW", "thetaH"]
if additionalWindowingFactor <= 1:
    sys.stderr.write("Not windowing our stats any further.\n")
else:
    sys.stderr.write("Plotting moving averages of our input windows (avg across %d windows)\n" %(additionalWindowingFactor))

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15, 7.5))
fig.subplots_adjust(right=0.75)
neutStatDir = "%s/%s/smallStats/%s/" %(baseDir, "realNeut", demogScenarioName)
neutExpectations = miscPlottingFuncs.readOverallStatMeansFromDir(neutStatDir, winSize)
bgsModels = os.listdir(baseDir)
panels = ['A', 'B', 'C', 'D']
simModels = {}
for i in range(len(bgsModels)):
    statDir = "%s/%s/smallStats/%s/" %(baseDir, bgsModels[i], demogScenarioName)
    statMeans, subWinMidpts = miscPlottingFuncs.readStatMeansFromDir(statDir, winSize)
    simModels[bgsModels[i]] = (statMeans, subWinMidpts)

simModels["hardSweep"] = miscPlottingFuncs.readStatMeansFromDir("%s/%s/testDataFwdpySmallStats/simHard_5/" %(sweepBaseDir, demogScenarioName), winSize)
simModels["softSweep"] = miscPlottingFuncs.readStatMeansFromDir("%s/%s/testDataFwdpySmallStats/simSoft_5/" %(sweepBaseDir, demogScenarioName), winSize)

scenariosToPlot = ["hardSweep", "softSweep", "absurdBgs", "realBgs"]
for i in range(len(scenariosToPlot)):
    statMeans, subWinMidpts = simModels[scenariosToPlot[i]]
    x, y = int(i / 2), i % 2
    plotStats(axes[x, y], statMeans, statsToPlot, subWinMidpts, neutExpectations, colors, additionalWindowingFactor, miscPlottingFuncs.scenarioName(scenariosToPlot[i]))
    axes[x, y].text(-0.05, 1.05, panels[i], transform=axes[x, y].transAxes, size=20)
fig.tight_layout()
fig.savefig(plotFileName)
