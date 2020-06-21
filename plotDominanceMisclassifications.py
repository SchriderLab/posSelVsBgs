import sys, os, gzip
import matplotlib
matplotlib.use('Agg')
from miscPlottingFuncs import *
import matplotlib.pyplot as plt
import numpy as np

"""usage example:
python plotDominanceMisclassifications.py simData/boykoDFE/dominanceTestResults/humanEquilib.out simData/boykoDFE/dominanceTestResults/tennessenEuro.out simData/flyHubDFE/dominanceTestResults/sheehanSong.out 1000 plots/domBarPlot.pdf
"""

humanEquilibFileName, tennessenEuroFileName, sheehanSongFileName, numReps, plotFileName = sys.argv[1:]
numReps = int(numReps)

def fnameToTitle(fname):
    pop = fname.split("/")[-1].split(".out")[0]
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

def binom_interval(success, total, confint=0.95):
    quantile = (1 - confint) / 2.
    lower = beta.ppf(quantile, success, total - success + 1)
    upper = beta.ppf(1 - quantile, success + 1, total - success)
    if success <= 0.0:
        lower = 0.0
    frac = success/total
    return (frac-lower, upper-frac)

def readDomTestFile(fname):
    first = True
    dominanceVals = []
    with open(fname) as f:
        for line in f:
            if first:
                header = line.strip().split()
                counts = {}
                for j in range(1, len(header)):
                    counts[header[j]] = []
                first = False
            else:
                line = line.strip().split()
                dominanceVals.append(float(line[0].split("=")[1])/2)
                for j in range(1, len(header)):
                    counts[header[j]].append(float(line[j]))
    return dominanceVals, counts['hard'], counts['soft']

fig, ax = plt.subplots(3,1, figsize=(4.5, 15))

i=0
for fname in [humanEquilibFileName, tennessenEuroFileName, sheehanSongFileName]:
    dominanceVals, hardCounts, softCounts = readDomTestFile(fname)

    counts = (np.array(hardCounts)+np.array(softCounts))*1000
    yerr = [binom_interval(int(x), numReps) for x in counts]
    lowers = [x[0] for x in yerr]
    uppers = [x[1] for x in yerr]
    p1 = ax[i].bar(np.arange(len(dominanceVals)), softCounts, width=0.6, color='violet')
    p2 = ax[i].bar(np.arange(len(dominanceVals)), hardCounts, width=0.6, bottom=softCounts, color='orange', yerr=(lowers, uppers))

    ax[i].set_title(fnameToTitle(fname), fontsize=14)
    ax[i].set_ylabel("Fraction of BGS simulations misclassified as sweep", fontsize=12)
    ax[i].legend((p1[0], p2[0]), ('Classified as soft', 'Classified as hard'), loc="upper left")
    ax[i].set_xlabel('Dominance', fontsize=14)
    ax[i].set_xticks(np.arange(len(dominanceVals)))
    ax[i].set_xticklabels(dominanceVals)
    ax[i].set_ylim((0, 0.16))
    i += 1

fig.tight_layout()
fig.savefig(plotFileName)
