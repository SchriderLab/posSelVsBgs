import os, gzip
import numpy as np

def scenarioName(selType):
    if selType == "absurdBgs":
        return "Central BGS"
    elif selType == "realBgs":
        return "Real BGS"
    elif selType == "hardSweep":
        return "Hard Sweep"
    elif selType == "softSweep":
        return "Soft Sweep"
    else:
        raise ValueError

def readStatsFromDir(statDir, winSize):
    firstHeader = True
    statVecs = {}
    for fileName in os.listdir(statDir):
        if fileName.endswith("gz"):
            fopen = gzip.open
        else:
            fopen = open
        with fopen(statDir + "/" + fileName) as f:
            first = True
            subWinMidpts = []
            for line in f:
                if first:
                    header = line.strip().split()
                    first = False
                    if firstHeader:
                        for i in range(1, len(header)):
                            statVecs[header[i]] = []
                        firstHeader = False
                    currStatVec = {}
                    for i in range(1, len(header)):
                        currStatVec[header[i]] = []
                else:
                    line = line.strip().split()
                    subWinIndex = int(line[0])
                    subWinMidpts.append(subWinIndex*winSize + winSize/2)
                    for i in range(1, len(header)):
                        currStatVec[header[i]].append(float(line[i]))
            for statName in statVecs:
                statVecs[statName].append(currStatVec[statName])
    return statVecs, subWinMidpts

def readStatMeansFromDir(statDir, winSize):
    statVecs, subWinMidpts = readStatsFromDir(statDir, winSize)

    statMeans = {}
    for statName in statVecs:
        statMeans[statName] = []
        for i in range(len(subWinMidpts)):
            statMeans[statName].append(np.nanmean([x[i] for x in statVecs[statName]]))
    return statMeans, subWinMidpts

def readOverallStatMeansFromDir(statDir, winSize):
    statVecs, subWinMidpts = readStatsFromDir(statDir, winSize)

    overallStatMeans = {}
    for statName in statVecs:
        currStats = []
        for i in range(len(subWinMidpts)):
            currStats += [x[i] for x in statVecs[statName]]
        overallStatMeans[statName] = np.nanmean(currStats)
    return overallStatMeans

def formatStatName(name):
    if name == "pi":
        return r"$\pi$"
    elif name == "thetaW":
        return r"$\theta_W$"
    elif name == "thetaH":
        return r"$\theta_H$"
    elif name == "tajD":
        return r"Tajima's $D$"
    elif name == "fayWuH":
        return r"Fay and Wu's $H$"
    elif name == "maxFDA":
        return "max(DAF)"
    elif name == "HapCount":
        return "# of distinct haplotpyes"
    elif name == "H12":
        return r"$H_{12}$"
    elif name == "H2/H1":
        return r"$H_{2}/H_{1}$"
    elif name == "ZnS":
        return r"$Z_{nS}$"
    elif name == "Omega":
        return r"Kim and Nielsen's $\omega$"
    elif name == "distVar":
        return r"Var[div]"
    elif name == "distSkew":
        return r"Skew[div]"
    elif name == "distKurt":
        return r"Kurt[div]"
    else:
        return name

def getStatYLimsDrosophila(name, minVal, maxVal):
    if name in ["thetaW", "thetaH", "H12"]:
        buf = maxVal
        return (0, min(1, maxVal+buf))
    elif name == "pi":
        return (0, 0.005)
    elif name == "tajD":
        return (min(-3, minVal), 1)
    elif name == "maxFDA":
        return (0.75, 1.0)
    elif name == "H2/H1":
        #buf = (maxVal-minVal)*0.2
        #return (minVal-buf, 1.0)
        return (0.05, 1.0)
    elif name == "HapCount":
        return (0, 100)
    elif name in ["fayWuH"]:
        return (-0.002, 0.006)
    elif name == "distSkew":
        return (-1, 4)
    elif name == "distKurt":
        return (-2, 15)
    elif name == "distVar":
        return (0, 1e-5)
    elif name == "ZnS":
        return (0, 0.25)
    elif name in ["Omega"]:
        buf = maxVal*0.5
        return (0, maxVal+buf)
    else:
        raise Exception

def getStatYLimsHuman(name, minVal, maxVal):
    if name in ["pi", "thetaW", "thetaH"]:
        buf = maxVal*0.2
        return (0, min(1, maxVal+buf))
    if name in ["H12"]:
        buf = maxVal*0.5
        return (0, min(1, maxVal+buf))
    elif name == "ZnS":
        return (0, 0.25)
    elif name == "distVar":
        return (0, 2.5e-7)
    elif name == "tajD":
        return (min(-2, minVal), max(2, maxVal))
    elif name == "maxFDA":
        return (0.5, 1.0)
    elif name == "H2/H1":
        return (0.15, 1.0)
    elif name == "HapCount":
        #buf = (maxVal-minVal)*1.75
        #return (minVal-buf, 100)
        return (30, 100)
    elif name in ["fayWuH"]:
        #buf = (maxVal-minVal)*0.5
        #return (minVal-buf, maxVal+buf)
        return (-0.0001, 0.00035)
    elif name in ["distSkew"]:
        return (-0.5, 2)
    elif name in ["distKurt"]:
        return (-1, 4.5)
        #buf = maxVal-minVal
        #return (minVal-buf*0.5, maxVal+buf*0.75)
    elif name in ["Omega"]:
        buf = maxVal*0.5
        return (0, max(maxVal+buf, 20))
    else:
        raise Exception

def running_mean(x, N): #nice code stolen from stackoverflow somewhere
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)

def trimMidpoints(mpts, N):
    newMpts = []
    for i in range(len(mpts)-N+1):
        newMpts.append((mpts[i]+mpts[i+N-1])/2)
    return newMpts

def directionOpp(d):
    if d == "left":
        return "right"
    elif d == "right":
        return "left"
    else:
        raise Exception

def overallMinAndMax(stats, statName):
    mins = []
    maxes = []
    for selType in stats:
        mins.append(min(stats[selType][statName]))
        maxes.append(max(stats[selType][statName]))
    return (min(mins), max(maxes))

def readStatsFromFile(fileName, winSize):
    if fileName.endswith("gz"):
        fopen = gzip.open
    else:
        fopen = open
    currStatVec = {}
    with fopen(fileName) as f:
        first = True
        subWinMidpts = []
        for line in f:
            if first:
                header = line.strip().split()
                first = False
                for i in range(1, len(header)):
                    currStatVec[header[i]] = []
            else:
                line = line.strip().split()
                subWinIndex = int(line[0])
                subWinMidpts.append(subWinIndex*winSize + winSize/2)
                for i in range(1, len(header)):
                    currStatVec[header[i]].append(float(line[i]))
    return currStatVec, subWinMidpts

def readStatMeansFromBaseDir(baseDir, demogScenarioName, winSize):
    statMeans = {}
    prevSubWinMidpts = None
    for subDir in os.listdir(baseDir):
        statVecs, subWinMidpts = readStatsFromDir(baseDir + "/" + subDir + "/stats/" + demogScenarioName, winSize)
        if prevSubWinMidpts:
            assert prevSubWinMidpts == subWinMidpts
        prevSubWinMidpts = subWinMidpts
        statMeans[subDir] = {}
        for statName in statVecs:
            statMeans[subDir][statName] = []
            for i in range(len(subWinMidpts)):
                statMeans[subDir][statName].append(np.mean([x[i] for x in statVecs[statName]]))
    return statMeans, subWinMidpts
