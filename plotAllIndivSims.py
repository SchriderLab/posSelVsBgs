import os
import runCmdAsJob

def demogScenarioToDfe(demogScenario):
    if demogScenario in ["tennessenEuro", "humanEquilib"]:
        return "boykoDFE"
    elif demogScenario == "sheehanSong":
        return "flyHubDFE"
    else:
        raise ValueError

baseDir = "simData"

demogScenarioLs = [("humanEquilib", 100000), ("tennessenEuro", 100000), ("sheehanSong", 10000)]
coalTypes = ["simHard_5", "simSoft_5"]
fwdTypeLs = ["absurdBgs", "realBgsWeakCNC", "realBgs", "realNeut"]

os.system("mkdir -p {}/indivPlotLogs".format(baseDir))

for demogScenario, winSize in demogScenarioLs:
    dfe = demogScenarioToDfe(demogScenario)
    for coalType in coalTypes:
        figOutDir = "{}/indivPlots/{}/{}".format(baseDir, demogScenario, coalType)
        cmd = "mkdir -p {}".format(figOutDir)
        print(cmd)
        os.system(cmd)

        subBaseDir = "{}/{}/randomExamples/".format(baseDir, dfe)
        statPrefix = "{}/sweeps/{}/testDataFwdpyStats/{}/{}_".format(baseDir, demogScenario, coalType, coalType)
        figPrefix = "{}.{}".format(demogScenario, coalType)
        cmd = "python plotBigStatsOneSim.py {} {} {} {}-{} None {}/{}".format(statPrefix, demogScenario, winSize, 0, 1000, figOutDir, figPrefix)
        print(cmd)
        logFileName = "{}/indivPlotLogs/{}.{}.log".format(baseDir, demogScenario, coalType)
        runCmdAsJob.runCmdAsJobWithoutWaitingWithLog(cmd, "pltIndv", "pltIndv.txt", "12:00:00", "general", "2G", logFileName)

    for fwdType in fwdTypeLs:
        if fwdType == "realBgs":
            domStrLs = ["0.0", "0.25", "0.5", "1.0"]
        else:
            domStrLs = ["0.25"]
        for domStr in domStrLs:
            figOutDir = "{}/indivPlots/{}/{}_{}_{}".format(baseDir, demogScenario, dfe, fwdType, domStr)
            cmd = "mkdir -p {}".format(figOutDir)
            print(cmd)
            os.system(cmd)

            if domStr == "0.25":
                subBaseDir = "{}/{}/randomExamples/".format(baseDir, dfe)
                statPrefix = "{}/{}/randomExamples/{}/stats/{}/{}_".format(baseDir, dfe, fwdType, demogScenario, demogScenario)
                logPrefix = "{}/{}/randomExamples/{}/simLogs/{}/{}_".format(baseDir, dfe, fwdType, demogScenario, demogScenario)
            else:
                subBaseDir = "{}/{}/dominanceRealBgsRandom/dominance_{}/shicFeatureVecs/".format(baseDir, dfe, float(domStr)*2)
                statPrefix = "{}/{}/dominanceRealBgsRandom/dominance_{}/stats/{}/{}_".format(baseDir, dfe, float(domStr)*2, demogScenario, demogScenario)
                logPrefix = "{}/{}/dominanceRealBgsRandom/dominance_{}/simLogs/{}/{}_".format(baseDir, dfe, float(domStr)*2, demogScenario, demogScenario)

            if not fwdType.startswith("realBgs"):
                logPrefix = "None"

            figPrefix = "{}.{}.{}.{}".format(demogScenario, dfe, fwdType, domStr)
            cmd = "python plotBigStatsOneSim.py {} {} {} {}-{} {} {}/{}".format(statPrefix, demogScenario, winSize, 0, 1000, logPrefix, figOutDir, figPrefix)
            print(cmd)
            logFileName = "{}/indivPlotLogs/{}.{}.{}.log".format(baseDir, demogScenario, fwdType, domStr)
            runCmdAsJob.runCmdAsJobWithoutWaitingWithLog(cmd, "pltIndv", "pltIndv.txt", "12:00:00", "general", "2G", logFileName)
