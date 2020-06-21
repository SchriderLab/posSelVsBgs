import sys, os, math

discoalCmdFileName, u, L, trainingOutDir, testOutDir = sys.argv[1:]
u = float(u)
L = int(L)
maxS = 0.05
minS = 0.0001

def readAndModifyDiscoalCmd(discoalCmdFileName, sampleNumber, maxNumSites=220000):
    with open(discoalCmdFileName, "rt") as discoalCmdFile:
        for line in discoalCmdFile:
            if line.startswith("discoal"):
                discoalCmd = ["discoal"] + line.strip().split()[1:]
                break
    for i in range(len(discoalCmd)):
        if discoalCmd[i] == "-t":
            theta = float(discoalCmd[i+1])
            break
        elif discoalCmd[i] == "-Pt":
            theta = (float(discoalCmd[i+1])+float(discoalCmd[i+2]))/2.0
            break
    discoalCmd[2] = str(sampleNumber)
    numSites = int(discoalCmd[3])
    if numSites > maxNumSites:
        discoalCmd[3] = str(maxNumSites)
    return " ".join(discoalCmd), theta

sampleNumber = 2000
discoalCmd, theta = readAndModifyDiscoalCmd(discoalCmdFileName, sampleNumber)

N0 = theta/(4*u*L)

alphaHigh = maxS*2*N0
alphaLow = minS*2*N0

tauLow = 0
tauHigh = tauHigh = 200/(4*N0)
selStr = " -ws 0 -Pa %f %f -Pu %f %f" %(alphaLow, alphaHigh, tauLow, tauHigh)
softStr = " -Pf 0 0.05"

sweepLocStr = " -x $x"

print("#!/bin/bash")
for outDir, simTitle in [(trainingOutDir, "training data"), (testOutDir, "test data")]:
    print("\n#generating %s\n" %(simTitle))
    print("mkdir -p %s" %(outDir))
    neutDiscoalCmd = "python runDiscoalIfNotComplete.py %s" %(discoalCmd)
    print("python runCmdAsJob.py \"%s %s/simNeut.msOut.gz\" neutSims neutSims.txt 48:00:00 general 128G %s/simNeut.msOut.log" %(neutDiscoalCmd, outDir, outDir))
    print("i=0")
    print("for x in 0.045454545454545456 0.13636363636363635 0.22727272727272727 0.3181818181818182 0.4090909090909091 0.5 0.5909090909090909 0.6818181818181818 0.7727272727272727 0.8636363636363636 0.9545454545454546;\ndo")
    hardDiscoalCmd = neutDiscoalCmd + selStr + sweepLocStr
    print("    python runCmdAsJob.py \"%s -i 40 %s/simHard_$i.msOut.gz\" hardSims hardSims.txt 48:00:00 general 128G %s/simHard_$i.msOut.log" %(hardDiscoalCmd, outDir, outDir))
    softDiscoalCmd = neutDiscoalCmd + selStr + softStr + sweepLocStr
    print("    python runCmdAsJob.py \"%s -i 40 %s/simSoft_$i.msOut.gz\" softSims softSims.txt 48:00:00 general 128G %s/simSoft_$i.msOut.log" %(softDiscoalCmd, outDir, outDir))
    print("    i=$((i + 1))\ndone")
    print("")
