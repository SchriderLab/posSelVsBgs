import sys

def parseDiscoalArgs(discoalCmd):
    relSizeChanges = []
    for argStr in discoalCmd.split(" -")[1:]:
        argStr = argStr.split()
        argFlag = argStr[0]
        if argFlag == "en":
            t, pop, size = argStr[1:]
            assert pop == "0"
            relSizeChanges.append((float(t), float(size)))
        elif argFlag == "Pt":
            thetaMin, thetaMax = [float(x) for x in argStr[1:]]
            thetaRange = [thetaMin, thetaMax]
        elif argFlag == "t":
            theta = float(argStr[1])
            thetaRange = [theta, theta]
        elif argFlag == "Pre":
            rhoMean = float(argStr[1])
            rhoMax = float(argStr[2])
        else:
            raise NotImplementedError
    return thetaRange, rhoMean, rhoMax, relSizeChanges

def parseDiscoalCmd(discoalCmd, u, r, L):
    splitCmd = discoalCmd.split()
    progName, sampleSize, numReps, numSites = splitCmd[:4]
    sampleSize = int(sampleSize)
    thetaRange, rhoMean, rhoMax, relSizeChanges = parseDiscoalArgs(discoalCmd)
    thetaMean = (thetaRange[0]+thetaRange[-1])/2.0
    N0 = thetaMean/(4*u*L)
    N02 = rhoMean/(4*r*L)
    assert abs(N0-N02) < 1
    prevT = 0
    popSizeChanges = []
    prevSize=int(round(N0))
    for t, size in sorted(relSizeChanges):
        size = int(round(size*N0))
        t = int(round(t*4*N0))
        interval = 1
        for gen in range(prevT, t, interval):
            popSizeChanges.append(prevSize)
        prevT = t
        prevSize = size
    return sampleSize, prevSize, popSizeChanges, thetaRange, rhoMean, rhoMax

def getParamBoundsAndPopSizeChangesFromDiscoalCmdFile(discoalCmdFileName, u, r, L, burnTime):
    with open(discoalCmdFileName) as discoalCmdFile:
        for line in discoalCmdFile:
            if line.startswith("python"):
                discoalCmd = line.strip().split('"')[1]
            elif line.startswith("discoal"):
                discoalCmd = line.strip()
        sampleSize, ancSize, popSizeChanges, thetaRange, rhoMean, rhoMax = parseDiscoalCmd(discoalCmd, u, r, L)
        nlist = [ancSize]*(burnTime*ancSize) + popSizeChanges[::-1]
        sys.stderr.write("ancestral population size: %d; burn length: %d; total num gens: %d; present-day pop size: %d \n" %(ancSize, burnTime*ancSize, len(nlist), nlist[-1]))
        return sampleSize, nlist, thetaRange, rhoMean, rhoMax

def main(): #for debugging only
    fileName = sys.argv[1]
    params = getParamBoundsAndPopSizeChangesFromDiscoalCmdFile(fileName, 1.2e-07, 1e-07, 1100000, 10)
    prevN = params[1][-1]
    count = 1
    for n in params[1]:
        print(n, count)
        count += 1
        prevN = n

if __name__ == "__main__":
    main()
