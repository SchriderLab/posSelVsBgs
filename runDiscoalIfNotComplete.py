import sys, os, gzip

def incompleteDiscoalOutput(outFileName, numReps, sampleSize):
    if os.path.isfile(outFileName):
        repsFound = 0
        linesFound = 0
        if outFileName.endswith(".gz"):
            fopen = gzip.open
        else:
            fopen = open
        neededLines = 2 + numReps*(4+sampleSize)
        try:
            with fopen(outFileName, "rt") as outFile:
                for line in outFile:
                    linesFound += 1
                    if line.startswith("//"):
                        repsFound += 1
            return repsFound != numReps or linesFound != neededLines
        except Exception:
            return True
    else:
        return True

def runDiscoalIfNotComplete(discoalCmd, outFileName):
    sampleSize, numReps = [int(x) for x in discoalCmd.split()[1:3]]
    if incompleteDiscoalOutput(outFileName, numReps, sampleSize):
        if outFileName.endswith(".gz"):
            cmd = "%s | gzip > %s" %(discoalCmd, outFileName)
        else:
            cmd = "%s > %s" %(discoalCmd, outFileName)
        print("simulation results in %s are incomplete; re-running" %(outFileName))
        os.system(cmd)
    else:
        print("simulation results in %s already complete; not running" %(outFileName))

def main():
    discoalCmd = " ".join(sys.argv[1:-1])
    outFileName = sys.argv[-1]
    runDiscoalIfNotComplete(discoalCmd, outFileName)

if __name__ == "__main__":
    main()
