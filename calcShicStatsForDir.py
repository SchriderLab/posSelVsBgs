import sys, os

msOutDir, numReps, outStatsDir, fvecFileName = sys.argv[1:]
numReps=int(numReps)

def readLines(fname, skip=0):
    with open(fname, "rt") as f:
        lines = f.readlines()[skip:]
    return lines

def writeAllMsOutInDirToFile(msOutDir, numReps, outFileName):
    msFileNames = os.listdir(msOutDir)
    msFileNames.sort(key=lambda x: int(x.split("_")[1].split(".msOut")[0]))
    with open(outFileName, "wt") as outFile:
        outFile.write("fwdpy 100 %s 110000 etc\nblah\n" %(numReps))
        for i in range(len(msFileNames)):
            for line in readLines(msOutDir + "/" + msFileNames[i], skip=2):
                outFile.write(line.strip() + "\n")

tempFileName = fvecFileName + ".tmp.msOut"
print("creating tmp ms out file")
writeAllMsOutInDirToFile(msOutDir, numReps, tempFileName)
cmd = "python ~/diploSHIC/diploSHIC.py fvecSim --totalPhysLen 110000 --outStatsDir %s haploid %s %s" %(outStatsDir, tempFileName, fvecFileName)
print("calculating stats")
os.system(cmd)
print("Done! Deleting tmp ms out file and peacing out.")
os.system("rm %s" %(tempFileName))
