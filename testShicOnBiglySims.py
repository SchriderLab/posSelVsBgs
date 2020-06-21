import sys, os, uuid
import numpy as np

classToIndex = {'hard':0, 'linkedHard':1, 'soft':2, 'linkedSoft':3, 'neutral':4}
classLabels = ['hard', 'linkedHard', 'soft', 'linkedSoft', 'neutral']

def generatePredFileName(tmpFilePrefix):
    return tmpFilePrefix + "." + uuid.uuid4().hex + ".tmp"

def readPredsIntoVector(tmpPredFileName):
    global classToIndex
    currPreds = np.loadtxt(tmpPredFileName, skiprows=1, dtype="str")
    predFracs = [0]*len(classToIndex)
    for predLine in currPreds:
        predFracs[classToIndex[predLine[4]]] += 1
    denom = sum(predFracs)
    for i in range(len(predFracs)):
        predFracs[i] /= denom
    return predFracs, currPreds

def writeTmpFvecInForPreds(headerLine, fvecLines, tmpFvecInFileName):
    with open(tmpFvecInFileName, "wt") as wFile:
        wFile.write("chrom\tclassifiedWinStart\tclassifiedWinEnd\tbigWinRange\t" + headerLine)
        for fvecLine in fvecLines:
            wFile.write("chrom\tclassifiedWinStart\tclassifiedWinEnd\tbigWinRange\t" + fvecLine)

def readAllFvecsForRepInOrder(fvecPrefixForRep):
    fvecPath = fvecPrefixForRep.split("/")
    prefix = fvecPath[-1]
    d = "/".join(fvecPath[:-1])

    fvecLs = []
    for fname in os.listdir(d):
        if fname.startswith(prefix):
            if (fname.endswith(".tmp") or fname.endswith(".tmp.input")) and not "central" in fname:
                sys.stderr.write("removing {}".format(fname) + "\n")
                os.system("rm " + d + "/" + fname)
            else:
                rep, winHeader, winIndex, extension = fname.split(".")
                winIndex = int(winIndex)

                currFvecLines = []
                with open(d + "/" + fname, "rt") as f:
                    for line in f:
                        currFvecLines.append(line)

                fvecLs.append((winIndex, currFvecLines))
    fvecLs.sort()

    headerLine = fvecLs[0][1][0]
    fvecLines = []
    for winIndex, currFvecLines in fvecLs:
        fvecLines.extend(currFvecLines[1:])
    return headerLine, fvecLines

def predictAllWindowsForRep(modelArch, modelWeights, fvecPrefixForRep):
    headerLine, fvecLines = readAllFvecsForRepInOrder(fvecPrefixForRep)

    tmpPredFileName = generatePredFileName(fvecPrefixForRep)
    tmpFvecInFileName = tmpPredFileName + ".input"
    writeTmpFvecInForPreds(headerLine, fvecLines, tmpFvecInFileName)
    cmd = "python ~/diploSHIC/diploSHIC.py predict %s %s %s %s" %(modelArch, modelWeights, tmpFvecInFileName, tmpPredFileName)
    os.system(cmd)
    predFracs, preds = readPredsIntoVector(tmpPredFileName)
    os.system("rm %s %s" %(tmpPredFileName, tmpFvecInFileName))
    return predFracs, preds

def main():
    global classLabels, classToIndex
    modelArch, modelWeights, fvecPrefix, firstRep, lastRep,  outFileName = sys.argv[1:]

    headerLine = "predClass\tprob(neutral)\tprob(likedSoft)\tprob(linkedHard)\tprob(soft)\tprob(hard)"
    with open(outFileName, "w") as outFile:
        outFile.write("repIndex\twinIndex" + "\t" + headerLine + "\n")

        for i in range(int(firstRep), int(lastRep)+1):
            fvecPrefixForRep = fvecPrefix + "_{}.win".format(i)
            predData = predictAllWindowsForRep(modelArch, modelWeights, fvecPrefixForRep)
            for j in range(len(predData[1])):
                outFile.write(str(i) + "\t" + str(j) + "\t" + "\t".join([str(x) for x in predData[1][j][4:]]) + "\n")

if __name__ == "__main__":
    main()
