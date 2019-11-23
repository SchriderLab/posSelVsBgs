import sys
import summStatFuncs, msTools

msOutFileName, L, bigStatFilePrefix, smallStatFilePrefix, fvecFilePrefix = sys.argv[1:]
L = int(L)

statNames = ["pi", "thetaW", "tajD", "thetaH", "fayWuH", "maxFDA", "HapCount", "H1", "H12", "H2/H1", "ZnS", "Omega", "distVar", "distSkew", "distKurt"]
smallStatNames = ["pi", "thetaW", "thetaH", "tajD", "fayWuH", "ZnS", "Omega"]

trainingDataFileObj, sampleSize, numInstances = msTools.openMsOutFileForSequentialReading(msOutFileName)

for repI in range(numInstances):
    gameteStrs, positionArray = msTools.readNextMsRepToGameteStrs(trainingDataFileObj, sampleSize, L, discretizePositions=False)
    positionArray = [x*L for x in positionArray]
    windowedStats = summStatFuncs.calculateWindowedStats(positionArray, gameteStrs, L, statNames)

    summStatFuncs.writeStats(windowedStats, statNames, bigStatFilePrefix + "_%d.stats" %(repI))
    smallWindowedStats = summStatFuncs.calculateWindowedStats(positionArray, gameteStrs, L, smallStatNames, numSubWins=1100)
    summStatFuncs.writeStats(smallWindowedStats, smallStatNames, smallStatFilePrefix + "_%d.stats" %(repI))
    if not fvecFilePrefix.lower() in ["none", "false"]:
        summStatFuncs.writeFvec(windowedStats, statNames, fvecFilePrefix + "_%d.fvec" %(repI))

    sys.stderr.write("Done %d reps\n" %(repI))

msTools.closeMsOutFile(trainingDataFileObj)
