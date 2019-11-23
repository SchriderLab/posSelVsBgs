#!/bin/bash

dataDir=simData
popLs=(humanEquilib sheehanSong tennessenEuro)
#popLs=(tennessenEuro)
for ((i=0; i<${#popLs[*]}; i++));
do
    pop=${popLs[$i]}
    for simType in trainingData testData;
    do
        simDir=$dataDir/sweeps/$pop/$simType
        for fileName in `ls $simDir | grep -v log`;
        do
            outStatsDir=$dataDir/sweeps/$pop/${simType}NiceStats
            fvecDir=$dataDir/sweeps/$pop/${simType}FeatureVecs
            fvecFileName=$fvecDir/$fileName
            logDir=$dataDir/sweeps/$pop/${simType}FeatureVecLogs
            logFileName=$logDir/$fileName

            fvecFileName=${fvecFileName/.msOut/.fvec}
            fvecFileName=${fvecFileName/.gz/}
            logFileName=${logFileName/.msOut/.log}
            logFileName=${logFileName/.gz/}

            mkdir -p $outStatsDir $fvecDir $logDir
            cmd="python ~/diploSHIC/diploSHIC.py fvecSim --totalPhysLen 110000 --outStatsDir $outStatsDir/ haploid $simDir/$fileName $fvecFileName"
            python runCmdAsJob.py "$cmd" fvecSim fvecSim.txt 12:00:00 general 4000 $logFileName
        done
    done
done
