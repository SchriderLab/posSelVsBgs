#!/bin/bash

dataDir=simData
popLs=(humanEquilib tennessenEuro sheehanSong)
LLs=(1100000 1100000 110000)
for ((i=0; i<${#popLs[*]}; i++));
do
    pop=${popLs[$i]}
    L=${LLs[$i]}
    for simType in trainingData testData;
    do
        simDir=$dataDir/sweeps/$pop/$simType
        bigStatsDir=$dataDir/sweeps/$pop/${simType}FwdpyStats
        smallStatsDir=$dataDir/sweeps/$pop/${simType}FwdpySmallStats
        fvecDir=$dataDir/sweeps/$pop/${simType}FwdpyFeatureVecs
        logDir=$dataDir/sweeps/$pop/${simType}FwdpyStatLogs
        mkdir -p $logDir
        for fileName in simNeut.msOut.gz simHard_5.msOut.gz simSoft_5.msOut.gz;
        do
            prefix=${fileName/.msOut/}
            prefix=${prefix/.gz/}
            logFileName=$logDir/$prefix.log

            mkdir -p $bigStatsDir/$prefix $smallStatsDir/$prefix $fvecDir/$prefix
            cmd="python calcStatsForMsOutFile.py $simDir/$fileName $L $bigStatsDir/$prefix/$prefix $smallStatsDir/$prefix/$prefix $fvecDir/$prefix/$prefix"
            python runCmdAsJob.py "$cmd" shicSbStats shicSbStats.txt 96:00:00 general 4000 $logFileName
        done
    done
done
