#!/bin/bash

dataDir=simData
popLs=(sheehanSong humanEquilib tennessenEuro)
dfeLs=(flyHubDFE boykoDFE boykoDFE)
#popLs=(tennessenEuro)
#dfeLs=(boykoDFE)
numReps=1000
for ((i=0; i<${#popLs[*]}; i++));
do
    pop=${popLs[$i]}
    dfe=${dfeLs[$i]}
    for dominance in 0.0 0.5 1.0;
    do
        fvecDir=$dataDir/$dfe/dominanceRealBgsRandom/dominance_$dominance/shicFeatureVecs
        logDir=$dataDir/$dfe/dominanceRealBgsRandom/dominance_$dominance/shicFeatureVecLogs
        msOutDir=$dataDir/$dfe/dominanceRealBgsRandom/dominance_$dominance/msOut/$pop
        outStatsDir=$dataDir/$dfe/dominanceRealBgsRandom/dominance_$dominance/shicStats/$pop
        fvecFileName=$fvecDir/$pop.fvec
        logFileName=$logDir/$pop.log

        mkdir -p $outStatsDir $fvecDir $logDir
        cmd="python calcShicStatsForDir.py $msOutDir/ $numReps $outStatsDir/ $fvecFileName"
        python runCmdAsJob.py "$cmd" domStats domStats.txt 48:00:00 general 4000 $logFileName
    done
done
