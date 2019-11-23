#!/bin/bash

dataDir=simData
selScenarioLs=(realBgsWeakCNC realBgs realNeut absurdBgs)
popLs=(sheehanSong humanEquilib tennessenEuro)
dfeLs=(flyHubDFE boykoDFE boykoDFE)
numReps=1000
for ((i=0; i<${#popLs[*]}; i++));
do
    pop=${popLs[$i]}
    dfe=${dfeLs[$i]}
    for ((j=0; j<${#selScenarioLs[*]}; j++));
    do
        selScenario=${selScenarioLs[$j]}
        fvecDir=$dataDir/$dfe/randomExamples/$selScenario/shicFeatureVecs
        logDir=$dataDir/$dfe/randomExamples/$selScenario/shicFeatureVecLogs
        msOutDir=$dataDir/$dfe/randomExamples/$selScenario/msOut/$pop
        outStatsDir=$dataDir/$dfe/randomExamples/$selScenario/shicStats/$pop
        fvecFileName=$fvecDir/$pop.fvec
        logFileName=$logDir/$pop.log

        mkdir -p $outStatsDir $fvecDir $logDir
        cmd="python calcShicStatsForDir.py $msOutDir/ $numReps $outStatsDir/ $fvecFileName"
        python runCmdAsJob.py "$cmd" bgsStats bgsStats.txt 120:00 general 4000 $logFileName
    done
done
