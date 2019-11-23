#!/bin/bash

dataDir=simData
popLs=(humanEquilib sheehanSong tennessenEuro)
dfeLs=(boykoDFE flyHubDFE boykoDFE)
for ((i=0; i<${#popLs[*]}; i++));
do
    pop=${popLs[$i]}
    dfe=${dfeLs[$i]}
    trainDir=$dataDir/sweeps/$pop/trainingDataSets/
    testDir=$dataDir/sweeps/$pop/testDataSets/
    outDir=$dataDir/sweeps/$pop/classifier/
    outputModel=$outDir/clf.mod
    logFileName=$outDir/training_testing.log
    figOutDir=plots/covfefe
    figFileName=$figOutDir/$pop.pdf
    mkdir -p $outDir $figOutDir
    cmd="python ~/diploSHIC/diploSHIC.py train $trainDir $testDir $outputModel"
    python runCmdAsJob.py "$cmd" trnShc trnShc.txt 120:00 general 4000 $logFileName
done
