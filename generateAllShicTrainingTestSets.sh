#!/bin/bash

dataDir=simData
popLs=(humanEquilib sheehanSong tennessenEuro)
for ((i=0; i<${#popLs[*]}; i++));
do
    pop=${popLs[$i]}
    for simType in trainingData testData;
    do
        fvecDir=$dataDir/sweeps/$pop/${simType}FeatureVecs
        outDir=$dataDir/sweeps/$pop/${simType}Sets/
        iogDir=$dataDir/sweeps/$pop/${simType}SetLogs/
        mkdir -p $logDir $outDir
        python ~/diploSHIC/diploSHIC.py makeTrainingSets $fvecDir/simNeut.fvec $fvecDir/simSoft $fvecDir/simHard 5 0,1,2,3,4,6,7,8,9,10 $outDir
    done
done
