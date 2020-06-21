#!/bin/bash

dataDir=simData
popLs=(humanEquilib sheehanSong tennessenEuro)
dfeLs=(boykoDFE flyHubDFE boykoDFE)
for ((i=0; i<${#popLs[*]}; i++));
do
    pop=${popLs[$i]}
    dfe=${dfeLs[$i]}
    dom00=$dataDir/$dfe/dominanceRealBgsRandom/dominance_0.0/shicFeatureVecs/$pop.fvec
    dom10=$dataDir/$dfe/dominanceRealBgsRandom/dominance_0.5/shicFeatureVecs/$pop.fvec
    dom20=$dataDir/$dfe/dominanceRealBgsRandom/dominance_1.0/shicFeatureVecs/$pop.fvec
    dom05=$dataDir/$dfe/randomExamples/realBgs/shicFeatureVecs/$pop.fvec
    outDir=$dataDir/$dfe/dominanceTestResults/
    mkdir -p $outDir

    modelDir=$dataDir/sweeps/$pop/classifier/
    outputModel=$modelDir/clf.mod

    cmd="python testShicOnDominanceSims.py $outputModel.json $outputModel.weights.hdf5 $dom00 $dom05 $dom10 $dom20 $outDir/$pop.out"
    python runCmdAsJob.py "$cmd" domTest domTest.txt 120:00 general 4000 $outDir/$pop.log
done
