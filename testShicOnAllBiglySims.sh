#!/bin/bash

dataDir=simData
popLs=(humanEquilib sheehanSong tennessenEuro)
dfeLs=(boykoDFE flyHubDFE boykoDFE)
for ((i=0; i<${#popLs[*]}; i++));
do
    pop=${popLs[$i]}
    dfe=${dfeLs[$i]}
    outDir=$dataDir/$dfe/biglyTestResults/
    mkdir -p $outDir

    modelDir=$dataDir/sweeps/$pop/classifier/
    outputModel=$modelDir/clf.mod

    cmd="python testShicOnBiglySims.py $outputModel.json $outputModel.weights.hdf5 $dataDir/$dfe/randomExamplesBigly/realBgs/featureVectors/$pop/$pop 0 99 $outDir/$pop.out"
    python runCmdAsJob.py "$cmd" biglyTest biglyTest.txt 300:00 general 4000 $outDir/$pop.log
    cmd="python testShicOnBiglySimsCentralWinOnly.py $outputModel.json $outputModel.weights.hdf5 $dataDir/$dfe/randomExamplesBigly/realBgs/featureVectors/$pop/$pop 0 99 $outDir/${pop}_central.out"
    python runCmdAsJob.py "$cmd" biglyTestCentral biglyTestCentral.txt 300:00 general 4000 $outDir/${pop}_central.log
done
