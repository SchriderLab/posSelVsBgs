#!/bin/bash

popLs=(humanEquilib sheehanSong tennessenEuro)
dfeLs=(boykoDFE flyHubDFE boykoDFE)
winSizesSmall=(1000 100 1000)
for ((i=0; i<${#popLs[*]}; i++));
do
    pop=${popLs[$i]}
    dfe=${dfeLs[$i]}
    winSizeSmall=${winSizesSmall[$i]}
    winSize=$((winSizeSmall * 100))
    echo $pop
    baseDir=simData/$dfe/randomExamples/
    plotDir=plots/standard/$dfe/$pop
    mkdir -p $plotDir

    python plotBigStats.py $baseDir simData/sweeps/ $pop $winSize $plotDir/bigStats.pdf
    python plotStatsFullFig.py $baseDir simData/sweeps/ $pop $winSizeSmall 10 $plotDir/smallStatsFullFig.pdf
done
