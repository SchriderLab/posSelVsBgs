#!/bin/bash

dataDir=simData
discoalCmdDir=discoalCmds/
simLaunchScriptDir=simLaunchScripts/
mkdir -p $dataDir/sweeps $simLaunchScriptDir
popLs=(sheehanSong humanEquilib tennessenEuro)
uLs=(5e-9 1.2e-8 1.2e-8)
LLs=(110000 1100000 1100000)
for ((i=0; i<${#popLs[*]}; i++));
do
    pop=${popLs[$i]}
    u=${uLs[$i]}
    L=${LLs[$i]}
    discoalCmdFileName=$discoalCmdDir/$pop.sh
    trainingOutDir=$dataDir/sweeps/$pop/trainingData
    testOutDir=$dataDir/sweeps/$pop/testData
    python generateSimLaunchScript.py $discoalCmdFileName $u $L $popFileName $trainingOutDir $testOutDir > $simLaunchScriptDir/$pop.sh
done
