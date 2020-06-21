#!/bin/bash

popSizeRescaleFactor=0.01
prespecifiedCoords=None
totalL=1210000
meanMu=5e-9
meanR=2.3e-8
dfe=flyHubDFE
dfeMean=-0.000133
dfeShape=0.35
dominance=0.25
delMutFracCoding=0.75
delMutFracCnc=0.75
annotDir=drosophilaAnnotations/
winSize=110000
stepSize=1000
for pop in sheehanSong
do
    discoalLaunchFile=discoalCmds/$pop.sh

    selScenarioLs=(realBgs)
    geneAnnotFileLs=($annotDir/refSeqAnnot.dm3.ucsc.12212018.gtf.gz)
    cncAnnotFileLs=($annotDir/phastConsElements15way_dm3.txt.gz)
    recRateFileLs=($annotDir/comeron.allMajArms.normalized.wig)
    chrLenFileLs=($annotDir/chromolens_dm3_autosOnly.txt)
    gapFileLs=($annotDir/gaps.bed)
    cncSelRatioLs=(1.0)

    for ((i=0; i<${#selScenarioLs[*]}; i++));
    do
        selScenario=${selScenarioLs[$i]}
        dataDir=simData/$dfe/randomExamplesBigly/$selScenario
        geneAnnotFile=${geneAnnotFileLs[$i]}
        cncAnnotFile=${cncAnnotFileLs[$i]}
        recRateFile=${recRateFileLs[$i]}
        chrLenFile=${chrLenFileLs[$i]}
        gapFile=${gapFileLs[$i]}
        cncSelRatio=${cncSelRatioLs[$i]}
        statsDir=$dataDir/stats
        smallStatsDir=$dataDir/smallStats
        msOutDir=$dataDir/msOut
        logDir=$dataDir/simLogs
        fvecDir=$dataDir/featureVectors
        mkdir -p $statsDir/$pop $smallStatsDir/$pop $msOutDir/$pop $logDir/$pop $fvecDir/$pop

        filePrefix="${pop}_\$SLURM_ARRAY_TASK_ID"
        logPrefix="${pop}_%a"
        cmd="python runBgsSimFromDiscoalCmdFileUsingAnnotAndRecMapFixedRatesSlidingWins.py $discoalLaunchFile $prespecifiedCoords $winSize $stepSize $chrLenFile $gapFile $geneAnnotFile $cncAnnotFile $recRateFile $totalL $dfeMean $dfeShape $dominance $cncSelRatio $meanMu $meanR $delMutFracCoding $delMutFracCnc $popSizeRescaleFactor $statsDir/$pop/${filePrefix} $fvecDir/$pop/${filePrefix} \$SLURM_ARRAY_TASK_ID"
        python runCmdAsJobArray.py "$cmd" bgsRepBigly bgsRepBigly.txt 10-00:00:00 general 24G $logDir/$pop/${logPrefix}.log 0-99
    done
done
