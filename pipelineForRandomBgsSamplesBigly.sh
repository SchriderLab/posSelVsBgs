#!/bin/bash

popSizeRescaleFactor=0.1
prespecifiedCoords=None
totalL=12100000
meanMu=1.2e-8
meanR=1e-8
dfe=boykoDFE
dfeMean=-0.030
dfeShape=0.206
dominance=0.25
delMutFracCoding=0.75
delMutFracCnc=0.75
annotDir=humanAnnotations/
winSize=1100000
stepSize=10000
for pop in humanEquilib tennessenEuro;
do
    discoalLaunchFile=discoalCmds/$pop.sh

    selScenarioLs=(realBgs)
    geneAnnotFileLs=($annotDir/refSeqAnnot.hg19.ucsc.12102018.gtf.gz $annotDir/refSeqAnnot.hg19.ucsc.12102018.gtf.gz fakeAnnotations/noSelRegion.gtf fakeAnnotations/centralSelRegion.1100kb.gtf)
    cncAnnotFileLs=($annotDir/phastConsElements100way_hg19.txt.gz $annotDir/phastConsElements100way_hg19.txt.gz fakeAnnotations/fakeEmptyCNCs.txt fakeAnnotations/fakeEmptyCNCs.txt)
    recRateFileLs=($annotDir/Kong2010_SexAveraged.wig $annotDir/Kong2010_SexAveraged.wig fakeAnnotations/constantRecRate.wig fakeAnnotations/constantRecRate.wig)
    chrLenFileLs=($annotDir/chromolens_hg19_autosOnly.txt $annotDir/chromolens_hg19_autosOnly.txt fakeAnnotations/fakeChrLens.1100kb.txt fakeAnnotations/fakeChrLens.1100kb.txt)
    gapFileLs=($annotDir/gaps.bed $annotDir/gaps.bed None None)
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
        msOutDir=$dataDir/msOut
        logDir=$dataDir/simLogs
        fvecDir=$dataDir/featureVectors
        mkdir -p $statsDir/$pop $msOutDir/$pop $logDir/$pop $fvecDir/$pop

        filePrefix="${pop}_\$SLURM_ARRAY_TASK_ID"
        logPrefix="${pop}_%a"
        cmd="python runBgsSimFromDiscoalCmdFileUsingAnnotAndRecMapFixedRatesSlidingWins.py $discoalLaunchFile $prespecifiedCoords $winSize $stepSize $chrLenFile $gapFile $geneAnnotFile $cncAnnotFile $recRateFile $totalL $dfeMean $dfeShape $dominance $cncSelRatio $meanMu $meanR $delMutFracCoding $delMutFracCnc $popSizeRescaleFactor $statsDir/$pop/${filePrefix} $fvecDir/$pop/${filePrefix} \$SLURM_ARRAY_TASK_ID"
        python runCmdAsJobArray.py "$cmd" bgsRepBigly bgsRepBigly.txt 96:00:00 general 24G $logDir/$pop/${logPrefix}.log 0-99
    done
done
