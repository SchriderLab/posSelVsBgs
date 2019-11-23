#!/bin/bash

popSizeRescaleFactor=0.1
prespecifiedCoords=None
L=1100000
meanMu=1.2e-8
meanR=1e-8
dfe=boykoDFE
dfeMean=-0.030
dfeShape=0.206
delMutFracCoding=0.75
delMutFracCnc=0.75
annotDir=humanAnnotations/
for pop in tennessenEuro humanEquilib;
do
    discoalLaunchFile=discoalCmds/$pop.sh

    selScenarioLs=(realBgsWeakCNC realBgs realNeut absurdBgs)
    geneAnnotFileLs=($annotDir/refSeqAnnot.hg19.ucsc.12102018.gtf.gz $annotDir/refSeqAnnot.hg19.ucsc.12102018.gtf.gz fakeAnnotations/noSelRegion.gtf fakeAnnotations/centralSelRegion.1100kb.gtf)
    cncAnnotFileLs=($annotDir/phastConsElements100way_hg19.txt.gz $annotDir/phastConsElements100way_hg19.txt.gz fakeAnnotations/fakeEmptyCNCs.txt fakeAnnotations/fakeEmptyCNCs.txt)
    recRateFileLs=($annotDir/Kong2010_SexAveraged.wig $annotDir/Kong2010_SexAveraged.wig fakeAnnotations/constantRecRate.wig fakeAnnotations/constantRecRate.wig)
    chrLenFileLs=($annotDir/chromolens_hg19_autosOnly.txt $annotDir/chromolens_hg19_autosOnly.txt fakeAnnotations/fakeChrLens.1100kb.txt fakeAnnotations/fakeChrLens.1100kb.txt)
    gapFileLs=($annotDir/gaps.bed $annotDir/gaps.bed None None)
    cncSelRatioLs=(0.1 1.0 1.0 1.0)

    for ((i=0; i<${#selScenarioLs[*]}; i++));
    do
        selScenario=${selScenarioLs[$i]}
        dataDir=simData/$dfe/randomExamples/$selScenario
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
        cmd="python runBgsSimFromDiscoalCmdFileUsingAnnotAndRecMap.py $discoalLaunchFile $prespecifiedCoords $chrLenFile $gapFile $geneAnnotFile $cncAnnotFile $recRateFile $L $dfeMean $dfeShape $cncSelRatio $meanMu $meanR $delMutFracCoding $delMutFracCnc $popSizeRescaleFactor $statsDir/$pop/${filePrefix}.stats $smallStatsDir/$pop/${filePrefix}.stats $fvecDir/$pop/${filePrefix}.fvec \$SLURM_ARRAY_TASK_ID > $msOutDir/$pop/${filePrefix}.msOut"
        python runCmdAsJobArray.py "$cmd" bgsRep bgsRep.txt 120:00 general 4000 $logDir/$pop/${logPrefix}.log 0-999
    done
done
