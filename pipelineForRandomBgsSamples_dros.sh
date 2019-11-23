#!/bin/bash

popSizeRescaleFactor=0.01
prespecifiedCoords=None
L=110000
meanMu=5e-9
meanR=2.3e-8
dfe=flyHubDFE
dfeMean=-0.000133
dfeShape=0.35
delMutFracCoding=0.75
delMutFracCnc=0.75
annotDir=drosophilaAnnotations/
for pop in sheehanSong
do
    discoalLaunchFile=discoalCmds/$pop.sh

    selScenarioLs=(realBgsWeakCNC realBgs realNeut absurdBgs)
    geneAnnotFileLs=($annotDir/refSeqAnnot.dm3.ucsc.12212018.gtf.gz $annotDir/refSeqAnnot.dm3.ucsc.12212018.gtf.gz fakeAnnotations/noSelRegion.gtf fakeAnnotations/centralSelRegion.110kb.gtf)
    cncAnnotFileLs=($annotDir/phastConsElements15way_dm3.txt.gz $annotDir/phastConsElements15way_dm3.txt.gz fakeAnnotations/fakeEmptyCNCs.txt fakeAnnotations/fakeEmptyCNCs.txt)
    recRateFileLs=($annotDir/comeron.allMajArms.normalized.wig $annotDir/comeron.allMajArms.normalized.wig fakeAnnotations/constantRecRate_110kb.wig fakeAnnotations/constantRecRate_110kb.wig)
    chrLenFileLs=($annotDir/chromolens_dm3_autosOnly.txt $annotDir/chromolens_dm3_autosOnly.txt fakeAnnotations/fakeChrLens.110kb.txt fakeAnnotations/fakeChrLens.110kb.txt)
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
        python runCmdAsJobArray.py "$cmd" bgsRep bgsRep.txt 12:00:00 general 4000 $logDir/$pop/${logPrefix}.log 0-999
    done
done
