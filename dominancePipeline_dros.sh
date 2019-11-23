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

    geneAnnotFile=$annotDir/refSeqAnnot.dm3.ucsc.12212018.gtf.gz
    cncAnnotFile=$annotDir/phastConsElements15way_dm3.txt.gz
    recRateFile=$annotDir/comeron.allMajArms.normalized.wig
    chrLenFile=$annotDir/chromolens_dm3_autosOnly.txt
    gapFile=$annotDir/gaps.bed
    cncSelRatio=1.0

    for dominance in 0.0 1.0 2.0;
    do
        dataDir=simData/$dfe/dominanceRealBgsRandom/dominance_$dominance
        statsDir=$dataDir/stats
        smallStatsDir=$dataDir/smallStats
        msOutDir=$dataDir/msOut
        logDir=$dataDir/simLogs
        fvecDir=$dataDir/featureVectors
        mkdir -p $statsDir/$pop $smallStatsDir/$pop $msOutDir/$pop $logDir/$pop $fvecDir/$pop

        filePrefix="${pop}_\$SLURM_ARRAY_TASK_ID"
        logPrefix="${pop}_%a"
        cmd="python runBgsSimFromDiscoalCmdFileUsingAnnotAndRecMap.py $discoalLaunchFile $prespecifiedCoords $chrLenFile $gapFile $geneAnnotFile $cncAnnotFile $recRateFile $L $dfeMean $dfeShape $dominance $cncSelRatio $meanMu $meanR $delMutFracCoding $delMutFracCnc $popSizeRescaleFactor $statsDir/$pop/${filePrefix}.stats $smallStatsDir/$pop/${filePrefix}.stats $fvecDir/$pop/${filePrefix}.fvec \$SLURM_ARRAY_TASK_ID > $msOutDir/$pop/${filePrefix}.msOut"
        python runCmdAsJobArray.py "$cmd" domRep domRep.txt 12:00:00 general 4000 $logDir/$pop/${logPrefix}.log 0-999
    done
done
