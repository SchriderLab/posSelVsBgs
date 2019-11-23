#!/bin/bash

popSizeRescaleFactor=0.1
L=1100000
meanMu=1.2e-8
meanR=1e-8
dfe=boykoDFE
dfeMean=-0.030
dfeShape=0.206
delMutFracCoding=0.75
delMutFracCnc=0.75
dominance=0.5

for pop in tennessenEuro humanEquilib;
do
    discoalLaunchFile=discoalCmds/$pop.sh
    annotDir=humanAnnotations/
    geneAnnotFile=$annotDir/refSeqAnnot.hg19.ucsc.12102018.gtf.gz
    cncAnnotFile=$annotDir/phastConsElements100way_hg19.txt.gz
    recRateFile=$annotDir/Kong2010_SexAveraged.wig
    chrLenFile=$annotDir/chromolens_hg19_autosOnly.txt
    gapFile=$annotDir/gaps.bed
    repeatedCoordsFile=repeatedCoords/human_hg19.txt
    declare -a simCoordsLs
    readarray -t simCoordsLs < $repeatedCoordsFile
    cncSelRatio=1.0
    for ((i=0; i<${#simCoordsLs[*]}; i++));
    do
        simCoords=${simCoordsLs[$i]}
        dataDir=simData/$dfe/repeatedExamples/$simCoords
        statsDir=$dataDir/stats
        smallStatsDir=$dataDir/smallStats
        msOutDir=$dataDir/msOut
        logDir=$dataDir/simLogs
        fvecDir=$dataDir/featureVectors
        mkdir -p $statsDir/$pop $smallStatsDir/$pop $msOutDir/$pop $logDir/$pop $fvecDir/$pop

        filePrefix="${pop}_\$SLURM_ARRAY_TASK_ID"
        logPrefix="${pop}_%a"
        cmd="python runBgsSimFromDiscoalCmdFileUsingAnnotAndRecMap.py $discoalLaunchFile $simCoords $chrLenFile $gapFile $geneAnnotFile $cncAnnotFile $recRateFile $L $dfeMean $dfeShape $dominance $cncSelRatio $meanMu $meanR $delMutFracCoding $delMutFracCnc $popSizeRescaleFactor $statsDir/$pop/${filePrefix}.stats $smallStatsDir/$pop/${filePrefix}.stats $fvecDir/$pop/${filePrefix}.fvec \$SLURM_ARRAY_TASK_ID > $msOutDir/$pop/${filePrefix}.msOut"
        python runCmdAsJobArray.py "$cmd" bgsRep bgsRep.txt 12:00:00 general 4000 $logDir/$pop/${logPrefix}.log 0-999
    done
done
