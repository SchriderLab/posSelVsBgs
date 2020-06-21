#!/bin/bash

popSizeRescaleFactor=0.01
L=110000
meanMu=5e-9
meanR=2.3e-8
#dfe=flyHubDFE
#dfeMean=-0.000133
#dfeShape=0.35
delMutFracCoding=0.75
delMutFracCnc=0.75
#dfe=flyAewDFE
#dfeMean=0.00305 #-1800/Ne; keightley and eyre-walker 2007; Ne is harmonic mean of terhorst model (590317.0315343805)
#dfeShape=0.38
dfe=flyHubDFE
dfeMean=-0.000133
dfeShape=0.35
dominance=0.25

pop=sheehanSong
discoalLaunchFile=discoalCmds/$pop.sh
annotDir=drosophilaAnnotations/
geneAnnotFile=$annotDir/refSeqAnnot.dm3.ucsc.12212018.gtf.gz
cncAnnotFile=$annotDir/phastConsElements15way_dm3.txt.gz
recRateFile=$annotDir/comeron.allMajArms.normalized.wig
chrLenFile=$annotDir/chromolens_dm3_autosOnly.txt
gapFile=$annotDir/gaps.bed
repeatedCoordsFile=repeatedCoords/drosophila_dm3.txt
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
    python runCmdAsJobArray.py "$cmd" bgsRep bgsRep.txt 24:00:00 general 4000 $logDir/$pop/${logPrefix}.log 0-999
done
