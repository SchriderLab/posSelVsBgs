#!/bin/bash

dataDir=simData
popLs=(humanEquilib sheehanSong tennessenEuro)
dfeLs=(boykoDFE flyHubDFE boykoDFE)

geneAnnotFileLs=(humanAnnotations/refSeqAnnot.hg19.ucsc.12102018.gtf.gz drosophilaAnnotations/refSeqAnnot.dm3.ucsc.12212018.gtf.gz humanAnnotations/refSeqAnnot.hg19.ucsc.12102018.gtf.gz)
cncAnnotFileLs=(humanAnnotations/phastConsElements100way_hg19.txt.gz drosophilaAnnotations/phastConsElements15way_dm3.txt.gz humanAnnotations/phastConsElements100way_hg19.txt.gz)
recRateFileLs=(humanAnnotations/Kong2010_SexAveraged.wig drosophilaAnnotations/comeron.allMajArms.normalized.wig humanAnnotations/Kong2010_SexAveraged.wig)
totalPhysLenLs=(1100000 110000 1100000)
rMeanLs=(1e-7 2.3e-6 1e-7)

for ((i=0; i<${#popLs[*]}; i++));
do
    rMean=${rMeanLs[$i]}
    geneAnnotFileName=${geneAnnotFileLs[$i]}
    cncAnnotFileName=${cncAnnotFileLs[$i]}
    recRateFileName=${recRateFileLs[$i]}
    pop=${popLs[$i]}
    dfe=${dfeLs[$i]}
    testDir=$dataDir/sweeps/$pop/testDataSets/
    absurdBgsFileName=$dataDir/$dfe/randomExamples/absurdBgs/shicFeatureVecs/$pop.fvec
    realBgsFileName=$dataDir/$dfe/randomExamples/realBgs/shicFeatureVecs/$pop.fvec
    realBgsWeakCNCFileName=$dataDir/$dfe/randomExamples/realBgsWeakCNC/shicFeatureVecs/$pop.fvec
    realNeutFileName=$dataDir/$dfe/randomExamples/realNeut/shicFeatureVecs/$pop.fvec
    realBgsLogDir=$dataDir/$dfe/randomExamples/realBgs/simLogs/$pop/
    realBgsWeakCNCLogDir=$dataDir/$dfe/randomExamples/realBgsWeakCNC/simLogs/$pop/
    realNeutLogDir=$dataDir/$dfe/randomExamples/realNeut/simLogs/$pop/
    outDir=$dataDir/sweeps/$pop/classifier/
    outputModel=$outDir/clf.mod
    summaryOutDir=$dataDir/sweeps/$pop/classificationSummaries/
    logFileName=$outDir/testing_summarizing.log
    figOutDir=plots/covfefe
    figFileName=$figOutDir/${pop}.pdf
    mkdir -p $outDir $figOutDir $summaryOutDir
    cmd="python buildHeatmapForClassifierAndSummarizeRegions.py $outputModel.json $outputModel.weights.hdf5 $testDir $absurdBgsFileName $realBgsFileName $realBgsWeakCNCFileName $realNeutFileName $realBgsLogDir $realBgsWeakCNCLogDir $realNeutLogDir $geneAnnotFileName $cncAnnotFileName $recRateFileName $rMean $summaryOutDir $figFileName"
    echo $cmd
    python runCmdAsJob.py "$cmd" tstShc tstShc.txt 24:00:00 general 16000 $logFileName
done
