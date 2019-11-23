#!/bin/bash

annotDir=humanAnnotations/
flyAnnotDir=drosophilaAnnotations/
popLs=(humanEquilib sheehanSong tennessenEuro)
cncAnnotFileLs=($annotDir/phastConsElements100way_hg19.txt.gz $flyAnnotDir/phastConsElements15way_dm3.txt.gz $annotDir/phastConsElements100way_hg19.txt.gz)
geneAnnotFileLs=($annotDir/refSeqAnnot.hg19.ucsc.12102018.gtf.gz $flyAnnotDir/refSeqAnnot.dm3.ucsc.12212018.gtf.gz $annotDir/refSeqAnnot.hg19.ucsc.12102018.gtf.gz)
dfeLs=(boykoDFE flyHubDFE boykoDFE)
cncDenomLs=(100000 10000 100000)
winSizeLs=(100000 10000 100000)
cncMaxLs=(0.5 1.5 0.5)
for ((i=0; i<${#popLs[*]}; i++));
do
    pop=${popLs[$i]}
    dfe=${dfeLs[$i]}
    geneAnnotFileName=${geneAnnotFileLs[$i]}
    cncAnnotFileName=${cncAnnotFileLs[$i]}
    cncDenom=${cncDenomLs[$i]}
    winSize=${winSizeLs[$i]}
    cncMax=${cncMaxLs[$i]}
    echo $pop
    baseDir=simData/$dfe/repeatedExamples/
    plotDir=plots/repeated/$dfe/$pop
    mkdir -p $plotDir

    for region in `ls $baseDir`;
    do
        python plotBigStatsOneSet.py $baseDir/$region/stats/$pop/ simData/$dfe/randomExamples/ $geneAnnotFileName $cncAnnotFileName $winSize $cncDenom $cncMax $plotDir/$region.pdf
    done
done
