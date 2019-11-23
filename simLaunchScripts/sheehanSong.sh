#!/bin/bash

#generating training data

mkdir -p /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/sheehanSong/trainingData
python ~/pytools/runCmdAsJob.py "python ~/pytools/runDiscoalIfNotComplete.py discoal 100 2000 110000 -Pt 254.36 2543.6 -Pre 6435.308 19305.924 -en 0.03931435760339676 0 0.2684384337159931 -en 0.3931435760339676 0 1.0169838024846674 /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/sheehanSong/trainingData/simNeut.msOut.gz" neutSims neutSims.txt 48:00:00 general 128G /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/sheehanSong/trainingData/simNeut.msOut.log
i=0
for x in 0.045454545454545456 0.13636363636363635 0.22727272727272727 0.3181818181818182 0.4090909090909091 0.5 0.5909090909090909 0.6818181818181818 0.7727272727272727 0.8636363636363636 0.9545454545454546;
do
    python ~/pytools/runCmdAsJob.py "python ~/pytools/runDiscoalIfNotComplete.py discoal 100 2000 110000 -Pt 254.36 2543.6 -Pre 6435.308 19305.924 -en 0.03931435760339676 0 0.2684384337159931 -en 0.3931435760339676 0 1.0169838024846674 -ws 0 -Pa 127.180000 63590.000000 -Pu 0.000000 0.000079 -x $x -i 40 /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/sheehanSong/trainingData/simHard_$i.msOut.gz" hardSims hardSims.txt 48:00:00 general 128G /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/sheehanSong/trainingData/simHard_$i.msOut.log
    python ~/pytools/runCmdAsJob.py "python ~/pytools/runDiscoalIfNotComplete.py discoal 100 2000 110000 -Pt 254.36 2543.6 -Pre 6435.308 19305.924 -en 0.03931435760339676 0 0.2684384337159931 -en 0.3931435760339676 0 1.0169838024846674 -ws 0 -Pa 127.180000 63590.000000 -Pu 0.000000 0.000079 -Pf 0 0.05 -x $x -i 40 /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/sheehanSong/trainingData/simSoft_$i.msOut.gz" softSims softSims.txt 48:00:00 general 128G /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/sheehanSong/trainingData/simSoft_$i.msOut.log
    i=$((i + 1))
done


#generating test data

mkdir -p /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/sheehanSong/testData
python ~/pytools/runCmdAsJob.py "python ~/pytools/runDiscoalIfNotComplete.py discoal 100 2000 110000 -Pt 254.36 2543.6 -Pre 6435.308 19305.924 -en 0.03931435760339676 0 0.2684384337159931 -en 0.3931435760339676 0 1.0169838024846674 /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/sheehanSong/testData/simNeut.msOut.gz" neutSims neutSims.txt 48:00:00 general 128G /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/sheehanSong/testData/simNeut.msOut.log
i=0
for x in 0.045454545454545456 0.13636363636363635 0.22727272727272727 0.3181818181818182 0.4090909090909091 0.5 0.5909090909090909 0.6818181818181818 0.7727272727272727 0.8636363636363636 0.9545454545454546;
do
    python ~/pytools/runCmdAsJob.py "python ~/pytools/runDiscoalIfNotComplete.py discoal 100 2000 110000 -Pt 254.36 2543.6 -Pre 6435.308 19305.924 -en 0.03931435760339676 0 0.2684384337159931 -en 0.3931435760339676 0 1.0169838024846674 -ws 0 -Pa 127.180000 63590.000000 -Pu 0.000000 0.000079 -x $x -i 40 /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/sheehanSong/testData/simHard_$i.msOut.gz" hardSims hardSims.txt 48:00:00 general 128G /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/sheehanSong/testData/simHard_$i.msOut.log
    python ~/pytools/runCmdAsJob.py "python ~/pytools/runDiscoalIfNotComplete.py discoal 100 2000 110000 -Pt 254.36 2543.6 -Pre 6435.308 19305.924 -en 0.03931435760339676 0 0.2684384337159931 -en 0.3931435760339676 0 1.0169838024846674 -ws 0 -Pa 127.180000 63590.000000 -Pu 0.000000 0.000079 -Pf 0 0.05 -x $x -i 40 /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/sheehanSong/testData/simSoft_$i.msOut.gz" softSims softSims.txt 48:00:00 general 128G /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/sheehanSong/testData/simSoft_$i.msOut.log
    i=$((i + 1))
done

