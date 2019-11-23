#!/bin/bash

#generating training data

mkdir -p /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/humanEquilib/trainingData
python ~/pytools/runCmdAsJob.py "python ~/pytools/runDiscoalIfNotComplete.py discoal 100 2000 220000 -Pt 96 960 -Pre 440 1320 /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/humanEquilib/trainingData/simNeut.msOut.gz" neutSims neutSims.txt 48:00:00 general 128G /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/humanEquilib/trainingData/simNeut.msOut.log
i=0
for x in 0.045454545454545456 0.13636363636363635 0.22727272727272727 0.3181818181818182 0.4090909090909091 0.5 0.5909090909090909 0.6818181818181818 0.7727272727272727 0.8636363636363636 0.9545454545454546;
do
    python ~/pytools/runCmdAsJob.py "python ~/pytools/runDiscoalIfNotComplete.py discoal 100 2000 220000 -Pt 96 960 -Pre 440 1320 -ws 0 -Pa 2.000000 1000.000000 -Pu 0.000000 0.005000 -x $x -i 40 /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/humanEquilib/trainingData/simHard_$i.msOut.gz" hardSims hardSims.txt 48:00:00 general 128G /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/humanEquilib/trainingData/simHard_$i.msOut.log
    python ~/pytools/runCmdAsJob.py "python ~/pytools/runDiscoalIfNotComplete.py discoal 100 2000 220000 -Pt 96 960 -Pre 440 1320 -ws 0 -Pa 2.000000 1000.000000 -Pu 0.000000 0.005000 -Pf 0 0.05 -x $x -i 40 /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/humanEquilib/trainingData/simSoft_$i.msOut.gz" softSims softSims.txt 48:00:00 general 128G /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/humanEquilib/trainingData/simSoft_$i.msOut.log
    i=$((i + 1))
done


#generating test data

mkdir -p /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/humanEquilib/testData
python ~/pytools/runCmdAsJob.py "python ~/pytools/runDiscoalIfNotComplete.py discoal 100 2000 220000 -Pt 96 960 -Pre 440 1320 /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/humanEquilib/testData/simNeut.msOut.gz" neutSims neutSims.txt 48:00:00 general 128G /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/humanEquilib/testData/simNeut.msOut.log
i=0
for x in 0.045454545454545456 0.13636363636363635 0.22727272727272727 0.3181818181818182 0.4090909090909091 0.5 0.5909090909090909 0.6818181818181818 0.7727272727272727 0.8636363636363636 0.9545454545454546;
do
    python ~/pytools/runCmdAsJob.py "python ~/pytools/runDiscoalIfNotComplete.py discoal 100 2000 220000 -Pt 96 960 -Pre 440 1320 -ws 0 -Pa 2.000000 1000.000000 -Pu 0.000000 0.005000 -x $x -i 40 /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/humanEquilib/testData/simHard_$i.msOut.gz" hardSims hardSims.txt 48:00:00 general 128G /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/humanEquilib/testData/simHard_$i.msOut.log
    python ~/pytools/runCmdAsJob.py "python ~/pytools/runDiscoalIfNotComplete.py discoal 100 2000 220000 -Pt 96 960 -Pre 440 1320 -ws 0 -Pa 2.000000 1000.000000 -Pu 0.000000 0.005000 -Pf 0 0.05 -x $x -i 40 /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/humanEquilib/testData/simSoft_$i.msOut.gz" softSims softSims.txt 48:00:00 general 128G /pine/scr/d/s/dschride/data/posSelVsBgs/sweeps/humanEquilib/testData/simSoft_$i.msOut.log
    i=$((i + 1))
done

