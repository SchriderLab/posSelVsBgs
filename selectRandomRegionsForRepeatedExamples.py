import sys, random
import numpy as np
import selRegionsFromAnnot

"""
This script was used to spit out 10 random regions for both human and drosophila.
For each of these regions we later simulated 1000 replicates.

usage eg:
python selectRandomRegionsForRepeatedExamples.py drosophilaAnnotations/chromolens_dm3_autosOnly.txt drosophilaAnnotations/gaps.bed drosophilaAnnotations/refSeqAnnot.dm3.ucsc.12212018.gtf.gz 110000 10
"""

chrLenFileName, gapFileName, geneAnnotFileName, L, numRegionsToDraw = sys.argv[1:]
numRegionsToDraw = int(numRegionsToDraw)

npSeed = random.randint(0,(2**31)-1)
np.random.seed(npSeed)

if gapFileName.lower() in ["none", "false"]:
    gapFileName = None

L=int(L)
sys.stderr.write("reading chromosome lengths, selected sites, and recombination map\n")
chrLens = selRegionsFromAnnot.readChrLens(chrLenFileName)
for i in range(numRegionsToDraw):
    winC, winS, winE = selRegionsFromAnnot.pickRandomWindow(L, chrLens, gapFileName=gapFileName, state=np.random.get_state())
    sys.stderr.write("modeling simulated region after %s:%d-%d\n" %(winC, winS, winE))
