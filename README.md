# Overview
This repository contains the code for my investigation of the separability of models of background selection and genetic hitchhiking in simulations designed to model these processes in humans and *Drosophila* based on these genomes' gene annotations, recombination maps, esimated distributions of fitness effects (DFEs), etc (see preprint at https://doi.org/10.1101/2019.12.13.876136). The goal is to create a more accurate model of how BGS behaves in actual genomes than has been done previously. Note that all annotations are from the GRCh37/hg19 and BDGP R5/dm3dm3 versions of the human and *Drosophila* assemblies, respectively.

This repo is not a cohesive software package, but rather a set of pipelines for conducting coalescent and forward simulations of hitchhiking and background selection, calculating population genetic summary statistics from these simulations, using classifiers to discriminate between the two models, and plotting summaries and results.

The code in this repo was used to run and process data from large numbers of simulations on the high-performance computing resources here at UNC. It assumes that it is being run on a Linux system using the SLURM scheduler and that has a partition called 'general'--otherwise the code will have to be edited accordingly (i.e. `runCmdAsJob.py` and `runCmdAsJobArray.py` and all calls to them). Also note that the full set of simulations created by running this code is fairly large, so you may with to alter the output path (which is simply set to `simData` in the `bash` and `python` scripts in this repo) to one on a volume with a fair amount of free space (at least 100 GB). I am viewing this repo primarily as a way to share methodological details for my manuscript that others may be interested in borrowing for their own analyses. However, if you are interesting in reproducing my analyses and encounter any issues when getting the pipelines here to run properly on your system please don't hesitate to contact me.

# Dependencies:

- `scikit-allel`: I used version 1.2.1. Installation instructions at https://scikit-allel.readthedocs.io/en/stable/
- `fwdpy11`: Must use version 0.1.4; code is not compatible with newer version that use tree-sequence recording.
- `diploSHIC`: The code in this repo assumes that `diploSHIC` is present in directory ~/diploSHIC/, that this directory is present in your environment's `PYTHONPATH` variable, and that its dependencies are also installed (all of which are available via conda or pip). For more information see https://github.com/kern-lab/diploSHIC
- `discoal`: This coalescent simulator is used for selective sweep simulations. It must be installed and the executable must be in a directory in your environment's `PATH` variable. See https://github.com/kern-lab/discoal

# Installation

Install of the dependencies above, then clone this repo and run `python setup.py install`.

# Pipeline

## Pipeline for simulating background selection (BGS) and BGS + sweeps via `fwdpy11`

To run BGS simulations and calculate some statistics summarize them, use the following `bash` scripts in any order:

- `./pipelineForRandomBgsSamples.sh`
- `./pipelineForRandomBgsSamples_dros.sh`
- `./dominancePipeline.sh`
- `./dominancePipeline_dros.sh`
- `./pipelineForRepeatedBgsExamples.sh`
- `./pipelineForRepeatedBgsExamples_dros.sh`
- `./pipelineForRandomBgsSamplesBigly.sh`
- `./pipelineForRandomBgsSamplesBigly_dros.sh`
- `./pipelineForRandomBgsSamplesSweep.sh`

These will each launch a series of SLURM array jobs. You will have to wait for them all to finish (which may take several days or longer depending on your computing resources) before proceeding with downstream analyses (though the coalescent simulations below could be launched in parallel).

## Pipeline for simulating selective sweeps

Simlpy run the following bash scripts in any order:

- `./simLaunchScripts/humanEquilib.sh`
- `./simLaunchScripts/tennessenEuro.sh`
- `./simLaunchScripts/sheehanSong.sh`

These scripts are generated via `./generateAllSweepSimLaunchFiles.sh`. Again, these will be launched to your SLURM queue and may take some time to complete.

## Calculating statistics

You man run the following commands in any order

- `./generateAllSmallAndBigStatsForShicSims.sh`
- `./generateAllShicSummaryStatsFwdpy.sh`
- `./generateAllShicSummaryStats.sh`
- `./generateAllShicSummaryStatsFwdpyDominance.sh`

These will calculate statistics on the coalescent simulations and also some additional statistics for our forward simulations. Again, these will be launched to your SLURM queue and may take some time to complete.

## Training S/HIC classifiers for discriminating between sweeps and neutrally evolving regions:

Once the calculations above have finished, run the following in this order:

1. `./generateAllShicTrainingTestSets.sh`
2. `./trainShic.sh`

Again, these will be launched to your SLURM queue and may take some time to complete.

## Running S/HIC classifiers on the BGS and sweep simulations:

Once training is complete, run the following in any order:

- `./testShicAndGetSummaries.sh`
- `./testShicOnAllDominanceSims.sh`
- `./testShicOnAllBiglySims.sh`
- `python sampleUniformlyFromSweepParamRange.py` (This will also generate a plot results for the BGS+sweep simulations in the form of a heatmap.)

These also go to the queue and are fairly slow (could be sped up substantially if needed).

## Plotting results

To plot the mean values of various summary statistics across simulations of random genomic regions (top command), mean values across replicates of a small number of randomly selected regions (middle command), and the values from each individual simulation replicate (bottom command), simply run the following in any order:

- `./plotAllStatsRandom.sh`
- `./plotAllStatsRepeated.sh`
- `python plotAllIndivSims.py`
- `python plotBiglyResults.py`
- `python plotDominanceMisclassifications.py` (See commented usage example at the top of the script.)
- `python plotRecRateMisclassifications.py`

All plots will appear in the `plots` directory and subdirectories therein.

# Additional contents

In addition to the code described above, this repo contains annotations from the human (hg19) and *Drosophila* (dm3) genomes that we used to model BGS. See manuscript for more information (https://doi.org/10.1101/2019.12.13.876136).
