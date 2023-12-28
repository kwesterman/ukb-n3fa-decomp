#!/bin/bash


#$ -l h_vmem=50G
#$ -l h_rt=5:00:00

#$ -cwd
#$ -j y

source /broad/software/scripts/useuse
use R-4.1

Rscript data_prep/prep_nmr.R
