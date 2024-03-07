#!/bin/bash


#$ -l h_vmem=70G
#$ -l h_rt=8:00:00

#$ -cwd
#$ -j y

source /broad/software/scripts/useuse
use R-4.1

Rscript data_prep/prep_ukb_phenos.R
