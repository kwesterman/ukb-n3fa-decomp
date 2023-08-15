#!/bin/sh


#$ -l h_vmem=160G
#$ -l h_rt=24:00:00
#$ -cwd
#$ -j y


source /broad/software/scripts/useuse
use Anaconda3


source activate ldsc


ldref_dir=../data/processed/ld_ref
plink19_dir=../opt/plink1.9
varmap_dir=../data/raw/1000GP_Phase3


# Sort plinkset using --make-bed and update cM map
#wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz  # From https://mathgen.stats.ox.ac.uk/impute/1000GP%20Phase%203%20haplotypes%206%20October%202014.html
#tar xzvf 1000GP_Phase3.tgz
#mv 1000GP_Phase3 ../data/raw/

cut -f 1 ../data/raw/ldsc/w_hm3.snplist > ${ldref_dir}/hm3_snps.txt
${plink19_dir}/plink \
        --bfile ${ldref_dir}/ukb_20k_hg19 \
	--extract ${ldref_dir}/hm3_snps.txt \
	--seed 123 \
	--thin-indiv-count 5000 \
        --make-bed \
	--memory 10000 \
        --out ${ldref_dir}/ukb_20k_hg19_withCM
	#--maf 0.05 \
${plink19_dir}/plink \
        --bfile ${ldref_dir}/ukb_20k_hg19_withCM \
        --cm-map ${varmap_dir}/genetic_map_chr@_combined_b37.txt \
        --make-just-bim \
	--memory 10000 \
        --out ${ldref_dir}/ukb_20k_hg19_withCM


# Estimate LD scores using reference panel
python ../opt/ldsc/ldsc.py \
	--bfile ${ldref_dir}/ukb_20k_hg19_withCM \
	--l2 \
	--ld-wind-cm 1 \
	--out ${ldref_dir}/ukb_20k_hg19_withCM
