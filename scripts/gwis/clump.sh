#!/bin/bash


#$ -l h_vmem=10G
#$ -l h_rt=2:00:00

#$ -cwd
#$ -j y


fn_prefix=$1


ldref_dir=../data/processed/ld_ref


plink2=../opt/plink2.0/plink2
plink19=../opt/plink1.9/plink


${plink19} \
	--bfile ${ldref_dir}/ukb_20k_hg19 \
	--clump ${fn_prefix}_merged_nom \
	--clump-p1 0.00001 \
	--clump-p2 0.001 \
	--clump-r2 0.2 \
	--clump-kb 5000\
	--clump-snp-field RSID \
	--clump-field robust_P_int \
	--out ${fn_prefix}

##${plink2} \
##	--bfile ${ldref_dir}/ukb_20k_hg19 \
##	--clump ${fn_prefix}_merged_nom \
##	--clump-p1 0.00001 \
##	--clump-p2 0.001 \
##	--clump-r2 0.2 \
##	--clump-kb 5000\
##	--clump-id-field RSID \
##	--clump-p-field robust_P_int \
##	--clump-a1-field EA \
##	--out ${fn_prefix}
