#!/bin/bash


#$ -l h_vmem=10G
#$ -l h_rt=3:00:00

#$ -cwd
#$ -j y


source /broad/software/scripts/useuse
use GCC-5.2


magma_dir=../opt/magma_v1.10
data_dir=../data/raw/magma
ldref_dir=../data/processed/ld_ref
output_dir=../data/processed/magma

geneloc_file=${data_dir}/NCBI37.3.gene.loc


# Create SNP-gene map
${magma_dir}/magma \
	--annotate window=0,0 \
	--snp-loc ${ldref_dir}/ukb_snp_loc_maf0001.txt \
	--gene-loc ${geneloc_file} \
	--out ${output_dir}/ukb_20k_hg19_0.0

${magma_dir}/magma \
	--annotate window=2,1 \
	--snp-loc ${ldref_dir}/ukb_snp_loc_maf0001.txt \
	--gene-loc ${geneloc_file} \
	--out ${output_dir}/ukb_20k_hg19_2.1

${magma_dir}/magma \
	--annotate window=10,2 \
	--snp-loc ${ldref_dir}/ukb_snp_loc_maf0001.txt \
	--gene-loc ${geneloc_file} \
	--out ${output_dir}/ukb_20k_hg19_10.2
