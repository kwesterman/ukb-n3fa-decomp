#!/bin/bash


#$ -l h_vmem=10G
#$ -l h_rt=6:00:00

#$ -cwd
#$ -j y


source /broad/software/scripts/useuse
use GCC-5.2
use R-4.0


sumstats_file=$1
prefix=$2


magma_dir=../opt/magma_v1.10
data_dir=../data/raw/magma
ldref_dir=../data/processed/ld_ref
output_dir=../data/processed/magma

geneloc_file=${data_dir}/NCBI37.3.gene.loc


#gunzip < ${sumstats_file} > ${prefix}_sumstats_tmp
#R --vanilla <<EOF
#library(tidyverse)
#read_tsv("${prefix}_sumstats_tmp") %>%
#  mutate(SNP = paste0(CHR, ":", POS)) %>%
#  write_tsv("${prefix}_sumstats_tmp")
#EOF

# Gene analysis (based on p-values)
${magma_dir}/magma \
	--bfile ${ldref_dir}/ukb_20k_hg19 \
	--pval ${sumstats_file} use=SNP,P ncol=N \
	--gene-model snp-wise=mean \
	--gene-annot ${output_dir}/ukb_20k_hg19_2.1.genes.annot \
	--genes-only \
	--out ${output_dir}/${prefix}
	#--pval ${prefix}_sumstats_tmp use=SNP,P ncol=N \
	#--gene-model snp-wise=multi \

R --vanilla <<EOF
library(tidyverse)
magma_res <- read_table("${output_dir}/${prefix}.genes.out", col_types=cols(CHR="c"))
gene_annot <- read_tsv("${geneloc_file}", col_names=c("ID", "CHR", "START", "END", "STRAND", "SYMBOL"),
		       col_types=cols(CHR="c"))
magma_res %>%
  left_join(select(gene_annot, ID, SYMBOL), by=c("GENE"="ID")) %>%
  write_tsv("${output_dir}/${prefix}.genes.out")
EOF

#rm ${prefix}_sumstats_tmp
