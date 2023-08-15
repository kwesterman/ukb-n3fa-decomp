#!/bin/sh


#$ -l os=RedHat7
#$ -l h_vmem=5G
#$ -l h_rt=12:00:00

#$ -pe smp 8
#$ -binding linear:8
#$ -R y

#$ -j y
#$ -cwd


exposure=$1
pheno=$2

chr=$SGE_TASK_ID


gPC_arr=(gPC1 gPC2 gPC3 gPC4 gPC5)
gPC_int_arr=( "${gPC_arr[@]/#/${exposure}By}" )
covars=$(cat ../data/processed/gwis_covariates.txt | tr '\n' ' ')
covars="${covars} ${gPC_int_arr[@]}"


source /broad/software/scripts/useuse
use R-4.0
R --no-save <<EOF
library(tidyverse)
read_csv("../data/processed/ukb_gwis_phenos.csv") %>%
  mutate(across(contains("gPC"), ~. * ${exposure}, .names="${exposure}By{.col}")) %>%
  write_csv("../data/processed/${exposure}_${pheno}_phenos_chr${chr}.tmp")
EOF

singularity_dir=~/kw/singularity
singularity exec \
	-B ../data/processed:/data \
	-B /broad/ukbb/imputed_v3:/bgendir \
	-B /humgen/florezlab/UKBB_app27892:/sampledir \
	${singularity_dir}/gem-v1.4.1a-workflow.simg \
	/bin/bash <<EOF

/GEM/GEM \
	--bgen /bgendir/ukb_imp_chr${chr}_v3.bgen \
	--sample /sampledir/ukb27892_imp_chrAUT_v3_s487395.sample \
	--pheno-file /data/${exposure}_${pheno}_phenos_chr${chr}.tmp \
	--sampleid-name id \
	--pheno-name ${pheno} \
	--exposure-names ${exposure} \
	--covar-names ${covars} \
	--delim , \
	--maf 0.01 \
	--missing-value NA \
	--robust 1 \
	--threads 8 \
	--output-style minimum \
	--out /data/gwis/${exposure}_${pheno}_chr${chr}

EOF

rm ../data/processed/${exposure}_${pheno}_phenos_chr${chr}.tmp
