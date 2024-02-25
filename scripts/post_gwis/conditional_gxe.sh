#!/bin/bash


pfile_prefix=$1
snp_covar_fn=$2
e=$3
pheno=$4
output_fn=$5

fu_dir=../data/processed/gene_followup
mkdir -p ${fu_dir}
geno_dir=../data/processed/genotypes
mkdir -p ${geno_dir}


source /broad/software/scripts/useuse
use R-4.1


gPC_arr=(gPC1 gPC2 gPC3 gPC4 gPC5 gPC6 gPC7 gPC8 gPC9 gPC10)
gPC_int_arr=( "${gPC_arr[@]/#/${e}By}" )
covars=$(cat ../data/processed/gwis_covariates.txt | tr '\n' ' ')
covars="${covars} ${gPC_int_arr[@]}"

R --no-save <<EOF
library(tidyverse)
phenos <- read_csv("../data/processed/ukb_gwis_phenos.csv") %>%
  mutate(across(contains("gPC"), ~ . * ${e}, .names="${e}By{.col}"))
phenos %>%
  write_csv("${output_fn}_phenos.tmp")
EOF
  
if [ -n "${snp_covar_fn}" ]; then

R --no-save <<EOF
library(tidyverse)
phenos <- read_csv("${output_fn}_phenos.tmp")
nominal_snp_df <- read_tsv("${snp_covar_fn}") %>% 
  select(id = IID, everything(), -c(FID, PAT, MAT, SEX, PHENOTYPE))
nominal_snps <- setdiff(names(nominal_snp_df), "id")
phenos %>%
  inner_join(nominal_snp_df, by = "id") %>%
  mutate(across(all_of(nominal_snps), ~ . * ${e}, .names = "{.col}_gxe")) %>%
  write_csv("${output_fn}_phenos.tmp")
write(c(nominal_snps, paste0(nominal_snps, "_gxe")), 
      "${output_fn}_snp_covars.tmp")
EOF

snp_covars=$(cat ${output_fn}_snp_covars.tmp | tr '\n' ' ')
rm ${output_fn}_snp_covars.tmp

fi

singularity_dir=~/kw/singularity
singularity exec \
	-B ${pfile_prefix}.pgen:/pgen \
	-B ${pfile_prefix}.pvar:/pvar \
	-B ${pfile_prefix}.psam:/psam \
	-B ${output_fn}_phenos.tmp:/phenofile \
	-B ../data/processed/gene_followup:/outputdir \
	${singularity_dir}/gem-v1.5.2-workflow.simg \
	/bin/bash <<EOF
	#-B ${bgen_fn}:/bgenfile.bgen \
	#-B ${sample_fn}:/samplefile \

/GEM/GEM \
	--pgen /pgen \
	--pvar /pvar \
	--psam /psam \
	--pheno-file /phenofile \
	--sampleid-name id \
	--pheno-name ${pheno} \
	--exposure-names ${e} \
	--covar-names ${covars} ${snp_covars} \
	--delim , \
	--missing-value NA \
	--cat-threshold 3 \
	--maf 0.0001 \
	--robust 1 \
	--output-style meta \
	--out /outputdir/$(basename ${output_fn}) \
	--threads 1
	#--bgen /bgenfile.bgen \
	#--sample /samplefile \

EOF

rm ${output_fn}_phenos.tmp
