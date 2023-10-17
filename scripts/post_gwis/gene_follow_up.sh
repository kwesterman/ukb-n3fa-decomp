#!/bin/bash


#$ -l h_vmem=20G
#$ -l h_rt=6:00:00

#$ -cwd
#$ -j y


chr=$1
range=$2
gene_symbol=$3
e=$4
pheno=$5

fu_dir=../data/processed/gene_followup
mkdir -p ${fu_dir}
geno_dir=../data/processed/genotypes
mkdir -p ${geno_dir}


source /broad/software/scripts/useuse
use R-4.1
use GCC-5.2


# Subset genome-wide summary stats to the gene region
R --no-save <<EOF
library(tidyverse)
pos_vec <- as.numeric(str_split(gsub(".*:", "", "${range}"), "-", simplify=TRUE))  # Already includes upstream/downstream padding
gene_ss <- read_tsv("../data/processed/gwis/${e}_${pheno}_chr${chr}") %>%
  filter(POS >= pos_vec[1], POS <= pos_vec[2]) 
gene_ss %>%
  write_tsv("${fu_dir}/${e}_${pheno}_${gene_symbol}_subset")
top_rsid_row <- gene_ss[which.min(gene_ss[["robust_P_Value_Interaction"]]), ]
write(top_rsid_row[["RSID"]], "${fu_dir}/${e}_${pheno}_${gene_symbol}_top_rsid.txt")
write(top_rsid_row[["Effect_Allele"]], "${fu_dir}/${e}_${pheno}_${gene_symbol}_top_rsid_allele.txt")
EOF

# Use qctool v2 to subset (by position) to variants in the gene
qctool="~/kw/opt/qctool_v2.0.6-CentOS\ Linux7.3.1611-x86_64/qctool"
eval "${qctool}" \
	-g /broad/ukbb/imputed_v3/ukb_imp_chr${chr}_v3.bgen \
	-s /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
	-incl-range "${range}" \
	-og ${geno_dir}/${gene_symbol}_genotypes.bgen \
	-os ${geno_dir}/${gene_symbol}_genotypes.sample

# Convert genotype file to PLINK2 format and export to raw and VCF
plink2=../opt/plink2.0/plink2
${plink2} \
	--bgen ${geno_dir}/${gene_symbol}_genotypes.bgen ref-first \
	--sample ${geno_dir}/${gene_symbol}_genotypes.sample \
	--make-pgen \
	--out ${geno_dir}/${gene_symbol}_genotypes
${plink2} \
	--pfile ${geno_dir}/${gene_symbol}_genotypes \
	--export vcf bgz vcf-dosage=DS-force id-paste=iid \
	--out ${geno_dir}/${gene_symbol}_genotypes 
${plink2} \
	--pfile ${geno_dir}/${gene_symbol}_genotypes \
	--export A \
	--out ${geno_dir}/${gene_symbol}_genotypes 
${plink2} \
	--pfile ${geno_dir}/${gene_symbol}_genotypes \
	--extract ${fu_dir}/${e}_${pheno}_${gene_symbol}_top_rsid.txt \
	--export A \
	--export-allele ${fu_dir}/${e}_${pheno}_${gene_symbol}_top_rsid_allele.txt \
	--out ${fu_dir}/${gene_symbol}_top_rsid

# Test GxE at all variants in the region (no MAF filter)

gPC_arr=(gPC1 gPC2 gPC3 gPC4 gPC5 gPC6 gPC7 gPC8 gPC9 gPC10)
gPC_int_arr=( "${gPC_arr[@]/#/${e}By}" )
covars=$(cat ../data/processed/gwis_covariates.txt | tr '\n' ' ')
covars="${covars} ${gPC_int_arr[@]}"

R --no-save <<EOF
library(tidyverse)
phenos <- read_csv("../data/processed/ukb_gwis_phenos.csv") %>%
  mutate(across(contains("gPC"), ~ . * ${e}, .names="${e}By{.col}"))
top_snp_df <- read_tsv("${fu_dir}/${gene_symbol}_top_rsid.raw") %>% 
  select(id = IID, matches("^rs")) %>%
  select(1, 2) %>%
  setNames(c("id", "top_snp"))
phenos %>%
  inner_join(top_snp_df, by = "id") %>%
  mutate(top_snp_gxe = top_snp * ${e}) %>%
  write_csv("../data/processed/${e}_${pheno}_phenos_${gene_symbol}.tmp")
EOF

singularity_dir=~/kw/singularity
singularity exec \
	-B ../data/processed:/data \
	-B /broad/ukbb/imputed_v3:/bgendir \
	-B /humgen/florezlab/UKBB_app27892:/sampledir \
	${singularity_dir}/gem-v1.5.2-workflow.simg \
	/bin/bash <<EOF

# Unconditional
/GEM/GEM \
	--bgen /data/genotypes/${gene_symbol}_genotypes.bgen \
	--sample /data/genotypes/${gene_symbol}_genotypes.sample \
	--pheno-file /data/${e}_${pheno}_phenos_${gene_symbol}.tmp \
	--sampleid-name id \
	--pheno-name ${pheno} \
	--exposure-names ${e} \
	--covar-names ${covars} \
	--delim , \
	--missing-value NA \
	--cat-threshold 3 \
	--maf 0.001 \
	--robust 1 \
	--output-style meta \
	--out /data/gene_followup/${e}_${pheno}_${gene_symbol}_regressions

# Conditional on top variant
/GEM/GEM \
	--bgen /data/genotypes/${gene_symbol}_genotypes.bgen \
	--sample /data/genotypes/${gene_symbol}_genotypes.sample \
	--pheno-file /data/${e}_${pheno}_phenos_${gene_symbol}.tmp \
	--sampleid-name id \
	--pheno-name ${pheno} \
	--exposure-names ${e} \
	--covar-names "${covars} top_snp top_snp_gxe" \
	--delim , \
	--missing-value NA \
	--cat-threshold 3 \
	--maf 0.001 \
	--robust 1 \
	--output-style meta \
	--out /data/gene_followup/${e}_${pheno}_${gene_symbol}_regressions_cond

EOF

rm ../data/processed/${e}_${pheno}_phenos_${gene_symbol}.tmp
