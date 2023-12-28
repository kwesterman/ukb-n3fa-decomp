#!/bin/bash


#$ -l h_vmem=20G
#$ -l h_rt=1:00:00

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
use Anaconda3


# Lift gene position bounds to GRCh38
echo ${range} | sed 's/[:-]/\t/g' > ${geno_dir}/${gene_symbol}_bounds_hg19.bed

liftover=../opt/liftover/liftOver
${liftover} \
	${geno_dir}/${gene_symbol}_bounds_hg19.bed \
        ../opt/liftover/hg19ToHg38.over.chain.gz \
	${geno_dir}/${gene_symbol}_bounds_hg38.bed \
	${geno_dir}/${gene_symbol}_bounds_hg38_unlifted.bed 
range_hg38=$(awk '{print $1":"$2"-"$3}' ${geno_dir}/${gene_symbol}_bounds_hg38.bed)


# Use bgenix to subset (by position) to variants in the gene
##qctool="../opt/qctool_v2.0.6-CentOS\ Linux7.3.1611-x86_64/qctool"
##eval "${qctool}" \
##	-g /humgen/florezlab/users/pschroeder/UKBB/TOPMedImputation/imputed_data/UKBB_UKBL_chr${chr}/region_merge_union/chr${chr}.dose.vcf.bgen \
##	-s /humgen/florezlab/users/pschroeder/UKBB/TOPMedImputation/GWAS_UKBB_UKBL/final.bgen.sample \
##	-incl-range "${range_hg38}" \
##	-og ${geno_dir}/${gene_symbol}_genotypes.bgen \
##	-os ${geno_dir}/${gene_symbol}_genotypes.sample
##	#-g /broad/ukbb/imputed_v3/ukb_imp_chr${chr}_v3.bgen \
##	#-s /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
source activate bgen
bgenix \
	-g /humgen/florezlab/users/pschroeder/UKBB/TOPMedImputation/imputed_data/UKBB_UKBL_chr${chr}/region_merge_union/chr${chr}.dose.vcf.bgen \
	-incl-range "${range_hg38}" \
	> ${geno_dir}/${gene_symbol}_genotypes.bgen
cp /humgen/florezlab/users/pschroeder/UKBB/TOPMedImputation/GWAS_UKBB_UKBL/final.bgen.sample ${geno_dir}/${gene_symbol}_genotypes.sample
conda deactivate


# Convert genotype file to PLINK2 format, lift back to hg19, and export to tabular format
plink2=../opt/plink2.0/plink2
${plink2} \
	--bgen ${geno_dir}/${gene_symbol}_genotypes.bgen ref-first \
	--sample ${geno_dir}/${gene_symbol}_genotypes.sample \
	--make-pgen \
	--out ${geno_dir}/${gene_symbol}_genotypes
mv ${geno_dir}/${gene_symbol}_genotypes.pvar ${geno_dir}/${gene_symbol}_genotypes_hg38.pvar
awk 'BEGIN {OFS="\t"} {print "chr" $1, $2, $2}' ${geno_dir}/${gene_symbol}_genotypes_hg38.pvar | tail -n +2 > ${geno_dir}/${gene_symbol}_genotypes_hg38.bed
${liftover} \
        ${geno_dir}/${gene_symbol}_genotypes_hg38.bed \
        ../opt/liftover/hg38ToHg19.over.chain.gz \
        ${geno_dir}/${gene_symbol}_genotypes_hg19.bed \
        ${geno_dir}/${gene_symbol}_genotypes_hg19.unlifted.bed 
paste <(tail -n +2 ${geno_dir}/${gene_symbol}_genotypes_hg38.pvar) ${geno_dir}/${gene_symbol}_genotypes_hg19.bed \
	| awk 'BEGIN {OFS="\t"} {print "chr" $1, $7, "chr"$1":"$7":"$4":"$5, $4, $5}' \
	> ${geno_dir}/${gene_symbol}_genotypes.pvar
cat <(echo "#CHROM  POS	ID	REF	ALT") <(cat ${geno_dir}/${gene_symbol}_genotypes.pvar) > ${geno_dir}/${gene_symbol}_genotypes.pvar_tmp
mv ${geno_dir}/${gene_symbol}_genotypes.pvar_tmp ${geno_dir}/${gene_symbol}_genotypes.pvar

##${plink2} \
##	--pfile ${geno_dir}/${gene_symbol}_genotypes \
##	--export vcf bgz vcf-dosage=DS-force id-paste=iid \
##	--out ${geno_dir}/${gene_symbol}_genotypes 
${plink2} \
	--pfile ${geno_dir}/${gene_symbol}_genotypes \
	--export A \
	--out ${geno_dir}/${gene_symbol}_genotypes \
	--memory 15000


# Subset genome-wide summary stats to the gene region

R --no-save <<EOF
library(tidyverse)

pos_vec <- as.numeric(str_split(gsub(".*:", "", "${range}"), "-", simplify=TRUE))  # Already includes upstream/downstream padding
gene_ss <- read_tsv("../data/processed/gwis/${e}_${pheno}_chr${chr}") %>%
  filter(POS >= pos_vec[1], POS <= pos_vec[2]) 
gene_ss %>%
  write_tsv("${fu_dir}/${e}_${pheno}_${gene_symbol}_subset")

missing_rsids <- scan("${fu_dir}/missing_rsids.txt", what = character())
top_snp_row <- gene_ss %>%
  filter(!(RSID %in% missing_rsids)) %>%
  slice(which.min(robust_P_Value_Interaction))
top_rsid_str <- top_snp_row[["RSID"]]
if (!grepl("^rs", top_rsid_str)) top_rsid_str <- paste0("chr", top_rsid_str)
write(top_rsid_str, "${fu_dir}/${e}_${pheno}_${gene_symbol}_top_rsid.txt")
top_snpid_str <- paste(paste0("chr", str_remove(top_snp_row[["CHR"]], "^0+")), top_snp_row[["POS"]], top_snp_row[["Non_Effect_Allele"]], top_snp_row[["Effect_Allele"]], sep = ":")
write(top_snpid_str, "${fu_dir}/${e}_${pheno}_${gene_symbol}_top_snpid.txt")
write(top_snp_row[["Effect_Allele"]], "${fu_dir}/${e}_${pheno}_${gene_symbol}_top_snp_allele.txt")
EOF


# Extract and export the top SNP

${plink2} \
	--pfile ${geno_dir}/${gene_symbol}_genotypes \
	--extract ${fu_dir}/${e}_${pheno}_${gene_symbol}_top_snpid.txt \
	--export A \
	--export-allele ${fu_dir}/${e}_${pheno}_${gene_symbol}_top_snp_allele.txt \
	--out ${fu_dir}/${e}_${pheno}_${gene_symbol}_top_snp


# Test GxE at all variants in the region (no MAF filter)

gPC_arr=(gPC1 gPC2 gPC3 gPC4 gPC5 gPC6 gPC7 gPC8 gPC9 gPC10)
gPC_int_arr=( "${gPC_arr[@]/#/${e}By}" )
covars=$(cat ../data/processed/gwis_covariates.txt | tr '\n' ' ')
covars="${covars} ${gPC_int_arr[@]}"

R --no-save <<EOF
library(tidyverse)
phenos <- read_csv("../data/processed/ukb_gwis_phenos.csv") %>%
  mutate(across(contains("gPC"), ~ . * ${e}, .names="${e}By{.col}"))
top_snp_df <- read_tsv("${fu_dir}/${e}_${pheno}_${gene_symbol}_top_snp.raw") %>% 
  select(id = IID, last_col()) %>%
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
	--maf 0.0001 \
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
	--covar-names ${covars} top_snp top_snp_gxe \
	--delim , \
	--missing-value NA \
	--cat-threshold 3 \
	--maf 0.0001 \
	--robust 1 \
	--output-style meta \
	--out /data/gene_followup/${e}_${pheno}_${gene_symbol}_regressions_cond

EOF

rm ../data/processed/${e}_${pheno}_phenos_${gene_symbol}.tmp


# Save top conditional variant

R --no-save <<EOF
library(tidyverse)

gene_ss_cond <- read_tsv("${fu_dir}/${e}_${pheno}_${gene_symbol}_regressions_cond") %>%
  filter(AF >= 0.01, AF <= 0.99)
top_snp_row <- gene_ss_cond[which.min(gene_ss_cond[["robust_P_Value_Interaction"]]), ]
top_snpid_str <- paste(paste0("chr", top_snp_row[["CHR"]]), top_snp_row[["POS"]], top_snp_row[["Non_Effect_Allele"]], top_snp_row[["Effect_Allele"]], sep = ":")
write(top_snpid_str, "${fu_dir}/${e}_${pheno}_${gene_symbol}_top_snpid2.txt")
EOF
