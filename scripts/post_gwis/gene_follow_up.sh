#!/bin/bash


#$ -l h_vmem=30G
#$ -l h_rt=2:00:00

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

plink2=../opt/plink2.0/plink2


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
source activate bgen
bgenix \
	-g /humgen/florezlab/users/pschroeder/UKBB/TOPMedImputation/imputed_data/UKBB_UKBL_chr${chr}/region_merge_union/chr${chr}.dose.vcf.bgen \
	-incl-range "${range_hg38}" \
	> ${geno_dir}/${gene_symbol}_genotypes.bgen
cp /humgen/florezlab/users/pschroeder/UKBB/TOPMedImputation/GWAS_UKBB_UKBL/final.bgen.sample ${geno_dir}/${gene_symbol}_genotypes.sample
conda deactivate


# Convert genotype file to PLINK2 format, lift back to hg19, and export to tabular format
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
cat <(echo "#CHROM	POS	ID	REF	ALT") <(cat ${geno_dir}/${gene_symbol}_genotypes.pvar) > ${geno_dir}/${gene_symbol}_genotypes.pvar_tmp
mv ${geno_dir}/${gene_symbol}_genotypes.pvar_tmp ${geno_dir}/${gene_symbol}_genotypes.pvar

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

topmed_imp_pvar <- read_tsv("${geno_dir}/${gene_symbol}_genotypes.pvar")
top_snp_row <- gene_ss %>%
  arrange(robust_P_Value_Interaction) %>%
  filter(POS %in% topmed_imp_pvar[["POS"]]) %>%  # Interested in top variant that is also present in TOPMed-imputed genotype dataset
  slice(1)
top_rsid_str <- top_snp_row[["RSID"]]
if (!grepl("^rs", top_rsid_str)) top_rsid_str <- paste0("chr", top_rsid_str)
write(top_rsid_str, "${fu_dir}/${e}_${pheno}_${gene_symbol}_top_rsid.txt")
top_snp_str <- paste(paste0("chr", str_remove(top_snp_row[["CHR"]], "^0+")), top_snp_row[["POS"]], top_snp_row[["Non_Effect_Allele"]], top_snp_row[["Effect_Allele"]], sep = ":")
write(top_snp_str, "${fu_dir}/${e}_${pheno}_${gene_symbol}_top_snp.txt")
tibble(top_snp_str, top_snp_row[["Effect_Allele"]]) %>%
  write_tsv("${fu_dir}/${e}_${pheno}_${gene_symbol}_top_snp_allele.txt", col_names = FALSE)

nominal_snps <- gene_ss %>%
  filter(robust_P_Value_Interaction < 0.05) %>%
  mutate(snpid = paste(paste0("chr", str_remove(CHR, "^0+")), POS, Non_Effect_Allele, Effect_Allele, sep = ":"))
write(nominal_snps[["snpid"]], "${fu_dir}/${e}_${pheno}_${gene_symbol}_nominal_snps.txt")
EOF


# Extract and export the top SNP and nominal SNP set

${plink2} \
	--pfile ${geno_dir}/${gene_symbol}_genotypes \
	--extract ${fu_dir}/${e}_${pheno}_${gene_symbol}_top_snp.txt \
	--export A \
	--export-allele ${fu_dir}/${e}_${pheno}_${gene_symbol}_top_snp_allele.txt \
	--out ${fu_dir}/${e}_${pheno}_${gene_symbol}_top_snp

${plink2} \
	--pfile ${geno_dir}/${gene_symbol}_genotypes \
	--extract ${fu_dir}/${e}_${pheno}_${gene_symbol}_nominal_snps.txt \
	--export A \
	--out ${fu_dir}/${e}_${pheno}_${gene_symbol}_nominal_snps


# GxE for all gene variants (low MAF filter), including conditional on nominal GWIS variants

post_gwis/conditional_gxe.sh \
	${geno_dir}/${gene_symbol}_genotypes \
	"" \
	${e} \
	${pheno} \
	${fu_dir}/${e}_${pheno}_${gene_symbol}_regressions

post_gwis/conditional_gxe.sh \
	${geno_dir}/${gene_symbol}_genotypes \
	${fu_dir}/${e}_${pheno}_${gene_symbol}_nominal_snps.raw \
	${e} \
	${pheno} \
	${fu_dir}/${e}_${pheno}_${gene_symbol}_regressions_cond


# Save top conditional variant

R --no-save <<EOF
library(tidyverse)
gene_ss_cond <- read_tsv("${fu_dir}/${e}_${pheno}_${gene_symbol}_regressions_cond") %>%
  filter(AF < 0.01 | AF > 0.99)
top_snp_row <- gene_ss_cond[which.min(gene_ss_cond[["robust_P_Value_Interaction"]]), ]
top_snp_str <- paste(paste0(top_snp_row[["CHR"]]), top_snp_row[["POS"]], top_snp_row[["Non_Effect_Allele"]], top_snp_row[["Effect_Allele"]], sep = ":")
write(top_snp_str, "${fu_dir}/${e}_${pheno}_${gene_symbol}_top_cond_snp.txt")
write(top_snp_row[["Effect_Allele"]], "${fu_dir}/${e}_${pheno}_${gene_symbol}_top_cond_snp_allele.txt")
EOF


# Extract and export the top conditional SNP

${plink2} \
	--pfile ${geno_dir}/${gene_symbol}_genotypes \
	--extract ${fu_dir}/${e}_${pheno}_${gene_symbol}_top_cond_snp.txt \
	--export A \
	--export-allele ${fu_dir}/${e}_${pheno}_${gene_symbol}_top_cond_snp_allele.txt \
	--out ${fu_dir}/${e}_${pheno}_${gene_symbol}_top_cond_snp

