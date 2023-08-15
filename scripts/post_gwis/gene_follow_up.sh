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
R --vanilla <<EOF
library(tidyverse)
pos_vec <- as.numeric(str_split(gsub(".*:", "", "${range}"), "-", simplify=TRUE))
gene_ss <- read_tsv("../data/processed/gwis/${e}_${pheno}_chr${chr}") %>%
  filter(POS >= pos_vec[1], POS <= pos_vec[2]) 
gene_ss %>%
  write_tsv("${fu_dir}/${e}_${pheno}_${gene_symbol}_subset")
top_rsid <- gene_ss[["RSID"]][which.min(gene_ss[["robust_P_Value_Interaction"]])]
write(top_rsid, "${fu_dir}/${e}_${pheno}_${gene_symbol}_top_rsid.txt")
EOF

# Use qctool v2 to subset (by position) to variants in the gene
qctool="~/kw/opt/qctool_v2.0.6-CentOS\ Linux7.3.1611-x86_64/qctool"
eval "${qctool}" \
	-g /broad/ukbb/imputed_v3/ukb_imp_chr${chr}_v3.bgen \
	-s /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
	-incl-range "${range}" \
	-og ${geno_dir}/${gene_symbol}_genotypes.bgen \
	-os ${geno_dir}/${gene_symbol}_genotypes.sample

# Convert genotype file to PLINK2 format and export to VCF
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

# Test GxE at all variants in the region (no MAF filter)
R --vanilla <<EOF
library(tidyverse)
library(vcfR)

run_gxe <- function(g, e, y, 
                    covars = c("age", "sex"), 
                    df = phenos,
                    std = FALSE) {
  if (std) df <- mutate(df, across(all_of(c(e, y)), ~ scale(.)))
  lm_str <- paste0(y, " ~ ", g, " * ", e, " + ", paste0(covars, collapse = " + "))
  lm_summ <- tryCatch({
    lm_fit <- lm(as.formula(lm_str), data = df) 
    lm_fit %>%
      broom::tidy() %>%
      filter(term == paste0(g, ":", e)) %>%
      mutate(residual_df = lm_fit[["df.residual"]])
  }, error = tibble(NA))
  lm_summ
}

minimal_covars <- c("age", "age_squared", "sex", "ageBySex")
ses_hl_covars <- c("ac", "income", "education", "smoking", "alcohol")
ffq_covars <- c("cooked_veg", "raw_veg", "fresh_fruit", "prmeat", "whole_bread")
covar_sets <- list(
  minimal = minimal_covars,
  adj = c(minimal_covars, ses_hl_covars),
  mdsAdj = c(minimal_covars, ses_hl_covars, "mds"),
  ffqAdj = c(minimal_covars, ses_hl_covars, ffq_covars)
)

gt_df <- read_tsv("${geno_dir}/${gene_symbol}_genotypes.raw", name_repair = "unique") %>%
  select(id = IID, matches("^rs")) %>%
  mutate(id = as.character(id))
  #rename_with(~ gsub("_.*", "", .), everything())

regression_df <- read_csv("../data/processed/ukb_gwis_phenos.csv") %>%
  mutate(id = as.character(id)) %>%
  left_join(gt_df, by = "id")

all_rsids <- setdiff(colnames(gt_df), "id")
gwis_covars <- scan("../data/processed/gwis_covariates.txt", what=character())
lm_res_df <- tibble(rsid = all_rsids) %>%
  rowwise() %>%
  mutate(model_res = list(run_gxe(rsid, "${e}", "${pheno}", 
				  covars = gwis_covars, df = regression_df))) %>%
  unnest(model_res)
lm_res_df %>%
  write_csv("${fu_dir}/${e}_${pheno}_${gene_symbol}_regressions")

#top_rsid <- scan("${fu_dir}/${e}_${pheno}_${gene_symbol}_top_rsid.txt", what = character())
#gt_df %>%
#  select(all_of(c("id", top_rsid))) %>%
#  write_csv("${fu_dir}/${e}_${pheno}_${gene_symbol}_top_rsid_genotypes.csv")
EOF

#vcfr_obj <- read.vcfR("${geno_dir}/${gene_symbol}_genotypes.vcf.gz")
#vcf_ids <- getID(vcfr_obj)
#gt_df <- extract.gt(vcfr_obj, element = "DS", mask = !duplicated(vcf_ids), as.numeric = TRUE) %>%
#  t() %>%
#  as_tibble(rownames = "id")

#lm_res_df <- tibble(rsid = all_rsids) %>%
#  slice(1:5) %>%
#  rowwise() %>%
#  mutate(model_res = list(run_gxe(rsid, "fish_oil", "Omega_3_pct", gwis_covars, regression_df))) %>%
#  unnest(model_res)
#
