#!/bin/bash


#$ -l h_vmem=100G
#$ -l h_rt=6:00:00

#$ -cwd
#$ -j y


source /broad/software/scripts/useuse
use GCC-5.2
use R-4.0


plink19_dir=../opt/plink1.9
plink2_dir=../opt/plink2.0
data_dir=../data/raw/magma
ldref_dir=../data/processed/ld_ref
varmap_dir=../data/raw/1000GP_Phase3

ukb_plinkset_dir=/humgen/florezlab/users/rmandla/ukbb-topmed/cojo/ukbb_geno


# Merge UKB plinksets to create LD reference
echo -n "" > ${ldref_dir}/ukb_20k_merge_list.txt
for chr in {1..22} X; do
	echo "${ukb_plinkset_dir}/chr${chr}-rename" >> ${ldref_dir}/ukb_20k_merge_list.txt
done
${plink2_dir}/plink2 \
	--pmerge-list ${ldref_dir}/ukb_20k_merge_list.txt bfile \
	--merge-max-allele-ct 2 \
	--maf 0.0001 \
	--make-bed \
	--out ${ldref_dir}/ukb_20k

# Lift LD reference plinkset coordinates from h38 to hg19
awk 'BEGIN {OFS="\t"} {print "chr" $1, $4, $4}' ${ldref_dir}/ukb_20k.bim > ${ldref_dir}/bim_hg38.bed
liftover=../opt/liftover/liftOver
${liftover} \
	${ldref_dir}/bim_hg38.bed \
	../opt/liftover/hg38ToHg19.over.chain.gz \
	${ldref_dir}/bim_hg19.bed \
	${ldref_dir}/ukb_20k_unlifted.bed
cat ${ldref_dir}/ukb_20k_unlifted.bed | grep -v "Deleted in new" | awk '{print $1 ":" $2}' | sed 's/chr//g' > ${ldref_dir}/ukb_20k_unlifted.exclude
${plink2_dir}/plink2 \
	--bfile ${ldref_dir}/ukb_20k \
	--exclude ${ldref_dir}/ukb_20k_unlifted.exclude \
	--make-bed \
	--out ${ldref_dir}/ukb_20k_hg19
cp ${ldref_dir}/ukb_20k_hg19.bim ${ldref_dir}/ukb_20k_hg19.bim_bkp
paste \
	<(cut -f 1,3 ${ldref_dir}/ukb_20k_hg19.bim_bkp) \
	<(cut -f 2 ${ldref_dir}/bim_hg19.bed) \
	<(cut -f 5,6 ${ldref_dir}/ukb_20k_hg19.bim_bkp) \
	| awk 'BEGIN {OFS="\t"} {print $1, $1 ":" $3, $2, $3, $4, $5}' \
	> ${ldref_dir}/ukb_20k_hg19.bim

# Merge UKB variant annotations to create SNP location file
echo -n "" > ${ldref_dir}/ukb_snp_loc_maf0001.txt
for chr in {1..22} X; do
	awk -v chr=${chr} 'BEGIN {OFS="\t"} $6 > 0.0001 {print $2, chr, $3}' /broad/ukbb/imputed_v3/ukb_mfi_chr${chr}_v3.txt >> ${ldref_dir}/ukb_snp_loc_maf0001.txt
done	

# Add rsIDs to LD reference .bim file
R --vanilla <<EOF
library(tidyverse)
bim <- read_tsv("${ldref_dir}/ukb_20k_hg19.bim", col_names=c("CHR", "ID", "POS", "BP", "A1", "A2"), col_types=cols(CHR="c", ID="c"))
ref <- read_tsv("${ldref_dir}/ukb_snp_loc_maf0001.txt", col_names=c("rsID", "CHR", "BP"), col_types=cols(CHR="c")) %>%
  distinct(CHR, BP, .keep_all=TRUE)  # Note: this arbitrarily uses the first SNP where there are multiple in the same location
left_join(bim, ref, by=c("CHR", "BP")) %>%
  select(CHR, rsID, POS, BP, A1, A2) %>%
  write_tsv("${ldref_dir}/ukb_20k_hg19.bim", col_names=FALSE)
EOF
