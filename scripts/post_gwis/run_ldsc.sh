#!/bin/sh


#$ -l h_vmem=10G
#$ -l h_rt=1:00:00
#$ -cwd
#$ -j y


source /broad/software/scripts/useuse
use Anaconda3

source activate ldsc


gwis_dir=../data/processed/gwis
ldref_dir=../data/processed/ld_ref
ldsc_datadir=../data/processed/ldsc


e1=$1
e2=$2
pheno=$3


working_dir=${e1}_${e2}_${pheno}_ldsc_dir


# Munge summary statistics to prep for LDSC
mkdir -p ${working_dir}
for e in ${e1} ${e2}; do 
cut  -f 2- ${gwis_dir}/${e}_${pheno}_merged > ${working_dir}/${e}_${pheno}_ss
../opt/ldsc/munge_sumstats.py \
        --sumstats ${working_dir}/${e}_${pheno}_ss \
        --merge-alleles ../data/raw/ldsc/w_hm3.snplist \
        --snp RSID \
        --N-col N_Samples \
        --a1 Effect_Allele \
        --a2 Non_Effect_Allele \
        --p robust_P_Value_Interaction \
        --frq AF \
        --signed-sumstats Beta_G-${e},0 \
        --chunksize 500000 \
        --out ${working_dir}/${e}_${pheno}
rm ${working_dir}/${e}_${pheno}_ss
done


# Run LDSC to calculate genetic correlations
../opt/ldsc/ldsc.py \
	--rg ${working_dir}/${e1}_${pheno}.sumstats.gz,${working_dir}/${e2}_${pheno}.sumstats.gz \
	--ref-ld ${ldref_dir}/ukb_20k_hg19_withCM \
	--w-ld ${ldref_dir}/ukb_20k_hg19_withCM \
	--out ${ldsc_datadir}/${e1}_${e2}_${pheno}

rm -r ${working_dir}

conda deactivate
