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
pheno1=$2
e2=$3
pheno2=$4


tag1=${e1}_${pheno1}
tag2=${e2}_${pheno2}

working_dir=${tag1}_${tag2}_ldsc_dir


# Munge summary statistics to prep for LDSC
mkdir -p ${working_dir}

cut  -f 2- ${gwis_dir}/${tag1}_merged > ${working_dir}/${tag1}_ss
../opt/ldsc/munge_sumstats.py \
        --sumstats ${working_dir}/${tag1}_ss \
        --merge-alleles ../data/raw/ldsc/w_hm3.snplist \
        --snp RSID \
        --N-col N_Samples \
        --a1 Effect_Allele \
        --a2 Non_Effect_Allele \
        --p robust_P_Value_Interaction \
        --frq AF \
        --signed-sumstats Beta_G-${e1},0 \
        --chunksize 500000 \
        --out ${working_dir}/${tag1}
rm ${working_dir}/${tag1}_ss

cut  -f 2- ${gwis_dir}/${tag2}_merged > ${working_dir}/${tag2}_ss
../opt/ldsc/munge_sumstats.py \
        --sumstats ${working_dir}/${tag2}_ss \
        --merge-alleles ../data/raw/ldsc/w_hm3.snplist \
        --snp RSID \
        --N-col N_Samples \
        --a1 Effect_Allele \
        --a2 Non_Effect_Allele \
        --p robust_P_Value_Interaction \
        --frq AF \
        --signed-sumstats Beta_G-${e2},0 \
        --chunksize 500000 \
        --out ${working_dir}/${tag2}
rm ${working_dir}/${tag2}_ss


# Run LDSC to calculate genetic correlations
../opt/ldsc/ldsc.py \
	--rg ${working_dir}/${tag1}.sumstats.gz,${working_dir}/${tag2}.sumstats.gz \
	--ref-ld ${ldref_dir}/ukb_20k_hg19_withCM \
	--w-ld ${ldref_dir}/ukb_20k_hg19_withCM \
	--out ${ldsc_datadir}/${tag1}_${tag2}

rm -r ${working_dir}

conda deactivate
