input_filename=$1
gene=$2


ANNOVARDIR=../opt/annovar
OUTPUTDIR=../data/processed/annovar
TARGETDIR=../data/processed/annovar


# whr_loci %>%
#   mutate(Locus=paste(start, end, sep="-"),
#          index_snp=whr_ss$rsID[match(index_snp, whr_ss$SNPID)],
#          p=format(p, digits=2, scientific=T),
#          marg_locus=ifelse(marg_locus, "Yes", "No"),
#          Pulit_sexDM=ifelse(Pulit_sexDM, "Yes", "No")) %>%
#   select(Chr=chr, Locus, `Index SNP`=index_snp, P=p, `Marginal locus`=marg_locus, `Pulit et al. interaction locus`=Pulit_sexDM) %>%
# write_tsv("../output/interaction_test_annovar_input.tsv")
awk '{print $3, $4, $4, $5, $6, $2}' $input_filename | tail -n +2 > $TARGETDIR/tmp_annovar_input.tsv

perl $ANNOVARDIR/annotate_variation.pl \
        -out $OUTPUTDIR/${gene}_variants \
        -build hg19 \
        $TARGETDIR/tmp_annovar_input.tsv \
        $ANNOVARDIR/humandb/
