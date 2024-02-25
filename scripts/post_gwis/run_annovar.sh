gene=$1


ANNOVARDIR=../opt/annovar
OUTPUTDIR=../data/processed/annovar
TARGETDIR=../data/processed/annovar


perl $ANNOVARDIR/annotate_variation.pl \
        -out $OUTPUTDIR/${gene}_variants \
        -build hg19 \
        $TARGETDIR/${gene}_annovar_input.tsv \
        $ANNOVARDIR/humandb/
