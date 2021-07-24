TRIOS=trios_ok_210723.txt

cat $TRIOS | parallel -j10 ./denovo_extract.sh
cut -d, -f3 $TRIOS | parallel -j20 ./denovo_filter.sh {}.mendelian.vcf.gz
cut -d, -f3 $TRIOS | parallel -j20 ./denovo_annotate.sh {}.mendelian.filt_indel.vcf.gz
cut -d, -f3 $TRIOS | parallel -j20 ./denovo_annotate.sh {}.mendelian.filt_snv.vcf.gz


