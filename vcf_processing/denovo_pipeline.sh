cat trios_ok_210622.txt | parallel -j10 ./denovo_extract.sh
ls *.mendelian.vcf.gz | parallel -j20 ./denovo_filter.sh
ls *.mendelian.filt.vcf.gz | parallel -j20 ./denovo_annotate.sh
