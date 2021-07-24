
SAMPLES2INCLUDE=$1
VCF=$2

cat $SAMPLES2INCLUDE | parallel -j 20 ./sv_get_af_distibution_per_sample.sh $VCF 

zcat *.vcf.gz | gzip -c > sv_multisample_210519.per_sample_stats.tsv.gz 

