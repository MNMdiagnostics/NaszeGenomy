TRIO=$1
PROBAND=`echo $TRIO | cut -d, -f3`
MULITSAMPLE_VCF=/tmp/multisample_20210716.dv.bcfnorm.filt.vcf.gz
echo $TRIO $PROBAND

bcftools view -f "PASS,." -s $TRIO -Ou $MULITSAMPLE_VCF \
	| bcftools view -e 'INFO/AC<1' -Ou \
	| bcftools +mendelian -t $TRIO -mx -Oz -o ${PROBAND}.mendelian.vcf.gz
