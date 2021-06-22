TRIO=$1
PROBAND=`echo $TRIO | cut -d, -f3`

echo $TRIO $PROBAND

bcftools view -f "PASS,." -s $TRIO -Ou /tmp/multisample_20210519.dv.bcfnorm.vcf.gz \
	| bcftools view -e 'INFO/AC<1' -Ou \
	| bcftools +mendelian -t $TRIO -mx -Oz -o ${PROBAND}.mendelian.vcf.gz
