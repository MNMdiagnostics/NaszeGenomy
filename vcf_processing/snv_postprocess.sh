#
# Takes approx. 8h
#
# filter PASS-variants, remove related individuals
# recalculate INFO/AF (AC and AN are updated automatically)
# remove variants with callrate<90%

date=210716
prefix=${1%".vcf.gz"}

bcftools view -f "PASS,." -S SAMPLES_TO_INCLUDE_${date}.txt -Ou ${1} \
	| bcftools view -e "INFO/AC < 1" -Ou \
	| bcftools +fill-tags -Ou -- -t AF \
	| bcftools view -e 'F_MISSING > 0.1' -Oz -o ${prefix}.unrelated.vcf.gz
tabix -p vcf ${prefix}.unrelated.vcf.gz

#
# Create full, and no-singleton AF databases
#

bcftools view --drop-genotypes -Oz -o ${prefix}.unrelated.nogt.vcf.gz ${prefix}.unrelated.vcf.gz
tabix -p vcf ${prefix}.unrelated.nogt.vcf.gz

bcftools view -e "INFO/AC < 2" -Oz -o ${prefix}.unrelated.nogt.ACgt1.vcf.gz ${prefix}.unrelated.nogt.vcf.gz
tabix -p vcf ${prefix}.unrelated.nogt.ACgt1.vcf.gz
md5sum ${prefix}.unrelated.nogt.ACgt1.vcf.gz > ${prefix}.unrelated.nogt.ACgt1.vcf.gz.md5
