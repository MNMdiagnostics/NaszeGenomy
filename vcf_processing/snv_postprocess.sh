#
# Takes approx. 8h
#
# filter PASS-variants, remove related individuals
# recalculate INFO/AF (AC and AN are updated automatically)
# remove variants with callrate<90%

prefix=${1%".vcf.gz"}

bcftools view -f "PASS,." -S SAMPLES_TO_INCLUDE_210519.txt -Ou ${1} \
	| bcftools view -e "INFO/AC < 1" -Ou \
	| bcftools +fill-tags -Ou -- -t AF \
	| bcftools view -e 'F_MISSING > 0.1' -Oz -o ${prefix}.filtered.vcf.gz

#
# Create full, and no-singleton AF databases
#

bcftools view --drop-genotypes -Oz -o ${prefix}.filtered.nogt.vcf.gz ${prefix}.filtered.vcf.gz

bcftools view -e "INFO/AC < 2" -Oz -o ${prefix}.filtered.nogt.ACgt1.vcf.gz ${prefix}.filtered.nogt.vcf.gz

