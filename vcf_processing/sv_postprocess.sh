# select samples
# drop genotypes
# drop annotation and header lines that expose sample IDs
# filter out AC=0 sites

prefix=${1%".vcf.gz"}

bcftools view -S ../SAMPLES_TO_INCLUDE_210716.txt -Ov $1 \
	| grep -v '^##smoove_count_stats=' \
	| grep -v '^##SAMPLE=' \
	| grep -v '^##bcftools_merge'  \
	| bcftools annotate -x INFO/SNAME,INFO/smoove_gene -Ou \
	| bcftools view -e "INFO/AC<1" -Ou \
        | bcftools +fill-tags -Ou -- -t AF,HWE \
	| bcftools filter -m+ -s HWE0 -e 'INFO/HWE<0.01' -Ou \
	| bcftools filter -m+ -s LOW_MSHQ -e 'INFO/MSHQ<3' -Ou \
	| bcftools view -e "INFO/AC<1" -Oz -o ${prefix}.unrelated.vcf.gz \

tabix -p vcf ${prefix}.unrelated.vcf.gz

bcftools view --drop-genotypes -Oz -o ${prefix}.unrelated.nogt.vcf.gz ${prefix}.unrelated.vcf.gz
tabix -p vcf ${prefix}.unrelated.nogt.vcf.gz

bcftools view -e "INFO/AC < 2" -Oz -o ${prefix}.unrelated.nogt.ACgt1.vcf.gz ${prefix}.unrelated.nogt.vcf.gz
tabix -p vcf ${prefix}.unrelated.nogt.ACgt1.vcf.gz
md5sum ${prefix}.unrelated.nogt.ACgt1.vcf.gz > ${prefix}.unrelated.nogt.ACgt1.vcf.gz.md5


