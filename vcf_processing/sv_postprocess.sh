# select samples
# drop genotypes
# drop annotation and header lines that expose sample IDs
# filter out AC=0 sites

prefix=${1%".vcf.gz"}

bcftools view -S ../SAMPLES_TO_INCLUDE_210519.txt --drop-genotypes -Ou $1 \
	| bcftools annotate -Ov -x INFO/SNAME \
	| grep -v '^##smoove_count_stats=' \
	| grep -v '^##SAMPLE=' \
	| grep -v '^##bcftools_merge'  \
	| bcftools view -e "INFO/AC<1" -Oz -o ${prefix}.unrelated.nogt.vcf.gz
