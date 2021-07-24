VCF=$1
s=$2
bcftools view -IUs $s -Ou $VCF \
	| bcftools view -e 'GT="ref"' \
	| bcftools query -f "$s\t%ID\t%CHROM\t%INFO/AC\t%INFO/AN\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/HWE\t%INFO/MSHQ\t%INFO/SU\t%INFO/SR\t%INFO/PE\t[%GT\t%GQ\t%DP\t%AB\t%DHFFC\t%SHQ]\n" \
	| bgzip -c > ${s}.af.tsv.gz
