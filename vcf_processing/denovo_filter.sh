VCF=$1
OUT1=${VCF%'vcf.gz'}"filt_snv.vcf.gz"
OUT2=${VCF%'vcf.gz'}"filt_indel.vcf.gz"

# require variant to have sufficient depth in all samples, including the alternate allele in the proband (AD[*:1]
# require variant to be present in the proband (3rd samples) with reasonable allelic balance (0.25), and not a single read in the parents

MUM_FILTER="(FORMAT/DP[0:0]>9 & FORMAT/AD[0:1]==0)"
DAD_FILTER="(FORMAT/DP[1:0]>9 & FORMAT/AD[1:1]==0)"
CHD_FILTER="(FORMAT/DP[2:0]>9 & FORMAT/AD[2:1]>4 & FORMAT/VAF[2:0]>0.25)"

SNV_FILTER='(INFO/AF<0.01 & TYPE="snp")'
INDEL_FILTER='(INFO/AF<0.01 & TYPE="indel")'

#echo "'${MUM_FILTER} && ${DAD_FILTER} && ${CHD_FILTER}'"

bcftools +fill-tags -Ou $VCF -- -t FORMAT/VAF \
       | bcftools filter -Ou -i '(FORMAT/DP[0:0]>9 & FORMAT/AD[0:1]==0) && (FORMAT/DP[1:0]>9 & FORMAT/AD[1:1]==0) && (FORMAT/DP[2:0]>9 & FORMAT/AD[2:1]>4 & FORMAT/VAF[2:0]>0.25)' \
       | bcftools filter -Oz -o $OUT1 -i '(INFO/AF<0.01 & TYPE="snp")'
    
bcftools +fill-tags -Ou $VCF -- -t FORMAT/VAF \
       | bcftools filter -Ou -i '(FORMAT/DP[0:0]>9 & FORMAT/AD[0:1]==0) && (FORMAT/DP[1:0]>9 & FORMAT/AD[1:1]==0) && (FORMAT/DP[2:0]>9 & FORMAT/AD[2:1]>4 & FORMAT/VAF[2:0]>0.25)' \
       | bcftools filter -Oz -o $OUT2 -i '(INFO/AF<0.01 & TYPE="indel")'


