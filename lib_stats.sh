### Genotypes
x=$1
y=${x##*/}
y=${y%.vcf.gz}

## Creates output dir if not existing
mkdir -p input
mkdir -p output

# filter_vep -i $x --filter "PICK is 1" | bgzip -c --threads 200 > input/$y.PICK.vcf.gz 

tabix -p vcf input/$y.PICK.vcf.gz &&
bcftools query -l input/$y.PICK.vcf.gz > input/id.txt &&
bcftools stats input/$y.PICK.vcf.gz -S input/id.txt > output/$y.stats &&

### AF filtering

bcftools filter -i 'AF>0.05' input/$y.PICK.vcf.gz -o output/$1.PICK.common.vcf.gz -O z --threads 200 &&
bcftools filter -i "AF<0.005" input/$y.PICK.vcf.gz -o output/$1.PICK.rare.vcf.gz -O z --threads 200 &&
bcftools filter -i "AF>0.005 & AF<0.05"  input/$y.PICK.vcf.gz -o output/$1.PICK.mediumrare.vcf.gz -O z --threads 200 &&

filter_vep -i $1.PICK.common.vcf.gz --filter "(IMPACT is HIGH)" --force_overwrite | bgzip -c --threads 200 > output/$y.PICK.common.HIGH.vcf.gz 
filter_vep -i $1.PICK.common.vcf.gz --filter "(IMPACT is MODERATE)" --force_overwrite | bgzip -c --threads 200 > output/$y.PICK.common.MODERATE.vcf.gz
filter_vep -i $1.PICK.common.vcf.gz --filter "(IMPACT is LOW)" --force_overwrite | bgzip -c --threads 200 > output/$y.PICK.common.LOW.vcf.gz 

filter_vep -i $1.PICK.rare.vcf.gz --filter "(IMPACT is HIGH)" --force_overwrite | bgzip -c --threads 200 > output/$y.PICK.rare.HIGH.vcf.gz 
filter_vep -i $1.PICK.rare.vcf.gz --filter "(IMPACT is MODERATE)" --force_overwrite | bgzip -c --threads 200 > output/$y.PICK.rare.MODERATE.vcf.gz
filter_vep -i $1.PICK.rare.vcf.gz --filter "(IMPACT is LOW)" --force_overwrite | bgzip -c --threads 200 > output/$y.PICK.rare.LOW.vcf.gz 

filter_vep -i $1.PICK.medium_rare.vcf.gz --filter "(IMPACT is HIGH)" --force_overwrite | bgzip -c --threads 200 > output/$y.PICK.medium_rare.HIGH.vcf.gz 
filter_vep -i $1.PICK.medium_rare.vcf.gz --filter "(IMPACT is MODERATE)" --force_overwrite | bgzip -c --threads 200 > output/$y.PICK.medium_rare.MODERATE.vcf.gz
filter_vep -i $1.PICK.medium_rare.vcf.gz --filter "(IMPACT is LOW)" --force_overwrite | bgzip -c --threads 200 > output/$y.PICK.medium_rare.LOW.vcf.gz 


### ROH
mkdir -p output/ROH
bcftools roh -S input/id.txt input/$y.PICK.vcf.gz -m input/genetic-map/genetic_map_chr{CHROM}_combined_b37.txt -o output/ROH/ROH.txt






