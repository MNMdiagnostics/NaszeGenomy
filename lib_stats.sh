### Genotypes
x=$1
y=${x##*/}
y=${y%.vcf.gz}

## Creates output dir if not existing
mkdir -p input
mkdir -p output

filter_vep -i $x --filter "PICK is 1" --force_overwrite | bgzip -c > input/$y.PICK.vcf.gz

tabix -p vcf input/$y.PICK.vcf.gz &&
bcftools query -l input/$y.PICK.vcf.gz > input/id.txt &&
bcftools stats input/$y.PICK.vcf.gz -S input/id.txt > output/$y.stats &&

### AF filtering
bcftools filter -i "AF<0.005" $y.PICK.vcf.gz  | bgzip -c > output/$y.PICK.rare.vcf.gz &&
bcftools filter -i "AF>0.005 & AF<0.05"  $y.PICK.vcf.gz | bgzip -c > output/$y.PICK.mediumrare.vcf.gz &&
bcftools filter -i "AF > 0.05"  $y.PICK.vcf.gz | bgzip -c > output/$y.common.vcf.gz  &&

## Parsing depth files
head -qn 1 /mnt/vault1/mnmorigin_pipeline_out/*/coverage/*wgstats | head -n 1 > input/depth_concat.txt && \
tail -qn +2 /mnt/vault1/mnmorigin_pipeline_out/*/coverage/*wgstats >> input/depth_concat.txt &&

## Parsing flagstat files
multiqc /mnt/covid4/mnt/vault1/mnmorigin_pipeline_out/*/*.flagstat -k tsv -f -o input

## Getting library stats
Rscript --vanilla library_stats.r
