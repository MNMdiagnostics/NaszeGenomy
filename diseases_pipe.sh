bcftools query -f "%CHROM %POS %END %AF\n" /mnt/vault1/mnmorigin_pipeline_out/multisample_20210519.dv.norm.vcf.gz > input/AFlist.txt

$1 = multisample_20210519.dv.norm.vep.tsv.gz

 filter_vep  -i /mnt/vault1/mnmorigin_pipeline_out/multisample_20210519.dv.norm.vep.tsv.gz \
 --filter "Gene in input/diseases/Ensembl_ACMG_v3.txt and IMPACT is HIGH" -o input/diseases/acmg_filtered.tsv --force_overwrite &

 filter_vep  -i /mnt/vault1/mnmorigin_pipeline_out/multisample_20210519.dv.norm.vep.tsv.gz \
   --filter "CLIN_SIG is pathogenic" -o input/diseases/clinsig_filtered.tsv --force_overwrite &

 filter_vep  -i /mnt/vault1/mnmorigin_pipeline_out/multisample_20210519.dv.norm.vep.tsv.gz \
  --filter "(IMPACT is HIGH or MODERATE) and AF < 0.005" -o input/diseases/putative_filtered.tsv --force_overwrite &

  wait

  Rscript disease_variant_prep.R

 filter_vep  -i /mnt/vault1/mnmorigin_pipeline_out/multisample_20210519.dv.norm.vep.tsv.gz \
 --filter "AF<0.005" -o input/diseases/ultra_rare.tsv --force_overwrite &

 filter_vep  -i /mnt/vault1/mnmorigin_pipeline_out/multisample_20210519.dv.norm.vep.tsv.gz \
 --filter " AF>0.001 and AF<0.005" -o input/diseases/medium_rare.tsv --force_overwrite &
 
 filter_vep  -i /mnt/vault1/mnmorigin_pipeline_out/multisample_20210519.dv.norm.vep.tsv.gz \
 --filter "AF>0.005" -o input/diseases/common.tsv --force_overwrite 




