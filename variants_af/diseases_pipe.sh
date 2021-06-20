bcftools query -f '%ID\t%AF\t%AN\t%AC\n' /genom_polaka/input/multisample_20210519.dv.bcfnorm.filtered.vcf.gz
 > input/multisample_20210519.dv.bcfnorm.filtered.AFlist.txt

Rscript disease_variant_prep.R


