bcftools query -f '%ID\t%AF\t%AN\t%AC\n' ../../input/multisample_20210716.dv.bcfnorm.filt.unrelated.vcf.gz \
 > ../input/multisample_20210519.dv.bcfnorm.filtered.AFlist.txt

Rscript disease_variant_prep.R


