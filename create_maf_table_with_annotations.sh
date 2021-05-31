#!/bin/bash


VEP_TABLE=multisample_${date}.dv.norm.tsv.gz

bcftools query -f "%INFO/AF\n" /mnt/vault1/mnmorigin_pipeline_out/multisample_${date}.dv.norm.ACANAF.vcf.gz | paste ${VEP_TABLE} - > ${VEP_TABLE}.MAF.tsv

