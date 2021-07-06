#!/bin/bash

singularity exec -B /mnt/:/mnt/ -B /tmp:/tmp -B /scratch:/scratch \
	/scratch/singularity_images/ensembl-vep_97_4.simg vep --dir_plugins /scratch/References/vep-cache/Plugins_98 \
	--offline --max_sv_size 1000000000 --cache --dir_cache /scratch/References/vep-cache/ \
	--assembly GRCh38 -i $1 -o $1_annotation.tsv.gz --tab \
	--no_stats --no_check_variants_order --force_overwrite --compress_output bgzip \
	--pick --everything --custom /scratch/References//gnomad.genomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz,gnomADg,vcf,exact,0,AF,AC,AN,AF_NFE \
	--custom /scratch/References//gnomad.genomes.r3.0.sites.noVEP.noAC0_AF_nfe.vcf.gz,gnomAD3g,vcf,exact,0,AF,AC,AN,AF_NFE --fork 64 --buffer_size 50000 \
	--custom /scratch/genom_polaka/input/clinvar210522.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN
