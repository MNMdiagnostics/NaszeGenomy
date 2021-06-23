prefix=${1%".vcf.gz"}

singularity exec -B /scratch:/scratch \
	/scratch/singularity_images/ensembl-vep_97_4.simg vep --dir_plugins /scratch/References/vep-cache/Plugins_98 \
	--offline --max_sv_size 1000000000 --cache --dir_cache /scratch/References/vep-cache/ \
	--assembly GRCh38 -i $1 -o ${prefix}.vep.tsv.gz --tab \
	--no_stats --no_check_variants_order --force_overwrite --compress_output bgzip \
	--pick --everything --buffer_size 50000
