#!/bin/bash
alias bcftools=/scratch/tools/bcftools-1.12/bcftools
export BCFTOOLS_PLUGINS=/scratch/tools/bcftools-1.12/plugins


### Input
x=$1
y=${x##*/}
y=${y%.vcf.gz}

### AC AN tags
bcftools +fill-tags $x -Ob  -- -t AN,AC | \
### Dropping genotypes
bcftools view  --drop-genotypes -Oz |  \
### Removing singletons
bcftools filter -e "INFO/AC[0]<2" | \
### Fixing variants ID
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' -o $y.ACAFAN.ACgt1.ID.vcf.gz

