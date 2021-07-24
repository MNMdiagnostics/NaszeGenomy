input=$1
prefix=$(basename $1)

bcftools query -f '%ID\t%AF\t%AN\t%AC\n' $input > input/$prefix.AFlist.txt



