bcftools query -f '%CHROM \t %POS \t %REF \t %ALT \t %TYPE  \t [%VF] \t %CSQ \n' U1.vep.ann.vcf | tr '|' '\t'


| tr '|' ' ' \

| cut -d' ' -f 1-4,7 \
| grep 'UTR' \
| head \
| column -t

 \t [%VF] \t %CSQ