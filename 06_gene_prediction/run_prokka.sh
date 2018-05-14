#!/bin/bash

pout=/media/harddrive/caseiGroup/annotation/prokka_out/
pin=/media/harddrive/caseiGroup/annotation/

export PATH=$PATH:/media/harddrive/tools/prokka/bin/
export PATH=$PATH:/media/harddrive/tools/barrnap/bin/

cd ${pin}

cut -f1 /media/harddrive/caseiGroup/quality_control/genomeTable.tsv > genomesUsed.tsv

while read genome; do
	prokka /media/harddrive/data/genomes/ncbi_caseiGroup_plus_isolates_fna/${genome}.fna --outdir ${pout}/${genome} --prefix ${genome} --compliant --genus Lactobacillus --usegenus --cpus 6
done < genomesUsed.tsv



