#!/bin/bash

genomes="
Lactobacillus\ssp.
"

pin=/media/harddrive/caseiGroup/03_add_Unclassified/
pout_genomes=/media/harddrive/data/genomes/lactobacillus_unclassified_fna/

cd ${pin}

for genome in $genomes; do
	grep ${genome} assembly_summary.txt >assembly_summary_${genome#*\\s}txt
done

cd ${pout_genomes}

for next in $(cat ${pin}assembly_summary_sp.txt | cut -f20); do
 	wget ${next}/${next##*\/}_genomic.fna.gz
done

rename 's/\..*\.fna/.fna/' *.fna.gz
