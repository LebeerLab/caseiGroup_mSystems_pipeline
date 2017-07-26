#!/bin/bash

genomes="
Lactobacillus\srhamnosus
Lactobacillus\scasei
Lactobacillus\sparacasei
Lactobacillus\szeae
"

pin=/media/harddrive/caseiGroup/02_download_genomes/ 
pout_genomes=/media/harddrive/data/genomes/ncbi_caseiGroup_plus_isolates_fna/

cp /dev/null assembly_summary_caseigroup.txt

for genome in $genomes; do
	grep ${genome} assembly_summary.txt >assembly_summary_${genome#*\\s}.txt
	cat assembly_summary_${genome#*\\s}.txt >>assembly_summary_caseigroup.txt
done

cd ${pout_genomes}

for next in $(cat ${pin}assembly_summary_caseigroup.txt | cut -f20); do
 	wget ${next}/${next##*\/}_genomic.fna.gz
done

rename 's/\..*\.fna/.fna/' *.fna.gz
