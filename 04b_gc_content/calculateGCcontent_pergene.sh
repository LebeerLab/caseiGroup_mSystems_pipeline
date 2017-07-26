#!/bin/bash


pin=/media/harddrive/caseiGroup/06_gene_prediction/prokka_out/
pout=/media/harddrive/caseiGroup/04b_gc_content/outfiles/

for ffn in $pin*/*.ffn
do
	infoseq -auto -only -name -length -pgc -description -columns N -delimiter '\t' $ffn > $pout${ffn##*/}.gccontent.tsv
done
