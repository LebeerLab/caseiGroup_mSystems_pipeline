#!/bin/bash

fin=/media/harddrive/caseiGroup/03_add_unclassified/blastagainstrdp/GenomestoUse.txt
pin_genomes=/media/harddrive/data/genomes/lactobacillus_unclassified_fna/
pout_genomes=/media/harddrive/data/genomes/ncbi_caseiGroup_plus_isolates_fna/

while read genome; do
      cp ${pin_genomes}${genome}.fna $pout_genomes
done < $fin

