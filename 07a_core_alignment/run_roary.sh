# pin needs to contain folders with .gff files
pin=/media/harddrive/caseiGroup/06_gene_prediction/prokka_out/
pout=/media/harddrive/caseiGroup/07a_core_alignment/roary_out/

echo STARTING ROARY WITH SPLIT PARALOGS ENABLED
echo CORE ALIGNMENT IS ALSO ENABLED, ELSE ROARY DOESN\'T OUTPUT PAN_GENOME_SEQUENCES 
echo 
roary -e --mafft -cd 96 -i 70 -p 8 -v -r -z -f $pout -o ${pout}clusters \
${pin}*/*.gff

