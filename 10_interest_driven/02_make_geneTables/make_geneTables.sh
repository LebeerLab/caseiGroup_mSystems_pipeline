pin_prokka=/media/harddrive/caseiGroup/annotation/prokka_out/
pout=/media/harddrive/caseiGroup/interest_driven/make_geneTables/

cd $pout

echo EXTRACT GENETABLE AND CONTIGTABLE FROM EACH GFF
echo
for file in ${pin_prokka}*/*.gff; do
	genome=${file%.gff}
	genome=${genome##*/}
	grep -v '#' $file | awk -v FS="\t" -v genome=$genome '$3 == "CDS" {print $1,genome,$4,$5,$7,$9}' | cut -d'|' -f3 | cut -d';' -f1 > ${genome}_geneTable.tsv
	grep '##sequence-region' $file | cut -d'|' -f3 > ${genome}_contigTable.tsv
done

echo CAT ALL FILES INTO ONE SUPERFILE
echo
cat *_geneTable.tsv > geneTable_allGenomes.tsv
cat *_contigTable.tsv > contigTable_allGenomes.tsv
echo 
