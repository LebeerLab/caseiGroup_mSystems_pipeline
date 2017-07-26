query=/media/harddrive/caseiGroup/blast_against_orthofinder_out/query.fasta
target=/media/harddrive/caseiGroup/orthofinder/ConcatedSequences.fa
pout=/media/harddrive/caseiGroup/blast_against_orthofinder_out/

cd $pout 

echo TRANSLATING QUERY
echo
transeq $query query.faa
echo

echo MAKING BLAST DATABASE ORTHOFINDER OUTPUT
echo
makeblastdb -in $target -dbtype prot -title target -out target -parse_seqids
echo

echo QUERYING BLAST DATABASE 
echo
blastp -task blastp -db target -query query.faa -out genes_hits.tsv -num_threads 6 -outfmt '6 qseqid sseqid pident length qcovs'
rm target.p*
echo

rm -f query.faa target.faa 


