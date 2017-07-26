#/bin/bash

pout=/media/harddrive/caseiGroup/interest_driven/map_orthogroups_to_eggnogFam/
pin=/media/harddrive/caseiGroup/interest_driven/extract_refSeqs_fromOrthofinder/

cd $pout

if false; then

echo DOWNLOADING HMMS OF EGGNOG BACILLUS COGS
echo
wget http://eggnogdb.embl.de/download/eggnog_4.5/data/bacNOG/bacNOG.hmm.tar.gz
tar xvzf bacNOG.hmm.tar.gz
echo

echo DOWNLOADING FUNCTIONAL ANNOTATIONS OF EGGNOG BACILLUS COGS
echo
wget http://eggnogdb.embl.de/download/eggnog_4.5/data/bacNOG/bacNOG.annotations.tsv.gz
gunzip bacNOG.annotations.tsv.gz
echo 

echo MAKING HMM DATABASE OF BACILLUS COGS
echo
cat bacNOG_hmm/*.hmm >bacDB.hmmer
hmmpress bacDB.hmmer
echo

fi 

echo SCANNING THE HMM DATABASE USING THE ORTHOFINDER REPRESENTATIVE SEQUENCES
echo
hmmscan --tblout hmmer_hits.tsv bacDB.hmmer ${pin}pan_genome_reference.fasta
tr -s ' ' < hmmer_hits.tsv | grep -v '^#' | tr ' ' '\t' | \
cut -f 1,3,5 > hmmer_hits_essential.tsv
echo

