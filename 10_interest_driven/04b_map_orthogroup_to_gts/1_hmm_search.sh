#!/bin/bash

pin=/media/harddrive/caseiGroup/interest_driven/extract_refSeqs_fromOrthofinder/
pout=/media/harddrive/caseiGroup/interest_driven/map_orthogroup_to_cazyFam

cd $pout

echo SCANNING THE HMM DATABASE USING THE ORTHOFINDER REPRESENTATIVE SEQUENCES
echo
hmmscan --tblout hmmer_hits.tsv HMMs_dbCAN_pfam.txt ${pin}pan_genome_reference.fasta
tr -s ' ' < hmmer_hits.tsv | grep -v '^#' | tr ' ' '\t' | \
cut -f 1,3,5 > hmmer_hits_essential.tsv
echo

