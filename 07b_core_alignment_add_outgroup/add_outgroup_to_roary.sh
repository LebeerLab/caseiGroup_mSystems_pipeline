pin_roary=/media/harddrive/caseiGroup/07a_core_alignment/roary_out/
url_outgroup=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/434/705/GCF_001434705.1_ASM143470v1/GCF_001434705.1_ASM143470v1_genomic.fna.gz
pin_scripts=/media/harddrive/caseiGroup/07b_core_alignment_add_outgroup/add_outgroup_to_roary_scripts/
pout=/media/harddrive/caseiGroup/07b_core_alignment_add_outgroup/add_outgroup_out/

threads=8

mkdir $pout
cd $pout

echo DOWNLOADING OUTGROUP GENOME
echo
wget $url_outgroup
rename 's/\..*\.fna/.fna/' *.fna.gz
echo

echo ANNOTATING OUTGROUP GENOME
echo
export PATH=$PATH:/media/harddrive/tools/prokka/bin/
export PATH=$PATH:/media/harddrive/tools/barrnap/bin/
gunzip *.fna.gz
prokka *.fna --outdir outgroup_annotated --prefix outgroup --compliant \
--genus Lactobacillus --usegenus --cpus $threads
echo

echo EXTRACTING CORE GENE SEQUENCES FROM ROARY OUTPUT
echo
grep locus_tag ${pin_roary}core_alignment_header.embl | cut -d'=' -f2 > core_geneGroups.txt
${pin_scripts}extract_genes_fromPanRefSeqs.py ${pin_roary}pan_genome_reference.fa \
core_geneGroups.txt > core_genes.ffn
rm core_geneGroups.txt
echo 

echo TRANSLATING CORE GENE SEQUENCES TO PROTEINS
echo
transeq -sequence core_genes.ffn -outseq core_genes.faa
echo 

echo MAKING BLAST DATABASE OF OUTGROUP GENOME PROTEINS
echo
makeblastdb -in outgroup_annotated/*.faa -dbtype prot -title outgroup -out outgroup \
-parse_seqids
echo 

echo QUERYING BLAST DATABASE OF OUTGROUP GENOME WITH CORE GENES
echo 
blastp -task blastp -db outgroup -query core_genes.faa -out core_genes_hits.tsv \
-max_target_seqs 1 -num_threads $threads -outfmt '6 qseqid sseqid pident length qcovs'
rm outgroup.p*
echo 

echo GENERATING BLAST LENGTH-SCORE PLOT
echo 
Rscript ${pin_scripts}createplot_addOutgroup.R $pout
echo 

echo RUNNING R SCRIPT THAT ADDS OUTGROUP GENES TO PRESENCE/ABSENCE TABLE
echo
Rscript ${pin_scripts}addOutgroup_to_gene_presence_absence_file.R ${pin_roary}gene_presence_absence.csv core_genes_hits.tsv \
gene_presence_absence_withOutgroup.csv geneGroup_and_outgroupSequence.csv
echo 

echo ADDING OUTGROUP SEQUENCES TO THE GENE GROUP ALIGNMENTS
echo
cp -rf ${pin_roary}pan_genome_sequences/ pan_genome_sequences_withOutgroup/
awk -v FS="," -v OFS="\t" -v outgroup=outgroup_annotated/*.ffn -v pin_roary=${pin_roary} -v pin_scripts=${pin_scripts} '{
         system(pin_scripts"extract_seq.py "outgroup" "$2" > seq.fasta");
         system("mafft --add seq.fasta "pin_roary"pan_genome_sequences/"$1".fa.aln > pan_genome_sequences_withOutgroup/"$1".fa.aln");
        system("rm seq.fasta")
 }' geneGroup_and_outgroupSequence.csv
rm geneGroup_and_outgroupSequence.csv 
for line in $(cat geneGroup_and_outgroupSequence.csv); do
	gene=${line%%","*}
	sequence=${line##*","}
	${pin_scripts}extract_seq.py outgroup_annotated/*.ffn ${sequence} > seq.fasta
	mafft --add seq.fasta ${pin_roary}pan_genome_sequences/${gene}.fa.aln > pan_genome_sequences_withOutgroup/${gene}.fa.aln
	rm seq.fasta
done
echo

echo RUNNING ROARY SCRIPT THAT MAKES CORE ALIGNMENT INCLUDING OUTGROUP
echo
pan_genome_core_alignment -o core_gene_alignment_withOutgroup.aln -cd 96 \
-m pan_genome_sequences_withOutgroup/ \
-s gene_presence_absence_withOutgroup.csv -z -v
echo 

echo COUNTING NUMBER OF N\'S PER GENOME
echo
${pin_scripts}number_of_Ns.py core_gene_alignment_withOutgroup.aln | sort -n -r -k2 > number_of_Ns_per_genome.txt
