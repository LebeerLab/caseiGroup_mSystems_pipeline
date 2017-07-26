pout=/media/harddrive/caseiGroup/03_add_unclassified/blastagainstrdp/ 
p_genomes=/media/harddrive/data/genomes/lactobacillus_unclassified_fna/

cd $pout 

echo MAKING BLAST DATABASE FROM RDP file
echo
makeblastdb -in rdp_20170219_Lactobacillus_nogaps.fasta -dbtype nucl -title rdp -out rdp -parse_seqids
echo

echo QUERYING BLAST DATABASE 
for next in $(cat ${pout}../assembly_summary_sp.txt | cut -f1); do
        echo PROCESSING ${next%.*}
	blastn -db rdp -query ${p_genomes}${next%.*}.fna -out out/${next%.*}.tsv -num_threads 8 -outfmt '6 qseqid sseqid pident length qcovs'

done

rm rdp.*

echo
echo EXTRACTING BEST HITS
echo 
head -n1 out/* > besthits.tsv
