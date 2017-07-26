genomeFolder=/media/harddrive/data/genomes/ncbi_caseiGroup_plus_isolates_fna/
quastPath=/media/harddrive/tools/quast-4.3/quast.py

python $quastPath -o qc_quast_allGenomes_results --threads 6 ${genomeFolder}*.fna*

