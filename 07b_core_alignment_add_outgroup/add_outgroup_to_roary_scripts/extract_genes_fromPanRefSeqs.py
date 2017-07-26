#!/usr/bin/python

from Bio import SeqIO
import sys

fastaFile = sys.argv[1]
wantedGenesFile = sys.argv[2]

wantedGenes = [line.strip() for line in open(wantedGenesFile).readlines()]

genes = list(SeqIO.parse(fastaFile, 'fasta'))
genes = [gene for gene in genes if gene.description.split(" ")[1] in wantedGenes]

SeqIO.write(genes, sys.stdout, 'fasta')
