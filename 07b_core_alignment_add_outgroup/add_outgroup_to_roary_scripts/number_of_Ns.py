#!/usr/bin/python

from Bio import SeqIO
import sys

seqs = list(SeqIO.parse(sys.argv[1], 'fasta'))

for seq in seqs:
	name, sequence = seq.id, str(seq.seq)
	print name, sequence.count('N')
