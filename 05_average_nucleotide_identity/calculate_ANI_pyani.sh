#!/bin/bash

pin=/media/harddrive/caseiGroup/05_average_nucleotide_identity_ANI/pyani_input
pout=/media/harddrive/caseiGroup/05_average_nucleotide_identity_ANI/pyani_output 

average_nucleotide_identity.py -i $pin -o $pout/ANIb -m ANIb --workers 8 -v
#average_nucleotide_identity.py -i $pin -o $pout/ANIm -m ANIm --workers 8 -v
average_nucleotide_identity.py -i $pin -o $pout/TETRA -m TETRA --workers 8 -v 


