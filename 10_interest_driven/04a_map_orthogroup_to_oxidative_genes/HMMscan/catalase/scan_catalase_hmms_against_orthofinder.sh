query=/media/harddrive/caseiGroup/orthofinder/ConcatedSequences.fa

targetA=/media/harddrive/caseiGroup/interest_driven/catalase_SOD_orthogroups/HMMscan/catalase/HMMs/catalase.hmm
targetB=/media/harddrive/caseiGroup/interest_driven/catalase_SOD_orthogroups/HMMscan/catalase/HMMs/Mn_catalase.hmm

pout=/media/harddrive/caseiGroup/interest_driven/catalase_SOD_orthogroups/HMMscan/catalase

cd $pout 


echo PERFORMING HMMSCAN ON FIRST TARGET 
echo
hmmscan --domtblout targetA.out $targetA $query 
echo

echo PERFORMING HMMSCAN ON SECOND TARGET
echo
hmmscan --domtblout targetB.out $targetB $query 
echo

echo PARSING OUTPUT
echo
./hmmscan-parser.sh targetA.out > targetA.parsed
./hmmscan-parser.sh targetB.out > targetB.parsed
cat targetA.parsed targetB.parsed > scan_catalase.out
rm -f query.fa targetB.out targetA.out targetA.parsed targetB.parsed



