query=/media/harddrive/caseiGroup/orthofinder/ConcatedSequences.fa

targetA=/media/harddrive/caseiGroup/interest_driven/catalase_SOD_orthogroups/HMMscan/SOD/HMMs/PF00080.hmm
targetB=/media/harddrive/caseiGroup/interest_driven/catalase_SOD_orthogroups/HMMscan/SOD/HMMs/PF00081.hmm
targetC=/media/harddrive/caseiGroup/interest_driven/catalase_SOD_orthogroups/HMMscan/SOD/HMMs/PF02777.hmm
targetD=/media/harddrive/caseiGroup/interest_driven/catalase_SOD_orthogroups/HMMscan/SOD/HMMs/PF09055.hmm

pout=/media/harddrive/caseiGroup/interest_driven/catalase_SOD_orthogroups/HMMscan/SOD/

cd $pout 

echo PERFORMING HMMSCAN ON FIRST TARGET 
echo
hmmscan --domtblout targetA.out $targetA $query 
echo

echo PERFORMING HMMSCAN ON SECOND TARGET
echo
hmmscan --domtblout targetB.out $targetB $query 
echo

echo PERFORMING HMMSCAN ON THIRD TARGET
echo
hmmscan --domtblout targetC.out $targetC $query 
echo

echo PERFORMING HMMSCAN ON FOURTH TARGET
echo
hmmscan --domtblout targetD.out $targetD $query 
echo


echo PARSING OUTPUT
echo
./hmmscan-parser.sh targetA.out > targetA.parsed
./hmmscan-parser.sh targetB.out > targetB.parsed
./hmmscan-parser.sh targetC.out > targetC.parsed
./hmmscan-parser.sh targetD.out > targetD.parsed

cat targetA.parsed targetB.parsed > scan.out
rm -f query.fa targetB.out targetA.out targetA.parsed targetB.parsed targetC.out targetC.parsed targetD.out targetD.parsed



