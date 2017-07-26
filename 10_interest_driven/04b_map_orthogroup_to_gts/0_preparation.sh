wget http://csbl.bmb.uga.edu/dbCAN/download/dbCAN-fam-HMMs.txt
wget http://csbl.bmb.uga.edu/dbCAN/download/hmmscan-parser.sh

hmmconvert -a dbCAN-fam-HMMs.txt > dbCAN-fam-HMMs_converted.txt
fromdos HMMs_pfam/*
cat HMMs_pfam/* dbCAN-fam-HMMs_converted.txt > HMMs_dbCAN_pfam.txt
hmmpress HMMs_dbCAN_pfam.txt
