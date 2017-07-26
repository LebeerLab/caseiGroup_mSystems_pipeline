
library(stringr)
library(ggplot2)

args = commandArgs(trailingOnly=T)
# Following files are required as input:
# - first argument: gene_presence_absence_splittedParalogs.csv
# - second argument: core_genes_hits.tsv
# Following output files need to be specified: 
# - third argument: output file (presence/absence file with outgroup added)
# - fourth argument: output file 2 (gene names and outgroup sequences of gene groups present in outgroup)

# read in the tables
presenceTable = read.table(args[1], header=T, sep=",", stringsAsFactors=F)
outgroupHits = read.table(args[2], sep="\t", stringsAsFactors=F)
names(outgroupHits) = c("query", "hit", "identity", "length", "coverage")

# put cut-off on blastp hits
ggplot(data=outgroupHits, aes(x=coverage, y=identity)) +
  geom_point()
outgroupHits = outgroupHits[outgroupHits$identity>50 & outgroupHits$coverage>75,]

# add outgroup genes to presence absence table 
nrow = nrow(presenceTable)
ncol = ncol(presenceTable)
outgroupHits$presenceRow = match(str_sub(outgroupHits$query, 1, 14), c(as.matrix(presenceTable[,15:ncol])))%%nrow
presenceTable$outgroup = character(nrow)
presenceTable[outgroupHits$presenceRow, "outgroup"] = outgroupHits$hit

# write updated presence absence table 
write.table(presenceTable, file=args[3], 
            row.names=F, quote=T, sep=",")

# write table with gene name and sequence name of gene groups present in outgroup
write.table(presenceTable[presenceTable$outgroup!="",c("Gene", "outgroup")], file=args[4], 
            row.names=F, col.names=F, quote=F, sep=",")
