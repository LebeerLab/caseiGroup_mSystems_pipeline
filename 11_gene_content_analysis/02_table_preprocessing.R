
library(tidyverse)
library(stringr)
library(ape)
source("function_order_contigs.R")

# load necessary tables and tree

load(file = "intermediate/T_gene.Robject")
load(file = "intermediate/T_contig.Robject")
load(file = "intermediate/T_genome.Robject")
tree <- read.tree('input/RAxML_bipartitions.caseigroup')

# remove "empty space" at beginning and endings of contigs

T_gene <- T_gene %>%
  group_by(contig) %>%
  mutate(contigStart = min(start)) %>%
  ungroup() %>%
  mutate(start = start - contigStart + 1000) %>%
  mutate(stop = stop - contigStart + 1000) 

T_contig <- T_gene %>%
  group_by(contig) %>%
  summarize(length = max(stop) + 1000)

# add "supercontig" variable to T_contig

T_contig <- T_contig %>%
  mutate(supercontig = contig) %>%
  mutate(positionInSupercontig = 1) 

# order contigs of clade A with ATCC 334 as reference genome (takes some minutes)

cladeAGenomes <- T_genome %>%
	filter(clade == "cladeA") %>%
	pull(genome)

for (genome in cladeAGenomes) {
	output <- order_contigs(T_gene, T_contig,
												 refGenome = "GCA_000014525", targetGenome = genome)
	T_gene <- output[[1]]
	T_contig <- output[[2]]
}

# order contigs of clade B with AMBR2 as reference genome (takes some minutes)

cladeBGenomes <- T_genome %>%
  filter(clade == "cladeB") %>%
  pull(genome)

for (genome in cladeBGenomes) {
  output <- order_contigs(T_gene, T_contig,
                           refGenome = "AMBR2", targetGenome = genome)
  T_gene <- output[[1]]
  T_contig <- output[[2]]
}

# order contigs of clade C with LGG as reference genome (takes some minutes)

cladeCGenomes <- T_genome %>%
  filter(clade == "cladeC") %>%
  pull(genome)

for (genome in cladeCGenomes) {
  output <- order_contigs(T_gene, T_contig,
                         refGenome = "GCA_000026505", targetGenome = genome)
  T_gene <- output[[1]]
  T_contig <- output[[2]]
}

# add coordinates to contig and gene tables

V_contig_withGenome <- T_gene %>%
  select(contig, genome) %>%
  distinct() %>%
  right_join(T_contig)

T_contig_withCo <- V_contig_withGenome %>%
  group_by(genome) %>%
  arrange(supercontig, positionInSupercontig) %>%
  mutate(contig_start = c(0, cumsum(length)[-length(length)]) + 1) %>%
  ungroup() %>%
  select(contig, length, supercontig, positionInSupercontig, contig_start)
save(T_contig_withCo, file = "intermediate/T_contig_withCo.Robject")

T_gene_withCo <- left_join(T_gene, T_contig_withCo) %>%
  mutate(start_inGenome = contig_start + start - 1) %>%
  mutate(stop_inGenome = contig_start + stop - 1) %>%
  select(-length, -contig_start)
save(T_gene_withCo, file = "intermediate/T_gene_withCo.Robject")

# add tree_position to T_genome

tree <- tree %>%
  ape::root(outgroup = "outgroup") %>%
  drop.tip("outgroup")

genomes_with_tree_position <- tibble(tip_label_index = tree$edge[,2]) %>%
  filter(tip_label_index <= length(tree$tip.label)) %>%
  mutate(
    tree_position = 1:n(),
    genome = tree$tip.label[tip_label_index]
  ) %>%
  select(- tip_label_index)

T_genome <- T_genome %>%
  left_join(genomes_with_tree_position)

# make tree-ordered versions of strain and genome

T_genome <- T_genome %>%
  mutate(
    strain_ordered = factor(strain, levels = strain[order(tree_position)]),
    genome_ordered = factor(genome, levels = genome[order(tree_position)])
  )

# make strain names unique

T_genome <- T_genome %>%
  split(f = .$strain) %>%
  lapply(FUN = function(group) within(group, {
    n <- length(strain)
    if (n > 1) {
      strain <- str_c(strain, 1:n, sep = " dupl ")
    }
    rm(n)
  })) %>%
  do.call(what = rbind) 

save(T_genome, file = "intermediate/T_genome.Robject")

# make "view" with number of genes of each orthogroup in each genome

V_genome_orthogroup = T_gene %>%
  select(gene, genome, orthogroup) %>%
  group_by(genome, orthogroup) %>%
  summarise(ngenes = n()) %>%
  ungroup()
save(V_genome_orthogroup, file = "intermediate/V_genome_orthogroup.Robject")
