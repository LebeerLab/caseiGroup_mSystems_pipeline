
library(tidyverse)
library(stringr)
source("function_read_hmmerHitsTable.R")

# construct T_gene

T_gene1 <- read_delim(delim = " ", file = "input/geneTable_allGenomes.tsv", col_names = F) %>%
  rename(contig = X1, genome = X2, start = X3, stop = X4, direction = X5, gene = X6) %>%
  mutate(gene = str_split_fixed(gene, "=", n = 2)[,2])

T_gene2_genome_orthogroup <- read_tsv(file = "input/Orthogroups.csv", col_names = T) %>%
  rename(orthogroup = X1) %>%
  gather(key = genome, value = genes, -orthogroup) %>%
  filter(! is.na(genes) ) %>%
  mutate(ngenes = str_count(genes, ",") + 1)

T_gene2_oneGene <- T_gene2_genome_orthogroup %>%
  filter(ngenes == 1) %>%
  select(-ngenes) %>%
  rename(gene = genes)

T_gene2_moreGenes <- T_gene2_genome_orthogroup %>%
  filter(ngenes > 1) %>%
  select(-ngenes) %>%
  split(f = 1:nrow(.)) %>%
  lapply(FUN = function(row) {
    gene = str_split(row$genes, pattern = ", ")[[1]]
    return(data.frame(orthogroup = row$orthogroup,
                      genome = row$genome,
                      gene = gene,
                      stringsAsFactors = F))
  }) %>% bind_rows()

T_gene2_genome_orthogroup <- rbind(T_gene2_oneGene, T_gene2_moreGenes)
T_gene <- left_join(T_gene1, T_gene2_genome_orthogroup)
save(T_gene, file = "intermediate/T_gene.Robject")

# construct T_contig

T_contig <- read_delim(delim = " ", file = "input/contigTable_allGenomes.tsv", col_names = F) %>%
  rename(contig = X1, unknown = X2, length = X3) %>%
  select(- unknown)
save(T_contig, file = "intermediate/T_contig.Robject")

# construct T_genome

T_genome <- read_tsv(file = "input/genomeTableWithClades.tsv", col_names = F) %>%
  rename(genome = X1, species = X2, clade = X3, strain = X4) %>%
  mutate(strain = ifelse(genome == "AMBR2", "AMBR2", strain))
save(T_genome, file = "intermediate/T_genome.Robject")

# construct T_orthogroup

T_orthogroup1_cazyFam <- read_hmmerHitsTable(file = "input/map_orthogroup_to_cazyFam/hmmer_hits_essential.tsv",
                                           parseTargetNames = F) %>%
  rename(cazyFam = target, orthogroup = query) %>%
  filter(e_value < 1e-18) %>% # threshold taken from dbCAN website
  select(-e_value)

T_orthogroup2_eggnogFam <- read_hmmerHitsTable(file = "input/map_orthogroup_to_eggnogFam/hmmer_hits_essential.tsv") %>%
  rename(eggnogFam = target, orthogroup = query) %>%
  select(-e_value)

T_orthogroup <- full_join(T_orthogroup1_cazyFam, T_orthogroup2_eggnogFam)
save(T_orthogroup, file = "intermediate/T_orthogroup.Robject")

# construct T_eggnogFam

T_eggnogFam <- read_tsv(file = "input/bacNOG.annotations.tsv", col_names = F) %>%
  select(X2, X6) %>%
  rename(eggnogFam = X2, eggnogFunction = X6)
save(T_eggnogFam, file = "intermediate/T_eggnogFam.Robject")

# construct T_eggnogFam_functCat

T_eggnogFam_functCat <- read_tsv(file = "input/bacNOG.annotations.tsv", col_names = F) %>%
  select(X2, X5) %>%
  rename(eggnogFam = X2, functCat = X5) %>%
  split(1:nrow(.)) %>%
  lapply(FUN = function(row) {
    if (str_length(row$functCat) == 1) return(row)
    data.frame(eggnogFam = row$eggnogFam,
               functCat = str_split(row$functCat, "")[[1]]) %>%
      return()
  }) %>% bind_rows() 
save(T_eggnogFam_functCat, file = "intermediate/T_eggnogFam_functCat.Robject")
