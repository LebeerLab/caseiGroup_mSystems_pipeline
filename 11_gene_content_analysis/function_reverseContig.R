
library(tidyverse)

reverse_contig = function(T_gene, T_contig, contigToReverse) {
  
  geneSplittingVar = rep("other", times = nrow(T_gene))
  geneSplittingVar[T_gene$contig == contigToReverse] = "target"
  T_gene = split(T_gene, f = geneSplittingVar)
  
  contigLength = T_contig %>%
    filter(contig == contigToReverse) %>%
    .$length
  
  T_gene$target = T_gene$target %>%
    mutate(length = contigLength) %>%
    reverse_genes() %>%
    select(- length)
  
  T_gene = bind_rows(T_gene)
  
  return(T_gene)
  
}

# assumes that T_gene_withLenth has a "length" variable with 
# the length of the contig
reverse_genes = function(T_gene_withLength) {
  
  T_gene_withLength = T_gene_withLength %>%
    mutate(direction = ifelse(direction == "+", "-", "+")) %>%
    mutate(startNew = length - stop + 1) %>%
    mutate(stopNew = length - start + 1) %>%
    mutate(start = startNew, stop = stopNew) %>%
    select(- startNew, - stopNew) %>%
    return()
  
}