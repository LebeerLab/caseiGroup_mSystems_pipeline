
library(tidyverse)
source("function_reverseContig.R")

make_clusterTable = function(V_gene_large, orthogroup, before = 4000, after = 60000) {
  
  V_gene_large %>%
    split(f = .$genome) %>%
    lapply(FUN = function(genomeTable) {
      isOrthogroup = ! is.na(genomeTable$orthogroup) & genomeTable$orthogroup == orthogroup
      if (sum(isOrthogroup) != 1) return(NULL)
      rightSupercontig = genomeTable$supercontig[isOrthogroup]
      plus = genomeTable$direction[isOrthogroup] == "+"
      if (plus) {
        clustBound = c(genomeTable$start_inGenome[isOrthogroup] - before, 
                       genomeTable$stop_inGenome[isOrthogroup] + after) 
      } else {
        clustBound = c(genomeTable$start_inGenome[isOrthogroup] - after, 
                       genomeTable$stop_inGenome[isOrthogroup] + before) 
      }
      genomeTable = genomeTable %>%
        filter(supercontig == rightSupercontig) %>%
        filter(pmin(start_inGenome, stop_inGenome) >= min(clustBound)) %>%
        filter(pmax(start_inGenome, stop_inGenome) <= max(clustBound))
      if (plus) {
        genomeTable = genomeTable %>%
          mutate(start = start_inGenome - min(clustBound)) %>%
          mutate(contig_start = contig_start - min(clustBound)) %>%
          mutate(contig_start = ifelse(contig_start < 0, NA, contig_start)) %>%
          mutate(stop = stop_inGenome - min(clustBound))
      } else {
        genomeTable = genomeTable %>%
          mutate(start = max(clustBound) - stop_inGenome) %>%
          mutate(contig_start = max(clustBound) - contig_start + length) %>%
          mutate(contig_start = ifelse(contig_start < 0, NA, contig_start)) %>%
          mutate(stop = max(clustBound) - start_inGenome) %>%
          mutate(direction = ifelse(direction == "+", "-", "+"))
      }
      return(genomeTable)
    }) %>% bind_rows() %>%
    return()
  
}