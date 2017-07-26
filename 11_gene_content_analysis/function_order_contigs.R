
library(tidyverse)
source("function_reverseContig.R")

order_contigs = function(T_gene, T_contig, refGenome, targetGenome) {
  
  if (refGenome == targetGenome) return(list(T_gene, T_contig))
  
  # PREPROCESSING
  
  geneSplittingVar = rep("other", times = nrow(T_gene))
  geneSplittingVar[T_gene$genome == refGenome] = "reference"
  geneSplittingVar[T_gene$genome == targetGenome] = "target"
  T_gene = split(T_gene, f = geneSplittingVar)
  
  referenceContigs = unique(T_gene$reference$contig)
  targetContigs = unique(T_gene$target$contig)
  
  contigSplittingVar = rep("other", times = nrow(T_contig))
  contigSplittingVar[T_contig$contig %in% referenceContigs] = "reference"
  contigSplittingVar[T_contig$contig %in% targetContigs] = "target"
  T_contig = split(T_contig, f = contigSplittingVar)
  
  # T_contig$reference = T_contig$reference %>%
  #   left_join(select(T_gene$reference, contig, gene)) %>%
  #   group_by(contig) %>%
  #   mutate(ngenes = n()) %>%
  #   ungroup() %>%
  #   select(- gene) %>%
  #   distinct()
  # 
  # T_contig$target = T_contig$target %>%
  #   left_join(select(T_gene$target, contig, gene)) %>%
  #   group_by(contig) %>%
  #   mutate(ngenes = n()) %>%
  #   ungroup() %>%
  #   select(- gene) %>%
  #   distinct()
  
  T_gene$reference = T_gene$reference %>%
    group_by(contig) %>%
    arrange(start) %>%
    mutate(position = 1:n()) %>%
    ungroup()
  
  # CONTIG-CONTIG COMPARISON
  
  contigTable = merge(T_contig$reference %>%
                        select(contig) %>%
                        rename(contig_reference = contig),
                      T_contig$target %>%
                        select(contig) %>%
                        rename(contig_target = contig)) %>%
    rowwise() %>%
    mutate(commonGenes = count_nCommonOrthogroups(T_gene, contig_reference, contig_target)) %>%
    ungroup()
  
  # SUPERCONTIG ASSIGNMENT
  
  T_mapping = contigTable %>%
    filter(commonGenes != 0) %>%
    group_by(contig_target) %>%
    filter(commonGenes == max(commonGenes)) %>%
    filter(n() == 1) %>%
    ungroup() %>%
    rename(contig = contig_target)
  
  T_contig$target = T_contig$target %>%
    left_join(T_mapping) %>%
    mutate(supercontig = ifelse(is.na(contig_reference), 
                                supercontig, contig_reference)) %>%
    select(- contig_reference)
  
  # POSITION IN SUPERCONTIG UPDATE
  
  T_contig$target = T_contig$target %>%
    group_by(supercontig) %>%
    mutate(nTargetsOnSupercontig = n()) %>%
    ungroup()
  
  for (row in 1:nrow(T_contig$target)) {
    
    nTargetsOnSupercontig = T_contig$target[row, "nTargetsOnSupercontig"]
    if(is.na(nTargetsOnSupercontig) | nTargetsOnSupercontig <= 1) next
    
    output = map_targetToReference(T_gene$target, T_gene$reference, T_contig$target,
                                   T_contig$target$contig[row],
                                   T_contig$target$supercontig[row])
    
    T_gene$target = output[[1]]
    
    T_contig$target[row, "positionInSupercontig"] = output[[2]]
    
  }
  
  # POSTPROCESSING
  
  T_gene$reference = select(T_gene$reference, - position)
  T_contig$target = select(T_contig$target, - commonGenes, - nTargetsOnSupercontig)
  
  T_gene = bind_rows(T_gene)
  T_contig = bind_rows(T_contig)
  
  return(list(T_gene, T_contig))
  
}

count_nCommonOrthogroups = function(T_gene, refContig, targetContig) {
  
  # rem.: not using "unique" makes this slightly faster
  referenceOrthogroups = T_gene$reference %>%
    filter(contig == refContig) %>%
    .$orthogroup
  targetOrthogroups = T_gene$target %>%
    filter(contig == targetContig) %>%
    .$orthogroup
  
  intersect(referenceOrthogroups, targetOrthogroups) %>%
    length() %>%
    return()
  
}

map_targetToReference = function(T_gene_target, T_gene_reference, T_contig_target, 
                                 targetContig, refContig) {
  
  T_gene_targetContig = T_gene_target %>%
    filter(contig == targetContig)
  targetOrthogroups = T_gene_targetContig$orthogroup
  
  T_gene_refContig = T_gene_reference %>%
    filter(contig == refContig)
  refOrthogroups = T_gene_refContig$orthogroup  
  
  commonOrthogroups = intersect(targetOrthogroups, refOrthogroups)
  
  T_gene_targetContig_matching = T_gene_targetContig %>%
    filter(orthogroup %in% commonOrthogroups)
  targetDirection = T_gene_targetContig_matching %>%
    filter(orthogroup == commonOrthogroups[1]) %>%
    .$direction %>%
    .[1]
  # targetDirectionBalance = T_gene_targetContig_matching %>%
  #   summarize(balance = sum(direction == "+") - sum(direction == "-")) %>%
  #   .$balance    
  
  T_gene_refContig_matching = T_gene_refContig %>%
    filter(orthogroup %in% commonOrthogroups)
  refDirection = T_gene_refContig_matching %>%
    filter(orthogroup == commonOrthogroups[1]) %>%
    .$direction %>%
    .[1]
  # refDirectionBalance = T_gene_refContig_matching %>%
  #   summarize(balance = sum(direction == "+") - sum(direction == "-")) %>%
  #   .$balance 
    
  positionInReference = T_gene_refContig_matching %>%
    .$position %>%
    median()
  
  if (targetDirection != refDirection) {
    T_gene_target = reverse_contig(T_gene_target, T_contig_target,
                                   contigToReverse = targetContig)
  }
  
  return(list(T_gene_target, positionInReference))
  
}