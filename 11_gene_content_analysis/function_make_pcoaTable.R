
make_pcoaTable = function(T_genome_orthogroup) {
  
  T_genome_orthogroup %>%
    select(orthogroup, genome, ngenes) %>%
    spread(key = orthogroup, value = ngenes, fill = 0) %>%
    remove_rownames() %>%
    column_to_rownames(var = "genome") %>%
    as.matrix() %>%
    vegdist(method = "bray") %>%
    cmdscale(k = 2) %>%
    as.data.frame() %>%
    rownames_to_column(var = "genome") %>%
    return()
  
}