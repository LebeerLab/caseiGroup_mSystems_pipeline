
read_hmmerHitsTable = function(file, parseTargetNames = T) {
  
  # Read table with hmmer hits
  hmmerHits = read_tsv(file, col_names = F)
  names(hmmerHits) = c("target", "query", "e_value")
  
  # Filter hmmer hits: keep only one row for each query (with max score)
  rowsToKeep = match(unique(hmmerHits$query), hmmerHits$query)
  offset = tapply(hmmerHits$e_value, hmmerHits$query, FUN=which.min) - 1
  rowsToKeep = rowsToKeep + offset
  hmmerHits = hmmerHits[rowsToKeep,]
  
  # Extract the useful part of the target names
  if (parseTargetNames) {
    hmmerHits$target = str_split_fixed(hmmerHits$target, "\\.", n=3)[,2]
  }
  
  return(hmmerHits)
  
}

