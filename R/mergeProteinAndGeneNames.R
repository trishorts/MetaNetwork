mergeAndWriteWGCNAWorkbook <- function(selectedDatabase, 
                                     allData, 
                                     dataList){
  userInputDatabaseSelectedColumns <- tibble(userInputDatabase$Entry, 
                                             userInputDatabase$`Gene names`, 
                                             userInputDatabase$`Protein names`)
  colnames(userInputDatabaseSelectedColumns) <- c("accession", "gene name", "protein name")
  dat.resMerged <- left_join(tibble(allData), 
                             userInputDatabaseSelectedColumns, 
                             by = "accession")
  dat.resMerged <- tibble(allData)
  list.cluster.datMerged <- list()
  for(i in seq_along(dataList)){
    list.cluster.datMerged[[i]] <- left_join(dataList[[i]], 
                                             userInputDatabaseSelectedColumns, 
                                             by = "accession")
  }
  names(list.cluster.datMerged) <- names(list.cluster.dat)
  createResultsWGCNAExcelWorkbook(dat.resMerged, list.cluster.datMerged)
}


