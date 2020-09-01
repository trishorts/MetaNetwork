createResultsWGCNAExcelWorkbook <- function(fullDataSet, modulesList){
  wb = createWorkbook()
  
  addWorksheet(wb, paste("WGCNA_workbook"))
  writeData(wb, sheet = 1, fullDataSet)
  
  # write modules only tabs
  for(ii in 1:length(modulesList)) {
    addWorksheet(wb, names(modulesList)[ii])
    writeData(wb, names(modulesList)[ii], modulesList[[ii]])
  }
  
  saveWorkbook(wb, file = paste("Results/", "ResultsWGCNA", ".xlsx", sep = ""), overwrite=TRUE)
}