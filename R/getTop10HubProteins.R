getTop10HubProteins <- function(df, file.name){
  ##get the module color:
  mod.color <- df$moduleColors[1]
  ## Get the number of unique proteins
  protein.count <- length(unique(df$accession))
  #get the column with the
  #get the GO terms from GOTerms, using the df.
  top10df <- df[1:10,]
  ggplot(data = top10df, mapping = aes_string(x = colnames(top10df)[1], y = colnames(top10df)[ncol(top10df)] ))+
    geom_col()+
    ggtitle(paste("Top KMEs for", mod.color, "module. ", sep = " "))+
    theme_classic()
  ggsave(filename = file.name, last_plot())
}

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