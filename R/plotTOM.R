plotTOM <- function(TOMData, 
                   proteinDendro, 
                   moduleColors, 
                   fileName = "Results/NetworkHeatmap.png"){
  tomplot <- TOMplot(TOMData, proteinDendro, moduleColors)
  
  dev.new()
  png(filename = fileName)
  tomplot
  dev.off()
}
