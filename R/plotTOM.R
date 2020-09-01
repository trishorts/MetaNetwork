plotTOM <- function(TOMData, 
                   proteinDendro, 
                   moduleColors, 
                   fileName = "Results/NetworkHeatmap.png"){
  message("Begin generating TOM plot.")
  message("Please be patient. This may take a while")
  tomplot <- TOMplot(TOMData, proteinDendro, moduleColors)
  
  dev.new()
  png(filename = fileName)
  tomplot
  dev.off()
  message("TOM plot successfully generated")
}
