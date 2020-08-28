plotEigenproteinClusteringPostMerge <- function(preMergeDendroData, 
                                                postMergeDendroData, 
                                                fileName = "Results/ModuleEigenproteinMergeDendrogram.pdf", 
                                                widthInches = 10){
  premergeDendro <- plot(preMergeDendroData, main = "Clustering of module eigenproteins, pre-merging", xlab = "", sub = "")
  postmergeDendro <- plot(postMergeDendroData, main = "Clustering of module eigenproteins, post-merging", xlab = "",
                          sub = "")
  
  pdf(file = fileName, width = widthInches)
  par(mfrow = c(2,1))
  par(cex = 0.6)
  premergeDendro
  abline(h = MCutHeight, col = "red")
  postmergeDendro
  dev.off()  
  
}

