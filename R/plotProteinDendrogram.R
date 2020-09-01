plotProteinDendrogram <- function(proteinDendrogram,
  dynamicColors,
  mergedColors){
  DendroColorMergedClust <- plotDendroAndColors(proteinDendrogram,
                                                colors = cbind(dynamicColors, mergedColors),
                                                groupLabels = c("Dynamic Tree Cut", "Merged dynamic"),
                                                dendroLabels = FALSE, hang = 0.03,
                                                addGuide = TRUE, guideHang = 0.05)
  pdf("Results/DendroColorMergedClust.pdf")
  DendroColorMergedClust
  dev.off()
  message("Protein dendrogram generated")
}
