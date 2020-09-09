plotProteinDendrogram <- function(proteinDendrogram,
  dynamicColors,
  mergedColors){

  pdf("Results/DendroColorMergedClust.pdf")
  plotDendroAndColors(proteinDendrogram,
                      colors = cbind(dynamicColors, mergedColors),
                      groupLabels = c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  message("Protein dendrogram generated")
}
