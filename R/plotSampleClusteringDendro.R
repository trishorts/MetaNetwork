plotSampleClusteringDendro <- function(dendro, fontSize = 15, 
                                       fileName = "SampleClustering.pdf"){
sampleClusteringQC  <- ggdendrogram(dendro)+
  ggtitle("Sample Clustering to Detect Outliers")+
  xlab("Sample")+
  ylab("Height")+
  theme_bw(base_size = fontSize)
ggsave(filename = fileName, path = "Results", plot = last_plot())
}
