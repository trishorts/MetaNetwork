# Server logic ----
# Helper functions required to make this work ----
# GetGOterms
getGOterms <- function(df, GOTerms, plot = TRUE, file.name){
  ##get the module color:
  mod.color <- df$moduleColors[1]
  ## Get the number of unique proteins
  protein.count <- length(unique(df$accession))
  #get the GO terms from GOTerms, using the df.
  cluster.uniprot <- df$accession
  index <- which(GOTerms$`Input Accession` %in% cluster.uniprot, useNames = TRUE, arr.ind = TRUE)
  
  GO.listModule <- GOTerms$`Input Accession`[index]
  GO.listID <- GOTerms$`GO:ID`[index]
  GO.listTerms <- GOTerms$`GO Term Name`[index]
  GO.listAspect <- GOTerms$Aspect[index]
  GO.df <- tibble(GO.listModule, GO.listID, GO.listTerms, GO.listAspect)
  ##Filtering the GO terms for P-type GO terms only
  GO.count <- GO.df%>%
    group_by(`GO.listID`)%>%
    filter(GO.listAspect == "F")%>% ## for both process and function: GO.listAspect == "F"|GO.listAspect == "P
    count(`GO.listTerms`)%>%
    arrange(-n)
  ##Plotting function
  if(plot == TRUE){
    p <- ggplot(data = GO.count[1:20,], mapping = aes(x = n, y = reorder(GO.listTerms, -n)))+
      geom_col()+
      ggtitle(paste("GO Terms for", mod.color, "module. ", protein.count, "proteins", sep = " "))+
      theme_classic()
    ggsave(filename = file.name, last_plot())
  }else{return(GO.count)}
}

# GetGOEnrichment
getGOEnrichment <- function(df, fullGOTermsList, plot = TRUE, file.name){
  
  #get the module color for plotting
  modcolor <- names(df)
  #
  fisher.vec <- vector()
  for(i in seq_along(df$`GO:ID`)){
    C_a <- df$`GO:ID`[[i]]
    index.val <- which(fullGOTermsList$`GO:ID` %in% C_a)
    A <- df$n[i]
    B <- sum(df$n)
    C <- fullGOTermsList$n[index.val] - A
    D <- sum(fullGOTermsList$n) - B
    fisher.tib <- tibble(rbind(A,C), rbind(B,D))
    inter1 <- fisher.test(fisher.tib)
    fisher.vec[i] <- inter1$p.value
  }
  
  
  GOterms.pvalue <- df %>%
    add_column(fisher.vec)%>%
    arrange(fisher.vec)
  top.20.pvalues <- GOterms.pvalue[1:20,]
  
  if(plot == TRUE){
    p <- ggplot(data = top.20.pvalues, mapping = aes(x = -log10(fisher.vec), y = reorder(GOTermName, -fisher.vec)))+
      geom_col()+
      geom_vline(xintercept = -log10(0.05/length(df$`GO:ID`)), col = "red")+
      scale_y_discrete()+
      ggtitle(file.name)+
      ylab("Top 20 GO Terms")+
      xlab("-log10(p-value)")+
      theme_bw()
    ggsave(file.name, plot = last_plot())
  }
}

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

fetchMart <- function(species){
  ensembl <- useMart("ensembl")
  if(species == "Human"){
  ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  } else if(species == "Mouse"){
    ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl)
    mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  }
}

read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

enrichAllez <- function(GeneSymbols, GeneUniverse, SpeciesLibrary = "org.Hs.eg", idtype = "SYMBOL", 
                        alter = TRUE, Lowersetsize = 5, Uppsersetsize = 500, outprefix, ...){  
  Scores <- rep(0, length(GeneUniverse))
  names(Scores) <- as.vector(GeneUniverse)
  Scores[intersect(names(Scores), GeneSymbols)] <- 1
  allezOutput <- allez(scores = Scores, 
                       lib = SpeciesLibrary, 
                       idtype = "SYMBOL", library.loc = GeneSymbols)
}

addPvaluesToAllezOutput <- function(outputAllez, Lowersetsize = 5, Uppersetsize = 500, side = "T"){
  if(side=="F")outputAllez$p.value <- pnorm(-abs(outputAllez$z.score))# two tailed
  if(side=="T"){
    prb <- pnorm(outputAllez$z.score)# one tailed
    outputAllez$p.value <- ifelse(1-prb>prb, prb, 1-prb)*2
  }
  outputAllez$p.adj <- p.adjust(outputAllez$p.value, method="BH")
  outputAllez <- outputAllez[which(outputAllez$set.size>Lowersetsize),]
  outputAllez <- outputAllez[which(outputAllez$set.size<Uppersetsize),]
  outputAllezOut <- outputAllez[order(outputAllez$p.value),c("Term","p.value","p.adj","z.score","set.size","set.mean","set.sd")]
  message("sets with size < ",Lowersetsize, " or > ", Uppersetsize, " are not considered" )
  return(outputAllezOut)
}

createAllezEnrichmentXLSX <- function(geneUniverse, AllezPvalues){
  wb = createWorkbook()
  
  addWorksheet(wb, paste("GeneUniverse"))
  writeData(wb, sheet = 1, geneUniverse)
  
  # write modules only tabs
  for(ii in 1:length(AllezPvalues)) {
    addWorksheet(wb, names(AllezPvalues)[ii])
    writeData(wb, names(AllezPvalues)[ii], AllezPvalues[[ii]])
  }
  saveWorkbook(wb, file = "Results/AllezEnrichment.xlsx", overwrite=TRUE)
}

server <- shinyServer(function(input, output) {
  # Generate the file preview ----
  output$preview <- renderTable({
    req(input$dataFile)
    req(input$groupsFile)
    req(input$dataFile)
    
    datafile.tibb <- read.csv(input$dataFile$datapath)
    groups.tibb <- read.csv(input$groupsFile$datapath)
    head(datafile.tibb)
  })
  # Action button reactive statement to activate the WGCNA workflow ----
  
  observeEvent(input$action, {
    # Action button for tab 1 ----
    #require a file input in order to do the workflow. ----
    req(input$dataFile)
    req(input$groupsFile)
    req(input$databaseFile)
    # User input defined/default setting for the WGCNA workflow ----
    RCutoff <- as.numeric(input$rcutoff)
    MCutHeight <- as.numeric(input$mcutheight)
    PowerUpper <- as.numeric(input$powerupper)
    minModuleSize <- as.numeric(input$minmodsize)
    IgnoreCols <- as.numeric(input$ignorecols)
    
    # Load the data files ----
    # Read the files in with read.csv, not read_csv
    if(input$ttestbool == TRUE){
      pre.allDataFile <- read.csv(input$dataFile$datapath)
      groupsFile <- read.csv(file = input$groupsFile$datapath)
      ttest_scores <- tibble(pre.allDataFile$accession, pre.allDataFile$TTEST)
      index_ttest <- which(colnames(pre.allDataFile) == "TTEST")
      allDataFile <- pre.allDataFile[-15]
    }else{
      allDataFile <- read.csv(file = input$dataFile$datapath)
      groupsFile <- read.csv(file = input$groupsFile$datapath)
    }
    # Transpose data so that genes are columns and samples are rows ----
    allDataFile_t <- as.data.frame(t(allDataFile[,-c(1:IgnoreCols)]) )
    # column names and row names are now switched.
    rownames(allDataFile_t) <- colnames(allDataFile)[-c(1:IgnoreCols)]
    colnames(allDataFile_t) <- allDataFile[,1]
    
    #WGCNA Workflow starts here ----
    allowWGCNAThreads()
    # Choose a set of soft-thresholding powers
    powers <- c(c(1:10), seq(from = 12, to=PowerUpper, by=2))
    
    # Call the network topology analysis function
    sft <- pickSoftThreshold(allDataFile_t, powerVector = powers, RsquaredCut = RCutoff, verbose = 5)
    # Output all files, so sft needs to output to the working directory. I will need to figure
    # out how to output all the files into a temporary directory with a bunch of subfolders
    #Dendrograms ----
    softPower <- 12
    # Build the adjacency table - use "signed" for proteomics data
    adjacency <- adjacency(allDataFile_t, power = softPower, type="signed")
    
    
    # Turn adjacency into topological overlap distance
    TOM <- TOMsimilarity(adjacency)
    dissTOM <- 1-TOM
    
    path <- getwd()
    dir.create(file.path(path, "Results"))
    # Clustering using TOM-based dissimilarity
    proTree <- hclust(as.dist(dissTOM), method = "average");
    
    
    # Module identification using dynamic tree cut
    dynamicMods <- cutreeDynamic(dendro = proTree, distM = dissTOM,
                                 deepSplit = 2, pamRespectsDendro = FALSE,
                                 minClusterSize = minModuleSize)
    dynamicColors <- labels2colors(dynamicMods)
    # Merge clusters
    ## The old way was that the modules were merged without first calculating the module eigenproteins, which I believe
    ## was incorrect
    # Old code remnant: mergedClust <- mergeCloseModules(allDataFile_t, dynamicColors, cutHeight = MCutHeight, verbose = 3)
    # New code starts here: 
    # Calculate the module eigenproteins
    #The threshold for the merge. Default is 0.25, corresponding to a correlation of 0.75
    allDataModuleEigenProteins <- moduleEigengenes(expr = allDataFile_t, colors = dynamicColors)
    preMergeModuleEigenProteins <- allDataModuleEigenProteins$eigengenes
    MEDiss <- 1 - cor(preMergeModuleEigenProteins)
    METree <- hclust(as.dist(MEDiss), method = "average")
    mergedClust <- mergeCloseModules(allDataFile_t, dynamicColors, cutHeight = MCutHeight, verbose = 3) 
    
    mergedColors <- mergedClust$colors
    mergedMEs <- mergedClust$newMEs
    moduleColors <- mergedColors
    
    MEs <- mergedMEs
    
    rownames(MEs) <- rownames(allDataFile_t)
    
    
    # reorder MEs by color names of modules
    MEs <- MEs[,order(colnames(MEs))]
    #Create the EigenGeneNetworks Plot
    DendrogramEigenproteins <- plotEigengeneNetworks(MEs, "Eigenprotein Network", marHeatmap = c(3,4,2,2), marDendro = c(3,4,2,5),
                                                     plotDendrograms = TRUE, xLabelsAngle = 90,heatmapColors=blueWhiteRed(50))
    # Plot module profiles with eigenproteins overlaid
    WGCNAClusterID <- moduleColors
    
    ## aggregate groups
    gp <- groupsFile
    noClusters <- nlevels(as.factor(moduleColors))
    
    ##tidy data format for ease of plotting
    
    col.length <- length(colnames(allDataFile_t)) - 1
    
    allDataFile.1 <- add_column(allDataFile_t, "Experiment" = rownames(allDataFile_t))
    
    dataFile_tidy <- gather(allDataFile.1, key = "Gene", value = "Expression", 1:col.length) ##Because of the add_column function one above this, the last column needs to be n-1
    
    genes_colors.df <- data.frame(colnames(allDataFile_t), moduleColors)
    colnames(genes_colors.df) <- c("Gene", "moduleColor")
    dataColTidy <- left_join(dataFile_tidy, genes_colors.df, "Gene")
    
    #convert to factor
    dataColTidy$moduleColor <- as.factor(dataColTidy$moduleColor)
    # maybe this epxression will be needed. We'll see.
    
    #make the groupsFile input match the format of dataColTidy so you can left_join them
    colnames(groupsFile) <- c("Experiment", "SampleID")
    groupsFile$Experiment <- str_replace_all(groupsFile$Experiment, pattern = " ", replacement = ".")
    
    allDataFinal <- left_join(dataColTidy, groupsFile, "Experiment")
    
    # Prepare the data for plotting ----
    # Get the eigenprotein means ----
    EigenproteinModuleMeans <- allDataFinal%>%
      group_by(SampleID, moduleColor)%>%
      summarize("Mean" = mean(Expression), "sd" = sd(Expression))

    #Create the TOM Plot
    #Create a new variable called TOMTransformationPower
    plotTOM <- dissTOM^7
    # Set diagonal to NA for a nicer plot
    diag(plotTOM) = NA
    # Get KME - module membership - correlation between proteins and eigenproteins ----
    
    kmes <- signedKME(allDataFile_t, MEs)
    
    
    # separate results by modules, order by kME, hub proteins on top, in preparation for
    
    dat.res <- data.frame(allDataFile, moduleColors , kmes)
    
    list.cluster.dat <- lapply(levels(as.factor(moduleColors)),
                               function(x) {dtemp = dat.res[dat.res$moduleColors == x,];
                               dtemp[order(dtemp[,paste0('kME',x)==colnames(dtemp)], decreasing=TRUE),
                                     -setdiff(grep("^kME", colnames(dtemp)), which(paste0('kME',x)==colnames(dtemp)))]} )
    
    names(list.cluster.dat) <- 	levels(as.factor(moduleColors))
    
    #Need to implement choices for uniprot database
   # userInputDatabase <- read_tsv(input$databaseFile$datapath)
    #userInputDatabaseSelectedColumns <- tibble(userInputDatabase$Entry, userInputDatabase$`Gene names`, userInputDatabase$`Protein names`)
    #colnames(userInputDatabaseSelectedColumns) <- c("accession", "gene name", "protein name")
# 
 #   dat.resMerged <- left_join(tibble(dat.res), userInputDatabaseSelectedColumns, by = "accession")
    dat.resMerged <- tibble(dat.res)
    #  list.cluster.datMerged <- list()
   # for(i in seq_along(list.cluster.dat)){
    #  list.cluster.datMerged[[i]] <- left_join(list.cluster.dat[[i]], userInputDatabaseSelectedColumns, by = "accession")
    #}
    list.cluster.datMerged <- list.cluster.dat
    names(list.cluster.datMerged) <- names(list.cluster.dat)
    
    
    #Calculate the module eigenproteins
    allDataModuleEigenProteins <- moduleEigengenes(expr =  allDataFile_t, colors = moduleColors)
    #Quality control on the sample clustering
    sampleTree <- hclust(dist(allDataFile_t), method = "average")
    
    #Get the Eigenproteins for each module ready to print
    EigenProteins <- tibble(rownames(allDataModuleEigenProteins$eigengenes), allDataModuleEigenProteins$eigengenes)
    EigenProteinsColumnNumber <- length(colnames(EigenProteins))

    EigenProteins_tidy <- EigenProteins%>%
      gather(key = "Samples", value = "ModuleEigenproteins", 2:EigenProteinsColumnNumber)
    colnames(EigenProteins_tidy) <- c("Experiment", "Module", "ModuleEigenprotein")
    EigenProteinsFin <- left_join(EigenProteins_tidy, tibble(groupsFile), "Experiment")
    
    ## Output the files
    path <- getwd()
    dir.create(file.path(path, "Results"))
    #Write the files
    # module memberships
    createResultsWGCNAExcelWorkbook(dat.resMerged, list.cluster.datMerged) #Need to update this with a filepath argument
    #Eigenproteins
    write_csv(EigenProteinsFin, path = "Results/Eigenproteins_by_sample.csv")
    
    #Create the .pdfs
    #Sample clustering quality control plot
    ggdendrogram(sampleTree)+
      ggtitle("Sample Clustering to Detect Outliers")+
      xlab("Sample")+
      ylab("Height")+
      theme_bw(base_family = "Arial", base_size = 15)
    
    # Scale free topology plot ----
    pdf("Results/ScaleFreeTopology.pdf", width = 10)
    par(mfrow = c(1,2));
    cex1 = 0.9;
    # Scale-free topology fit index as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
         main = paste("Scale independence"));
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=powers,cex=cex1,col="red");
    # this line corresponds to using an R^2 cut-off of h
    h = RCutoff
    abline(h=RCutoff,col="red")
    # Mean connectivity as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
         main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
    dev.off()

    # Plot the pre-merge and the merged module eigenprotein clustering
    pdf("Results/Module Eigenprotein PrePost MergeDendrogram.pdf", width = 10)
    par(mfrow = c(2,1))
    par(cex = 0.6)
    plot(METree, main = "Clustering of module eigenproteins, pre-merging", xlab = "", sub = "")
    abline(h = MCutHeight, col = "red")
    plot(mergedClust$dendro, main = "Clustering of module eigenproteins, post-merging", xlab = "",
         sub = "")
    dev.off()


    # Plot the dendrogram following cluster merging ----
    pdf("Results/DendroColorMergedClust.pdf")
    plotDendroAndColors(proTree, cbind(dynamicColors, mergedColors),
                        c("Dynamic Tree Cut", "Merged dynamic"),
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
    dev.off()

    #Eigenprotein Violin Plots
    EigenProteinPlot <- ggplot(data = EigenproteinModuleMeans,
                               mapping = aes(x = SampleID, y = Mean))+
      geom_violin(data = allDataFinal, mapping = aes(x = SampleID, y = Expression))+
      geom_point()+
      geom_line(data = EigenproteinModuleMeans, mapping = aes(x = SampleID, y = Mean, group = 1))+
      ylab("Protein Expression")+
      facet_wrap(~moduleColor)+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90))
    ggsave(plot = EigenProteinPlot, filename = "eigenprotein_cluster_profiles.pdf", path = "Results")

    #Eigeneproteins dendrogram
    pdf("Results/Dendrogram_eigenproteins.pdf", width = 10)
    plotEigengeneNetworks(MEs, "EigenproteinNetwork", marHeatmap = c(3,4,2,2), marDendro = c(3,4,2,5),
                          plotDendrograms = TRUE, xLabelsAngle = 90,heatmapColors=blueWhiteRed(50))
    dev.off()

    #Heatmap of the Module Eigenproteins
    pdf(file = "Results/Module Eigenproteins.pdf", width = 10)
    heatmap3(MEs, distfun = function(x) dist(x, method="euclidean"),
             main = "Module Eigenproteins",
             cexRow = 0.6, cexCol = 0.6)
    dev.off()

    #pdf("Results/Network_heatmap.pdf", width = 10)
    TOMplot <- TOMplot(plotTOM, proTree, moduleColors, main = "Network heatmap plot, all proteins")
    #as.pdf(TOMplot)
    #dev.off()

    #need to add the adjacency heatmap
    MET <- orderMEs(MEs = MEs)
    pdf(file = "Results/Eigenprotein_adjacency heatmap.pdf", width = 10)
    par(cex = 1.0)
    plotEigengeneNetworks(MET, "Eigenprotein Dendrogram",
                          marDendro = c(0,4,2,0), plotHeatmaps = FALSE)

    plotEigengeneNetworks(MET, "Eigenprotein adjacency heatmap",
                          marDendro = c(3,4,2,2), xLabelsAngle = 90)
    dev.off()
    
    #Output code
    output$workflowOutput <- renderText({
      print("WGCNA workflow completed. Results exported to folder.")
    })
  }) # this one is the close of the observeEvent function that wraps the WGCNA analysis workflow
  
  observeEvent(input$action2, {
    #require GO File to be uploaded----
    req(input$WGCNAResults)
    
    ## Load the WGCNA Results File
    wgcnaResults <- read_excel_allsheets(input$WGCNAResults$datapath)
    sheetNumber <- length(wgcnaResults)
    ## Load biomaRt (update this to make it faster -AVC 08/04/20)
    BiocManager::install("biomaRt", update = FALSE)
    library(biomaRt)
    ## Get the correct database
    if(input$organismID == "Human"){
      mart <- fetchMart("Human")
      BiocManager::install('org.Hs.eg.db')
      library(org.Hs.eg.db)
    } else if(input$organismID == "Mouse"){
      mart <- fetchMart("Mouse")
      BiocManager::install('org.Mm.eg.db')
      library(org.Mm.eg.db)
    }
    
    ## The next three things will probably end up as a single function
    ## Convert the Uniprot Acessions to Gene Symbols for each sample 
    sheetNumber <- length(wgcnaResults)
    UniprotAcessions <- list()
    for(i in seq_len(sheetNumber)){
      UniprotAcessions[[i]] <- wgcnaResults[[i]][,1]
    }
    ModuleNames <- names(wgcnaResults)
    ModuleNames2 <- names(wgcnaResults)[2:length(wgcnaResults)]
    names(UniprotAcessions) <- ModuleNames
    
    ConvertedGeneSymbols <- list()
    for(i in seq_len(sheetNumber)){
      ConvertedGeneSymbols[[i]] <- getBM(attributes = "external_gene_name", 
                                         filters = "uniprot_gn_id", mart = mart, 
                                         values = UniprotAcessions[[i]])$external_gene_name
    }
    
    geneUniverse <- unique(ConvertedGeneSymbols[[1]])
    ConvertedGeneSymbolsWithoutUniverse <- ConvertedGeneSymbols[2:length(ConvertedGeneSymbols)]
    names(ConvertedGeneSymbolsWithoutUniverse) <- ModuleNames[2:length(ConvertedGeneSymbols)]
    
    AllezEnriched <- list()
    #if(input$organismID == "Human"){
    for(i in 1:length(ConvertedGeneSymbolsWithoutUniverse)){
      AllezEnriched[[i]] <- enrichAllez(ConvertedGeneSymbolsWithoutUniverse[[i]], 
                                   GeneUniverse = geneUniverse)
    }
    #} else if(input$organismID == "Mouse"){
    for(i in 1:length(ConvertedGeneSymbolsWithoutUniverse)){
      AllezEnriched[[i]] <- enrichAllez(ConvertedGeneSymbolsWithoutUniverse[[i]], 
                                        GeneUniverse = geneUniverse, SpeciesLibrary = "org.Mm.eg")
    }
    #}
    
    AllezPvalues <- list()
    for(i in 1:length(AllezEnriched)){
      AllezPvalues[[i]] <- addPvaluesToAllezOutput(outputAllez = AllezEnriched[[i]][[1]])
    }
    names(AllezPvalues) <- ModuleNames2
    
    
    ## Create the workbook with the Allez sheets 
    createAllezEnrichmentXLSX(geneUniverse, AllezPvalues)
    message("Allez Enrichment workflow completed")
  })
}) #this one is the server's end
