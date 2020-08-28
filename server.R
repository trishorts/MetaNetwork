## Need to load all the R functions required to make this thing run
source(file = "./R/installFunctions.R")
installFunctions()

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
    allDataFile <- read.csv(file = input$dataFile$datapath)
    groupsFile <- read.csv(file = input$groupsFile$datapath)

    # Transpose data so that genes are columns and samples are rows ----
    allDataFile_t <- as.data.frame(t(allDataFile[,-c(1:IgnoreCols)]) )
    # column names and row names are now switched.
    rownames(allDataFile_t) <- colnames(allDataFile)[-c(1:IgnoreCols)]
    colnames(allDataFile_t) <- allDataFile[,1]
    
    #WGCNA Workflow starts here ----
    allowWGCNAThreads()
    # Choose a set of soft-thresholding powers
    powers <- c(c(1:10), seq(from = 12, to=PowerUpper, by=2))
    sft <- pickSoftThreshold(allDataFile_t, 
                               powerVector = powers, 
                               RsquaredCut = RCutoff, 
                               verbose = 5)
    softPower <- 12 ## Need to implement 
    # Build the adjacency table - use "signed" for proteomics data
    # create control flow for signed versus unsigned data
    adjacency <- adjacency(allDataFile_t, power = softPower, type="signed")
    TOM <- TOMsimilarity(adjacency)
    dissTOM <- 1-TOM
    # Clustering using TOM-based dissimilarity
    proTree <- hclust(as.dist(dissTOM), method = "average");
    
    path <- getwd()
    dir.create(file.path(path, "Results"))
    
    
    
    # Module identification using dynamic tree cut
    dynamicMods <- cutreeDynamic(dendro = proTree, distM = dissTOM,
                                 deepSplit = 2, pamRespectsDendro = FALSE,
                                 minClusterSize = minModuleSize)
    dynamicColors <- labels2colors(dynamicMods)
    # Merge clusters
    # Calculate the module eigenproteins
    #The threshold for the merge. Default is 0.25, corresponding to a correlation of 0.75
     
    mergedClust <- mergeModuleEigenproteins(allDataFile_t, moduleColors = dynamicColors, 
                             cutHeight = MCutHeight)
    mergedColors <- mergedClust$colors
    mergedMEs <- mergedClust$newMEs
    moduleColors <- mergedColors
    MEs <- mergedMEs
    rownames(MEs) <- rownames(allDataFile_t)
    
    # reorder MEs by color names of modules
    MEs <- MEs[,order(colnames(MEs))]
    #Create the EigenGeneNetworks Plot
    DendrogramEigenproteins <- plotEigengeneNetworks(MEs, "Eigenprotein Network", 
                                                     marHeatmap = c(3,4,2,2), 
                                                     marDendro = c(3,4,2,5),
                                                     plotDendrograms = TRUE, 
                                                     xLabelsAngle = 90,
                                                     heatmapColors=blueWhiteRed(50))
    # Plot module profiles with eigenproteins overlaid
    WGCNAClusterID <- moduleColors
    
    ##tidy data format for ease of plotting
    
    tidyDataFinal <- tidyModuleDataForOutput(dataInput = allDataFile_t, moduleColors = moduleColors, 
                            groupsFile = groupsFile)
    
    #Create the TOM Plot
    plotTOM <- dissTOM^7
    diag(plotTOM) = NA # Makes a nicer plot
    # Get KME - module membership - correlation between proteins and eigenproteins ----
    kmes <- signedKME(allDataFile_t, MEs)
    
    # separate results by modules, order by kME, hub proteins on top, in preparation for
    dat.res <- data.frame(allDataFile, moduleColors , kmes)
    
    list.cluster.dat <- lapply(levels(as.factor(moduleColors)),
                               function(x) {dtemp = dat.res[dat.res$moduleColors == x,];
                               dtemp[order(dtemp[,paste0('kME',x)==colnames(dtemp)], decreasing=TRUE),
                                     -setdiff(grep("^kME", colnames(dtemp)), which(paste0('kME',x)==colnames(dtemp)))]} )
    
    names(list.cluster.dat) <- 	levels(as.factor(moduleColors))
    
    #Calculate the module eigenproteins
    allDataModuleEigenProteins <- moduleEigengenes(expr =  allDataFile_t, colors = moduleColors)
    #Quality control on the sample clustering
    sampleTree <- hclust(dist(allDataFile_t), method = "average")
    
    #Get the Eigenproteins for each module ready to print
    EigenProteins <- tibble(rownames(allDataModuleEigenProteins$eigengenes), 
                            allDataModuleEigenProteins$eigengenes)
    EigenProteinsColumnNumber <- length(colnames(EigenProteins))

      EigenProteins_tidy <- EigenProteins%>%
        gather(key = "Samples", value = "ModuleEigenproteins", 2:EigenProteinsColumnNumber)
      colnames(EigenProteins_tidy) <- c("Experiment", "Module", "ModuleEigenprotein")
      EigenProteinsFin <- left_join(EigenProteins_tidy, tibble(groupsFile), "Experiment")
    
      ## Output the files
      path <- getwd()
      dir.create(file.path(path, "Results"))
      # Write the files
      # module memberships
      userInputDatabase <- read_tsv(input$databaseFile$datapath)
      mergeAndWriteWGCNAWorkbook(selectedDatabase = userInputDatabase, 
                                 allData = dat.res, 
                                 dataList = list.cluster.dat)
      
      #Eigenproteins
      write_csv(EigenProteinsFin, path = "Results/Eigenproteins_by_sample.csv")
      EigenproteinsWide <- spread(data = EigenProteinsFin, 
                                  key = `Module`, 
                                  value = ModuleEigenprotein)
      write_csv(EigenProteins, path = "Results/EigenproteinsWide.csv")
      # Proportion of variance explained
      writeVarianceExplained(datExpr = allDataFile_t, 
                             colors = mergedColors, 
                             MEs = mergedMEs)
    
     # Create the .pdfs
     # Sample clustering quality control plot
      plotSampleClusteringDendro(sampleTree)
    
      # Scale free topology plot ----
      plotScaleFreeTopology(scaleFreeThreshold = Rcutoff)
     
      # Plot the pre-merge and the merged module eigenprotein clustering
      plotEigenproteinClusteringPostMerge(preMergeDendroData = METree, 
                                        postMergeDendroData = mergedClust$dendro) 
    
      # Plot the dendrogram following cluster merging ----
      plotProteinDendrogram(proTree)
     
      #Eigeneproteins dendrogram
      plotEigenproteinsNetwork(moduleEigenproteins = MEs)
     
      #Heatmap of the Module Eigenproteins
      plotmoduleEigenproteinsHeatmap(moduleEigenproteins = MEs)
    
      # Plot network heatmap 
      plotTOM(TOMData = plotTOM, 
              proteinDendro = proTree, 
              moduleColors = moduleColors)
      #need to add the adjacency heatmap
      plotAdjacencyHeatmap(moduleEigenproteins = MEs)
    
      #Output code
      output$workflowOutput <- renderText({
      print("WGCNA workflow completed. Results exported to folder.")
      #Grid arrange the plots in preparation for output
      output$PLOTOUTPUT <-  renderPlot({
      grid.arrange(tomplotplot, sampleClusteringQC, 
      DendrogramEigenproteins, DendroColorMergedClust, ncol = 2)
      })
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
      speciesLibrary <- 'org.Hs.eg'
      library(org.Hs.eg.db)
    } else if(input$organismID == "Mouse"){
      mart <- fetchMart("Mouse")
      BiocManager::install('org.Mm.eg.db')
      library(org.Mm.eg.db)
      speciesLibrary <- 'org.Mm.eg'
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
    for(i in 1:length(ConvertedGeneSymbolsWithoutUniverse)){
      AllezEnriched[[i]] <- enrichAllez(ConvertedGeneSymbolsWithoutUniverse[[i]], 
                                   GeneUniverse = geneUniverse, 
                                   SpeciesLibrary = speciesLibrary) #species library 
      # is defined in the control flow in lines 462-473
    }
    
    
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
