## test-mergeAndWriteWGCNAWorkbook.R 
source("./R/mergeAndWriteWGCNAWorkbook.R")

testDatabase <- read_tsv(file = "TestDataSet/filtered_Hsapiens.tab")

testData <- rnorm(n = 3000, mean = 16, sd = 3.2)
set.seed(9876)
testUniprotAccessions <- base::sample(testDatabase$Entry, 
                                      size = length(testData), 
                                      replace = FALSE)
testAllDataTibble  <- tibble(testUniprotAccessions, testData)
colnames(testAllDataTibble) <- c("accession", "expression")

testModulesMembership <- list(testAllDataTibble[1:750,], 
                              testAllDataTibble[751:1500,], 
                              testAllDataTibble[1501:2250,], 
                              testAllDataTibble[2251:3000,]) 
names(testModulesMembership) <- c("blue", "red", "green", "purple")

mergeAndWriteWGCNAWorkbook(selectedDatabase = testDatabase, 
                           allData = testAllDataTibble, 
                           dataList = testModulesMembership)

realTestData <- read.csv(file = "TestDataSet/Unsupervised_data.csv")
