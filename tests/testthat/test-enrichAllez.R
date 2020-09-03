## Test data
require(org.Hs.eg.db)
require(org.Mm.eg.db)
testUniverse <- read_csv(file.path("tests", "TestDataSet", "Unsupervised_data.csv"))[[1]]
LightCyanTestAccessions <- read_table(file.path("tests", "TestDataSet", "LightCyanTestAccessions.txt"), 
                                      col_names = FALSE)[[1]]


testAllezEnrichment <- enrichAllez(GeneSymbols = LightCyanTestAccessions, 
            GeneUniverse = testUniverse)

