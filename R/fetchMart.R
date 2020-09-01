fetchMart <- function(species){
  ensembl <- useMart("ensembl")
  if(species == "Human"){
    ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
    mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  } else if(species == "Mouse"){
    ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl)
    mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  }
  message("Successfully retrieved database")
}