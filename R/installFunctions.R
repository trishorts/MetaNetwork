installFunctions <- function(){
 Rfiles <- list.files(path = "./R")
 RfilesPaths <- paste("R/", Rfiles)
 for(i in 1:length(RfilesPaths)){
   source(RfilesPaths)
 }
}