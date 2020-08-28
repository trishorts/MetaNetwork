enrichAllez <- function(GeneSymbols, GeneUniverse, SpeciesLibrary = "org.Hs.eg", idtype = "SYMBOL", 
                        alter = TRUE, Lowersetsize = 5, Uppsersetsize = 500, outprefix, ...){  
  Scores <- rep(0, length(GeneUniverse))
  names(Scores) <- as.vector(GeneUniverse)
  Scores[base::intersect(names(Scores), GeneSymbols)] <- 1
  allezOutput <- allez(scores = Scores, 
                       lib = SpeciesLibrary, 
                       idtype = "SYMBOL", library.loc = GeneSymbols)
}