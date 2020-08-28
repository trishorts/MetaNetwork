## Variance explained by eigenproteins
writeVarianceExplained <- function(datExpr, 
         colors, 
         MEs){
  varianceExplained <- propVarExplained(datExpr, 
                   colors, 
                   MEs, 
                   corFnc = "cor", 
                   corOptions = "use = 'p'")
  write_csv(x = varianceExplained, path = "Results/varianceExplained.csv")
}
