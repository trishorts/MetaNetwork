addPvaluesToAllezOutput <- function(outputAllez, Lowersetsize = 5, Uppersetsize = 500, side = "T"){
  if(side=="F")outputAllez$p.value <- pnorm(-abs(outputAllez$z.score))# two tailed
  if(side=="T"){
    prb <- pnorm(outputAllez$z.score)# one tailed
    outputAllez$p.value <- ifelse(1-prb>prb, prb, 1-prb)*2
  }
  outputAllez$p.adj <- p.adjust(outputAllez$p.value, method="BH")
  outputAllez <- outputAllez[which(outputAllez$set.size>Lowersetsize),]
  outputAllez <- outputAllez[which(outputAllez$set.size<Uppersetsize),]
  outputAllezOut <- outputAllez[order(outputAllez$p.value),
                                c("Term","p.value","p.adj","z.score",
                                  "set.size","set.mean","set.sd")]
  message("sets with size < ",Lowersetsize, " or > ", Uppersetsize, " are not considered" )
  LocalOut <- outputAl
  lezOut
  message("Successfully added p-values to Allez Output")
  return(outputAllezOut)
}