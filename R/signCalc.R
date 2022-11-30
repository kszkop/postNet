#Function calucate input for feature integration but without motif, just based on signature
signCalc  <- function(signature, #vector of genes
                      ads=NULL,
                      customBg=NULL,
                      geneList=NULL
){
  #Subset annot for only expressed genes
  if(!is.null(ads) | !is.null(customBg)){
    bg <- row.names(ads@dataP)
  } else {
    if(!is.null(customBg)){
      bg <- customBg
    } else {
      bg <- as.character(unlist(geneList))
    }
  }
  #
  Outvec <- rep(0,length(bg))
  names(Outvec) <- bg
  #
  Outvec[names(Outvec) %in% as.character(signature)] <- 1
  
  #
  return(Outvec)
}
