signCalc <- function(ptn,
                     signatures) {
  #
  outSign <- list()
  #
  for (i in 1:length(signatures)) {
    #
    Outvec <- rep(0, length(ptn_geneID(ptn, 'CDS')))
    names(Outvec) <- ptn_geneID(ptn, 'CDS')
    #
    tmpSignature <- signatures[[i]]
    #
    Outvec[names(Outvec) %in% as.character(tmpSignature)] <- 1
    #
    outSign[[names(signatures)[i]]] <- Outvec
  }
  #
  return(outSign)
}
