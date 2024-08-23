signCalc <- function(a2sU,
                     signatures) {
  #
  outSign <- list()
  #
  Outvec <- rep(0, length(a2sU_geneID(a2sU, 'CDS')))
  names(Outvec) <- a2sU_geneID(a2sU, 'CDS')
  for (i in 1:length(signatures)) {
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
