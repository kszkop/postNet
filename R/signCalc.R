signCalc <- function(ptn, signatures) {
  #
  outSign <- list()
  #
  for (i in seq_along(signatures)) {
    #
    Outvec <- rep(0, length(ptn_geneID(ptn, "CDS")))
    names(Outvec) <- ptn_geneID(ptn, "CDS")
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
