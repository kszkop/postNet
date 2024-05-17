signCalc <- function(a2sU,
                     addSign) {
  #
  outSign <- list()
  #
  Outvec <- rep(0, length(anota2seqUtilsGetCDSgeneID(a2sU)))
  names(Outvec) <- anota2seqUtilsGetCDSgeneID(a2sU)
  for (i in 1:length(addSign)) {
    #
    signature <- addSign[[i]]
    #
    Outvec[names(Outvec) %in% as.character(signature)] <- 1
    #
    outSign[[names(addSign)[i]]] <- Outvec
  }
  #
  return(outSign)
}
