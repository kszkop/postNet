signCalc <- function(addSign,
                     annot,
                     ads = NULL,
                     customBg = NULL,
                     geneList = NULL) {
  annotBg <- gSel(annot = annot, ads = ads, customBg = customBg, geneList = geneList)
  #
  outSign <- list()
  for (i in 1:length(addSign)) {
    Outvec <- rep(0, length(unique(annotBg$geneID)))
    names(Outvec) <- unique(annotBg$geneID)
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
