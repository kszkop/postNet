adjustSeq <- function(annot,
                      adjObj,
                      region_adj,
                      excl = FALSE,
                      keepAll = FALSE) {
  #
  colnames(adjObj) <- c("id", "sequence")

  annotTmp <- annot
  adjObj <- adjObj[adjObj$id %in% annotTmp$id, ]
  #
  if (isTRUE(excl)) {
    annotTmp <- annotTmp[annotTmp$id %in% adjObj$id, ]
  }
  #
  if (region_adj == "UTR5") {
    annotTmp[match(adjObj$id, annotTmp$id), "UTR5_seq"] <- adjObj$sequence
    if (!isTRUE(keepAll)) {
      #
      annotTmp <- annotTmp[(annotTmp$geneID %in% unique(annotTmp[annotTmp$id %in% adjObj$id, ]$geneID) & annotTmp$id %in% unique(annotTmp[annotTmp$id %in% adjObj$id, ]$id)) | (!annotTmp$geneID %in% unique(annotTmp[annotTmp$id %in% adjObj$id, ]$geneID)), ]
    }
  } else if (region_adj == "UTR3") {
    annotTmp[annotTmp$id %in% adjObj$id, "UTR3_seq"] <- adjObj$sequence
    if (!isTRUE(keepAll)) {
      annotTmp <- annotTmp[(annotTmp$geneID %in% unique(annotTmp[annotTmp$id %in% adjObj$id, ]$geneID) & annotTmp$id %in% unique(annotTmp[annotTmp$id %in% adjObj$id, ]$id)) | (!annotTmp$geneID %in% unique(annotTmp[annotTmp$id %in% adjObj$id, ]$geneID)), ]
    }
  } else {
    stop("Please provide correct region")
  }
  return(annotTmp)
}
