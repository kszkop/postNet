gseaAnalysis <- function(ptn,
                         genesSlopeFiltOut = NULL,
                         collection = NULL,
                         subcollection = NULL,
                         subsetNames = NULL,
                         geneSet = NULL,
                         maxSize = 500,
                         minSize = 10,
                         name = NULL) {
  check_ptn(ptn)
  if (is.null(geneSet) && is.null(collection)) {
    stop("Please provide an input for either 'geneSet' or 'collection'.")
  }
  #
  if (!check_number(maxSize) | !check_number(minSize)) {
    stop("The inputs for 'maxSize' and 'minSize' must be integers.")
  }
  if (minSize <= 0 | maxSize <= 0) {
    stop("The inputs for 'maxSize' and 'minSize' must be positive.")
  }
  if (maxSize <= minSize) {
    stop("'maxSize' must be greater than 'minSize'.")
  }
  #
  #
  # if(!is.null(ads)){
  #  tmpAds <- anota2seq::anota2seqGetOutput(ads,
  #                                          analysis = regulationGen,
  #                                          output = "full",
  #                                          selContrast = contrastSel,
  #                                          getRVM = TRUE)
  #  #
  #  if (!is.null(genesSlopeFiltOut)) {
  #    tmpAdsFilt <- tmpAds[!row.names(tmpAds) %in% genesSlopeFiltOut, ]
  #  }  else {
  #    tmpAdsFilt <- tmpAds
  #  }
  #  #
  #  tmpP <- tmpAdsFilt[, "apvRvmP"]
  #  tmpEff <- tmpAdsFilt[, "apvEff"]
  # rankedRVMP <- rank(-log10(tmpP) * sign(tmpEff))
  #  rankIn <- tmpEff[order(tmpEff,decreasing = T)]
  # } else if (!is.null(rankIn)){
  #  rankIn <- rankIn[order(rankIn,decreasing = T)]
  # } else {
  #  stop("No anota2seq object or ranks provided")
  # }
  #
  effTmp <- ptn_effect(ptn)
  if (!is.null(genesSlopeFiltOut)) {
    effIn <- effTmp[!names(effTmp) %in% genesSlopeFiltOut]
  } else {
    effIn <- effTmp
  }
  #
  rankIn <- effIn[order(effIn, decreasing = TRUE)]
  if (is.null(geneSet)) {
    species <- ptn_species(ptn)
    if (!species %in% c("human", "mouse")) {
      stop("This option is currently only available for human or mouse.")
    }
    checkCollection(collection)
    eh <- ExperimentHub::ExperimentHub()
    AnnotationHub::query(eh, "msigdb")

    versionTmp <- as.character(sort(as.numeric(msigdb::getMsigdbVersions()), decreasing = TRUE))[1]
    msigdbOut <- msigdb::getMsigdb(org = ifelse(species == "human", "hs", "mm"), id = "SYM", version = versionTmp)
    msigdbOut <- msigdb::appendKEGG(msigdbOut, version = versionTmp)
    #

    # Start with hallmark
    collectionTmp <- msigdb::subsetCollection(msigdbOut, collection = collection, subcollection = subcollection)
    if (!is.null(subsetNames)) {
      collectionTmp <- collectionTmp[names(collectionTmp) %in% subsetNames]
    }
    geneSet_ids <- GSEABase::geneIds(collectionTmp)
  } else {
    check_geneList(geneSet)
    geneSet_ids <- geneSet
  }
  resOut <- fgsea::fgsea(pathways = geneSet_ids, stat = rankIn, minSize = minSize, maxSize = maxSize)

  # format output
  resOut$Count <- unlist(lapply(resOut$leadingEdge, length))
  colnames(resOut) <- c("Term", "pvalue", "adjusted_pvalue", "log2err", "ES", "NES", "Size", "Genes", "Count")
  resOut <- resOut[, c(1, 5, 6, 4, 9, 7, 2, 3, 8)]
  gseaOut <- resOut[order(resOut$adjusted_pvalue), ]
  gseaOut$Genes <- sapply(gseaOut$Genes, function(x) paste(x, collapse = ":"))

  nameTmp <- ifelse(!is.null(name), paste(name, "gseaAnalysis", sep = "_"), "gseaAnalysis")
  data.table::fwrite(gseaOut, file = paste(nameTmp, ".txt", sep = ""), sep = "\t") # , sep2=c("", ":", ""))
  #
  ptn@analysis@GSEA <- gseaOut
  #
  return(ptn)
}
