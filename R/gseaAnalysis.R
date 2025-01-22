gseaAnalysis <- function(ptn,
                         genesSlopeFiltOut = NULL,
                         collection = NULL,
                         subcollection = NULL,
                         subsetNames = NULL,
                         geneSet = NULL,
                         maxSize = 500,
                         minSize = 10,
                         name = NULL){
  if (!checkPtn(ptn)) {
    stop("ptn is not a valid 'postNetData' object.")
  }
  if(is.null(geneSet) && is.null(collection)){
    stop("please provide geneSet or collection")
  }
  #
  #if(!is.null(ads)){
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
    #rankedRVMP <- rank(-log10(tmpP) * sign(tmpEff))
  #  rankIn <- tmpEff[order(tmpEff,decreasing = T)]
  #} else if (!is.null(rankIn)){
  #  rankIn <- rankIn[order(rankIn,decreasing = T)]
  #} else {
  #  stop("No anota2seq object or ranks provided")
  #}
  #
  effTmp <- ptn_eff(ptn)
  if (!is.null(genesSlopeFiltOut)) {
    effIn <- effTmp[!names(effTmp) %in% genesSlopeFiltOut ]
  }  else {
    effIn <- effTmp
  }
  #
  rankIn <- effIn[order(effIn,decreasing = T)]
  if(is.null(geneSet)){
    species <- ptn_species(ptn)
    if (!species %in% c("human","mouse")) {
      stop("This option is only  available for species: human and mouse at the moment")
    }
    checkCollection(collection)
    eh <- ExperimentHub::ExperimentHub()
    AnnotationHub::query(eh , 'msigdb')
  
    versionTmp <- as.character(sort(as.numeric(msigdb::getMsigdbVersions()),decreasing = T))[1]
    msigdbOut <- msigdb::getMsigdb(org = ifelse(species=="human",'hs', 'mm'), id = 'SYM', version = versionTmp)
    msigdbOut <- msigdb::appendKEGG(msigdbOut, version = versionTmp)
    #
  
    #Start with hallmark
    collectionTmp <- msigdb::subsetCollection(msigdbOut, collection = collection , subcollection = subcollection)
    if(!is.null(subsetNames)){
      collectionTmp <- collectionTmp[names(collectionTmp) %in% subsetNames]
    }
    geneSet_ids <- GSEABase::geneIds(collectionTmp)
  } else {
    checkGeneList(geneSet)
    geneSet_ids <- geneSet
  }
  resOut <- fgsea::fgsea(pathways = geneSet_ids, stat = rankIn, minSize  = minSize, maxSize  = maxSize)
  
  #format output
  resOut$Count <- unlist(lapply(resOut$leadingEdge, length))
  colnames(resOut) <- c("Term","pvalue",'adjusted_pvalue','log2err',"ES",'NES','Size',"Genes",'Count')
  resOut <- resOut[,c(1,5,6,4,9,7,2,3,8)]
  gseaOut <- resOut[order(resOut$adjusted_pvalue),]
  gseaOut$Genes <- sapply(gseaOut$Genes, function(x) paste(x, collapse=':'))

  nameTmp <- ifelse(!is.null(name), paste(name, "gseaAnalysis", sep='_'), "gseaAnalysis")
  data.table::fwrite(gseaOut, file=paste(nameTmp,".txt",sep=''), sep="\t")#, sep2=c("", ":", ""))
  #
  a2sU@analysis@GSEA <- gseaOut
  #
  return(a2sU)
}

####
gseaPlot <- function(ptn,
                     termNames,
                     genesSlopeFiltOut = NULL,
                     gseaParam = 1,
                     ticksSize = 0.3,
                     pdfName = NULL ){
  #
  if (!checkPtn(ptn)) {
    stop("ptn is not a valid 'postNetData' object.")
  }
  if(is.null(a2sU_gsea(a2sU))){
    stop("Please run gseaAnalysis first ")
  } else {
    gseaOut <- a2sU_gsea(a2sU)
  }
  if(!is_number(gseaParam) | !is_number(ticksSize)){
    stop("please provide numeric value")
  }
  #
  effTmp <- a2sU_eff(a2sU)
  if (!is.null(genesSlopeFiltOut)) {
    effIn <- effTmp[!names(effTmp) %in% genesSlopeFiltOut ]
  }  else {
    effIn <- effTmp
  }
  #
  rankIn <- effIn[order(effIn,decreasing = T)]
  #if(!is.null(ads)){
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
  #  #rankedRVMP <- rank(-log10(tmpP) * sign(tmpEff))
  #  rankIn <- tmpEff[order(tmpEff,decreasing = T)]
  #} else if (!is.null(rankIn)){
  #  rankIn <- rankIn[order(rankIn,decreasing = T)]
  #} else {
  #  stop("No anota2seq object or ranks provided")
  #}
  #
  #
  rnk <- rank(-rankIn)
  ord <- order(rnk)
  statsAdj <- rankIn[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  #
  for(tname in termNames){
    termTmp <- tname
    pathGenes <- unlist(strsplit(as.character(gseaOut[gseaOut$Term %in% termTmp,]$Genes),':'))
    #pathGenes<- unname(as.vector(na.omit(match(pathGenes, names(statsAdj)))))
    pathGenes <- match(pathGenes, names(statsAdj))
    pathGenes <- sort(pathGenes)
    gseaRes <- fgsea::calcGseaStat(statsAdj, selectedStats = pathGenes, returnAllExtremes = TRUE)
    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops
    n <- length(statsAdj)
    xs <- as.vector(rbind(pathGenes - 1, pathGenes))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
    diff <- (max(tops) - min(bottoms))/8
    x = y = NULL
    
    nameOut <- ifelse(!is.null(pdfName),paste(pdfName,"gsea", termTmp, sep="_"), paste("gsea", termTmp, sep="_"))
    
    pdf(paste(nameOut,".pdf",sep=''), width = 8, height = 4, useDingbats = F)
    par(mar = c(5, 5, 5, 10), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 0.9, cex.main = 0.7, cex.lab = 0.7)
    peOut <- ggplot2::ggplot(toPlot, ggplot2::aes(x = x, y = y)) + ggplot2::geom_line(color = "grey75", linetype=1, size = 0.75)  + ggplot2::geom_line(color = "grey75") + ggplot2::theme_bw() + 
      ggplot2::theme(panel.border = ggplot2::element_blank(),panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank()) + 
      ggplot2::theme(axis.line.x = ggplot2::element_line(color="black", linewidth = 0.5),axis.line.y = ggplot2::element_line(color="black", linewidth = 0.5)) +
      ggplot2::labs(title=termTmp, x = "rank", y = "enrichment score") +
      ggplot2::geom_segment(data = data.frame(x = pathGenes), mapping = ggplot2::aes(x = x, y = -diff/4, xend = x, yend = diff/4), size = ticksSize, color="firebrick1")
    print(peOut)
    dev.off()
  }
}

