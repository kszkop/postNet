gseaPlot <- function(ptn,
                     termNames,
                     genesSlopeFiltOut = NULL,
                     gseaParam = 1,
                     ticksSize = 0.3,
                     pdfName = NULL ){
  #
  check_ptn(ptn)
  if(is.null(ptn_GSEA(ptn))){
    stop("Please run gseaAnalysis first ")
  } else {
    gseaOut <- ptn_GSEA(ptn)
  }
  if(!check_number(gseaParam) | !check_number(ticksSize)){
    stop("please provide numeric value")
  }
  #
  effTmp <- ptn_effect(ptn)
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
    peOut <- ggplot2::ggplot(toPlot, ggplot2::aes(x = x, y = y)) + ggplot2::geom_line(color = "grey75", linetype=1, linewidth = 0.75)  + ggplot2::geom_line(color = "grey75") + ggplot2::theme_bw() + 
      ggplot2::theme(panel.border = ggplot2::element_blank(),panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank()) + 
      ggplot2::theme(axis.line.x = ggplot2::element_line(color="black", linewidth = 0.5),axis.line.y = ggplot2::element_line(color="black", linewidth = 0.5)) +
      ggplot2::labs(title=termTmp, x = "rank", y = "enrichment score") +
      ggplot2::geom_segment(data = data.frame(x = pathGenes), mapping = ggplot2::aes(x = x, y = -diff/4, xend = x, yend = diff/4), linewidth = ticksSize, color="firebrick1")
    print(peOut)
    dev.off()
  }
}