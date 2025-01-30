miRNAanalysis <- function(ptn,
                          miRNATargetScanFile,
                          genesSlopeFiltOut = NULL,
                          maxSize = 500,
                          minSize = 10) {
  #
  check_ptn(ptn)
  if(!check_number(maxSize) | !check_number(minSize)){
    stop("please provide numeric value")
  }
  if(minSize <= 0 | maxSize <= 0) {
    stop("size parameters must be positive")
  }
  if(maxSize <= minSize) {
    stop("maxSize must be greater than minSize")
  }
  miRNATargetScan <- checkFileColumns(miRNATargetScanFile)
  if(length(intersect(unique(miRNATargetScan$Gene.Symbol), ptn_background(ptn)))<2){
    stop('no overlap between genes the dataset and in the miRNATargetScanFile')
  }
  #
  miRNAOut <- list()
  #
  effTmp <- ptn_effect(ptn)
  if (!is.null(genesSlopeFiltOut)) {
    effIn <- effTmp[!names(effTmp) %in% genesSlopeFiltOut ]
  }  else {
    effIn <- effTmp
  }
  #
  rankIn <- effIn[order(effIn,decreasing = T)]
  #
  miRNATargetScanFilt <- miRNATargetScan[miRNATargetScan$weighted.context...score < -0.2 & miRNATargetScan$Site.Type %in% c(1, 2), ]
  miRNAList <- rep(list(NA), length(unique(miRNATargetScanFilt$miRNA)))
  names(miRNAList) <- unique(miRNATargetScanFilt$miRNA)
  for (miRNA in 1:length(miRNAList)) {
    miRNAList[[miRNA]] <- miRNATargetScanFilt$Gene.Symbol[miRNATargetScanFilt$miRNA %in% names(miRNAList)[miRNA]]
  }
  #
  miRNA_enrichment_targetScan <- gage::gage(exprs = rankIn, 
                                            gsets = miRNAList,
                                            set.size= c(minSize, maxSize), 
                                            rank.test = TRUE, 
                                            use.fold = FALSE)
  #
  miRNA<- new("postNetmiRNA",
                    miRNA_analysis = miRNA_enrichment_targetScan,
                    miRNA_to_gene = miRNAList)
  ptn@analysis@miRNA <-  miRNA
  #
  return(ptn)
}


