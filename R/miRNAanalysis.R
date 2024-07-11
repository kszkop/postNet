miRNAanalysis <- function(a2sU,
                          miRNATargetScanFile,
                          genesSlopeFiltOut = NULL,
                          maxSize = 500,
                          minSize = 10) {
  #
  if (!checkUtils(a2sU)) {
    stop("a2sU is not a valid 'anota2seqUtilsData' object.")
  }
  if(!is_number(maxSize) | !is_number(minSize) | !is_number(counts) |!is_number(FDR)){
    stop("please provide numeric value")
  }
  miRNATargetScan <- checkFileColumns(miRNATargetScanFile)
  if(length(intersect(unique(miRNATargetScan$Gene.Symbol), a2sU_bg(a2sU)))<2){
    stop('no overlap between genes the dataset and in the miRNATargetScanFile')
  }
  #
  miRNAOut <- list()
  #
  effTmp <- a2sU_eff(a2sU)
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
  miRNA<- new("anota2seqUtilsmiRNA",
                    miRNA_analysis = miRNA_enrichment_targetScan,
                    miRNA_to_gene = miRNAList)
  a2sU@analysis@miRNA <-  miRNA
  #
  return(a2sU)
}


