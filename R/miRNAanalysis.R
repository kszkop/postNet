miRNAanalysis <- function(a2sU,
                          miRNATargetScanFile,
                          genesSlopeFiltOut = NULL) {
  #
  if (!checkUtils(a2sU)) {
    stop("a2sU is not a valid 'anota2seqUtilsData' object.")
  }
  miRNATargetScan <- checkFileColumns(miRNATargetScanFile)
  if(length(intersect(unique(miRNATargetScan$Gene.Symbol), anota2seqUtilsGetBg(a2sU)))<2){
    stop('no overlap between genes the dataset and in the miRNATargetScanFile')
  }
  #
  miRNAOut <- list()
  #
  effTmp <- anota2seqUtilsGetEff(a2sU)
  if (!is.null(genesSlopeFiltOut)) {
    effIn <- effTmp[!names(effTmp) %in% genesSlopeFiltOut ]
  }  else {
    effIn <- effTmp
  }
  #
  rankIn <- effIn[order(effIn,decreasing = T)]
  #
  #if(!is.null(ads)){
  #  tmpAds <- anota2seq::anota2seqGetOutput(ads,
  #                                          analysis = regulationGen,
  #                                          output = "full",
  #                                          selContrast = contrastSel,
  #                                          getRVM = TRUE)
  #  #
  #  if (!is.null(genesSlopeFiltOut)) {
  ##    tmpAdsFilt <- tmpAds[row.names(tmpAds) %in% genesSlopeFiltOut, ]
  #  }  else {
  #    tmpAdsFilt <- tmpAds
  #  }
  #  #
  #  tmpP <- tmpAdsFilt[, "apvRvmP"]
  #  tmpEff <- tmpAdsFilt[, "apvEff"]
  #  rankedRVMP <- rank(-log10(tmpP) * sign(tmpEff))
  #} else if (!is.null(customBg)){
  #  tmpAds <- customBg
  #  #
  #  rankedRVMP <- rankIn
  #  names(rankedRVMP) <- tmpAds
  #} else {
  #  stop("No background provided")
  #}
  #miRNATargetScan$Gene.Tax.ID
  #weighted.context...score
  #Site.Type 
  #Gene.Symbol
  #miRNA
  #miRNATargetScan <- read.delim(miRNATargetScanFile)
  #miRNATargetScan <- miRNATargetScan[miRNATargetScan$Gene.Tax.ID == "9606",]
  #
  miRNATargetScanFilt <- miRNATargetScan[miRNATargetScan$weighted.context...score < -0.2 & miRNATargetScan$Site.Type %in% c(1, 2), ]
  miRNAList <- rep(list(NA), length(unique(miRNATargetScanFilt$miRNA)))
  names(miRNAList) <- unique(miRNATargetScanFilt$miRNA)
  for (miRNA in 1:length(miRNAList)) {
    miRNAList[[miRNA]] <- miRNATargetScanFilt$Gene.Symbol[miRNATargetScanFilt$miRNA %in% names(miRNAList)[miRNA]]
  }
  #
  miRNA_RVMP_enrichment_targetScan <- gage::gage(exprs = rankIn, 
                                                 gsets = miRNAList,
                                                 set.size = c(5, 500), 
                                                 same.dir = TRUE, 
                                                 compare = "paired",
                                                 FDR.adj = TRUE,
                                                 saaPrep = gage::gagePrep, 
                                                 saaTest = gage::gs.tTest,
                                                 rank.test = TRUE, 
                                                 use.fold = FALSE, 
                                                 saaSum = gage::gageSum,
                                                 use.stouffer = TRUE
  )
  #
  sigmiRNA_targetScan <- gage::sigGeneSet(miRNA_RVMP_enrichment_targetScan)
  idmiRNA <- rownames(sigmiRNA_targetScan$greater)
  
  #
  miSign <- list()
  if(length(idmiRNA)>0){
    for(i in 1:length(idmiRNA)){
      idTmp <- idmiRNA[i]
      miSign[[idTmp]] <- unique(miRNATargetScanFilt[miRNATargetScanFilt$miRNA==idTmp,]$Gene.Symbol)
    }
    outSign <- signCalc(a2sU, addSign = miSign)
    #
    miRNAOut <- append(list(analysisOut=miRNA_RVMP_enrichment_targetScan), outSign)
  } else {
    miRNAOut <- list(analysisOut=miRNA_RVMP_enrichment_targetScan)
  }
  a2sU@analysis@miRNA <- miRNAOut
  #
  return(miRNAOut)
}
