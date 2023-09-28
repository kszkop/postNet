motifAnalysis <- function(annot,
                          stremeThreshold = 0.05,
                          minwidth = 6,
                          memePath = NULL,
                          ads=NULL,
                          regulation = NULL,
                          contrast = NULL,
                          geneList = NULL,
                          customBg = NULL,
                          selection = "random",
                          region,
                          subregion = NULL,
                          subregionSel = NULL) {
  ####
  if (is.null(ads) & is.null(customBg)) {
    stop("No background provided")
  }
  ##
  annotBg <- gSel(annot = annot, ads = ads, customBg = customBg, geneList = geneList)
  #
  motifsStemeOutRegion <- list()
  for(reg in region){
    annotTmp <- regSel(annot = annotBg, region = reg)
    annotBgSel <- isoSel(annot = annotTmp, method = selection)
    #
    if (!is.null(subregion)) {
      #
      subSeq <- as.character(sapply(annotBgSel$seqTmp, function(x) subset_seq(x, pos = subregion, subregionSel = subregionSel)))
      #
      annotBgSel$seqTmp <- subSeq
    }
    annotBgSel <- annotBgSel[!is.na(annotBgSel$seqTmp), ]
    #
    seqForAnalysis <- annotBgSel$seqTmp
    names(seqForAnalysis) <- annotBgSel$geneID
    #
    resOut <- resSel(vIn = seqForAnalysis, ads = ads, regulation = regulation, contrast = contrast, customBg = customBg, geneList = geneList)
    
    #motifsStemeOutRegulation <- list()
    #for(regul in regulation){
      #
    if (!is.null(ads)) {
      motifsIn <- regulation
      if (!is.null(geneList)) {
        motifsIn <- append(motifsIn, names(geneList))
      }
    } else if (!is.null(geneList)) {
      motifsIn <- names(geneList)
    }
    #
    motifsTmpOut <- list()
    for (j in 1:length(motifsIn)) {
      # select regulated genes
      motifsSel <- motifsIn[j]

      if (isTRUE(motifsSel %in% regulation)) {
        regIn <- paste(motifsSel, "_c", contrast[j], sep = "")
      } else {
        regIn <- motifsSel
      }
      #
      regSeq <- resOut[[regIn]]
      controlSeq <- seqForAnalysis
      if (j == 1) {
        seqinr::write.fasta(sequences = as.list(as.character(controlSeq)), names = names(controlSeq), file.out = paste(paste("Control", reg, sep = "_"), ".fa", sep = ""))
      }
      seqinr::write.fasta(sequences = as.list(as.character(regSeq)), names = names(regSeq), file.out = paste(paste("Regulated", reg, regIn, sep = "_"), ".fa", sep = ""))
      #
      outdirTmp <- paste("stremeOut", reg, regIn, sep = "_")
      streme_out <- memes::runStreme(input = paste(paste("Regulated", reg, regIn, sep = "_"), ".fa", sep = ""), control = paste(paste("Control", reg, sep = "_"), ".fa", sep = ""), meme_path = memePath, alph = "dna", outdir = outdirTmp, minw = minwidth)
      #
      streme_out <- streme_out[streme_out$pval < stremeThreshold, ]
      motifsTmpOut[[regIn]] <- streme_out
    }
    motifsStemeOut <- append(list(motifsOut = as.character(unlist(lapply(motifsTmpOut, function(x) x$consensus)))), motifsTmpOut)
    motifsStemeOutRegion[[reg]] <- motifsStemeOut
  }
  return(motifsStemeOutRegion)
}
