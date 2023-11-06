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
  #
  if(!is.null(ads) && !is.null(geneList)){
    stop("please provide anota2seq object or genelist, not both.")
  }
  checkAnnot(annot)
  checkRegion(region)
  checkSelection(selection)

  if(!is.null(ads)){
    if (!checkAds(ads)) {
      stop("ads is not a valid 'Anota2seqDataSet' object.")
    }
    if (!is.null(regulation) && !is.character(regulation) && !regulation %in% c("translationUp","translationDown","translatedmRNAUp","translatedmRNADown","bufferingmRNAUp","bufferingmRNADown","mRNAAbundanceUp","mRNAAbundanceDown","totalmRNAUp","totalmRNADown")) {
      stop("'regulation' should be a character vector chosen from translationUp,translationDown,translatedmRNAUp,translatedmRNADown,bufferingmRNAUp,bufferingmRNADown,mRNAAbundanceUp,mRNAAbundanceDown,totalmRNAUp,totalmRNADown")
    }
    if (!is.null(regulation)){
      if(!is.null(contrast) && !is.numeric(contrast) && !length(contrast) == length(regulation) && !contrast %in% seq(1,ncol(ads@contrasts),1)){
        stop("'contrast' should be a numeric vector chosen from each regulation mode")
      }
    }
  } 
  if(is.null(ads)){
    if(is.null(geneList)){
      stop('Either anota2seq object of gene list must be provided')
    } else {
      if(!checkGeneList(geneList)){
        stop("'geneList' is empty or not named")
      }
    }
  }
  if(!is.null(customBg)){
    if(!is.character(customBg)){
      stop("'customBg' is not character vector")
    }
    if(!length(setdiff(unlist(geneList), customBg)==0)){
      stop("There are entries in geneList that are not in 'customBg'")
    }
  }
  
  if(!is.null(subregion) && (!is.numeric(subregion) || !length(subregion)==1)){
    stop("'subregion' must be a numeric and just number")
  }
  if (!is.null(subregionSel) && is.character(subregionSel) && length(subregionSel) == 1) {
    if (!subregionSel %in% c("select", "exclude")) {
      stop("'subregionSel' must be a character and only 'select' or 'exclude'")
    }
  }
  if(is.null(stremeThreshold) || !is_number(stremeThreshold)){
    stop("please provide numeric p-value threshold for motif selection")
  }
  if(is.null(minwidth) || !is_number(minwidth)){
    stop("please provide numeric minimal width for motif selection")
  }
  if(is.null(memePath)){
    stop("please provide path to meme suit")
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
      #if (!is.null(geneList)) {
      #  motifsIn <- append(motifsIn, names(geneList))
      #}
    } else if (!is.null(geneList)) {
      motifsIn <- names(geneList)
    } else {
      stop("'ads' or 'geneList' must be provided")
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
