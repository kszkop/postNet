contentMotifs <- function(annot,
                          motifsIn,
                          seqType = "dna",
                          dist = 1,
                          min_score = 47,
                          unitOut = "number",
                          resid = FALSE,
                          ads = NULL,
                          regulation = NULL,
                          contrast = NULL,
                          geneList = NULL,
                          geneListcolours = NULL,
                          customBg = NULL,
                          selection,
                          region,
                          subregion = NULL,
                          subregionSel=NULL,
                          comparisons = NULL,
                          pdfName = NULL,
                          plotOut = TRUE) {
  
  #
  if(!is.null(ads) && !is.null(geneList)){
    stop("please provide anota2seq object or genelist, not both.")
  }
  
  checkAnnot(annot)
  checkRegion(region)
  checkSelection(selection)
  
  if(!checkLogicalArgument(plotOut)){
    stop("'plotOut' can only be only be logical: TRUE of FALSE ")
  } 
  
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
      if (!is.null(geneListcolours) && !is.character(geneListcolours) && !length(geneListcolours)== length(geneList)) {
        stop("'geneListcolours' should be a character vector of the same length as geneList.")
      }
    }
  }
  if(!is.null(customBg)){
    if(!is.character(customBg)){
      stop("'customBg' is not character vector")
    }
    if(!length(setdiff(unlist(geneList), customBg))==0){
      stop("There are entries in geneList that are not in 'customBg'")
    }
  }
  if(!is.null(comparisons)){
    if(!checkComparisons(comparisons)){
      stop("'comparisons' must be a list of numeric vector for paired comparisons example: list(c(0,2),c(0,1)). 0 is always a background.")
    }
    if(length(which(unique(unlist(list(c(0,2),c(0,1))))==0)>0) && is.null(customBg) && is.null(ads)){
      stop(" 0 is always a background, but no background provided")
    }
  }
  if(!is.null(subregion) && (!is.numeric(subregion) || !length(subregion)==1)){
    stop("'subregion' must be a numeric and just number")
  }
  if (!is.null(subregionSel) && !subregionSel %in% c("select", "exclude")) {
    stop("'subregionSel' must be a character and only 'select' or 'exclude'")
  } 
  if(!is_number(dist)){
    stop("please provide numeric minimal distance between motifs")
  }
  if(!isUnitOut(unitOut)){
    stop("'unitOut' must be one from these: 'numeric' or 'position'")
  }
  if(!checkLogicalArgument(resid)){
    stop("'resid', i.e whether the values should be normalised for the length, can only be only be logical: TRUE of FALSE ")
  } 
  if(!is_valid_seq_type(seqType)){
    stop("'seqType' sequence type must be selected from one of these: 'dna', 'rna' or 'protein' ")
  } 
  if(!is_motifs(motifsIn)){
    stop("'motifsIn' should be not null character vector of sequence motifs ")
  } 

  ###
  annotBg <- gSel(annot = annot, ads = ads, customBg = customBg, geneList = geneList)
  #
  motifFinalRegion <- list()
  for(reg in region){
    annotTmp <- regSel(annot = annotBg, region = reg)
    annotBgSel <- isoSel(annot = annotTmp, method = selection)
    #
    if (tolower(seqType) == "protein") {
      proseqtmp <- as.character(sapply(annotBgSel$seqTmp, function(x) seqinr::c2s(seqinr::translate(seqinr::s2c(x)))))
      #
      annotBgSel$seqTmp <- proseqtmp
      #
      annotBgSel$lenTmp <- annotBgSel$lenTmp / 3
    }
    #
    if (!is.null(subregion)) {
      #
      subSeq <- as.character(sapply(annotBgSel$seqTmp, function(x) subset_seq(x, pos = subregion, subregionSel = subregionSel)))
      #
      annotBgSel$seqTmp <- subSeq
    }
    annotBgSel <- annotBgSel[!is.na(annotBgSel$seqTmp), ]
    #
    motifsFinal <- list()
    for (i in 1:length(motifsIn)) {
      motif <- motifsIn[i]
      #
      if (motif == "G4" & !tolower(seqType)=='protein') {
        if(!is_number(min_score)){
          stop("please provide numeric minimal score for g-quadruplexes selection")
        }
        motifOutTmp <- as.numeric(sapply(annotBgSel$seqTmp, calc_g4, min_score = min_score))
        names(motifOutTmp) <- as.character(annotBgSel$geneID)
      } else {
        motif <- toupper(motif)
        #
        if (tolower(seqType) == "dna" | tolower(seqType) == "rna") {
          motifTmp <- convertIUPAC(motif)
        } else {
          motifTmp <- replaceProtAmbig(motif)
        }
        #
        motifOutTmp <- lapply(annotBgSel$seqTmp, function(x) calc_motif(x, motifIn = motifTmp, dist = dist, unit = unitOut))
        names(motifOutTmp) <- as.character(annotBgSel$geneID)
        if(unitOut == "number"){
          motifOutTmp <- unlist(motifOutTmp)
        } 
      } 
      #
      if (isTRUE(resid) & tolower(unitOut) == 'number') {
        motifOut <- lm(as.numeric(motifOutTmp) ~ log2(as.numeric(annotBgSel$lenTmp)))$residuals
        names(motifOut) <- names(motifOutTmp)
      } else {
        motifOut <- motifOutTmp
      }
      #
      if (tolower(unitOut) == "number" & isTRUE(plotOut)) {
        nameTmp <- ifelse(is.null(pdfName), paste(region, motif, "content.pdf", sep = "_"), paste(pdfName, reg, motif, "content.pdf", sep = "_"))
        nameOut <- nameTmp
        #
        resOut <- resSel(vIn = motifOut, ads = ads, regulation = regulation, contrast = contrast, customBg = customBg, geneList = geneList)
        coloursOut <- coloursSel(resOut=resOut, geneList = geneList, geneListcolours = geneListcolours)
        # Plot
        pdf(nameOut, width = 8, height = 8, useDingbats = F)

        plotEcdf(resOut, coloursOut, comparisons)
        dev.off()
      }
      motifsFinal[[paste(reg, motif, sep = "_")]] <- motifOut
    }
    motifFinalRegion <-  append(motifFinalRegion,motifsFinal)
  }
  return(motifFinalRegion)
}
