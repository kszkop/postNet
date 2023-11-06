contentAnalysis <- function(annot,
                            contentIn,
                            ads = NULL,
                            regulation = NULL,
                            contrast = NULL,
                            geneList = NULL,
                            geneListcolours = NULL,
                            customBg = NULL,
                            selection = "random",
                            region,
                            subregion = NULL, 
                            subregionSel = NULL,
                            comparisons = NULL,
                            plotOut = TRUE,
                            plotType = "boxplot",
                            pdfName = NULL) {
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
  if(isTRUE(plotOut)){
    if(!is.null(plotType)){
      checkPlotType(plotType)
    }
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
    if(!length(setdiff(unlist(geneList), customBg)==0)){
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

  if(!is.null(contentIn) && !isDNAsequence(contentIn)){
    stop("'contentIn' must be a character vector with DNA sequences")
  }
  if(!is.null(subregion) && (!is.numeric(subregion) || !length(subregion)==1)){
    stop("'subregion' must be a numeric and just number")
  }
  if (!is.null(subregionSel) && is.character(subregionSel) && length(subregionSel) == 1) {
    if (!subregionSel %in% c("select", "exclude")) {
      stop("'subregionSel' must be a character and only 'select' or 'exclude'")
    }
  } 
  
  ####
  annotBg <- gSel(annot = annot, ads = ads, customBg = customBg, geneList = geneList)
  #
  contentFinal <- list()
  for(reg in toupper(region)){
    annotTmp <- regSel(annot = annotBg, region = reg)
    annotBgSel <- isoSel(annot = annotTmp, method = selection)
    #
    if (!is.null(subregion)) {
      if(is.null(subregionSel)){
        stop("You have chosen option to select subset of the sequence. Please provide parameter 'subregionSel' to 'select' or 'exclude'")
      }
      #
      subSeq <- as.character(sapply(annotBgSel$seqTmp, function(x) subset_seq(x, pos = subregion, subregionSel = subregionSel)))
      #
      annotBgSel$seqTmp <- subSeq
    }
    annotBgSel <- annotBgSel[!is.na(annotBgSel$seqTmp), ]
    #
    contentFinal <- list()
    for (i in 1:length(contentIn)) {
      content <- contentIn[i]
    
      contentOut <- as.numeric()
      for (i in 1:nrow(annotBgSel)) {
        tmpSeq <- annotBgSel$seqTmp[i]
        tmpCont <- sapply(seqinr::s2c(toupper(content)), function(x) calc_content(tmpSeq, x))
        contentOut[i] <- sum(tmpCont)
      }
      names(contentOut) <- annotBgSel$geneID
      #
      if (isTRUE(plotOut)) {
        #
        resOut <- resSel(vIn = contentOut, ads = ads, regulation = regulation, contrast = contrast, customBg = customBg, geneList = geneList)
        if(length(resOut)==0){
          stop('There are no regulated genes. Check the input or run without indicating regulation and comparisons')
        }
        coloursOut <- coloursSel(resOut=resOut, geneList = geneList, geneListcolours = geneListcolours)

        # Plot
        pdf(ifelse(is.null(pdfName), paste(reg, content, "Ncontent.pdf", sep = "_"), paste(pdfName, reg, content, "Ncontent.pdf", sep = "_")), width = 8, height = 8, useDingbats = F)
        
        if (plotType == "boxplot"){
          plotBoxplots(resOut, coloursOut, comparisons)
        } else if(plotType == "violin") {
          plotViolin(resOut, coloursOut, comparisons)
        } else if (plotType == "ecdf") {
          plotEcdf(resOut, coloursOut, comparisons)
        }
        dev.off()
      }
      contentFinal[[paste(reg, content, sep = "_")]] <- contentOut
    }
  }
  #
  return(contentFinal)
}

