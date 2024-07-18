contentAnalysis <- function(a2sU,
                            contentIn,
                            region,
                            subregion = NULL, 
                            subregionSel = NULL,
                            comparisons = NULL,
                            plotOut = TRUE,
                            plotType = "boxplot",
                            pdfName = NULL) {
  #
  if (!checkUtils(a2sU)) {
    stop("a2sU is not a valid 'anota2seqUtilsData' object.")
  }
  checkRegion(region)

  if(!is_logical(plotOut)){
    stop("'plotOut' can only be only be logical: TRUE of FALSE ")
  }
  if(isTRUE(plotOut)){
    if(!is.null(plotType)){
      checkPlotType(plotType)
    } else {
      stop("Please provide 'plotType' to select option for plotting, from: 'boxplot','violin ,'ecdf'. ")
    }
  }
  if(!is.null(comparisons)){
    if(!checkComparisons(comparisons)){
      stop("'comparisons' must be a list of numeric vector for paired comparisons example: list(c(0,2),c(0,1)). 0 is always a background.")
    }
    #
    if(length(which(unique(unlist(comparisons))==0))>0 && is.null(a2sU_bg(a2sU))){
      stop(" 0 is always a background, but no background provided")
    }
  }
  if(!isDNAsequence(contentIn)){
    stop("'contentIn' must be a character vector with DNA sequences")
  }
  if(!is.null(subregion) && (!is.numeric(subregion) || !length(subregion)==1)){
    stop("'subregion' must be a numeric and just number")
  }
  if (!is.null(subregionSel) && !subregionSel %in% c("select", "exclude")) {
      stop("'subregionSel' must be a character and only 'select' or 'exclude'")
  } 
  
  ####
  contentFinal <- list()
  for(reg in toupper(region)){
    #
    seqTmp <- a2sU_sequences(a2sU, region = reg)
    #
    if (!is.null(subregion)) {
      if(is.null(subregionSel)){
        stop("You have chosen option to select subset of the sequence. Please provide parameter 'subregionSel' to 'select' or 'exclude'")
      }
      #
      subSeq <- sapply(seqTmp, function(x) subset_seq(x, pos = subregion, subregionSel = subregionSel))
    }
    #
    for (i in 1:length(contentIn)) {
      content <- contentIn[i]
    
      contentOut <- as.numeric()
      for (i in 1:length(seqTmp)) {
        tmpSeq <- seqTmp[i]
        tmpCont <- sapply(seqinr::s2c(toupper(content)), function(x) calc_content(tmpSeq, x))
        contentOut[i] <- sum(tmpCont)
      }
      names(contentOut) <- a2sU_geneID(a2sU, region=reg)
      #
      if (isTRUE(plotOut)) {
        resOut <- resQuant(qvec = contentOut, a2sU = a2sU)
        if(length(resOut)==0){
          stop('There are no regulated genes. Check the input or run without indicating regulation and comparisons')
        }
        colOut <- colPlot(a2sU)
        # Plot
        pdf(ifelse(is.null(pdfName), paste(reg, content, "content.pdf", sep = "_"), paste(pdfName, reg, content, "content.pdf", sep = "_")), width = 8, height = 8, useDingbats = F)
        ylabel = paste(reg, '%', paste0(content, "content"), sep = "_")
        plotUtils(resOut, colOut, comparisons, ylabel = ylabel ,plotType = plotType)
        dev.off()
      }
      contentFinal[[paste(reg, content, sep = "_")]] <- contentOut
    }
  }
  #
  return(contentFinal)
}

