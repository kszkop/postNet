lengthAnalysis <- function(ptn,
                           region,
                           comparisons = NULL,
                           plotOut = TRUE,
                           plotType = "boxplot",
                           pdfName = NULL) {
  #
  check_ptn(ptn)
  check_region(region)
  #
  if(!check_logical(plotOut)){
    stop("'plotOut' can only be only be logical: TRUE of FALSE ")
  } 
  if(isTRUE(plotOut)){
    if(!is.null(plotType)){
      check_plotType(plotType)
    } else {
      stop("Please provide 'plotType' to select option for plotting, from: 'boxplot','violin ,'ecdf'. ")
    }
  }
  #
  if(!is.null(comparisons)){
    if(!check_comparisons(comparisons)){
      stop("'comparisons' must be a list of numeric vector for paired comparisons example: list(c(0,2),c(0,1)). 0 is always a background.")
    }
    #
    if(length(which(unique(unlist(comparisons))==0))>0 && is.null(ptn_background(ptn))){
      stop(" 0 is always a background, but no background provided")
    }
  }
  #
  lengthFinal <- list()
  for(reg in region){
    #
    seqTmp <- ptn_sequences(ptn,region = reg)
    #
    lenTmp <- as.numeric(sapply(seqTmp, function(x) length(seqinr::s2c(x))))
    #
    lenForAnalysis <- log2(as.numeric(lenTmp))
    names(lenForAnalysis) <- ptn_geneID(ptn, region=reg)
    #
    if (isTRUE(plotOut)) {
      resOut <- resQuant(qvec = lenForAnalysis, ptn = ptn)
      if(length(resOut)==0){
        stop('There are no regulated genes. Check the input or run without indicating regulation and comparisons')
      }
      colOut <- colPlot(ptn)
      # Plot
      pdf(ifelse(is.null(pdfName), paste(reg, plotType, "lengthAnalysis.pdf", sep = "_"), paste(pdfName, reg, plotType, "lengthAnalysis.pdf", sep = "_")), width = 8, height = 8, useDingbats = F)
      ylabel = 'Log2 length'
      plotUtils(resOut, colOut, comparisons, ylabel = ylabel ,plotType = plotType)
      dev.off()
    }
    lengthFinal[[paste(reg, 'length', sep = "_")]] <- lenForAnalysis
  }
  #
  return(lengthFinal)
}
