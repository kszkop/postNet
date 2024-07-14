lengthAnalysis <- function(a2sU,
                           region,
                           comparisons = NULL,
                           plotOut = TRUE,
                           plotType = "boxplot",
                           pdfName = NULL) {
  #
  if (!checkUtils(a2sU)) {
    stop("a2sU is not a valid 'anota2seqUtilsData' object.")
  }
  checkRegion(region)
  #
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
  #
  if(!is.null(comparisons)){
    if(!checkComparisons(comparisons)){
      stop("'comparisons' must be a list of numeric vector for paired comparisons example: list(c(0,2),c(0,1)). 0 is always a background.")
    }
    #
    if(length(which(unique(unlist(comparisons))==0))>0 && is.null(a2sU_bg(a2sU))){
      stop(" 0 is always a background, but no background provided")
    }
  }
  #
  lengthFinal <- list()
  for(reg in region){
    #
    seqTmp <- a2sU_sequences(a2sU,region = reg)
    #
    lenTmp <- as.numeric(sapply(seqTmp, function(x) length(seqinr::s2c(x))))
    #
    lenForAnalysis <- log2(as.numeric(lenTmp))
    names(lenForAnalysis) <- a2sU_geneID(a2sU, region=reg)
    #
    if (isTRUE(plotOut)) {
      resOut <- resQuant(qvec = lenForAnalysis, a2sU = a2sU)
      if(length(resOut)==0){
        stop('There are no regulated genes. Check the input or run without indicating regulation and comparisons')
      }
      colOut <- colPlot(a2sU)
      # Plot
      pdf(ifelse(is.null(pdfName), paste(reg, plotType, "lengthAnalysis.pdf", sep = "_"), paste(pdfName, reg, plotType, "lengthAnalysis.pdf", sep = "_")), width = 8, height = 8, useDingbats = F)
      if (tolower(plotType) == "boxplot"){
        plotBoxplots(resOut, colOut, comparisons = comparisons, ylabel = 'Log2 length')
      } else if (tolower(plotType) == "violin") {
        plotViolin(qvec = lenForAnalysis, a2sU = a2sU, comparisons = comparisons, ylabel = 'Log2 length')
      } else if (tolower(plotType) == "ecdf") {
        plotEcdf(qvec = lenForAnalysis, a2sU = a2sU, comparisons = comparisons, ylabel = 'Log2 length')
      }
      dev.off()
    }
    lengthFinal[[paste(reg, 'length', sep = "_")]] <- lenForAnalysis
  }
  #
  return(lengthFinal)
}
