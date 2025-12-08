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
  if (!check_logical(plotOut)) {
    stop("The input for 'plotOut' must be logical: TRUE or FALSE.")
  }
  if (isTRUE(plotOut)) {
    if (!is.null(plotType)) {
      check_plotType(plotType)
    } else {
      stop("Please provide an input for 'plotType'. The options are: 'boxplot', 'violin', or 'ecdf'.")
    }
  }
  #
  if (!is.null(comparisons)) {
    if (!check_comparisons(comparisons)) {
      stop("The input for 'comparisons' must be a list of numeric vectors of paired comparisons. For example: list(c(0,2),c(0,1)). 0 always \ denotes the background gene set.")
    }
    #
    if (length(which(unique(unlist(comparisons)) == 0)) > 0 && is.null(ptn_background(ptn))) {
      stop("0 always denotes the background, but no background has been provided.")
    }
  }
  #
  lengthFinal <- list()
  for (reg in region) {
    #
    seqTmp <- ptn_sequences(ptn, region = reg)
    #
    lenTmp <- as.numeric(sapply(seqTmp, function(x) length(seqinr::s2c(x))))
    #
    lenForAnalysis <- log2(as.numeric(lenTmp))
    names(lenForAnalysis) <- ptn_geneID(ptn, region = reg)
    #
    if (isTRUE(plotOut)) {
      resOut <- resQuant(qvec = lenForAnalysis, ptn = ptn)
      if (length(resOut) == 0) {
        stop("There are no regulated genes in your input. Please check the input or run without indicating 'regulation' and 'comparisons'.")
      }
      colOut <- colPlot(ptn)
      # Plot
      pdf(ifelse(is.null(pdfName), paste(reg, plotType, "lengthAnalysis.pdf", sep = "_"), paste(pdfName, reg, plotType, "lengthAnalysis.pdf", sep = "_")), width = 8, height = 8, useDingbats = F)
      ylabel <- "Length (Log2 scale)"
      plotPostNet(resOut, colOut, comparisons, ylabel = ylabel, plotType = plotType)
      dev.off()
    }
    lengthFinal[[paste(reg, "length", sep = "_")]] <- lenForAnalysis
  }
  #
  return(lengthFinal)
}
