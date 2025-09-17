codonCalc <- function(ptn,
                      featsel,
                      analysis = "codon",
                      unit = "count",
                      comparisons = NULL,
                      plotOut = TRUE,
                      plotType = "ecdf",
                      pdfName = NULL) {
  #
  check_ptn(ptn)

  if (!is.null(ptn_codonAnalysis(ptn))) {
    codonsAll <- ptn_codonAnalysis(ptn)
    check_codonIn(codonsAll)
  } else {
    stop("Codons analysis is NULL in the 'postNetData' object provided. Please run the codonUsage analysis function prior to running CodonCalc.")
  }
  #
  if (!check_logical(plotOut)) {
    stop("'plotOut' can only be logical: TRUE of FALSE ")
  }
  #
  if (isTRUE(plotOut)) {
    if (!is.null(plotType)) {
      check_plotType(plotType)
    } else {
      stop("Please provide a selection for 'plotType' to specify the method for plotting. Options include: 'boxplot','violin , or 'ecdf'. ")
    }
  }
  #
  if (!is.null(comparisons)) {
    if (!check_comparisons(comparisons)) {
      stop("The input for 'comparisons' must be a list of numeric vectors of paired comparisons. For example: list(c(0,2),c(0,1)). 0 always denotes the background gene set.")
    }
    #
    if (length(which(unique(unlist(comparisons)) == 0)) > 0 && is.null(ptn_background(ptn))) {
      stop("0 always denotes the background, but no background has been provided.")
    }
  }
  #
  if (!is_valid_analysis(analysis)) {
    stop("Please provide a selection for 'analysis'. The options are: 'codon' or 'AA' (amino acid).")
  }
  #
  if (analysis == "codon") {
    if (!check_codons(featsel)) {
      stop("Your input for 'featsel' does not contain valid codons.")
    }
  }
  if (!isUnit(unit)) {
    stop("The 'unit' parameter must be either 'count' or 'freq'.")
  }
  #
  codonCalcOut <- list()
  for (i in 1:length(featsel)) {
    #
    featTmp <- featsel[[i]]
    featNameTmp <- names(featsel)[i]
    if (featNameTmp == 0) {
      stop("The input for 'featSel' should be a named list.")
    }
    # if(is.null(featselName)){
    #  featNameTmp <- paste("codon",names(featsel)[i],sep="_")
    # }
    #
    regName <- names(featsel[i])
    nameTmp <- ifelse(is.null(pdfName), paste("features", regName, "codonCalc.pdf", sep = "_"), paste(pdfName, "features", regName, "codonCalc.pdf", sep = "_"))
    # nameOut <- paste(dirTmp,nameTmp, sep='/')
    nameOut <- nameTmp
    #
    if (analysis == "codon") {
      codonTmp <- codonsAll[codonsAll$codon %in% featTmp, ]
    } else if (analysis == "AA") {
      codonTmp <- codonsAll[codonsAll$AA %in% featTmp, ]
    }
    if (unit == "count") {
      tmp <- codonTmp %>%
        group_by(geneID) %>%
        summarise(count = sum(count))
      codonCalcOutTmp <- tmp$count
    } else if (unit == "freq") {
      tmp <- codonTmp %>%
        group_by(geneID) %>%
        summarise(freq = sum(frequency))
      codonCalcOutTmp <- tmp$freq
    }
    names(codonCalcOutTmp) <- tmp$geneID
    #
    codonCalcOut[[featNameTmp]] <- codonCalcOutTmp
    ##
    if (isTRUE(plotOut)) {
      #
      resOut <- resQuant(qvec = codonCalcOutTmp, ptn = ptn)
      if (length(resOut) == 0) {
        stop("There are no regulated genes in your input. Check the input or run without indicating regulation and comparisons.")
      }
      colOut <- colPlot(ptn)
      pdf(nameOut, width = 8, height = 8, useDingbats = F)
      ylabel <- paste("codon usage(", unit, ")", sep = "")
      plotPostNet(resOut, colOut, comparisons, ylabel = ylabel, plotType = plotType)
      # paste("codon usage(", unit, ")","\n", featNameTmp, ":" ,paste(featTmp,collapse=','),sep='')
      dev.off()
    }
  }
  #
  return(codonCalcOut)
}
