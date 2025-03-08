codonCalc <- function(ptn,
                      featsel,
                      analysis="codon",
                      unit = "count",
                      comparisons = NULL,
                      plotOut = TRUE,
                      plotType = 'ecdf',
                      pdfName = NULL) {
  #
  check_ptn(ptn)

  if(!is.null(ptn_codonAnalysis(ptn))){
    codonsAll <- ptn_codonAnalysis(ptn)
    check_codonIn(codonsAll)
  } else {
    stop("codons analysis is null, please provide valid 'postNetData' object and run codonUsage analysis")
  }
  #
  if(!check_logical(plotOut)){
    stop("'plotOut' can only be only be logical: TRUE of FALSE ")
  } 
  #
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
    if(length(which(unique(unlist(comparisons))==0))>0 && is.null(ptn_bg(ptn))){
      stop(" 0 is always a background, but no background provided")
    }
  }
  #
  if(!is_valid_analysis(analysis)){
    stop("Please provide an 'analysis'. It can only be 'codon' or 'AA'")
  } 
  #
  if(analysis=='codon'){
    if(!check_codons(featsel)){
      stop("featsel does not contain valid codons")
    } 
  }
  if(!isUnit(unit)){
    stop("'unit' must be one from these: 'count' or 'freq'")
  }
  #
  codonCalcOut <- list()
  for(i in 1:length(featsel)){
    #
    featTmp <- featsel[[i]]
    featNameTmp <- names(featsel)[i]
    if(featNameTmp==0){
      stop("'featSel should be a named list")
    }
    #if(is.null(featselName)){
    #  featNameTmp <- paste("codon",names(featsel)[i],sep="_")
    #} 
    #
    nameTmp <- ifelse(is.null(pdfName),paste("features",i,"codonCalc.pdf", sep = "_"), paste(pdfName,"features",i,"codonCalc.pdf", sep = "_"))
    #nameOut <- paste(dirTmp,nameTmp, sep='/')
    nameOut <- nameTmp
    #
    if (analysis == "codon") {
      codonTmp <- codonsAll[codonsAll$codon %in% featTmp, ]
    } else if (analysis == "AA") {
      codonTmp <- codonsAll[codonsAll$AA %in% featTmp, ]
    }
    if (unit == "count") {
      tmp <- codonTmp %>% group_by(geneID) %>% summarise(count = sum(count))
      codonCalcOutTmp <- tmp$count
    } else if (unit == "freq") {
      tmp <- codonTmp %>%  group_by(geneID) %>% summarise(freq = sum(frequency))
      codonCalcOutTmp <- tmp$freq
    }
    names(codonCalcOutTmp) <- tmp$geneID
    #
    codonCalcOut[[featNameTmp]] <- codonCalcOutTmp
    ## 
    if (isTRUE(plotOut)) {
      #
      resOut <- resQuant(qvec = codonCalcOutTmp, ptn = ptn)
      if(length(resOut)==0){
        stop('There are no regulated genes. Check the input or run without indicating regulation and comparisons')
      }
      colOut <- colPlot(ptn)
      pdf(nameOut, width = 8, height = 8, useDingbats = F)
      ylabel <- paste("codon usage(", unit, ")",sep = "")
      plotPostNet(resOut, colOut, comparisons, ylabel = ylabel ,plotType = plotType)
      #paste("codon usage(", unit, ")","\n", featNameTmp, ":" ,paste(featTmp,collapse=','),sep='')
      dev.off()
    }
  }
  #
  return(codonCalcOut)
}
