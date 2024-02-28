codonCalc <- function(codonIn,
                      featsel,
                      featselName=NULL,
                      analysis="codon",
                      unit = "count",
                      ads = NULL,
                      regulation = NULL,
                      contrast = NULL,
                      geneList = NULL,
                      geneListcolours = NULL,
                      customBg = NULL,
                      comparisons = NULL,
                      plotOut = TRUE,
                      plotType = NULL,
                      pdfName = NULL) {
  #
  if(!is.null(ads) && !is.null(geneList)){
    stop("please provide anota2seq object or genelist, not both.")
  }
  #
  if(!checkLogicalArgument(plotOut)){
    stop("'plotOut' can only be only be logical: TRUE of FALSE ")
  } 
  #
  if(isTRUE(plotOut)){
    if(!is.null(plotType)){
      checkPlotType(plotType)
    } else {
      stop("Please provide 'plotType' to select option for plotting, from: 'boxplot','violin ,'ecdf'. ")
    }
  }
  #
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
  #
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
  #
  if(!is.null(customBg)){
    if(!is.character(customBg)){
      stop("'customBg' is not character vector")
    }
    if(!length(setdiff(unlist(geneList), customBg))==0){
      stop("There are entries in geneList that are not in 'customBg'")
    }
  }
  #
  if(!is.null(comparisons)){
    if(!checkComparisons(comparisons)){
      stop("'comparisons' must be a list of numeric vector for paired comparisons example: list(c(0,2),c(0,1)). 0 is always a background.")
    }
    if(length(which(unique(unlist(comparisons))==0))>0 && is.null(customBg) && is.null(ads)){
      stop(" 0 is always a background, but no background provided")
    }
  }
  if(!is_valid_analysis(analysis)){
    stop("Please provide an 'analysis'. It can only be 'codon' or 'AA'")
  } 
  
  #
  codonCalcOut <- list()
  #
  for(i in 1:length(featsel)){
    #
    featTmp <- featsel[[i]]
    if(is.null(featselName)){
      featNameTmp <- paste("codon",names(featsel)[i],sep="_")
    } 
    #
    nameTmp <- ifelse(is.null(pdfName),paste("features",i,"codonCalc.pdf", sep = "_"), paste(pdfName,"features",i,"codonCalc.pdf", sep = "_"))
    #nameOut <- paste(dirTmp,nameTmp, sep='/')
    nameOut <- nameTmp
    #
    if (analysis == "codon") {
      codonTmp <- codonIn[codonIn$codon %in% featTmp, ]
    } else if (analysis == "AA") {
      codonTmp <- codonIn[codonIn$AA %in% featTmp, ]
    }
    if (unit == "count") {
      tmp <- codonTmp %>% group_by(geneID) %>% summarise(count = sum(codonCount))
      codonCalcOutTmp <- tmp$count
    } else if (unit == "freq") {
      tmp <- codonTmp %>%  group_by(geneID) %>% summarise(freq = sum(codonFreq))
      codonCalcOutTmp <- tmp$freq
    }
    names(codonCalcOutTmp) <- tmp$geneID
    #
    codonCalcOut[[featNameTmp]] <- codonCalcOutTmp
    ## 
    if (isTRUE(plotOut)) {
      #
      resOut <- resSel(vIn = codonCalcOutTmp, ads = ads, regulation = regulation, contrast = contrast, customBg = customBg, geneList = geneList)
      #
      coloursOut <- coloursSel(resOut, geneList = geneList, geneListcolours = geneListcolours)
      #
      pdf(nameOut, width = 8, height = 8, useDingbats = F)
      par(mar = c(5, 5, 8, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
      xlim_min <- as.numeric(quantile(as.numeric(unlist(resOut)), 0.01))
      xlim_max <- as.numeric(quantile(as.numeric(unlist(resOut)), 0.99))
      #
      plot(ecdf(resOut[[1]]), col = coloursOut[1], main = "", xlab = "", ylab = "", verticals = TRUE, do.p = FALSE, lwd = 3, bty = "n", yaxt = "none", font = 2, xlim = c(xlim_min, xlim_max), xaxt = "none")

      for (i in 2:length(resOut)) {
        lines(ecdf(resOut[[i]]), col = coloursOut[i], main = "", xlab = "", verticals = TRUE, do.p = FALSE, lwd = 4)
      }

      mtext(side = 1, line = 4, paste("codon usage \n", unit, sep = ""), col = "black", font = 2, cex = 1.2)
      mtext(side = 2, line = 3, "Fn(x)", col = "black", font = 2, cex = 1.2)

      axis(side = 1, seq(floor(xlim_min), ceiling(xlim_max), 1), font = 2, lwd = 2)
      axis(side = 2, seq(0, 1, 0.2), font = 2, las = 2, lwd = 2)

      # 
      if (!is.null(comparisons)) {
        addStats(comparisons, plotType, resOut, coloursOut)
      }
      dev.off()
    }
  }
  #
  return(codonCalcOut)
}
