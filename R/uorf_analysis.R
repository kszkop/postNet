uorf_analysis <- function(annot,
                          startCodon = "ATG",
                          KozakContext = "strong",
                          onlyUTR5 = FALSE,
                          unitOut = "number",
                          ads = NULL,
                          regulation = NULL,
                          contrast = NULL,
                          geneList = NULL,
                          geneListcolours = NULL,
                          customBg = NULL,
                          selection = "random",
                          comparisons = NULL,
                          plotOut = TRUE,
                          pdfName = NULL) {
  #
  if(!is.null(ads) && !is.null(geneList)){
    stop("please provide anota2seq object or genelist, not both.")
  }
  
  checkAnnot(annot)
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
  if(!isStartCodon(startCodon)){
    stop("'startCodon' must be a character vector of length one, and contain only 3 nucleotide sequence, ex. 'ATG'")
  }
  if(!isKozakContext(KozakContext)){
    stop("'KozakContext' must be one from these: 'strong','adequate1','adequate2','weak','any'")
  }
  if(!checkLogicalArgument(onlyUTR5)){
    stop("'onlyUTR5' can only be only be logical: TRUE of FALSE ")
  }
  if(!isUnitOut(unitOut)){
    stop("'unitOut' must be one from these: 'numeric' or 'position'")
  }
  
  #
  KozakContext <- tolower(KozakContext)
  if (KozakContext == "strong") {
    context <-  paste("[AG][ATGC][ATGC]", toupper(startCodon), "G", sep = "")
  } else if (KozakContext == "adequate1") {
    context <- paste("[AG][ATGC][ATGC]", toupper(startCodon), "[ATC]", sep = "")
  } else if (KozakContext == "adequate2") {
    context <- paste("[TC][ATGC][ATGC]", toupper(startCodon), "G", sep = "")
  } else if (KozakContext == "weak") {
    context <- paste("[TC][ATGC][ATGC]", toupper(startCodon), "[ATC]", sep = "")
  } else if (KozakContext == "any") {
    context <- paste("[ATGC][ATGC][ATGC]", toupper(startCodon), "[ATGC]", sep = "")
  } else {
    stop("Please provide correct Kozak context")
  }
  #
  uORFFinal <- list()
  #
  annotBg <- gSel(annot = annot, ads = ads, customBg = customBg, geneList = geneList)
  annotTmp <- regSel(annot = annotBg, region = "UTR5", ext = ifelse(isTRUE(onlyUTR5), FALSE, TRUE))
  annotBgSel <- isoSel(annot = annotTmp, method = selection)
  
  #
  if (!isTRUE(onlyUTR5)) {
    uorfOut <- mapply(calc_uORF, seqTmp=annotBgSel$seqTmp, ext = annotBgSel$extSeq, context = context, unit = tolower(unitOut), USE.NAMES=FALSE)
  } else {
    uorfOut <- sapply(annotBgSel$seqTmp, function(x) calc_uORF(x, ext=NULL, context = context, unit = tolower(unitOut)), USE.NAMES=FALSE)
  }
  #
  names(uorfOut) <- annotBgSel$geneID
  #
  if (tolower(unitOut) == "number" & isTRUE(plotOut)) {
    #
    resOut <- resSel(vIn = uorfOut, ads = ads, regulation = regulation, contrast = contrast, customBg = customBg, geneList = geneList)
    if(length(resOut)==0){
      stop('There are no regulated genes. Check the input or run without indicating regulation and comparisons')
    }
    #
    coloursOut <- coloursSel(resOut=resOut, geneList = geneList, geneListcolours = geneListcolours)
    #
    resProp <- as.numeric()
    for (i in 1:length(resOut)) {
      resProp[i] <- length(resOut[[i]][resOut[[i]] > 0]) / length(resOut[[i]])
    }
    # Plot
    pdf(ifelse(is.null(pdfName), paste("uORFs_", KozakContext, ".pdf", sep = ""), paste(pdfName, "_uORFs_", KozakContext, ".pdf", sep = "")), width = 8, height = 8, useDingbats = F)
    par(mar = c(8, 12, 5, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
    barplot(resProp, col = coloursOut, xaxt = "n", xlab = "", ylab = "", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, ylim = c(0, ifelse(is.null(comparisons), 1, 1 + (length(comparisons) * 0.1))), space = 0)

    axis(side = 2, font = 2, las = 2, lwd = 2, at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2))
    text(seq(0.5, length(resOut), 1), par("usr")[3] - 0.05, labels = names(resOut), xpd = NA, cex = 0.9, srt = 45, adj = 1.1)

    mtext(side = 2, line = 6, paste("proportion of uORFs in \n", KozakContext, " Kozak context", sep = ""), col = "black", font = 2, cex = 1.7, at = 0.5)

    # Plot stats
    if (!is.null(comparisons)) {
      for (j in 1:length(comparisons)) {
        if (names(resOut)[1] == 'background') {
          compTmp <- comparisons[[j]] + 1
        } else {
          compTmp <- comparisons[[j]]
        }
        # stats
        pvalTmp <- format(as.numeric(wilcox.test(resOut[[compTmp[1]]], resOut[[compTmp[2]]], alternative = "two.sided")[3]), scientific = TRUE, digits = 2)
        #
        yposTmp <- (max(as.numeric(resProp)) + 0.05) + (j * 0.1)
        rect(xleft = compTmp[1] - 0.5, xright = compTmp[2] - 0.5, ybottom = yposTmp, ytop = yposTmp, lwd = 2)
        #
        text((sum(compTmp) / 2) - 0.5, yposTmp + 0.05, pvalTmp, cex = 0.75)
      }
    }
    dev.off()
  }
  uORFFinal[[paste('uORFs',startCodon,KozakContext,sep='_')]] <- uorfOut
  #
  return(uORFFinal)
}
