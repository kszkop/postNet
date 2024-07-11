uorf_analysis <- function(a2sU,
                          startCodon = "ATG",
                          KozakContext = "strong",
                          onlyUTR5 = FALSE,
                          unitOut = "number",
                          comparisons = NULL,
                          plotOut = TRUE,
                          pdfName = NULL) {
  #
  if (!checkUtils(a2sU)) {
    stop("a2sU is not a valid 'anota2seqUtilsData' object.")
  }
  if(!checkLogicalArgument(plotOut)){
    stop("'plotOut' can only be only be logical: TRUE of FALSE ")
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
  
  seqTmp <- a2sU_sequences(a2sU,"UTR5")
  #
  if (!isTRUE(onlyUTR5)) {
    seq <- list()
    seq[[1]] <- a2sU_sequences(a2sU, 'CDS')
    seq[[2]] <- a2sU_sequences(a2sU, 'UTR3')
    extSeq <- combSeq(seqIn = seq)
    extSeq <- unlist(extSeq)
    #
    uorfOut <- mapply(calc_uORF, seqTmp=seqTmp, ext = extSeq, context = context, unit = tolower(unitOut), USE.NAMES=FALSE)
  } else {
    uorfOut <- sapply(seqTmp, function(x) calc_uORF(x, ext=NULL, context = context, unit = tolower(unitOut)), USE.NAMES=FALSE)
  }
  #
  names(uorfOut) <- names(seqTmp)
  #
  if (tolower(unitOut) == "number" & isTRUE(plotOut)) {
    #
    resOut <- resQuant(qvec = uorfOut, a2sU = a2sU)
    
    if(length(resOut)==0){
      stop('There are no regulated genes. Check the input or run without indicating regulation and comparisons')
    }
    colOut <- colPlot(a2sU)
    #
    resProp <- as.numeric()
    for (i in 1:length(resOut)) {
      resProp[i] <- length(resOut[[i]][resOut[[i]] > 0]) / length(resOut[[i]])
    }
    # Plot
    pdf(ifelse(is.null(pdfName), paste("uORFs_", KozakContext, ".pdf", sep = ""), paste(pdfName, "_uORFs_", KozakContext, ".pdf", sep = "")), width = 8, height = 8, useDingbats = F)
    par(mar = c(8, 12, 5, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
    barplot(resProp_any, col = colOut, xaxt = "n", xlab = "", ylab = "", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, ylim = c(0, ifelse(is.null(comparisons), 1, 1 + (length(comparisons) * 0.1))), space = 0)
    barplot(resProp_strong, add=T, col = colorspace::darken(colOut,0.25),xaxt = "n", xlab = "", ylab = "", main = "", lwd = 1, bty = "n", yaxt = "n", ylim = c(0, ifelse(is.null(comparisons), 1, 1 + (length(comparisons) * 0.1))), space = 0)
    
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
