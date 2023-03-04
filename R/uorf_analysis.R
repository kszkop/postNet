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
  if (KozakContext == "strong") {
    context <-  paste("[AG][ATGC][ATGC]", startCodon, "G", sep = "")
  } else if (KozakContext == "adequate1") {
    context <- paste("[AG][ATGC][ATGC]", startCodon, "[ATC]", sep = "")
  } else if (KozakContext == "adequate2") {
    context <- paste("[TC][ATGC][ATGC]", startCodon, "G", sep = "")
  } else if (KozakContext == "weak") {
    context <- paste("[TC][ATGC][ATGC]", startCodon, "[ATC]", sep = "")
  } else if (KozakContext == "any") {
    context <- paste("[ATGC][ATGC][ATGC]", startCodon, "[ATGC]", sep = "")
  } else {
    stop("Please provide correct Kozak context")
  }
  #
  annotBg <- gSel(annot = annot, ads = ads, customBg = customBg, geneList = geneList)
  annotTmp <- regSel(annot = annotBg, region = "UTR5", ext = ifelse(isTRUE(onlyUTR5), FALSE, TRUE))
  annotBgSel <- isoSel(annot = annotTmp, method = selection)
  
  #
  if (!isTRUE(onlyUTR5)) {
    uorfOut <- mapply(calc_uORF, seqTmp=annotBgSel$seqTmp, ext = annotBgSel$extSeq, context = context, unit = unitOut, USE.NAMES=FALSE)
  } else {
    uorfOut <- sapply(annotBgSel$seqTmp, function(x) calc_uORF(x, ext=NULL, context = context, unit = unitOut), USE.NAMES=FALSE)
  }
  #
  names(uorfOut) <- annotBgSel$geneID
  #
  if (unitOut == "number" & isTRUE(plotOut)) {
    #
    resOut <- resSel(vIn = uorfOut, ads = ads, regulation = regulation, contrast = contrast, customBg = customBg, geneList = geneList)
    #
    coloursOut <- coloursSel(ads = ads, regulation = regulation, geneList = geneList, geneListcolours = geneListcolours, customBg = customBg)
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
        if (!is.null(ads) | !is.null(customBg)) {
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
  #
  return(uorfOut)
}
