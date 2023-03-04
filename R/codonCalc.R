codonCalc <- function(codonIn,
                      featsel,
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
                      plotType=NULL,
                      pdfName = NULL) {
  if (analysis == "codon") {
    codonTmp <- codonIn[codonIn$codon %in% featsel, ]
  } else if (analysis == "AA") {
    codonTmp <- codonIn[codonIn$AA %in% featsel, ]
  }
  if (unit == "count") {
    tmp <- codonTmp %>% group_by(geneID) %>% summarise(count = sum(codonCount))
    codonCalcOut <- tmp$count
  } else if (unit == "freq") {
    tmp <- codonTmp %>%  group_by(geneID) %>% summarise(freq = sum(codonFreq))
    codonCalcOut <- tmp$freq
  }
  names(codonCalcOut) <- tmp$geneID

  ## 
  if (isTRUE(plotOut)) {
    #
    resOut <- resSel(vIn = codonCalcOut, ads = ads, regulation = regulation, contrast = contrast, customBg = customBg, geneList = geneList)
    #
    coloursOut <- coloursSel(ads = ads, regulation = regulation, geneList = geneList, geneListcolours = geneListcolours, customBg = customBg)
    #
    pdf(ifelse(is.null(pdfName), "codonCalc.pdf", paste(pdfName, "codonCalc.pdf", sep = "_")), width = 8, height = 8, useDingbats = F)
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
      addStats(comparisons, ads, customBg, plotType, resOut, coloursOut)
    }
    dev.off()
  }
  #
  return(codonCalcOut)
}
