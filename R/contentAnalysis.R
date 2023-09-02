contentAnalysis <- function(annot,
                            contentIn,
                            ads = NULL,
                            regulation = NULL,
                            contrast = NULL,
                            geneList = NULL,
                            geneListcolours = NULL,
                            customBg = NULL,
                            selection = "random",
                            region,
                            subregion = NULL, 
                            subregionSel,
                            comparisons = NULL,
                            plotOut = TRUE,
                            plotType = "boxplot",
                            pdfName = NULL) {
    
  ####
  annotBg <- gSel(annot = annot, ads = ads, customBg = customBg, geneList = geneList)
  #
  contentFinal <- list()
  for(reg in region){
    annotTmp <- regSel(annot = annotBg, region = reg)
    annotBgSel <- isoSel(annot = annotTmp, method = selection)
    #
    if (!is.null(subregion)) {
      #
      subSeq <- as.character(sapply(annotBgSel$seqTmp, function(x) subset_seq(x, pos = subregion, subregionSel = subregionSel)))
      #
      annotBgSel$seqTmp <- subSeq
    }
    annotBgSel <- annotBgSel[!is.na(annotBgSel$seqTmp), ]
    #
    contentFinal <- list()
    for (i in 1:length(contentIn)) {
      content <- contentIn[[i]]
    
      contentOut <- as.numeric()
      for (i in 1:nrow(annotBgSel)) {
        tmpSeq <- annotBgSel$seqTmp[i]
        tmpCont <- sapply(seqinr::s2c(toupper(content)), function(x) calc_content(tmpSeq, x))
        contentOut[i] <- sum(tmpCont)
      }
      names(contentOut) <- annotBgSel$geneID
      #
      if (isTRUE(plotOut)) {
        #
        resOut <- resSel(vIn = contentOut, ads = ads, regulation = regulation, contrast = contrast, customBg = customBg, geneList = geneList)
        coloursOut <- coloursSel(ads = ads, regulation = regulation, geneList = geneList, geneListcolours = geneListcolours, customBg = customBg)
        # Plot
        pdf(ifelse(is.null(pdfName), paste(reg, content, "contentAnalysis.pdf", sep = "_"), paste(pdfName, reg, content, "contentAnalysis.pdf", sep = "_")), width = 8, height = 8, useDingbats = F)
        if (plotType == "boxplot" | plotType == "violin") {
          par(mar = c(8, 12, 5, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
          # 
          if (!is.null(regulation)) {
            xlimIn <- c(0.5, length(regulation) + ifelse(!is.null(geneList), length(geneList), 0) + 1.5)
          } else {
            xlimIn <- c(0.5, length(geneList) + 1.5)
          }
          plot(1, 1, xlim = xlimIn, ylim = c(0, 100 + (length(comparisons) * 5)), xaxt = "n", xlab = "", ylab = "", type = "n", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, frame.plot = FALSE)

          axis(side = 2, font = 2, las = 2, lwd = 2, at = seq(0, 100, 25), labels = seq(0, 100, 25))
          text(1:length(resOut), par("usr")[3] - 0.45, labels = names(resOut), xpd = NA, cex = 0.9, srt = 45, adj = 1)
          mtext(side = 2, line = 6, paste(reg, toupper(content), "content (%)", sep = " "), col = "black", font = 2, cex = 1.7, at = 50)
          #
          if (!is.null(ads) | !is.null(customBg)) {
            abline(lty = 5, h = median(resOut[[1]]))
          }
          for (i in 1:length(resOut)) {
            if (plotType == "violin") {
              vioplot::vioplot(resOut[[i]], add = TRUE, at = i, col = coloursOut[i], xaxt = "n", xlab = "", ylab = "", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, frame.plot = FALSE)
            } else if (plotType == "boxplot") {
              boxplot(resOut[[i]], add = TRUE, at = i, col = coloursOut[i], xaxt = "n", xlab = "", ylab = "", type = "n", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, frame.plot = FALSE, outcol = "grey65", whiskcol = "grey65", outline = FALSE, medcol = "black", staplelty = 0, whisklty = 1)
            }
            text(i, 0, round(mean(resOut[[i]], 0)), font = 2)
          }
        } else if (plotType == "ecdf") {
          #
          xlim_min <- as.numeric(quantile(as.numeric(unlist(resOut)), 0.01))
          xlim_max <- as.numeric(quantile(as.numeric(unlist(resOut)), 0.99))
          #
          par(mar = c(5, 5, 8, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
          plot(ecdf(resOut[[1]]), col = coloursOut[1], main = "", xlab = "", ylab = "", verticals = TRUE, do.p = FALSE, lwd = 3, bty = "n", yaxt = "none", font = 2, xlim = c(xlim_min, xlim_max), xaxt = "none")

          for (i in 2:length(resOut)) {
            lines(ecdf(resOut[[i]]), col = coloursOut[i], main = "", xlab = "", verticals = TRUE, do.p = FALSE, lwd = 4)
          }

          mtext(side = 1, line = 4, paste(paste(reg, toupper(content), "content (%)", sep = " "), "\n", region, sep = ""), col = "black", font = 2, cex = 1.2)
          mtext(side = 2, line = 3, "Fn(x)", col = "black", font = 2, cex = 1.2)

          axis(side = 1, seq(floor(xlim_min), ceiling(xlim_max), 10), font = 2, lwd = 2)
          axis(side = 2, seq(0, 1, 0.2), font = 2, las = 2, lwd = 2)
        }
        # Plot stats
        if (!is.null(comparisons)) {
          addStats(comparisons, ads, customBg, plotType, resOut, coloursOut)
        }
        dev.off()
      }
      contentFinal[[paste(reg, content, sep = "_")]] <- contentOut
    }
  }
  #
  return(contentFinal)
}

