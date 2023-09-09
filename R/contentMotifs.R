contentMotifs <- function(annot,
                          motifsIn,
                          seqType = "dna",
                          dist = 1,
                          min_score = 47,
                          unitOut = "number",
                          resid = FALSE,
                          ads = NULL,
                          regulation = NULL,
                          contrast = NULL,
                          geneList = NULL,
                          geneListcolours = NULL,
                          customBg = NULL,
                          selection,
                          region,
                          subregion = NULL,
                          subregionSel,
                          comparisons = NULL,
                          outDir=NULL,
                          pdfName = NULL,
                          plotOut = TRUE) {
  ####
  if(isTRUE(plotOut)){
    if (!is.null(outDir)){
      dirTmp <- outDir
    } else {
      dirTmp <- paste('motifAnalysis', format(Sys.time(), "%Y%m%e_%X"),sep='_')
    }
    dir.create(dirTmp)
 
    nameTmp <- ifelse(is.null(pdfName), paste(region, motif, "content.pdf", sep = "_"), paste(pdfName, reg, motif, "content.pdf", sep = "_"))
    nameOut <- paste(dirTmp,nameTmp, sep='/')
  }
  #
  annotBg <- gSel(annot = annot, ads = ads, customBg = customBg, geneList = geneList)
  #
  motifFinalRegion <- list()
  for(reg in region){
    annotTmp <- regSel(annot = annotBg, region = reg)
    annotBgSel <- isoSel(annot = annotTmp, method = selection)
    #
    if (seqType == "protein") {
      proseqtmp <- as.character(sapply(annotBgSel$seqTmp, function(x) seqinr::c2s(seqinr::translate(seqinr::s2c(x)))))
      #
      annotBgSel$seqTmp <- proseqtmp
      #
      annotBgSel$lenTmp <- annotBgSel$lenTmp / 3
    }
    #
    if (!is.null(subregion)) {
      #
      subSeq <- as.character(sapply(annotBgSel$seqTmp, function(x) subset_seq(x, pos = subregion, subregionSel = subregionSel)))
      #
      annotBgSel$seqTmp <- subSeq
    }
    annotBgSel <- annotBgSel[!is.na(annotBgSel$seqTmp), ]

    #
    motifsFinal <- list()
    for (i in 1:length(motifsIn)) {
      motif <- motifsIn[[i]]
      #
      if (motif == "G4") {
        motifOutTmp <- as.numeric(sapply(annotBgSel$seqTmp, calc_g4, min_score = min_score))
        names(motifOutTmp) <- as.character(annotBgSel$geneID)
      } else {
        motif <- toupper(motif)
        #
        if (seqType == "dna" | seqType == "rna") {
          motifTmp <- convertIUPAC(motif)
        } else {
          motifTmp <- replaceProtAmbig(motif)
        }
        #
        motifOutTmp <- lapply(annotBgSel$seqTmp, function(x) calc_motif(x, motifIn = motifTmp, dist = dist, unit = unitOut))
        names(motifOutTmp) <- as.character(annotBgSel$geneID)
        if(unitOut == "number"){
          motifOutTmp <- unlist(motifOutTmp)
        } 
      } 
      #
      if (isTRUE(resid) & unitOut == 'number') {
        motifOut <- lm(as.numeric(motifOutTmp) ~ log2(as.numeric(annotBgSel$lenTmp)))$residuals
        names(motifOut) <- names(motifOutTmp)
      } else {
        motifOut <- motifOutTmp
      }
      #
      if (unitOut == "number" & isTRUE(plotOut)) {
        #
        resOut <- resSel(vIn = motifOut, ads = ads, regulation = regulation, contrast = contrast, customBg = customBg, geneList = geneList)
        coloursOut <- coloursSel(ads = ads, regulation = regulation, geneList = geneList, geneListcolours = geneListcolours, customBg = customBg)
        #
        # Plot
        pdf(nameOut, width = 8, height = 8, useDingbats = F)
        par(mar = c(5, 5, 8, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
        #
        xlim_min <- ifelse(isTRUE(resid), floor(quantile(as.numeric(unlist(resOut)), 0.01)), 0)
        xlim_max <- roundUpNice(abs(as.numeric(quantile(as.numeric(unlist(resOut)), 0.99))))

        plot(ecdf(resOut[[1]]), col = coloursOut[1], main = "", xlab = "", ylab = "", verticals = TRUE, do.p = FALSE, lwd = 3, bty = "n", yaxt = "none", font = 2, xlim = c(xlim_min, xlim_max), xaxt = "none")
        #
        for (i in 2:length(resOut)) {
          lines(ecdf(resOut[[i]]), col = coloursOut[i], main = "", xlab = "", verticals = TRUE, do.p = FALSE, lwd = 4)
        }

        mtext(side = 1, line = 4, paste("Number of motifs \n", motif, sep = ""), col = "black", font = 2, cex = 1.2)
        mtext(side = 2, line = 3, "Fn(x)", col = "black", font = 2, cex = 1.2)

        axis(side = 1, seq(xlim_min, xlim_max, 1), font = 2, lwd = 2)
        axis(side = 2, seq(0, 1, 0.2), font = 2, las = 2, lwd = 2)
        #
        if (!is.null(comparisons)) {
          addStats(comparisons, ads, customBg, plotType = "ecdf", resOut, coloursOut)
        }
        dev.off()
      }
      motifsFinal[[paste(reg, motif, sep = "_")]] <- motifOut
    }
    motifFinalRegion <-  append(motifFinalRegion,motifsFinal)
  }
  return(motifFinalRegion)
}
