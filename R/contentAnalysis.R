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
                            subregionSel = NULL,
                            comparisons = NULL,
                            plotOut = TRUE,
                            plotType = "boxplot",
                            pdfName = NULL) {
  #
  checkParameters(annot, ads, regulation, contrast, geneList, geneListcolours, customBg, selection, region, comparisons, plotOut, plotType, contentIn, subregion, subregionSel)

  ####
  annotBg <- gSel(annot = annot, ads = ads, customBg = customBg, geneList = geneList)
  #
  contentFinal <- list()
  for(reg in toupper(region)){
    annotTmp <- regSel(annot = annotBg, region = reg)
    annotBgSel <- isoSel(annot = annotTmp, method = selection)
    #
    if (!is.null(subregion)) {
      if(is.null(subregionSel)){
        stop("You have chosen option to select subset of the sequence. Please provide parameter 'subregionSel' to 'select' or 'exclude'")
      }
      #
      subSeq <- as.character(sapply(annotBgSel$seqTmp, function(x) subset_seq(x, pos = subregion, subregionSel = subregionSel)))
      #
      annotBgSel$seqTmp <- subSeq
    }
    annotBgSel <- annotBgSel[!is.na(annotBgSel$seqTmp), ]
    #
    contentFinal <- list()
    for (i in 1:length(contentIn)) {
      content <- contentIn[i]
    
      contentOut <- as.numeric()
      for (i in 1:nrow(annotBgSel)) {
        tmpSeq <- annotBgSel$seqTmp[i]
        tmpCont <- sapply(seqinr::s2c(toupper(content)), function(x) calc_content(tmpSeq, x))
        contentOut[i] <- sum(tmpCont)
      }
      names(contentOut) <- annotBgSel$geneID
      #
      if (isTRUE(plotOut)) {
        if(is.null(plotType)){
          stop("Please provide 'plotType', to choose from 'boxplot','violin','ecdf'")
        }
        #
        resOut <- resSel(vIn = contentOut, ads = ads, regulation = regulation, contrast = contrast, customBg = customBg, geneList = geneList)
        if(length(resOut)==0){
          stop('There are no regulated genes. Check the input or run without indicating regulation and comparisons')
        }
        coloursOut <- coloursSel(resOut=resOut, geneList = geneList, geneListcolours = geneListcolours)

        # Plot
        pdf(ifelse(is.null(pdfName), paste(reg, content, "Ncontent.pdf", sep = "_"), paste(pdfName, reg, content, "Ncontent.pdf", sep = "_")), width = 8, height = 8, useDingbats = F)
        
        if (plotType == "boxplot"){
          plotBoxplots(resOut, coloursOut, comparisons)
        } else if(plotType == "violin") {
          plotViolin(resOut, coloursOut, comparisons)
        } else if (plotType == "ecdf") {
          plotEcdf(resOut, coloursOut, comparisons)
        }
        dev.off()
      }
      contentFinal[[paste(reg, content, sep = "_")]] <- contentOut
    }
  }
  #
  return(contentFinal)
}

