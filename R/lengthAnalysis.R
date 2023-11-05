lengthAnalysis <- function(annot,
                           ads = NULL,
                           regulation = NULL,
                           contrast = NULL,
                           geneList = NULL,
                           geneListcolours = NULL,
                           customBg = NULL,
                           selection = "random",
                           region,
                           comparisons = NULL,
                           plotOut = TRUE,
                           plotType = "boxplot",
                           pdfName = NULL) {
  #
  checkParameters(annot, ads, regulation, contrast, geneList, geneListcolours, customBg, selection, region, comparisons, plotOut, plotType)
  #
  annotBg <- gSel(annot = annot, ads = ads, customBg = customBg, geneList = geneList)
  #
  lengthFinal <- list()
  for(reg in region){
    annotTmp <- regSel(annot = annotBg, region = reg)
    annotBgSel <- isoSel(annot = annotTmp, method = selection)
    #
    lenForAnalysis <- log2(as.numeric(annotBgSel$lenTmp))
    names(lenForAnalysis) <- annotBgSel$geneID
    #
    if (isTRUE(plotOut)) {
      if(is.null(plotType)){
        stop("Please provide 'plotType', to choose from 'boxplot','violin','ecdf'")
      }
      #
      resOut <- resSel(vIn = lenForAnalysis, ads = ads, regulation = regulation, contrast = contrast, customBg = customBg, geneList = geneList)
      if(length(resOut)==0){
        stop('There are no regulated genes. Check the input or run without indicating regulation and comparisons')
      }
      coloursOut <- coloursSel(resOut=resOut, geneList = geneList, geneListcolours = geneListcolours)
      
      # Plot
      pdf(ifelse(is.null(pdfName), paste(reg, plotType, "lengthAnalysis.pdf", sep = "_"), paste(pdfName, reg, plotType, "lengthAnalysis.pdf", sep = "_")), width = 8, height = 8, useDingbats = F)
      if (plotType == "boxplot"){
        plotBoxplots(resOut, coloursOut, comparisons)
      } else if(plotType == "violin") {
        plotViolin(resOut, coloursOut, comparisons)
      } else if (plotType == "ecdf") {
        plotEcdf(resOut, coloursOut, comparisons)
      }
      dev.off()
    }
    lengthFinal[[paste(reg, 'length', sep = "_")]] <- lenForAnalysis
  }
  #
  return(lengthFinal)
}
