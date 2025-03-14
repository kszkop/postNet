foldingEnergyAnalysis <- function(ptn,
                                  sourceFE = "load",
                                  customFileFE = NULL,
                                  residFE = FALSE,
                                  region = NULL,
                                  comparisons = NULL,
                                  plotOut = TRUE,
                                  plotType = "ecdf",
                                  pdfName = NULL) {
  #
  species <- ptn_species(ptn)
  version <- ptn_version(ptn)
  checkSourceFE(sourceFE)
  
  if(!check_logical(plotOut)){
    stop("'plotOut' can only be only be logical: TRUE of FALSE ")
  } 
  if(isTRUE(plotOut)){
    if(!is.null(plotType)){
      check_plotType(plotType)
    } else {
      stop("Please provide 'plotType' to select option for plotting, from: 'boxplot','violin ,'ecdf'. ")
    }
  }
 
  if(!is.null(comparisons)){
    if(!check_comparisons(comparisons)){
      stop("'comparisons' must be a list of numeric vector for paired comparisons example: list(c(0,2),c(0,1)). 0 is always a background.")
    }
    #
    if(length(which(unique(unlist(comparisons))==0))>0 && is.null(ptn_background(ptn))){
      stop(" 0 is always a background, but no background provided")
    }
  }
  if(!check_logical(residFE)){
      stop("'residFE', i.e whether the values should be normalised for the length, can only be only be logical: TRUE of FALSE ")
  }
  #
  if(sourceFE=='custom'){
    if(is.null(customFileFE)){
      stop("Please provide a custom file.")
    }
  }
  if(sourceFE=='load'){
    if(!is_valid_species(species)){
      stop("Please specify a species, at the moment only 'human' or 'mouse' are available).") 
    }
  }
  if (sourceFE == "custom") {
    #
    feOut <- list()
    #
    energyIn <- read.delim(customFileFE, stringsAsFactors = FALSE)
    energyIn$fold_energy <- as.numeric(energyIn$fold_energy)
    #
    feOutTmp <- runFE(energyIn = energyIn, ptn = ptn, residFE = residFE)
    feOut[['custom']] <- feOutTmp
    #
    if (isTRUE(plotOut)) {
      resOut <- resQuant(qvec = feOutTmp, ptn = ptn)
      if(length(resOut)==0){
        stop('There are no regulated genes. Check the input or run without indicating regulation and comparisons')
      }
      colOut <- colPlot(ptn)
      # Plot
      pdf(ifelse(is.null(pdfName), paste("custom", plotType, "foldEnergyAnalysis.pdf", sep = "_"), paste(pdfName, "custom", plotType, "foldEnergyAnalysis.pdf", sep = "_")), width = 8, height = 8, useDingbats = F)
      ylabel <- ifelse(isTRUE(residFE), 'residuals (fe ~ length)', 'folding energy')
      plotPostNet(resOut, colOut, comparisons, ylabel = ylabel ,plotType = plotType)
      dev.off()
    }
  } else if (sourceFE == "load") {
    #
    if(is.null(region)){
      stop("Please provide region")
    }
    check_region(region)
    #
    feOut <- list()
    # list existing species
    currTmp <- list.files(system.file("extdata/annotation/refseq/", package = "postNet"))
    
    if (!species %in% currTmp) {
      stop("This option is only  available for species: human and mouse at the moment. Please use option createFromFile")
    }
    #
    for(reg in region){
      if (species == "human") {
        energyIn <- read.delim(system.file(paste("extdata/annotation/refseq/human", version, sep = "/"), paste("humanDB_", reg, "_foldEnergy", ".txt.gz", sep = ""), package = "postNet"), stringsAsFactors = FALSE)
        energyIn$fold_energy <- as.numeric(energyIn$fold_energy)
      }
      if (species == "mouse") {
        energyIn <- read.delim(system.file(paste("extdata/annotation/refseq/mouse", version, sep = "/"), paste("mouseDB_", reg, "_foldEnergy", ".txt.gz", sep = ""), package = "postNet"), stringsAsFactors = FALSE)
        energyIn$fold_energy <- as.numeric(energyIn$fold_energy)
      }
      feOutTmp <- runFE(energyIn = energyIn, ptn = ptn, residFE = residFE)
      feOut[[paste(reg, "foldingEnergy",sep='_')]] <- feOutTmp
      #
      if (isTRUE(plotOut)) {
        resOut <- resQuant(qvec = feOutTmp, ptn = ptn)
        if(length(resOut)==0){
          stop('There are no regulated genes. Check the input or run without indicating regulation and comparisons')
        }
        colOut <- colPlot(ptn)
        # Plot
        pdf(ifelse(is.null(pdfName), paste(reg, plotType, "foldEnergyAnalysis.pdf", sep = "_"), paste(pdfName, reg, plotType, "foldEnergyAnalysis.pdf", sep = "_")), width = 8, height = 8, useDingbats = F)
        ylabel <- ifelse(isTRUE(residFE), 'residuals (fe ~ length)', 'folding energy')
        plotPostNet(resOut, colOut, comparisons, ylabel = ylabel ,plotType = plotType)
        dev.off()
      }
    }
    return(feOut)
  } else {
    stop("No correct option for source file provided")
  }
}

runFE <- function(energyIn,
                  residFE = FALSE,
                  ptn){
  #
  colnames(energyIn) <- c("id", "fold_energy", "length")
  energyIn$geneID <- ptn_geneID(ptn,region='CDS')[match(energyIn$id,ptn_id(ptn, region='CDS'))]
  energyIn <- na.omit(energyIn)
  #
  if (isTRUE(residFE)) {
    feForAnalysis <- lm(as.numeric(energyIn$fold_energy) ~ as.numeric(energyIn$length))$residuals
    names(feForAnalysis) <- energyIn$geneID
  } else {
    feForAnalysis <- as.numeric(energyIn$fold_energy)
    names(feForAnalysis) <- energyIn$geneID
  }
  return(feForAnalysis)
}
