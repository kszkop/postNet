foldingEnergyAnalysis <- function(a2sU,
                                  sourceFE = "load",
                                  fromFasta = FALSE,
                                  customFileFE = NULL,
                                  residFE = FALSE,
                                  region = NULL,
                                  comparisons = NULL,
                                  plotOut = TRUE,
                                  plotType = "ecdf",
                                  pdfName = NULL) {
  #
  checkRegion(region)
  
  species <- a2sU_species(a2sU)
  version <- a2sU_version(a2sU)
  
  # Validate the source input
  tryCatch({
    # Code that may throw an error
    checkSourceFE(sourceFE)
  }, error = function(e) {
    stop("Source check failed: ", e$message)
  })

  if(!is_logical(plotOut)){
    stop("'plotOut' can only be only be logical: TRUE of FALSE ")
  }
  if(isTRUE(plotOut)){
    if(!is.null(plotType)){
      checkPlotType(plotType)
    } else {
      stop("Please provide 'plotType' to select option for plotting, from: 'boxplot','violin ,'ecdf'. ")
    }
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
  if(!is_logical(residFE)){
      stop("'residFE', i.e whether the values should be normalised for the length, can only be only be logical: TRUE of FALSE ")
  }
  if(sourceFE=='create'){
    if(!is_logical(fromFasta)){
      stop("'fromFasta' can only be only be logical: TRUE of FALSE ")
    }
    if (isTRUE(fromFasta)) {
      if(is.null(customFileFE)){
        stop("Please provide a fasta file.")
      }
    } else {
      if(!is_valid_species(species)){
        stop("Please specify a species, at the moment only 'human' or 'mouse' are available).") 
      }
    }
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
  #
  if (sourceFE == "create") {
    ##
    message('The calculations can take a very long time, particulariy if there are many long sequences. Using pre-calculation is recommended (option: "load")')
    message('Furthermore, this option will only calculate folding energies and output txt file. Susequently, please use it with option "load" to use them for analysis')
    if (isTRUE(fromFasta)) {
      #
      runMfold(customFileFE)
      #
    } else {
      if(is.null(region)){
        stop("Please provide region")
      }
      #
      currTmp <- list.files(system.file("extdata/annotation/refseq/", package = "anota2seqUtils"))
      if (!species %in% currTmp) {
        stop("This option is only  available for species: human and mouse at the moment. Please use option fromFasta = TRUE")
      }
      #
      for(reg in region){
        #
        seqTmp <- a2sU_sequences(a2sU,region = reg)
        # Write out sequences
        seqinr::write.fasta(sequences = as.list(as.character(seqTmp)), names = names(seqTmp), file.out = paste(reg, ".fa", sep = ""))
        #
        runMfold(paste(reg, ".fa", sep = ""))
      }
    }
    message('Folding energy calculations finished. ')
  } else if (sourceFE == "custom") {
    #
    feOut <- list()
    #
    energyIn <- read.delim(customFileFE, stringsAsFactors = FALSE)
    energyIn$fold_energy <- as.numeric(energyIn$fold_energy)
    #
    feOutTmp <- runFE(energyIn = energyIn, a2sU = a2sU, residFE = residFE)
    #
    if (isTRUE(plotOut)) {
      resOut <- resQuant(qvec = feOutTmp, a2sU = a2sU)
      if(length(resOut)==0){
        stop('There are no regulated genes. Check the input or run without indicating regulation and comparisons')
      }
      colOut <- colPlot(a2sU)
      # Plot
      pdf(ifelse(is.null(pdfName), paste(reg, plotType, "foldEnergyAnalysis.pdf", sep = "_"), paste(pdfName, reg, plotType, "foldEnergyAnalysis.pdf", sep = "_")), width = 8, height = 8, useDingbats = F)
      ylabel <- ifelse(isTRUE(residFE), 'residuals (fe ~ length)', 'folding energy')
      if (tolower(plotType) == "boxplot"){
        plotBoxplots(resOut, colOut, comparisons = comparisons, ylabel = ylabel)
      } else if (tolower(plotType) == "violin") {
        plotViolin(qvec = lenForAnalysis, a2sU = a2sU, comparisons = comparisons, ylabel = ylabel)
      } else if (tolower(plotType) == "ecdf") {
        plotEcdf(qvec = lenForAnalysis, a2sU = a2sU, comparisons = comparisons, ylabel = ylabel)
      }
      dev.off()
    }
    feOut[['custom']] <- feOutTmp
    return(feOut)
    #
  } else if (sourceFE == "load") {
    #
    if(is.null(region)){
      stop("Please provide region")
    }
    #
    feOut <- list()
    # list existing species
    currTmp <- list.files(system.file("extdata/annotation/refseq/", package = "anota2seqUtils"))
    
    if (!species %in% currTmp) {
      stop("This option is only  available for species: human and mouse at the moment. Please use option createFromFile")
    }
    #
    for(reg in region){
      if (species == "human") {
        energyIn <- read.delim(system.file(paste("extdata/annotation/refseq/human", version, sep = "/"), paste("humanDB_", reg, "_foldEnergy", ".txt.gz", sep = ""), package = "anota2seqUtils"), stringsAsFactors = FALSE)
        energyIn$fold_energy <- as.numeric(energyIn$fold_energy)
      }
      if (species == "mouse") {
        energyIn <- read.delim(system.file(paste("extdata/annotation/refseq/mouse", version, sep = "/"), paste("mouseDB_", reg, "_foldEnergy", ".txt.gz", sep = ""), package = "anota2seqUtils"), stringsAsFactors = FALSE)
        energyIn$fold_energy <- as.numeric(energyIn$fold_energy)
      }
      feOutTmp <- runFE(energyIn = energyIn, a2sU = a2sU, residFE = residFE)
      feOut[[paste(reg, "foldingEnergy",sep='_')]] <- feOutTmp
      #
      if (isTRUE(plotOut)) {
        resOut <- resQuant(qvec = feOutTmp, a2sU = a2sU)
        if(length(resOut)==0){
          stop('There are no regulated genes. Check the input or run without indicating regulation and comparisons')
        }
        colOut <- colPlot(a2sU)
        # Plot
        pdf(ifelse(is.null(pdfName), paste(reg, plotType, "foldEnergyAnalysis.pdf", sep = "_"), paste(pdfName, reg, plotType, "foldEnergyAnalysis.pdf", sep = "_")), width = 8, height = 8, useDingbats = F)
        ylabel <- ifelse(isTRUE(residFE), 'residuals (fe ~ length)', 'folding energy')
        if (tolower(plotType) == "boxplot"){
          plotBoxplots(resOut, colOut, comparisons = comparisons, ylabel = ylabel)
        } else if (tolower(plotType) == "violin") {
          plotViolin(qvec = lenForAnalysis, a2sU = a2sU, comparisons = comparisons, ylabel = ylabel)
        } else if (tolower(plotType) == "ecdf") {
          plotEcdf(qvec = lenForAnalysis, a2sU = a2sU, comparisons = comparisons, ylabel = ylabel)
        }
        dev.off()
      }
    }
    return(feOut)
  } else {
    stop("No correct option for source file provided")
  }
}

#
runFE <- function(energyIn,
                  residFE = FALSE,
                  a2sU){
  #
  colnames(energyIn) <- c("id", "fold_energy", "length")
  energyIn$geneID <- a2sU_geneID(a2sU)[match(energyIn$id,a2sU_id(a2sU))]
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

