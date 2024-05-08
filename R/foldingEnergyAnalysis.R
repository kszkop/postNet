foldingEnergyAnalysis <- function(a2sU,
                                  sourceFE = "load",
                                  #version = NULL,
                                  #species = NULL,
                                  fromFasta = FALSE,
                                  customFileFE = NULL,
                                  #onlyRun = FALSE,
                                  residFE = FALSE,
                                  region = NULL,
                                  comparisons = NULL,
                                  plotOut = TRUE,
                                  plotType = "boxplot",
                                  pdfName = NULL) {
  #
  checkRegion(region)
  
  species <- anota2seqUtilsGetSpecies(a2sU)
  version <- anota2seqUtilsGetVersion(a2sU)
  
  # Validate the source input
  tryCatch({
    # Code that may throw an error
    checkSourceFE(sourceFE)
  }, error = function(e) {
    stop("Source check failed: ", e$message)
  })

  if(!checkLogicalArgument(plotOut)){
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
    if(length(which(unique(unlist(comparisons))==0))>0 && is.null(anota2seqUtilsGetBg(a2sU))){
      stop(" 0 is always a background, but no background provided")
    }
  }
  if(!checkLogicalArgument(residFE)){
    stop("'residFE', i.e whether the values should be normalised for the length, can only be only be logical: TRUE of FALSE ")
  } 
  
  if(sourceFE=='create'){
    if(!checkLogicalArgument(fromFasta)){
      stop("'fromFasta', can only be only be logical: TRUE of FALSE ")
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
    if (isTRUE(fromFasta)) {
      #
      runMfold(customFileFE)
      #
      #if (!isTRUE(onlyRun)) {
        energyIn <- read.delim(gsub(".fa", "_foldEnergy.txt", customFileFE), stringsAsFactors = FALSE)
        energyIn <- energyIn[!grepl("NM_", energyIn$fold_energy), ]
        energyIn$fold_energy <- as.numeric(energyIn$fold_energy)
      #}
    } else {
      currTmp <- list.files(system.file("extdata/annotation/refseq/", package = "anota2seqUtils"))
      
      if (!species %in% currTmp) {
        stop("This option is only  available for species: human and mouse at the moment. Please use option createFromFile")
      }
      #
      if (is.null(version)) {
        version <- checkAvailableVersions(species = species)
        # extract the latest
        versionInd <- sub("^[^.]*.", "", version)
        versionInd <- sort(versionInd, decreasing = T)[1]
        version <- version[grep(versionInd, version)]
      }
      
      annot <- retrieveFormatData(source='load', species='human', version = version)
      
      for(reg in region){
        # Select region of interest
        if (reg  == "UTR5") {
          # Write out sequences
          seqinr::write.fasta(sequences = as.list(annot$UTR5_seq), names = annot$id, file.out = paste(reg, ".fa", sep = ""))
          #
          runMfold(paste(reg, ".fa", sep = ""))
        }
        if (reg == "UTR3") {
          # Write out sequences
          seqinr::write.fasta(sequences = as.list(annot$UTR3_seq), names = annot$id, file.out = paste(reg, ".fa", sep = ""))
          #
          runMfold(paste(reg, ".fa", sep = ""))
        }
        if (region == "CDS") {
          # Write out sequences
          seqinr::write.fasta(sequences = as.list(annot$CDS_seq), names = annot$id, file.out = paste(reg, ".fa", sep = ""))
          #
          runMfold(paste(reg, ".fa", sep = ""))
        #
        }
      }
    }
  } else if (sourceFE == "custom") {
    #
    feOut <- list()
    #
    energyIn <- read.delim(customFileFE, stringsAsFactors = FALSE)
    energyIn$fold_energy <- as.numeric(energyIn$fold_energy)
    #
    feOutTmp <- runFE(energyIn = energyIn, a2sU, comparisons=comparisons, residFE = residFE,  plotOut=plotOut, plotType = plotType, pdfName=pdfName)
    feOut[['custom']] <- feOutTmp
    #
  } else if (sourceFE == "load") {
    #
    feOut <- list()
    # list existing species
    currTmp <- list.files(system.file("extdata/annotation/refseq/", package = "anota2seqUtils"))
    
    if (!species %in% currTmp) {
      stop("This option is only  available for species: human and mouse at the moment. Please use option createFromFile")
    }
    #
    if (is.null(version)) {
      version <- checkAvailableVersions(species = species)
      # extract the latest
      versionInd <- sub("^[^.]*.", "", version)
      versionInd <- sort(versionInd, decreasing = T)[1]
      version <- version[grep(versionInd, version)]
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
      feOutTmp <- runFE(energyIn = energyIn, annot = annot,  ads = ads, region = reg, regulation = regulation, contrast = contrast, customBg = customBg, geneList = geneList, geneListcolours = geneListcolours, comparisons=comparisons, selection = selection, residFE = residFE,  plotOut=plotOut, plotType = plotType, pdfName=pdfName)
      feOut[[paste(reg, "foldingEnergy",sep='_')]] <- feOutTmp
    }
  } else {
    stop("No correct option for source file provided")
  }
  #
  return(feOut)
}

#
runFE <- function(energyIn,
                  a2sU,
                  region = NULL,
                  comparisons = NULL,
                  residFE = FALSE, 
                  plotOut = TRUE,
                  plotType = "boxplot",
                  pdfName = NULL){
  #
  energyInGene <- merge(energyIn, annot[, c(1, 2)], by = "id", all.x = T)
  energyInGene <- na.omit(energyInGene)
  colnames(energyInGene) <- c("id", "fold_energy", "lenTmp", "geneID")
  
  energyInGeneBg <- gSel(annot = energyInGene, ads = ads, customBg = customBg, geneList = geneList)
  energyInGeneSel <- isoSel(annot = energyInGeneBg, method = selection)
  
  #
  if (isTRUE(residFE)) {
    feForAnalysis <- lm(as.numeric(energyInGeneSel$fold_energy) ~ as.numeric(energyInGeneSel$lenTmp))$residuals
    names(feForAnalysis) <- energyInGeneSel$geneID
  } else {
    feForAnalysis <- as.numeric(energyInGeneSel$fold_energy)
    names(feForAnalysis) <- energyInGeneSel$geneID
  }
  
  #
  if (isTRUE(plotOut)) {
    resOut <- resSel(vIn = feForAnalysis, ads = ads, regulation = regulation, contrast = contrast, customBg = customBg, geneList = geneList)
    coloursOut <- coloursSel(resOut=resOut, geneList = geneList, geneListcolours = geneListcolours)
    
    pdf(ifelse(is.null(pdfName), paste(ifelse(!is.null(region), region, ""), plotType, "foldEnergyAnalysis.pdf", sep = "_"), paste(pdfName, ifelse(!is.null(region), region, ""), plotType, "foldenergyAnalysis.pdf", sep = "_")), width = 8, height = 8, useDingbats = F)
    #
    if (tolower(plotType) == "boxplot"){
      plotBoxplots(resOut, coloursOut, comparisons)
    } else if (tolower(plotType) == "violin") {
      plotViolin(resOut, coloursOut, comparisons)
    } else if (tolower(plotType) == "ecdf") {
      plotEcdf(resOut, coloursOut, comparisons)
    }
    #
    dev.off()
  }
  return(feForAnalysis)
}
