foldingEnergyAnalysis <- function(annot,
                                  sourceFE = "load",
                                  version = NULL,
                                  species = NULL,
                                  fromFasta = FALSE,
                                  customFileFE = NULL,
                                  #onlyRun = FALSE,
                                  residFE = FALSE,
                                  ads = NULL,
                                  regulation = NULL,
                                  contrast = NULL,
                                  geneList = NULL,
                                  geneListcolours = NULL,
                                  customBg = NULL,
                                  selection,
                                  region = NULL,
                                  comparisons = NULL,
                                  plotOut = TRUE,
                                  plotType = "boxplot",
                                  pdfName = NULL) {
  #
  if(!is.null(ads) && !is.null(geneList)){
    stop("please provide anota2seq object or genelist, not both.")
  }
  
  checkAnnot(annot)
  checkRegion(region)
  checkSelection(selection)
  
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
  if(!is.null(ads)){
    if (!checkAds(ads)) {
      stop("ads is not a valid 'Anota2seqDataSet' object.")
    }
    if (!is.null(regulation) && !is.character(regulation) && !regulation %in% c("translationUp","translationDown","translatedmRNAUp","translatedmRNADown","bufferingmRNAUp","bufferingmRNADown","mRNAAbundanceUp","mRNAAbundanceDown","totalmRNAUp","totalmRNADown")) {
      stop("'regulation' should be a character vector chosen from translationUp,translationDown,translatedmRNAUp,translatedmRNADown,bufferingmRNAUp,bufferingmRNADown,mRNAAbundanceUp,mRNAAbundanceDown,totalmRNAUp,totalmRNADown")
    }
    if (!is.null(regulation)){
      if(!is.null(contrast) && !is.numeric(contrast) && !length(contrast) == length(regulation) && !contrast %in% seq(1,ncol(ads@contrasts),1)){
        stop("'contrast' should be a numeric vector chosen from each regulation mode")
      }
    }
  } 
  if(is.null(ads)){
    if(is.null(geneList)){
      stop('Either anota2seq object of gene list must be provided')
    } else {
      if(!checkGeneList(geneList)){
        stop("'geneList' is empty or not named")
      }
      if (!is.null(geneListcolours) && !is.character(geneListcolours) && !length(geneListcolours)== length(geneList)) {
        stop("'geneListcolours' should be a character vector of the same length as geneList.")
      }
    }
  }
  if(!is.null(customBg)){
    if(!is.character(customBg)){
      stop("'customBg' is not character vector")
    }
    if(!length(setdiff(unlist(geneList), customBg)==0)){
      stop("There are entries in geneList that are not in 'customBg'")
    }
  }
  if(!is.null(comparisons)){
    if(!checkComparisons(comparisons)){
      stop("'comparisons' must be a list of numeric vector for paired comparisons example: list(c(0,2),c(0,1)). 0 is always a background.")
    }
    if(length(which(unique(unlist(list(c(0,2),c(0,1))))==0)>0) && is.null(customBg) && is.null(ads)){
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
      if(is.null(customFileFE){
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
    if(is.null(customFileFE){
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
    feOutTmp <- runFE(energyIn = energyIn, annot = annot,  ads = ads, region = 'custom', regulation = regulation, contrast = contrast, customBg = customBg, geneList = geneList, geneListcolours = geneListcolours, comparisons=comparisons, selection = selection, residFE = residFE,  plotOut=plotOut, plotType = plotType, pdfName=pdfName)
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
runFE <- function(energyIn = energyIn, 
                  annot = annot, 
                  ads = ads,
                  region = reg,
                  regulation = regulation, 
                  contrast = contrast,
                  customBg = customBg, 
                  geneList = geneList, 
                  geneListcolours = geneListcolours, 
                  comparisons=comparisons,
                  selection = selection, 
                  residFE = residFE, 
                  plotOut=plotOut,
                  plotType = plotType,
                  pdfName=pdfName){
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
    coloursOut <- coloursSel(ads = ads, regulation = regulation, geneList = geneList, geneListcolours = geneListcolours, customBg = customBg)
    
    pdf(ifelse(is.null(pdfName), paste(ifelse(!is.null(region), region, ""), plotType, "foldEnergyAnalysis.pdf", sep = "_"), paste(pdfName, ifelse(!is.null(region), region, ""), plotType, "foldenergyAnalysis.pdf", sep = "_")), width = 8, height = 8, useDingbats = F)
    #
    xlim_set <- roundUpNice(abs(max(abs(as.numeric(quantile(as.numeric(unlist(resOut)), 0.01))), abs(as.numeric(quantile(as.numeric(unlist(resOut)), 0.99))))))
    #
    if (plotType == "boxplot" | plotType == "violin") {
      #
      if (!is.null(regulation)) {
        xlimIn <- c(0.5, length(regulation) + ifelse(!is.null(geneList), length(geneList), 0) + 1.5)
      } else {
        xlimIn <- c(0.5, length(geneList) + 1.5)
      }
      par(mar = c(8, 5, 8, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
      plot(1, 1, xlim = xlimIn, ylim = c(-xlim_set, ifelse(isTRUE(residFE), xlim_set, 50)), xaxt = "n", xlab = "", ylab = "", type = "n", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, frame.plot = FALSE)
      for (i in 1:length(resOut)) {
        if (plotType == "violin") {
          vioplot::vioplot(resOut[[i]], add = TRUE, at = i, col = coloursOut[i], xaxt = "n", xlab = "", ylab = "", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, frame.plot = FALSE)
        } else if (plotType == "boxplot") {
          boxplot(resOut[[i]], add = TRUE, at = i, col = coloursOut[i], xaxt = "n", xlab = "", ylab = "", type = "n", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, frame.plot = FALSE, outcol = "grey65", whiskcol = "grey65", outline = FALSE, medcol = "black", staplelty = 0, whisklty = 1)
        }
        text(i, -xlim_set, round(mean(resOut[[i]], 0)), font = 2)
      }
      axis(side = 2, font = 2, las = 2, lwd = 2, at = seq(-xlim_set, ifelse(isTRUE(residFE), xlim_set, 50), 20), labels = seq(-xlim_set, ifelse(isTRUE(residFE), xlim_set, 50), 20))
      
      mtext(side = 2, line = 6, "folding energy", col = "black", font = 2, cex = 1.7, at = 0)
      if (!is.null(ads) | !is.null(customBg)) {
        abline(lty = 5, h = median(resOut[[1]]))
      }
      text(1:length(resOut), par("usr")[3] - 0.45, labels = names(resOut), xpd = NA, cex = 0.9, srt = 45, adj = 1)
    } else if (plotType == "ecdf") {
      par(mar = c(5, 5, 8, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
      plot(ecdf(resOut[[1]]), col = coloursOut[1], main = "", xlab = "", ylab = "", verticals = TRUE, do.p = FALSE, lwd = 3, bty = "n", yaxt = "none", font = 2, xlim = c(-xlim_set, ifelse(isTRUE(residFE), xlim_set, 0)), xaxt = "none")
      
      mtext(side = 1, line = 4, paste("folding energy", "\n", paste(ifelse(!is.null(region), region, ""), sep = "")), col = "black", font = 2, cex = 1.2)
      mtext(side = 2, line = 3, "Fn(x)", col = "black", font = 2, cex = 1.2)
      
      axis(side = 1, seq(-xlim_set, ifelse(isTRUE(residFE), xlim_set, 0), 20), font = 2, lwd = 2)
      axis(side = 2, seq(0, 1, 0.2), font = 2, las = 2, lwd = 2)
      
      for (i in 2:length(resOut)) {
        lines(ecdf(resOut[[i]]), col = coloursOut[i], main = "", xlab = "", verticals = TRUE, do.p = FALSE, lwd = 4)
      }
      if (!is.null(comparisons)) {
        addStats(comparisons, ads, customBg, plotType, resOut, coloursOut)
      }
    }
    #
    dev.off()
  }
  return(feForAnalysis)
}
