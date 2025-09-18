foldingEnergyAnalysis <- function(ptn,
                                  sourceFE = "load",
                                  customFileFE = NULL,
                                  residFE = FALSE,
                                  region,
                                  comparisons = NULL,
                                  plotOut = TRUE,
                                  plotType = "ecdf",
                                  pdfName = NULL) {
  #
  species <- ptn_species(ptn)
  version <- ptn_version(ptn)
  checkSourceFE(sourceFE)

  if (!check_logical(plotOut)) {
    stop("The input for 'plotOut' must be logical: TRUE or FALSE.")
  }
  if (isTRUE(plotOut)) {
    if (!is.null(plotType)) {
      check_plotType(plotType)
    } else {
      stop("Please provide an input for 'plotType'. The options are: 'boxplot', 'violin', or 'ecdf'.")
    }
  }

  if (!is.null(comparisons)) {
    if (!check_comparisons(comparisons)) {
      stop("The input for 'comparisons' must be a list of numeric vectors of paired comparisons. For example: list(c(0,2),c(0,1)). 0 always \
           denotes the background gene set.")
    }
    #
    if (length(which(unique(unlist(comparisons)) == 0)) > 0 && is.null(ptn_background(ptn))) {
      stop("0 always denotes the background, but no background has been provided.")
    }
  }
  if (!check_logical(residFE)) {
    stop("The input for 'residFE' (specifying whether the values should be normalised for the sequence length) must be logical: TRUE of FALSE.")
  }
  #
  if (sourceFE == "custom") {
    if (is.null(customFileFE)) {
      stop("Please provide a custom file with folding energies.")
    }
  }
  if (sourceFE == "load") {
    if (!is_valid_species(species)) {
      stop("Please provide an input for 'species'. Currently 'human' or 'mouse' are available.")
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
    feOut[["custom"]] <- feOutTmp
    #
    if (isTRUE(plotOut)) {
      resOut <- resQuant(qvec = feOutTmp, ptn = ptn)
      if (length(resOut) == 0) {
        stop("There are no regulated genes in your input. Please check the input or run without indicating 'regulation' and 'comparisons'.")
      }
      colOut <- colPlot(ptn)
      # Plot
      pdf(ifelse(is.null(pdfName), paste("custom", plotType, "foldEnergyAnalysis.pdf", sep = "_"), paste(pdfName, "custom", plotType, "foldEnergyAnalysis.pdf", sep = "_")), width = 8, height = 8, useDingbats = F)
      ylabel <- ifelse(isTRUE(residFE), "Residuals (FE ~ Length)", "Folding Energy")
      plotPostNet(resOut, colOut, comparisons, ylabel = ylabel, plotType = plotType)
      dev.off()
    }
  } else if (sourceFE == "load") {
    #
    if (is.null(region)) {
      stop("Please provide an input for 'region'.")
    }
    check_region(region)
    #
    feOut <- list()
    # list existing species
    currTmp <- list.files(system.file("extdata/annotation/refseq/", package = "postNet"))

    if (!species %in% currTmp) {
      stop("The 'load' option for 'sourceFE' is currently only available for human and mouse. Custom folding energies can be provided by specifying the \
           'custom' option, and using the 'customFileFE' parameter.")
    }
    #
    for (reg in region) {
      if (species == "human") {
        energyIn <- read.delim(system.file(paste("extdata/annotation/refseq/human", version, sep = "/"), paste("humanDB_", reg, "_foldEnergy", ".txt.gz", sep = ""), package = "postNet"), stringsAsFactors = FALSE)
        energyIn$fold_energy <- as.numeric(energyIn$fold_energy)
      }
      if (species == "mouse") {
        energyIn <- read.delim(system.file(paste("extdata/annotation/refseq/mouse", version, sep = "/"), paste("mouseDB_", reg, "_foldEnergy", ".txt.gz", sep = ""), package = "postNet"), stringsAsFactors = FALSE)
        energyIn$fold_energy <- as.numeric(energyIn$fold_energy)
      }
      feOutTmp <- runFE(energyIn = energyIn, ptn = ptn, residFE = residFE)
      feOut[[paste(reg, "foldingEnergy", sep = "_")]] <- feOutTmp
      #
      if (isTRUE(plotOut)) {
        resOut <- resQuant(qvec = feOutTmp, ptn = ptn)
        if (length(resOut) == 0) {
          stop("There are no regulated genes in your input. Please check the input or run without indicating 'regulation' and 'comparisons'.")
        }
        colOut <- colPlot(ptn)
        # Plot
        pdf(ifelse(is.null(pdfName), paste(reg, plotType, "foldEnergyAnalysis.pdf", sep = "_"), paste(pdfName, reg, plotType, "foldEnergyAnalysis.pdf", sep = "_")), width = 8, height = 8, useDingbats = F)
        ylabel <- ifelse(isTRUE(residFE), "Residuals (FE ~ Length)", "Folding Energy")
        plotPostNet(resOut, colOut, comparisons, ylabel = ylabel, plotType = plotType)
        dev.off()
      }
    }
    return(feOut)
  } else {
    stop("The selection of the folding energy source file is not valid. Please check the inputs and refer to the help manual for details.")
  }
}

runFE <- function(energyIn,
                  residFE = FALSE,
                  ptn) {
  #
  colnames(energyIn) <- c("id", "fold_energy", "length")
  energyIn$geneID <- ptn_geneID(ptn, region = "CDS")[match(energyIn$id, ptn_id(ptn, region = "CDS"))]
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
