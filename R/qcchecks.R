checkParameters <- function(annot, 
                            ads, 
                            regulation, 
                            contrast, 
                            geneList, 
                            geneListcolours, 
                            customBg, 
                            selection, 
                            region=NULL, 
                            comparisons, 
                            plotOut,
                            plotType = NULL,
                            contentIn= NULL,
                            subregion= NULL,
                            subregionSel= NULL,
                            startCodon = NULL,
                            KozakContext = NULL,
                            onlyUTR5 = NULL,
                            unitOut = NULL
                            ){
  ####
  if(!is.null(ads) && !is.null(geneList)){
    stop("please provide anota2seq object or genelist, not both.")
  }
  # Check for other input conditions as needed
  if(!is.null(annot)){
    checkAnnot(annot)
  }
  if(!is.null(region)){
    checkRegion(region)
  }
  if(!is.null(selection)){
    checkSelection(selection)
  }
  
  if(!checkLogicalArgument(plotOut)){
    stop("'plotOut' can only be only be logical: TRUE of FALSE ")
  } 
  if(isTRUE(plotOut)){
    if(!is.null(plotType)){
      plotType <- checkPlotType(plotType)
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
  if(!is.null(contentIn) && !isDNAsequence(contentIn)){
    stop("'contentIn' must be a character vector with DNA sequences")
  }
  if(!is.null(subregion) && (!is.numeric(subregion) || !length(subregion)==1)){
    stop("'subregion' must be a numeric and just number")
  }
  if (!is.null(subregionSel) && is.character(subregionSel) && length(subregionSel) == 1) {
    if (!subregionSel %in% c("select", "exclude")) {
      stop("'subregionSel' must be a character and only 'select' or 'exclude'")
    }
  } 
  if(!is.null(startCodon) && !isStartCodon(startCodon)){
    stop("'startCodon' must be a character vector of length one, and contain only 3 nucleotide sequence, ex. 'ATG'")
  }
  if(!is.null(KozakContext) && !isKozakContext(KozakContext)){
    stop("'KozakContext' must be one from these: 'strong','adequate1','adequate2','weak','any'")
  }
  if(!is.null(onlyUTR5) && !checkLogicalArgument(onlyUTR5)){
    stop("'onlyUTR5' can only be only be logical: TRUE of FALSE ")
  }
  if(!is.null(unitOut) && !isUnitOut(unitOut)){
    stop("'unitOut' must be one from these: 'numeric' or 'position'")
  }
}

checkRegion <- function(region, convertToUppercase = TRUE) {
  valid_regions <- c('UTR3', 'CDS', 'UTR5')
  
  if (!is.character(region) || length(region) == 0) {
    stop("'region' must be a non-empty character vector with valid values, to choose from 'UTR3', 'CDS', 'UTR5'.")
  }
  
  if (convertToUppercase) {
    region <- toupper(region)
  }
  
  if (!all(region %in% valid_regions)) {
    stop("'region' must contain valid values: 'UTR3', 'CDS', 'UTR5'.")
  }
}

checkSelection <- function(selection, convertToLowercase = TRUE) {
  valid_selection <- c('random', 'longest', 'shortest')
  
  if (!is.character(selection) || length(selection) == 0) {
    stop("'selection' must be one of: 'random', 'longest', or 'shortest'.")
  }
  
  if (convertToLowercase) {
    selection <- tolower(selection)
  }
  
  if (!selection %in% tolower(valid_selection)) {
    stop("'selection' must be one of: 'random', 'longest', or 'shortest'.")
  }
}

checkPlotType <- function(plotType, convertToLowercase = TRUE) {
  valid_plottypes <- c('boxplot', 'violin', 'ecdf')
  
  if (!is.character(plotType) || length(plotType) == 0) {
    stop("'plotType' must be one of: 'boxplot', 'violin', or 'ecdf'.")
  }
  
  if (convertToLowercase) {
    plotType <- tolower(plotType)
  }
  
  if (!plotType %in% tolower(valid_plottypes)) {
    stop("'plotType' must be one of: 'boxplot', 'violin', or 'ecdf'.")
  }
}

checkAnnot <- function(annot, expectedCols = c("id", "geneID", "UTR5_seq", "CDS_seq", "UTR3_seq")) {
  if (!is.data.frame(annot)) {
    stop("'annot' should be a data frame.")
  }
  
  if (!all(expectedCols %in% colnames(annot))) {
    stop("The following columns are missing in 'annot': ", paste(expectedCols[!expectedCols %in% colnames(annot)], collapse = ", "))
  }
}

# Check if 'ads' is a valid 'Anota2seqDataSet' object
checkAds <- function(obj) {
  if (!inherits(obj, "Anota2seqDataSet")) {
    return(FALSE)
  }
  return(TRUE)
}

checkGeneList <- function(obj) {
  if (is.list(obj) && length(obj) > 0 && any(names(obj) != "")) {
    return(TRUE)
  }
  return(FALSE)
}

checkComparisons <- function(obj) {
  if (!is.list(obj)) {
    return(FALSE)
  }
  
  all(sapply(obj, function(x) is.numeric(x) && length(x) == 2))
}

checkLogicalArgument <- function(arg) {
  if (is.logical(arg) && length(arg) == 1) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Function to validate the source input
checkSource <- function(source) {
  valid_sources <- c("create", "createFromSourceFiles", "load", "custom", "createFromFiles")
  if (!(source %in% valid_sources)) {
    stop("Invalid source. Please provide a valid source option.")
  }
}

# Function to validate the species input
checkSpecies <- function(source, species) {
  if (source %in% c("create", "load") && is.null(species)) {
    stop("Please specify a species (e.g., 'human' or 'mouse').")
  }
  if (source %in% c("create", "load") && !(species %in% c("human", "mouse"))) {
    stop("This option is only available for species: human and mouse at the moment. Please use option createFromFile.")
  }
}

isDNAsequence <- function(contentIn) {
  if (!is.character(contentIn)) {
    return(FALSE)
  }
  
  # Define a regular expression pattern for valid DNA characters
  pattern <- "^[ACGTacgt]+$"
  
  # Use grepl to check if contentIn matches the pattern
  if (all(grepl(pattern, contentIn))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

isStartCodon <- function(startCodon) {
  if (is.character(startCodon) && nchar(startCodon) == 1) {
    start_codon <- toupper(seqinr::s2c(startCodon))
    valid_bases <- c("A", "C", "G", "T")
    if (all(start_codon %in% valid_bases) &&
        length(start_codon) == 3) {
      return(TRUE)
    }
  }
  return(FALSE)
}

isKozakContext <- function(KozakContext) {
  KozakContext <- tolower(KozakContext)
  valid_values <- c("strong", "adequate1", "adequate2", "weak", "any")
  if (is.character(KozakContext) && length(KozakContext) == 1 &&
      KozakContext %in% valid_values) {
    return(TRUE)
  }
  return(FALSE)
}

isUnitOut <- function(unitOut) {
  unitOut <- tolower(unitOut)
  valid_values <- c("number", "position")
  if (is.character(unitOut) && length(unitOut) == 1 && unitOut %in% valid_values) {
    return(TRUE)
  }
  return(FALSE)
}

# Function to validate specific input parameters
checkInput <- function(source, customFile, rna_gbff_file, rna_fa_file, genomic_gff_file, posFile) {
  if (source == "createFromSourceFiles") {
    if (is.null(rna_gbff_file)) {
      stop("Please provide an rna_gbff_file.")
    }
    if (is.null(rna_fa_file)) {
      stop("Please provide an rna_fa_file.")
    }
    if (is.null(genomic_gff_file)) {
      stop("Please provide a genomic_gff_file.")
    }
  } else if (source == "custom") {
    if (is.null(customFile)) {
      stop("Please provide a customFile.")
    }
  } else if (source == "createFromFiles") {
    if (is.null(posFile)) {
      stop("Please provide a posFile in the format: id, UTR5_len, CDS_stop, Total_len.")
    }
  }
}