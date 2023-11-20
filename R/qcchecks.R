checkRegion <- function(region, convertToUppercase = TRUE) {
  valid_regions <- c('UTR3', 'CDS', 'UTR5')
  
  if (is.null(region) || !is.character(region) || length(region) == 0) {
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
  
  if (is.null(selection) || !is.character(selection) || length(selection) == 0) {
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
  
  if (is.null(plotType) || !is.character(plotType) || length(plotType) == 0) {
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
  if (is.null(annot) || !is.data.frame(annot)) {
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

is_number <- function(x) {
  is.numeric(x) && !is.na(x)
}

# Function to validate the source input
checkSource <- function(source) {
  valid_sources <- c("create", "createFromSourceFiles", "load", "custom", "createFromFiles")
  if (!(source %in% valid_sources)) {
    stop("Invalid source. Please provide a valid source option.")
  }
}

# Function to validate the source input
checkSourceFE <- function(sourceFE) {
  valid_sourcesFE <- c("create", "load", "custom")
  if (!(sourceFE %in% valid_sourcesFE)) {
    stop("Invalid sourceFE. Please provide a valid sourceFE option.")
  }
}

is_valid_species <- function(species) {
  species <- tolower(species)
  if (!is.null(species) && (species == "human" || species == "mouse")) {
    return(TRUE)
  }
  return(FALSE)
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

is_valid_seq_type <- function(seqType) {
  valid_types <- c('dna', 'rna', 'protein')
  if (!is.null(seqType) && tolower(seqType) %in% valid_types) {
    return(TRUE)
  }
  return(FALSE)
}

isStartCodon <- function(startCodon) {
  if (is.character(startCodon) && length(startCodon) == 1) {
    start_codon <- toupper(seqinr::s2c(startCodon))
    valid_bases <- c("A", "C", "G", "T")
    if (all(start_codon %in% valid_bases) && length(start_codon) == 3) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    return(FALSE)
  }
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

is_motifs <- function(motifsIn) {
  if (!is.null(motifsIn) && is.character(motifsIn) && length(motifsIn) > 0) {
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

is_annotType <- function(annotType) {
  valid_types <- c('refseq', 'ccds', 'custom')
  
  if (!is.null(annotType) && tolower(annotType) %in% valid_types) {
    return(TRUE)
  }
  return(FALSE)
}

is_valid_sourceSeq <- function(sourceSeq) {
  if (is.null(sourceSeq)) {
    return(FALSE)
  }
  sourceSeq <- tolower(sourceSeq)
  if (sourceSeq %in% c("load", "create")) {
    return(TRUE)
  }
  return(FALSE)
}

is_valid_analysis <- function(analysis) {
  if (is.null(analysis)) {
    return(FALSE)
  }
  if (analysis %in% c("codon", "AA")) {
    return(TRUE)
  }
  return(FALSE)
}

checkDirectory <- function(path) {
  if (!dir.exists(path)) {
    stop("Directory does not exist.")
  }
}

checkcodSource <- function(codSource) {
  valid_codSource <- c("sequence", "riboseq")
  if (is.null(codSource)) {
    stop("'codSource' cannot be null.")
  } else {
    # Convert to lowercase
    codSource <- tolower(codSource)
    
    if (!(codSource %in% valid_codSource)) {
      stop("Invalid codSource. Allowed 'codSource' are 'sequence' or 'riboseq'.")
    } 
  }
}
