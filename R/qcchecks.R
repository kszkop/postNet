checkRegion <- function(region, convertToUppercase = TRUE) {
  valid_regions <- c('UTR3', 'CDS', 'UTR5', 'CCDS')
  
  if (is.null(region) || !is.character(region) || length(region) == 0) {
    stop("'region' must be a non-empty character vector with valid values, to choose from 'UTR3', 'CDS', 'UTR5', 'CCDS'")
  }
  
  if (convertToUppercase) {
    region <- toupper(region)
  }
  
  if (!all(region %in% valid_regions)) {
    stop("'region' must contain valid values: 'UTR3', 'CDS', 'UTR5','CCDS'.")
  }
}

check_adjObj <- function(adjObj) {
  if (!is.list(adjObj)) {
    stop("'adjObj' is not a list")
  }
  valid_names <- c('UTR5', 'UTR3')
  if (!all(names(adjObj) %in% valid_names)) {
    stop("names of the entries in the list should be only 'UTR3' or 'UTR5'")
  }
  for (name in names(adjObj)) {
    entry <- adjObj[[name]]
    if (!is.character(entry) || !all(nchar(entry) > 0)) {
      stop("the entries in list should be a character vector with nucleotide sequences")
    }
    if (!isDNAsequence(entry[1])) {
      stop("It looks like the sequences are not DNA sequences")
    }
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
    stop("'annot' should be a data frame with columns: 'id', 'geneID', 'UTR5_seq', 'CDS_seq', 'UTR3_seq'")
  }
  if (!all(expectedCols %in% colnames(annot))) {
    stop("The following columns are missing in 'annot': ", paste(expectedCols[!expectedCols %in% colnames(annot)], collapse = ", "))
  }
}

checkAnnotCod <- function(annot, expectedCols = c("id", "geneID", "CDS_seq")) {
  if (is.null(annot) || !is.data.frame(annot)) {
    stop("'customFileCod' should be a file in format data frame with columns: 'id', 'geneID', 'CDS_seq'")
  }
  if (!all(expectedCols %in% colnames(annot))) {
    stop("The following columns are missing in 'customFileCod': ", paste(expectedCols[!expectedCols %in% colnames(annot)], collapse = ", "))
  }
}

# Check if 'ads' is a valid 'Anota2seqDataSet' object
checkAds <- function(obj) {
  if (!inherits(obj, "Anota2seqDataSet")) {
    return(FALSE)
  }
  return(TRUE)
}

checkUtils <- function(obj) {
  if (!inherits(obj, "anota2seqUtilsData")) {
    return(FALSE)
  }
  return(TRUE)
}

checkComparisons <- function(obj) {
  if (!is.list(obj)) {
    return(FALSE)
  }
  all(sapply(obj, function(x) is.numeric(x) && length(x) == 2))
}

is_logical <- function(x) {
  is.logical(x) && length(x) == 1
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

is_by_3 <- function(seqs) {
  all(sapply(seqs, function(x) length(seqinr::s2c(x)) %% 3 == 0))
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

isUnit <- function(unit) {
  unitOut <- tolower(unit)
  valid_values <- c("count", "freq")
  if (is.character(unit) && length(unit) == 1 && unit %in% valid_values) {
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

check_codonIn<- function(codonIn) {
  #is it a list
  if (!is.list(codonIn)) {
    stop('"codonIn" is not a list, and so probably not an output of codonUsage function')
  }
  #'codonAll' exists and is a list
  if (!"codonAll" %in% names(codonIn) || !is.list(codonIn$codonAll)) {
    stop('"codonIn" does not contain codonAll element, and so probably not an output of codonUsage function')
  }
  
  #'codonAll' list contains required elements
  required_elements <- c("geneID", "codon", "AA", "codonCount", "codonFreq", "AACountPerGene")
  if (!all(required_elements %in% names(codonIn$codonAll))) {
    stop('"codonIn$codonAll" element does not contain all required elements, and so probably "codonIn" is not an output of codonUsage function')
  }
  
  # elements are of correct types
  if (!all(sapply(codonIn$codonAll[grep("geneID|codon|AA", names(codonIn$codonAll))], is.character)) ||
      !all(sapply(codonIn$codonAll[grep("codonCount|codonFreq|AACountPerGene", names(codonIn$codonAll))], is.double))) {
    stop('"codonIn$codonAll" elements are not of correct types, and so probably "codonIn" is not an output of codonUsage function')
  }
}

check_codons <- function(featSel) {
  # Define a vector of all possible codons
  all_codons <- c("AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC",
                  "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT",
                  "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC",
                  "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT",
                  "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAA", "TAC",
                  "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT",
                  "TTA", "TTC", "TTG", "TTT")
  
  # Check if all elements of featSel are either single codons or combinations of codons
  all(sapply(featSel, function(x) {
    if (nchar(x) %% 3 == 0 && all(strsplit(x, "")[[1]] %in% c("A", "C", "G", "T"))) {
      all(substring(x, seq(1, nchar(x), by = 3), seq(3, nchar(x), by = 3)) %in% all_codons)
    } else {
      x %in% all_codons
    }
  }))
}

isValidSlope <- function(slope) {
  return(!is.null(slope) && is.numeric(slope) && !is.na(slope))
}

checkSlopes <- function(minSlope, maxSlope) {
  if (!isValidSlope(minSlope)) {
    stop("minSlope is either NULL, not numeric, or NA.")
  }
  if (!isValidSlope(maxSlope)) {
    stop("maxSlope is either NULL, not numeric, or NA.")
  }
  return(TRUE)
}


checkFileColumns <- function(filePath) {
  if (is.null(filePath)) {
    stop("File path is NULL.")
  }

  if (!file.exists(filePath)) {
    stop("File does not exist.")
  }

  fileData <- read.delim(filePath)
  
  requiredColumns <- c('Gene.Tax.ID','weighted.context...score','Site.Type','Gene.Symbol','miRNA')
  missingColumns <- setdiff(requiredColumns, colnames(fileData))
  if (length(missingColumns) > 0) {
    stop(paste("The following required columns are missing:", paste(missingColumns, collapse = ", ")))
  }
  
  noS <- length(unique(fileData$Gene.Tax.ID))
  if(noS>1){
    stop('Please subset the file for only desired specie')
  }
  return(fileData)
}

checkCollection <- function(collection) {
  if (is.null(collection)) {
    stop('Please provide collection or geneSet')
  }
  
  collections <- c('c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'h')

  for (choice in collections) {
    if (!choice %in% collections) {
      stop('Please provide valid collections (to choose from c1,c2,c3,c4,c5,c6,h)')
    }
  }
}

checkGeneList <- function(obj) {
  if (!is.list(obj)) {
    stop("The input is not a list.")
  }
  
  if (length(obj) == 0) {
    stop("The list is empty.")
  }
  
  if (all(names(obj) == "")) {
    stop("The list is not a named list.")
  }
}

checkDirection <- function(direction) {
  if (is.null(direction)) {
    stop("The direction cannot be NULL.")
  }
  if (!direction %in% c("greater", "less")) {
    stop('The direction must be either "greater" or "less".')
  }
}

checkCategory <- function(category) {
  if (is.null(category)) {
    stop("The category cannot be NULL.")
  }
  selCat <- c("BP", "CC", "MF", "KEGG")
  if (!all(category %in% selCat)) {
    stop('The "category" must be a combination of "BP", "CC", "MF", and "KEGG" ')
  }
}

check_size <- function(size) {
  if (is.null(size) || !(size == "Count" || size == "geneRatio")) {
    stop("The 'size' must be not null and only can be 'Count' or 'geneRatio'" )
  }
}
