check_region <- function(region, convertToUppercase = TRUE) {
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
    stop("The names of the entries in the adjObj list should be only 'UTR3' or 'UTR5'.")
  }
  for (name in names(adjObj)) {
    entry <- adjObj[[name]]
    if (!is.character(entry) || !all(nchar(entry) > 0)) {
      stop("The entries in the adjObj list should be character vectors with DNA nucleotide sequences.")
    }
  }
    check_DNAsequence(entry) 
}

check_selection <- function(selection, convertToLowercase = TRUE) {
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

check_plotType <- function(plotType, convertToLowercase = TRUE) {
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

check_ptn <- function(obj) {
  if (!inherits(obj, "postNetData")) {
    stop("ptn is not a valid 'postNetData' object.")
  }
}

check_ads <- function(obj) {
  if (!inherits(obj, "Anota2seqDataSet")) {
    stop("ads is not a valid 'Anota2seqDataSet' object.")
  }
}

check_comparisons <- function(obj) {
  if (!is.list(obj)) {
    return(FALSE)
  }
  all(sapply(obj, function(x) is.numeric(x) && length(x) == 2))
}

is_valid_named_list <- function(obj) {
  # Check if the object is NULL
  if (is.null(obj)) {
    return(FALSE)
  }
  
  # Check if the object is a list
  if (!is.list(obj)) {
    return(FALSE)
  }
  
  # Check if the list has names
  if (is.null(names(obj)) || any(names(obj) == "")) {
    return(FALSE)
  }
  
  # Check if each element in the list is a numeric vector
  for (item in obj) {
    if (!is.numeric(item) || !is.vector(item)) {
      return(FALSE)
    }
  }
  return(TRUE)
}

is_numeric_vector <- function(obj) {
  return(is.numeric(obj) && is.vector(obj))
}

check_logical <- function(x) {
  is.logical(x) && !is.na(x) && length(x) == 1
}

check_number <- function(x) {
  is.numeric(x) && !is.na(x) && length(x) == 1
}


is_named_list_of_named_numeric_vectors <- function(x) {
  # Check if the input is a list
  if (!is.list(x)) {
    return(FALSE)
  }
  
  # Check if the list is named (all elements have names)
  if (is.null(names(x)) || any(names(x) == "")) {
    return(FALSE)
  }
  
  # Check if each element is a named numeric vector
  for (element in x) {
    if (!is.numeric(element) || is.null(names(element)) || any(names(element) == "")) {
      return(FALSE)
    }
  }
  
  # If all checks pass, return TRUE
  return(TRUE)
}

# Function to validate the source input
check_source <- function(source) {
  valid_sources <- c("create", "createFromSourceFiles", "load", "custom", "createFromFasta")
  if (!(source %in% valid_sources)) {
    stop("Invalid source. Please provide a valid source option.")
  }
}

# Function to validate the source input
checkSourceFE <- function(sourceFE) {
  valid_sourcesFE <- c("load", "custom")
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

check_DNAsequence <- function(contentIn) {
  if (!is.character(contentIn)) {
    stop("The adjObj must be a named list of character vectors with DNA sequences.")
  }
  # Define a regular expression pattern for valid DNA characters
  pattern <- "^[ACGTacgt]+$"
  
  # Use grepl to check if contentIn matches the pattern
  if (!all(grepl(pattern, contentIn))) {
    stop("The entries provided in adjObj do not appear to all be DNA sequences. Please check that the sequences are corect.")
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
check_input <- function(source, customFile, rna_gbff_file, rna_fa_file, genomic_gff_file, posFile, fastaFile) {
  check_source(source)
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
  } else if (source == "createFromFasta") {
    if (is.null(posFile)) {
      stop("Please provide a posFile in the format: id, UTR5_len, CDS_stop, Total_len.")
    }
    if (is.null(fastaFile)) {
      stop("Please provide a fastaFile.")
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

#checkcodSource <- function(codSource) {
#  valid_codSource <- c("sequence", "riboseq")
#  if (is.null(codSource)) {
#    stop("'codSource' cannot be null.")
#  } else {
#    # Convert to lowercase
#    codSource <- tolower(codSource)
#    
#    if (!(codSource %in% valid_codSource)) {
#      stop("Invalid codSource. Allowed 'codSource' are 'sequence' or 'riboseq'.")
#    } 
#  }
#}

check_codonIn<- function(codonIn) {
  #
  required_elements <- c("geneID", "codon", "AA", "count", "frequency", "AACountPerGene", "relative_frequency")
  if (!all(required_elements %in% colnames(codonIn)) || is.null(codonIn)) {
    stop('"codonsAll" element does not contain all required elements, and so probably it is not an output of codonUsage function')
  }
}


check_codons <- function(featsel) {
  # Check if input is a named list
  if (!is.list(featsel) || is.null(names(featsel)) || any(names(featsel) == "")) {
    stop("Input must be a named list.")
  }
  
  # Define a vector of all possible codons
  all_codons <- c("AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC",
                  "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT",
                  "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC",
                  "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT",
                  "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAA", "TAC",
                  "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT",
                  "TTA", "TTC", "TTG", "TTT")
  
  # Check if all elements of featSel (values in the named list) are valid codons
  all(sapply(featsel, function(codons) {
    # Ensure that each codon in the vector is valid
    all(codons %in% all_codons)
  }))
}

check_AA <- function(featSel) {
  if (!is.list(featSel) || is.null(names(featSel)) || any(names(featSel) == "")) {
    stop("Input must be a named list.")
  }
  #
  single_to_three <- c(A = "Ala", R = "Arg", N = "Asn", D = "Asp", C = "Cys", 
                       Q = "Gln", E = "Glu", G = "Gly", H = "His", I = "Ile", 
                       L = "Leu", K = "Lys", M = "Met", F = "Phe", P = "Pro", 
                       S = "Ser", T = "Thr", W = "Trp", Y = "Tyr", V = "Val")
  
  # 
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
    stop('Please subset the file the include only the desired species')
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

check_geneList <- function(obj) {
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

check_direction <- function(direction) {
  if (is.null(direction)) {
    stop("The direction cannot be NULL.")
  }
  if(length(direction) != 1){
    stop("Please provide only one: greater or less")
  }
  if (!direction %in% c("greater", "less")) {
    stop('The direction must be either "greater" or "less".')
  }
}

check_category <- function(category) {
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

check_analysis_type <- function(analysis_type){
  if (is.null(analysis_type)) {
    stop("Please provide 'analysis_type' argument. It should be'lm' for linear model or 'rf' for random forest")
  }
  if (!analysis_type %in% c("lm", "rf")) {
    stop("'analysis_type' can be only 'lm' for linear model or 'rf' for random forest")
  }
}

is_valid_NetModelSel <- function(NetModelSel) {
  if (is.null(NetModelSel)) {
    return(FALSE)
  }
  if (NetModelSel %in% c("omnibus", "adjusted", "univariate")) {
    return(TRUE)
  }
  return(FALSE)
}

check_model <- function(model, analysis_type) {
  check_analysis_type(analysis_type)
  if (is.null(model)) {
    stop("Please provide correct 'model' for a analysis type")
  }
  if(analysis_type == "lm"){
    if (!model %in% c("univariateModel", "stepwiseModel", "finalModel")) {
      stop("Please provide correct 'model'. For 'lm', choose one of these: 'univariateModel', 'stepwiseModel', 'finalModel'")
    }
  }
  if(analysis_type == "rf"){
    if(!model %in% c("preModel", "borutaModel", "finalModel")) {
      stop("Please provide correct 'model'. For 'rf', choose one of these: 'preModel', 'borutaModel', 'finalModel'")
    }
  }
}

check_features <- function(features) {
  #
  if (is.null(features)) {
    stop("Error: 'features' cannot be NULL.")
  }
  
  # 
  if (!is.list(features)) {
    stop("Error: 'features' must be a list.")
  }
  
  #
  if (is.null(names(features)) || any(names(features) == "")) {
    stop("Error: 'features' must be a named list.")
  }
  
  # 
  if (!all(sapply(features, is.vector))) {
    stop("Error: Each element in 'features' must be a vector.")
  }
  if(length(features) < 2){
    stop(" Please provide at least two features")
  }
}

checklmfeatGroup <- function(lmfeatGroup, numfeatures) {
  if (!is.null(lmfeatGroup)) {
    #
    if (!is.vector(lmfeatGroup)) {
      stop("Error: 'lmfeatGroup' must be a vector when provided.")
    }
    
    # 
    if (length(lmfeatGroup) != numfeatures) {
      stop("Error: Length of 'lmfeatGroup' must match the number of 'features' when 'lmfeatGroup' is not NULL.")
    }
  }
}

checklmfeatGroupColour <- function(lmfeatGroupColour, lmfeatGroup) {
  #
  if (!is.null(lmfeatGroupColour)) {
    #
    if (is.null(lmfeatGroup)) {
      stop("Error: 'lmfeatGroup' cannot be NULL when 'lmfeatGroupColour' is provided.")
    }
    # 
    unique_lmfeatGroup <- length(unique(lmfeatGroup))
    
    #
    if (length(lmfeatGroupColour) != unique_lmfeatGroup) {
      stop("Error: Length of 'lmfeatGroupColour' must match the number of unique values in 'lmfeatGroup'.")
    }
  }
}

check_shiftUnit <- function(unit) {
  if (is.null(unit)) {
    return(FALSE)
  }
  if (unit == "FDR") {
    return(TRUE)
  }
  if (grepl("^p[1-9][0-9]?$", unit)) {
    return(TRUE)
  }
  return(FALSE)
}

check_featSel <- function(featSel, features) {
  if (!is.null(featSel) && is.character(featSel) && length(featSel) >= 2) {
    if (all(featSel %in% colnames(features))) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    return(FALSE)
  }
}

check_predFeat <- function(predFeat) {
  is_list<- is.list(predFeat)
  has_rownames <- !is.null(rownames(predFeat)) && all(rownames(predFeat) != "")
  is_numeric <- all(sapply(predFeat, is.numeric))
  has_no_nas <- all(complete.cases(predFeat))
  
  return(is_dataframe && has_rownames && is_numeric && has_no_nas)
}
