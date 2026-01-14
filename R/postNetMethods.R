setGeneric(
  "ptn_sequences",
  function(ptn, region) standardGeneric("ptn_sequences")
)
setMethod(
  "ptn_sequences", "postNetData",
  function(ptn, region) {
    check_region(region)
    check_ptn(ptn)
    if (length(region) > 1) {
      stop("The input for 'region' can only be one of: 'UTR3', 'CDS', 'UTR5', or alternatively 'CCDS' if codon analysis was performed using the CCDS annotation.")
    }
    tmpReg <- slot(ptn@annot, region)
    seqOut <- tmpReg@sequences
    return(seqOut)
  }
)

setGeneric(
  "ptn_id",
  function(ptn, region) standardGeneric("ptn_id")
)
setMethod(
  "ptn_id", "postNetData",
  function(ptn, region) {
    check_region(region)
    check_ptn(ptn)
    if (length(region) > 1) {
      stop("The input for 'region' can only be one of: 'UTR3', 'CDS', 'UTR5', or alternatively 'CCDS' if codon analysis was performed using the CCDS annotation.")
    }
    tmpReg <- slot(ptn@annot, region)
    idOut <- tmpReg@id
    return(idOut)
  }
)

setGeneric(
  "ptn_geneID",
  function(ptn, region) standardGeneric("ptn_geneID")
)
setMethod(
  "ptn_geneID", "postNetData",
  function(ptn, region) {
    check_region(region)
    check_ptn(ptn)
    if (length(region) > 1) {
      stop("The input for 'region' can only be one of: 'UTR3', 'CDS', 'UTR5', or alternatively 'CCDS' if codon analysis was performed using the CCDS annotation.")
    }
    tmpReg <- slot(ptn@annot, region)
    geneIDOut <- tmpReg@geneID
    return(geneIDOut)
  }
)

setGeneric(
  "ptn_dataIn",
  function(ptn) standardGeneric("ptn_dataIn")
)
setMethod(
  "ptn_dataIn", "postNetData",
  function(ptn) {
    check_ptn(ptn)
    ptn@dataIn
  }
)

setGeneric(
  "ptn_geneList",
  function(ptn) standardGeneric("ptn_geneList")
)
setMethod(
  "ptn_geneList", "postNetData",
  function(ptn) {
    check_ptn(ptn)
    ptn@dataIn@geneList
  }
)

setGeneric(
  "ptn_background",
  function(ptn) standardGeneric("ptn_background")
)
setMethod(
  "ptn_background", "postNetData",
  function(ptn) {
    check_ptn(ptn)
    ptn@dataIn@background
  }
)

setGeneric(
  "ptn_effect",
  function(ptn) standardGeneric("ptn_effect")
)
setMethod(
  "ptn_effect", "postNetData",
  function(ptn) {
    check_ptn(ptn)
    ptn@dataIn@effect
  }
)

setGeneric(
  "ptn_colours",
  function(ptn) standardGeneric("ptn_colours")
)
setMethod(
  "ptn_colours", "postNetData",
  function(ptn) {
    check_ptn(ptn)
    ptn@dataIn@colours
  }
)

setGeneric(
  "ptn_species",
  function(ptn) standardGeneric("ptn_species")
)
setMethod(
  "ptn_species", "postNetData",
  function(ptn) {
    check_ptn(ptn)
    ptn@species
  }
)

setGeneric(
  "ptn_version",
  function(ptn) standardGeneric("ptn_version")
)
setMethod(
  "ptn_version", "postNetData",
  function(ptn) {
    check_ptn(ptn)
    ptn@version
  }
)

setGeneric(
  "ptn_selection",
  function(ptn) standardGeneric("ptn_selection")
)
setMethod(
  "ptn_selection", "postNetData",
  function(ptn) {
    check_ptn(ptn)
    ptn@selection
  }
)

setGeneric(
  "ptn_motifSelection",
  function(ptn, region) standardGeneric("ptn_motifSelection")
)
setMethod(
  "ptn_motifSelection", "postNetData",
  function(ptn, region) {
    check_ptn(ptn)
    check_region(region)
    if (length(region) > 1) {
      stop("The input for 'region' can only be one of: 'UTR3', 'CDS', 'UTR5'.")
    }
    tmpReg <- slot(ptn@analysis@motifs, region)
    motifsOut <- tmpReg$motifSelection
    return(motifsOut)
  }
)

setGeneric(
  "ptn_motifGeneList",
  function(ptn, region, geneList) standardGeneric("ptn_motifGeneList")
)
setMethod(
  "ptn_motifGeneList", "postNetData",
  function(ptn, region, geneList) {
    check_ptn(ptn)
    check_region(region)
    if (length(region) > 1) {
      stop("The input for 'region' can only be one of: 'UTR3', 'CDS', 'UTR5'.")
    }
    if (!any(geneList %in% names(ptn_geneList(ptn)))) {
      stop("The list of regulated genes 'geneList' is not stored in the postNetData object.")
    }
    tmpReg <- slot(ptn@analysis@motifs, region)
    motifsAnalysisOut <- tmpReg[[geneList]]
    return(motifsAnalysisOut)
  }
)

setGeneric(
  "ptn_codonAnalysis",
  function(ptn) standardGeneric("ptn_codonAnalysis")
)
setMethod(
  "ptn_codonAnalysis", "postNetData",
  function(ptn) {
    check_ptn(ptn)

    tmpOut <- ptn@analysis@codons@codonAnalysis
    out <- s4_to_dataframe(tmpOut)
    return(out)
  }
)

setGeneric(
  "ptn_codonSelection",
  function(ptn, comparison) standardGeneric("ptn_codonSelection")
)
setMethod(
  "ptn_codonSelection", "postNetData",
  function(ptn, comparison) {
    check_ptn(ptn)

    out <- ptn@analysis@codons@codonSelection[[comparison]]
    return(out)
  }
)

setGeneric(
  "ptn_features",
  function(ptn) standardGeneric("ptn_features")
)
setMethod(
  "ptn_features", "postNetData",
  function(ptn) {
    check_ptn(ptn)

    out <- ptn@features
    return(out)
  }
)


ptn_miRNA_analysis <- function(ptn,
                               direction,
                               threshold) {
  #
  check_direction(tolower(direction))
  check_ptn(ptn)
  if (!check_number(threshold)) {
    stop(paste("Please provide a single numeric value for ", threshold, sep = ""))
  }

  if (is.null(slot(ptn@analysis, "miRNA"))) {
    stop("Please run the miRNAanalysis() function first.")
  } else {
    miRNAres <- ptn@analysis@miRNA@miRNA_analysis
  }

  if (tolower(direction) == "greater") {
    resOut <- miRNAres$greater
  } else if (tolower(direction) == "less") {
    resOut <- miRNAres$less
  }
  #
  resOut <- resOut[which(resOut[, 4] < threshold), ]
  #
  if (nrow(resOut) > 0) {
    resOut <- resOut[, c(1, 2, 5, 3, 4)]
    resOut <- data.frame(id = row.names(resOut), resOut, row.names = NULL)
  } else {
    message("There are no enriched miRNAs to output.")
  }
  return(resOut)
}


ptn_miRNA_to_gene <- function(ptn,
                              miRNAs) {
  #
  check_ptn(ptn)

  if (is.null(slot(ptn@analysis, "miRNA"))) {
    stop("Please run the miRNAanalysis() function first.")
  } else {
    miRNATmp <- ptn@analysis@miRNA@miRNA_to_gene
  }
  miRNAsOut <- miRNATmp[which(names(miRNATmp) %in% miRNAs)]

  return(miRNAsOut)
}

###
ptn_GO <- function(ptn,
                   category,
                   geneList,
                   threshold) {
  #
  check_ptn(ptn)
  check_category(category)
  if (length(category) != 1) {
    stop("Please provide only one GO category.")
  }
  if (!check_number(threshold)) {
    stop(paste("Please provide a single numeric value for ", threshold, sep = ""))
  }
  #
  if (!any(geneList %in% names(ptn_geneList(ptn)))) {
    stop("None of the regulated genes included in 'geneList' are included in the postNetData object.")
  }
  #
  if (is.null(slot(ptn@analysis, "GO"))) {
    stop("Please run the goAnalysis() function first.")
  } else {
    GOres <- slot(ptn@analysis@GO, category)
  }
  GOresOut <- GOres[[which(geneList == names(GOres))]]@result
  GOresOut <- GOresOut[which(GOresOut$p.adjust < threshold), ]

  if (nrow(GOresOut) > 0) {
    GOresOut <- data.frame(GOresOut, row.names = NULL)
  } else {
    message("There are no enriched GO categories to output.")
  }
  return(GOresOut)
}

###
ptn_GSEA <- function(ptn,
                     threshold = NULL) {
  check_ptn(ptn)

  if (is.null(slot(ptn@analysis, "GSEA"))) {
    stop("Please run the gseaAnalysis() function first.")
  } else {
    gseaOut <- slot(ptn@analysis, "GSEA")
  }
  if (!is.null(threshold)) {
    if (!check_number(threshold)) {
      stop(paste("Please provide a single numeric value for ", threshold, sep = ""))
    }
    gseaOut <- gseaOut[which(gseaOut[, 8] < threshold), ]
  }
  if (nrow(gseaOut) == 0) {
    message("There are no enriched gene sets to output from GSEA.")
  }
  return(gseaOut)
}

ptn_GAGE <- function(ptn,
                     category,
                     direction,
                     threshold) {
  #
  check_direction(tolower(direction))
  check_ptn(ptn)
  check_category(category)
  if (!check_number(threshold)) {
    stop(paste("Please provide a single numeric value for ", threshold, sep = ""))
  }

  if (length(category) != 1) {
    stop("Please provide only one category.")
  }
  if (is.null(slot(ptn@analysis, "GAGE"))) {
    stop("Please run the gageAnalysis() function first.")
  } else {
    GAGEres <- slot(ptn@analysis@GAGE, category)
  }

  if (tolower(direction) == "greater") {
    resOut <- GAGEres$greater
  } else if (tolower(direction) == "less") {
    resOut <- GAGEres$less
  }
  #
  resOut <- resOut[which(resOut[, 6] < threshold), ]
  #
  if (nrow(resOut) == 0) {
    message("There are no enriched terms to output from GAGE.")
  } else {
    return(resOut)
  }
}

# ptn_check_comparisons <- function(ptn, analysis_type){
# check_ptn(ptn)
#  check_analysis_type(analysis_type)
#
#  tmpIn <- ptn@analysis@featureIntegration[[analysis_type]]
#  if(is.null(tmpIn)){
#    stop(paste('Please run ', analysis_type, ' analysis first', sep=''))
#  } else {
#    print(names(tmpIn))
#  }
# }

ptn_check_models <- function(ptn, analysis_type) {
  check_ptn(ptn)
  check_analysis_type(analysis_type)

  tmpIn <- ptn@analysis@featureIntegration[[analysis_type]]
  if (is.null(tmpIn)) {
    stop(paste("Please run ", analysis_type, " analysis first", sep = ""))
  } else {
    print(names(tmpIn))
  }
}

ptn_model <- function(ptn, analysis_type, model, comparison) {
  check_ptn(ptn)
  check_analysis_type(analysis_type)
  check_model(model, analysis_type = analysis_type)

  if (!check_number(comparison)) {
    stop(paste("Please provide a single numeric value for ", comparison, sep = ""))
  }

  tmpIn <- ptn@analysis@featureIntegration[[analysis_type]]
  if (comparison > length(tmpIn)) {
    stop(paste("There are only ", length(tmpIn), " comparisons", sep = ""))
  }
  tmpIn <- tmpIn[[comparison]]
  #
  tmpOut <- slot(tmpIn, model)
  #
  return(tmpOut)
}

ptn_selectedFeatures <- function(ptn, analysis_type, comparison) {
  check_ptn(ptn)
  check_analysis_type(analysis_type)
  if (!check_number(comparison)) {
    stop(paste("Please provide a single numeric value for ", comparison, sep = ""))
  }
  #
  tmpIn <- ptn@analysis@featureIntegration[[analysis_type]]
  if (comparison > length(tmpIn)) {
    stop(paste("There are only ", length(tmpIn), " comparisons", sep = ""))
  }
  tmpIn <- tmpIn[[comparison]]
  tmpOut <- slot(tmpIn, "selectedFeatures")
  #
  return(tmpOut)
}

ptn_networkGraph <- function(ptn, comparison) {
  check_ptn(ptn)
  if (!check_number(comparison)) {
    stop(paste("Please provide a single numeric value for ", comparison, sep = ""))
  }

  tmpIn <- ptn@analysis@featureIntegration$lm[[comparison]]

  tmpOut <- tmpIn@networkGraph
  return(tmpOut)
}
