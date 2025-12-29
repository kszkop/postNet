contentMotifs <- function(ptn,
                          motifsIn,
                          seqType = "dna",
                          dist = 1,
                          min_score = 47,
                          unitOut = "number",
                          resid = FALSE,
                          region,
                          subregion = NULL,
                          subregionSel = NULL,
                          comparisons = NULL,
                          pdfName = NULL,
                          plotOut = TRUE) {
  #
  check_ptn(ptn)
  check_region(region)

  if (!check_logical(plotOut)) {
    stop("The input for 'plotOut' must be logical: TRUE or FALSE.")
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

  if (!is.null(subregion) && (!is.numeric(subregion) || !length(subregion) == 1)) {
    stop("The input for 'subregion' must be an integer.")
  }
  if (!is.null(subregionSel) && !subregionSel %in% c("select", "exclude")) {
    stop("The input for 'subregionSel' must be either 'select' or 'exclude'.")
  }
  if (!check_number(dist)) {
    stop("For 'dist', please provide a numeric value specifying the minimal distance between motifs.")
  }
  if (!isUnitOut(unitOut)) {
    stop("The input for 'unitOut' must be either 'number' or 'position'.")
  }
  if (!check_logical(resid)) {
    stop("The input for 'resid' (specifying whether the values should be normalised for the sequence length) must be logical: TRUE or FALSE.")
  }
  if (!is_valid_seq_type(seqType)) {
    stop("The input for 'seqType' must be either 'dna', 'rna', or 'protein'.")
  }
  if (!is_motifs(motifsIn)) {
    stop("The input for 'motifsIn' should be a character vector of sequence motifs to be detected and quantified. Ambiguities can be \
    specified using  IUPAC codes or [ ] (bracket) annotations. Optionally, G-quadruplexes can be specified as 'G4'.")
  }
  #
  motifFinalRegion <- list()
  for (reg in region) {
    seqTmp <- ptn_sequences(ptn, region = reg)
    names(seqTmp) <- ptn_geneID(ptn, region = reg)

    #
    if (tolower(seqType) == "protein") {
      if (!is_by_3(seqTmp)) {
        stop("Not all sequences provided can be divided into codons (are multiples of 3) so cannot be translated into protein sequences.")
      }
      proseqtmp <- sapply(seqTmp, function(x) seqinr::c2s(seqinr::translate(seqinr::s2c(x))))
      seqTmp <- proseqtmp
    }
    #
    if (!is.null(subregion)) {
      #
      subSeq <- sapply(seqTmp, function(x) subset_seq(x, pos = subregion, subregionSel = subregionSel))
      if (length(which(is.na(subSeq))) > 0) {
        message("For some sequences, the selected subregion is longer than the sequence region. These sequences will be removed from the analysis.")
      }
      seqTmp <- subSeq
    }
    seqTmp <- seqTmp[!is.na(seqTmp)]

    motifsFinal <- list()
    for (i in 1:length(motifsIn)) {
      motif <- motifsIn[i]
      #
      if (motif == "G4" & !tolower(seqType) == "protein") {
        if (!check_number(min_score)) {
          stop("Please provide a numeric value specifying the minimal score for G-quadruplexes selection ('min_score' parameter).")
        }
        motifOutTmp <- sapply(seqTmp, calc_g4, min_score = min_score, unit = unitOut)
      } else {
        motif <- toupper(motif)
        #
        if (tolower(seqType) == "dna" | tolower(seqType) == "rna") {
          motifTmp <- convertIUPAC(motif)
        } else {
          motifTmp <- replaceProtAmbig(motif)
        }
        #
        motifOutTmp <- lapply(seqTmp, function(x) calc_motif(x, motifIn = motifTmp, dist = dist, unit = unitOut))

        if (unitOut == "number") {
          motifOutTmp <- unlist(motifOutTmp)
        }
      }
      #
      if (tolower(unitOut) == "number") {
        if (length(which(motifOutTmp > 0))) {
          if (isTRUE(resid)) {
            #
            lenTmp <- sapply(seqTmp, function(x) length(seqinr::s2c(x)))
            #
            motifOut <- lm(as.numeric(motifOutTmp) ~ log2(as.numeric(lenTmp)))$residuals
          } else {
            motifOut <- motifOutTmp
          }
          names(motifOut) <- names(seqTmp)

          if (tolower(unitOut) == "number" & isTRUE(plotOut)) {
            nameTmp <- ifelse(is.null(pdfName), paste(reg, motif, "content.pdf", sep = "_"), paste(pdfName, reg, motif, "content.pdf", sep = "_"))
            nameOut <- nameTmp
            #
            resOut <- resQuant(qvec = motifOut, ptn = ptn)

            colOut <- colPlot(ptn)
            # Plot
            pdf(nameOut, width = 8, height = 8, useDingbats = FALSE)
            ylabel <- paste(reg, motif, sep = "_")
            plotPostNet(resOut, colOut, comparisons, ylabel = ylabel, plotType = "ecdf")
            dev.off()
          }
          motifsFinal[[paste(reg, motif, sep = "_")]] <- motifOut
        } else {
          message(paste(motifTmp, "does not have any sites", sep = " "))
        }
      } else {
        if (!any(sapply(motifOutTmp, function(x) !all(is.na(x))))) {
          message(paste(motifTmp, "does not have any sites", sep = " "))
        } else {
          motifsFinal[[paste(reg, motif, sep = "_")]] <- motifOutTmp
        }
      }
    }
    #
    motifFinalRegion <- append(motifFinalRegion, motifsFinal)
  }
  return(motifFinalRegion)
}
