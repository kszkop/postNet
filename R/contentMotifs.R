contentMotifs <- function(a2sU,
                          motifsIn,
                          seqType = "dna",
                          dist = 1,
                          min_score = 47,
                          unitOut = "number",
                          resid = FALSE,
                          region,
                          subregion = NULL,
                          subregionSel=NULL,
                          comparisons = NULL,
                          pdfName = NULL,
                          plotType = "ecdf",
                          num_threads = 1,
                          plotOut = TRUE) {
  
  #
  if (!checkUtils(a2sU)) {
    stop("a2sU is not a valid 'anota2seqUtilsData' object.")
  }
  checkRegion(region)

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
  
  if(!is.null(subregion) && (!is.numeric(subregion) || !length(subregion)==1)){
    stop("'subregion' must be a numeric and just number")
  }
  if (!is.null(subregionSel) && !subregionSel %in% c("select", "exclude")) {
    stop("'subregionSel' must be a character and only 'select' or 'exclude'")
  } 
  if(!is_number(dist)){
    stop("please provide numeric minimal distance between motifs")
  }
  if(!is_number(num_threads)){
    stop("please provide numeric number of threads you wish to use 'num_threads' ")
  } else {
    num_avail <- parallel::detectCores()
    if(num_threads > 0.75 *num_avail){
      num_threads <- round(0.75 *num_avail, digits = 0)
      message(paste('Number of threads were reduced to: ',num_threads,', i.e around 75% of your available resources', sep=''))
    }
  }
  if(!isUnitOut(unitOut)){
    stop("'unitOut' must be one from these: 'numeric' or 'position'")
  }
  if(!is_logical(resid)){
      stop("'resid', i.e whether the values should be normalised for the length, can only be only be logical: TRUE of FALSE ")
  }
  if(!is_valid_seq_type(seqType)){
    stop("'seqType' sequence type must be selected from one of these: 'dna', 'rna' or 'protein' ")
  } 
  if(!is_motifs(motifsIn)){
    stop("'motifsIn' should be not null character vector of sequence motifs ")
  } 
  #
  motifFinalRegion <- list()
  for(reg in region){
    
    seqTmp <- a2sU_sequences(a2sU,region=reg)
    #
    if (tolower(seqType) == "protein") {
      if(!is_by_3(seqTmp)){
        stop(" Not all sequences provided are multipliers of 3 so cannot be translated into proteins ")
      }
      proseqtmp <- sapply(seqTmp, function(x) seqinr::c2s(seqinr::translate(seqinr::s2c(x))))
      seqTmp <- proseqtmp
    }
    #
    if (!is.null(subregion)) {
      #
      subSeq <- sapply(seqTmp, function(x) subset_seq(x, pos = subregion, subregionSel = subregionSel))
      seqTmp <- subSeq
    }
    
    cl <- parallel::makeCluster(num_threads)
    doParallel::registerDoParallel(cl)

    motifsFinal <- foreach::foreach(i = 1:length(motifsIn), .combine = 'c', .packages = c('anota2seqUtils','foreach')) %dopar% {
      
    #motifsFinal <- list()
    #for (i in 1:length(motifsIn)) {
      motif <- motifsIn[i]
      #
      if (motif == "G4" & !tolower(seqType)=='protein') {
        if(!is_number(min_score)){
          stop("please provide numeric minimal score for g-quadruplexes selection")
        }
        motifOutTmp <- sapply(seqTmp, calc_g4, min_score = min_score)
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

        #cl2 <- parallel::makeCluster(num_threads)
        #doParallel::registerDoParallel(cl2)
        #motifOutTmp <- foreach(x = seqTmp, .combine = 'c', .packages = c('anota2seqUtils')) %dopar% {
        #  calc_motif(x, motifIn = motifTmp, dist = dist, unit = unitOut)
        #}        
        #parallel::stopCluster(cl2)
        if(unitOut == "number"){
          motifOutTmp <- unlist(motifOutTmp)
        } 
      } 
      #
      if (isTRUE(resid) & tolower(unitOut) == 'number') {
        #
        lenTmp <- sapply(seqTmp, function(x) length(seqinr::s2c(x)))
        #
        motifOut <- lm(as.numeric(motifOutTmp) ~ log2(as.numeric(lenTmp)))$residuals
      } else {
        motifOut <- motifOutTmp
      }
      names(motifOut) <- a2sU_geneID(a2sU, region=reg)
      #
      if (tolower(unitOut) == "number" & isTRUE(plotOut)) {
        nameTmp <- ifelse(is.null(pdfName), paste(region, motif, "content.pdf", sep = "_"), paste(pdfName, reg, motif, "content.pdf", sep = "_"))
        nameOut <- nameTmp
        #
        resOut <- resQuant(qvec = motifOut, a2sU = a2sU)
        
        colOut <- colPlot(a2sU)
        # Plot
        pdf(nameOut, width = 8, height = 8, useDingbats = F)
        ylabel <- paste(region, motif, sep = "_")
        plotUtils(resOut, colOut, comparisons, ylabel = ylabel ,plotType = plotType)
        dev.off()
      }
      #motifsFinal[[paste(reg, motif, sep = "_")]] <- motifOut
      setNames(list(motifOut), paste(reg, motif, sep = "_"))
      #list(paste(reg, motif, sep = "_") = motifOut)
    }
    parallel::stopCluster(cl)
    #
    motifFinalRegion <-  append(motifFinalRegion,motifsFinal)
  }
  return(motifFinalRegion)
}
