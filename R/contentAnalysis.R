contentAnalysis <- function(ptn,
                            contentIn,
                            region,
                            subregion = NULL, 
                            subregionSel = NULL,
                            comparisons = NULL,
                            plotOut = TRUE,
                            plotType = "boxplot",
                            pdfName = NULL) {
  #
  check_ptn(ptn)
  check_region(region)

  if(!check_logical(plotOut)){
    stop("'plotOut' can only be only be logical: TRUE of FALSE ")
  } 
  if(isTRUE(plotOut)){
    if(!is.null(plotType)){
      check_plotType(plotType)
    } else {
      stop("Please provide 'plotType' to select option for plotting, from: 'boxplot','violin ,'ecdf'. ")
    }
  }
  if(!is.null(comparisons)){
    if(!check_comparisons(comparisons)){
      stop("'comparisons' must be a list of numeric vector for paired comparisons example: list(c(0,2),c(0,1)). 0 is always a background.")
    }
    #
    if(length(which(unique(unlist(comparisons))==0))>0 && is.null(ptn_background(ptn))){
      stop(" 0 is always a background, but no background provided")
    }
  }
  check_DNAsequence(contentIn)

  if(!is.null(subregion) && (!is.numeric(subregion) || !length(subregion)==1)){
    stop("'subregion' must be a numeric and just number")
  }
  if (!is.null(subregionSel) && !subregionSel %in% c("select", "exclude")) {
      stop("'subregionSel' must be a character and only 'select' or 'exclude'")
  } 
  
  ####
  contentFinal <- list()
  for(reg in toupper(region)){
    #
    seqTmp <- ptn_sequences(ptn, region = reg)
    names(seqTmp) <- ptn_geneID(ptn, region=reg)
    #
    if (!is.null(subregion)) {
      if(is.null(subregionSel)){
        stop("You have chosen option to select subset of the sequence. Please provide parameter 'subregionSel' to 'select' or 'exclude'")
      }
      #
      subSeq <- sapply(seqTmp, function(x) subset_seq(x, pos = subregion, subregionSel = subregionSel))
      
      if(length(which(is.na(subSeq)))>0){
        message('For some of the sequences the selected subregion is longer than the sequence region and these sequences will be removed')
      }
      seqTmp <- subSeq
    } 
    #
    for (i in 1:length(contentIn)) {
      content <- contentIn[i]
    
      ##
      contentTmp <- nPos_extract(content)
      
      if(!is.na(contentTmp$positions[1]) & reg != 'CDS'){
        next
      }

      contentOut <- as.numeric()
      for (i in 1:length(seqTmp)) {
        tmpSeq <- seqTmp[i]
        
        if (!is.na(contentTmp$positions[1])) {
          nPos <- contentTmp$positions
          tmpCont <- sapply(seqinr::s2c(toupper(contentTmp$nucleotide)), function(x) calc_content_pos(tmpSeq, nPos, x))
        } else {
          tmpCont <- sapply(seqinr::s2c(toupper(content)), function(x) calc_content(tmpSeq, x))
        } 
        contentOut[i] <- sum(tmpCont)
      }
      names(contentOut) <- names(seqTmp)
      contentOut <- contentOut[!is.na(contentOut)]
      
      #
      if (isTRUE(plotOut)) {
        resOut <- resQuant(qvec = contentOut, ptn = ptn)
        if(length(resOut)==0){
          stop('There are no regulated genes. Check the input or run without indicating regulation and comparisons')
        }
        if(diff(range(as.numeric(unlist(resOut)))) < .Machine$double.eps ^ 0.5){
          message(paste('No plot will be produced as all values are the same, (equal ', as.numeric(names(table(unlist(resOut)))), ') for', content, sep=' '))
        } else { 
          colOut <- colPlot(ptn)
          # Plot
          pdf(ifelse(is.null(pdfName), paste(reg, content, "content.pdf", sep = "_"), paste(pdfName, reg, content, "content.pdf", sep = "_")), width = 8, height = 8, useDingbats = F)
          ylabel = paste(paste0(content, " content"), 'in ', reg,  '(%)', sep = " ")
          plotPostNet(resOut, colOut, comparisons, ylabel = ylabel ,plotType = plotType)
          dev.off()
        }
      }
      contentFinal[[paste(reg, content, sep = "_")]] <- contentOut
    }
  }
  #
  return(contentFinal)
}

