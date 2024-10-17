motifAnalysis <- function(a2sU,
                          stremeThreshold = 0.05,
                          minwidth = 6,
                          memePath = NULL,
                          seqType = "dna",
                          region,
                          subregion = NULL,
                          subregionSel = NULL) {
  ####
  if (!checkUtils(a2sU)) {
    stop("a2sU is not a valid 'anota2seqUtilsData' object.")
  }
  checkRegion(region)

  if (is.null(a2sU_bg(a2sU))){
    stop("Background must be provided for motif analysis")
  }
  if(!is.null(subregion) && (!is.numeric(subregion) || !length(subregion)==1)){
    stop("'subregion' must be a numeric and just number")
  }
  if (!is.null(subregionSel) && !subregionSel %in% c("select", "exclude")) {
    stop("'subregionSel' must be a character and only 'select' or 'exclude'")
  } 
  if(!is_number(stremeThreshold)){
    stop("please provide numeric p-value threshold for motif selection")
  }
  if(!is_number(minwidth)){
    stop("please provide numeric minimal width for motif selection")
  }
  if(is.null(memePath)){
    stop("please provide path to meme suit")
  }
  if(!is_valid_seq_type(seqType)){
    stop("'seqType' sequence type must be selected from one of these: 'dna', 'rna' or 'protein' ")
  } 
  
  motifsOut <- new("anota2seqUtilsMotifs",
                  UTR5 = NULL,
                  CDS = NULL,
                  UTR3 = NULL)
  #
  motifsStemeOutRegion <- list()
  for(reg in region){
    #
    seqTmp <- a2sU_sequences(a2sU, region = reg)
    
    if (tolower(seqType) == "protein") {
      if(!is_by_3(seqTmp)){
        stop(" Not all sequences provided are multipliers of 3 so cannot be translated into proteins ")
      }
      proseqtmp <- sapply(seqTmp, function(x) seqinr::c2s(seqinr::translate(seqinr::s2c(x))))
      seqTmp <- proseqtmp
    }
    
    is_by_3 <- function(seqs) {
      all(sapply(seqs, function(x) length(seqinr::c2s(x)) %% 3 == 0))
    }
    
    #
    if (!is.null(subregion)) {
      #
      subSeq <- sapply(seqTmp, function(x) subset_seq(x, pos = subregion, subregionSel = subregionSel))
      seqTmp <- subSeq
    }
    #
    seqForAnalysis <- seqTmp
    names(seqForAnalysis) <- a2sU_geneID(a2sU, region=reg)
    #
    resOut <- resQuant(qvec = seqForAnalysis, a2sU = a2sU)
    
    controlSeq <- resOut[[1]]
    seqinr::write.fasta(sequences = as.list(as.character(controlSeq)), names = names(controlSeq), file.out = paste(paste("Control", reg, sep = "_"), ".fa", sep = ""))
    
    motifsTmpOut <- list()
    for (j in 2:length(resOut)) {
      #
      regSeq <- resOut[[j]]
      seqinr::write.fasta(sequences = as.list(as.character(regSeq)), names = names(regSeq), file.out = paste(paste("Regulated", reg, names(resOut)[j], sep = "_"), ".fa", sep = ""))
      #
      outdirTmp <- paste("stremeOut", reg, names(resOut)[j], sep = "_")
      streme_out <- memes::runStreme(input = paste(paste("Regulated", reg, names(resOut)[j], sep = "_"), ".fa", sep = ""), control = paste(paste("Control", reg, sep = "_"), ".fa", sep = ""), meme_path = memePath, alph = tolower(seqType), outdir = outdirTmp, minw = minwidth)
      if(nrow(streme_out)==0){
        message(paste('No motifs found in: ',paste(reg, names(resOut)[j], sep = "_"),sep=''))
      }
      #
      streme_out <- streme_out[streme_out$pval < stremeThreshold, ]
      if(nrow(streme_out)==0){
        message(paste('No motifs passed thresholds in: ',paste(reg, names(resOut)[j], sep = "_"),sep=''))
      }
      motifsTmpOut[[names(resOut)[j]]] <- streme_out
    }
    motifsStemeOut <- append(list(motifsOut = as.character(unlist(lapply(motifsTmpOut, function(x) x$consensus)))), motifsTmpOut)
    if(reg == 'UTR5'){
      motifsOut@UTR5 <- motifsStemeOut
    } else if (reg == 'CDS'){
      motifsOut@CDS <- motifsStemeOut
    } else if (reg == 'UTR3'){
      motifsOut@UTR3 <- motifsStemeOut
    }
  }
  a2sU@analysis@motifs <- motifsOut
  
  return(a2sU)
}
