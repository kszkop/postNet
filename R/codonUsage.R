codonUsage <- function(a2sU,
                       annotType = NULL,#"ccds", # option: 'refseq' or 'ccds'
                       sourceSeq = "load", # option to 'load' available or create new 'create',
                       analysis,
                       codonN = 1,
                       codSource = "sequence", # option: 'sequence' or if ribosomal profiling 'riboseq'
                       dpn_path = NULL, # path to dpn files if riboseq
                       cds_filt = TRUE,
                       pAdj = 0.01,
                       rem5=TRUE,
                       plotHeatmap = TRUE,
                       thresX1 = 0.3,
                       thresY1 = 0.3,
                       thresX2 = 0.3,
                       thresY2 = 0.3,
                       subregion = NULL, # number of nucleotides from start if positive or end if negative.
                       subregionSel= NULL, # select or exclude , required if subregion is not null.
                       inputTable=NULL,
                       comparisons,
                       plotType_index = 'boxplot',
                       pdfName = NULL) {
  #
  if (!checkUtils(a2sU)) {
    stop("a2sU is not a valid 'anota2seqUtilsData' object.")
  }
  if(!is_annotType(annotType)){
    stop("Please provide 'annotType', i.e source of annotation'")
  }
  if(annotType=="ccds"){
    species <- a2sU_species(a2sU)
    if(!is_valid_species(species)){
      stop("Please specify a species, at the moment only 'human' or 'mouse' are available).") 
    }
    selectionTmp <- slot(a2sU, 'selection')
    checkSelection(selectionTmp)
    if(!is_valid_sourceSeq(sourceSeq)){
      stop("Please provide 'sourceSeq', i.e. 'load' or 'create'")
    } 
    if (sourceSeq == "create") {
      if (species == "human") {
        #
        download.file("ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS.current.txt", destfile = "CCDS_human.txt")
        ccds <- read.delim(file = "CCDS_human.txt", as.is = c(F, T, T, F, T, F, F, T, T, T, F))
        unlink("CCDS_human.txt")
        #
        ccds$ccds_id_cl <- gsub("\\..*", "", ccds$ccds_id)
        
        #
        download.file("ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS_nucleotide.current.fna.gz",
                      destfile = "CCDS_nucleotide_human.fna.gz"
        )
        R.utils::gunzip("CCDS_nucleotide_human.fna.gz")
        ccdsSeq <- seqinr::read.fasta("CCDS_nucleotide_human.fna", seqtype = "AA")
        unlink("CCDS_nucleotide_human.fna")
        
        # Below, I just exclude these duplicated sequences
        ccdsSeq <- ccdsSeq[!duplicated(sapply(strsplit(names(ccdsSeq), "\\|"), function(x) x[1]))]
        
        ccdsSeq <- data.frame(id = sub("\\..*", "", names(ccdsSeq)), CDS_seq = t(as.data.frame(lapply(ccdsSeq, function(x) paste(x, collapse = "")))), row.names = NULL, stringsAsFactors = FALSE)
        
        annotTmp <- merge(unique(ccds[, c("ccds_id_cl", "gene")]), ccdsSeq, by.x = "ccds_id_cl", by.y = "id", all.y = TRUE, all.x = FALSE)
        #
        colnames(annotTmp)[2] <- 'geneID'
        lenTmp <- as.numeric(sapply(annotTmp$CDS_seq, function(x) length(seqinr::s2c(x))))
        annotTmp$lenTmp <- lenTmp
        #
        annotSel <- isoSel(annot = annotTmp, method = a2sU_selection(a2sU))
        colnames(annotSel)[1:3] <- c('id','geneID','CDS_seq')
        #
        annot <- new("anota2seqUtilsRegion",
                     id = annotSel$id,
                     geneID = annotSel$geneID,
                     seq = annotSel$CDS_seq)
        
        a2sU@annot@CCDS <- annot
      } else if (species == "mouse") {
        #
        download.file("ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_mouse/CCDS.current.txt", destfile = "CCDS_mouse.txt")
        ccds <- read.delim(file = "CCDS_mouse.txt", as.is = c(F, T, T, F, T, F, F, T, T, T, F))
        unlink("CCDS_mouse.txt")
        #
        ccds$ccds_id_cl <- gsub("\\..*", "", ccds$ccds_id)
        
        #
        download.file("ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_mouse/CCDS_nucleotide.current.fna.gz",
                      destfile = "CCDS_nucleotide_mouse.fna.gz"
        )
        R.utils::gunzip("CCDS_nucleotide_mouse.fna.gz")
        ccdsSeq <- seqinr::read.fasta("CCDS_nucleotide_mouse.fna", seqtype = "AA")
        unlink("CCDS_nucleotide_mouse.fna")
        
        # Below, I just exclude these duplicated sequences
        ccdsSeq <- ccdsSeq[!duplicated(sapply(strsplit(names(ccdsSeq), "\\|"), function(x) x[1]))]
        
        ccdsSeq <- data.frame(id = sub("\\..*", "", names(ccdsSeq)), CDS_seq = t(as.data.frame(lapply(ccdsSeq, function(x) paste(x, collapse = "")))), row.names = NULL, stringsAsFactors = FALSE)
        
        annotTmp <- merge(unique(ccds[, c("ccds_id_cl", "gene")]), ccdsSeq, by.x = "ccds_id_cl", by.y = "id", all.y = TRUE, all.x = FALSE)
        #
        colnames(annotTmp)[2] <- 'geneID'
        lenTmp <- as.numeric(sapply(annotTmp$CDS_seq, function(x) length(seqinr::s2c(x))))
        annotTmp$lenTmp <- lenTmp
        #
        annotSel <- isoSel(annot = annotTmp, method = a2sU_selection(a2sU))
        colnames(annotSel)[1:3] <- c('id','geneID','CDS_seq')
        #
        annot <- new("anota2seqUtilsRegion",
                     id = annotSel$id,
                     geneID = annotSel$geneID,
                     seq = annotSel$CDS_seq)
        
        a2sU@annot@CCDS <- annot
      }
    } else if (sourceSeq == "load") {
      # list existing species
      currTmp <- list.files(system.file("extdata/annotation/ccds", package = "anota2seqUtils"))
      #
      if (!species %in% currTmp) {
        stop("This option is only  available for species: human and mouse at the moment. Please use option createFromFile")
      }
      if (species == "human") {
        annotTmp <- read.delim(system.file(paste("extdata/annotation/ccds/human", sep = "/"), "humanDB_ccds.txt.gz", package = "anota2seqUtils"), stringsAsFactors = FALSE)
      }
      if (species == "mouse") {
        annotTmp <- read.delim(system.file(paste("extdata/annotation/ccds/mouse", sep = "/"), "mouseDB_ccds.txt.gz", package = "anota2seqUtils"), stringsAsFactors = FALSE) # }
      }
      lenTmp <- as.numeric(sapply(annotTmp$CDS_seq, function(x) length(seqinr::s2c(x))))
      annotTmp$lenTmp <- lenTmp
      #
      annotSel <- isoSel(annot = annotTmp, method = a2sU_selection(a2sU))
      colnames(annotSel)[1:3] <- c('id','geneID','CDS_seq')
      #
      annot <- new("anota2seqUtilsRegion",
                   id = annotSel$id,
                   geneID = annotSel$geneID,
                   seq = annotSel$CDS_seq)
      
      a2sU@annot@CCDS <- annot
    }
  } #else if(annotType=="refseq" | annotType=='riboseq'){
  #  
  #  annot <- a2sU_sequences(a2sU,region='CDS')
  #  
  #} 
  ####
  if(!is_valid_analysis(analysis)){
    stop("Please provide an 'analysis'. It can only be 'codon' or 'AA'")
  }
  if(!is_number(codonN)){
    stop("Please provide numerical value for 'codonN'")
  }
  if(!is_logical(rem5)){
    stop("rem5 can be only TRUE or FALSE")
  }
  checkcodSource(codSource)
  codSourse <- tolower(codSource)
  #
  if(!is.null(subregion) && (!is.numeric(subregion) || !length(subregion)==1)){
    stop("'subregion' must be a numeric and just number")
  }
  if (!is.null(subregionSel) && !subregionSel %in% c("select", "exclude")) {
    stop("'subregionSel' must be a character and only 'select' or 'exclude'")
  }
  
  if(!is.null(comparisons)){
    if(!checkComparisons(comparisons)){
      stop("'comparisons' must be a list of numeric vector for paired comparisons example: list(c(0,2),c(0,1)). 0 is always a background.")
    }
    #
    if(length(which(unique(unlist(comparisons))==0))>0 && is.null(a2sU_bg(a2sU))){
      stop(" 0 is always a background, but no background provided")
    }
  } else {
    stop(
      "For further steps pairs of comparison of codon composition between regulations should be provided with direction for each.",
      call. = FALSE
    )
  }
  #
  nameTmp <- ifelse(is.null(pdfName), analysis, paste(pdfName, analysis, sep = "_"))
  nameOut <- nameTmp
  #
  if(annotType=="ccds"){
    seqTmp <- a2sU_sequences(a2sU, region = 'CCDS')
    names(seqTmp) <- a2sU_geneID(a2sU, region = 'CCDS')
  } else {
    seqTmp <- a2sU_sequences(a2sU, region = 'CDS')
    names(seqTmp) <- a2sU_geneID(a2sU, region = 'CDS')
  }
  #remove last codon (stop codon)
  seqTmp <- remove_last3(seqTmp)
  #
  if (!is.null(subregion)) {
    #
    subSeq <- sapply(seqTmp, function(x) subset_seq(x, pos = subregion, subregionSel = subregionSel))
    #
    seqTmp <- subSeq
  }
  #
  if (codSource == "sequence") {
    geneIDs <- names(seqTmp)
    #ptm <- proc.time()
    codonTmp <- list()
    for (i in 1:length(geneIDs)) {
      codonTmp[[i]] <- codonCount(gene = geneIDs[i], seq = seqTmp[i], codonN = codonN)
    }
    #proc.time() - ptm
    codonAll <- data.table::rbindlist(codonTmp, use.names = TRUE, fill = TRUE)
    #proc.time() - ptm
  
    if (codonN == 1) {
      codonAll <- codonAll[!codonAll$AA %in% c("Stop", "W", "M", "O", "U"), ]
    }
    codonsAllOut <- new("anota2seqUtilsCodonsAll",
                       geneID = codonAll$geneID,
                       codon  = codonAll$codon,
                       AA = codonAll$AA,
                       count = codonAll$count,
                       frequency = codonAll$frequency,
                       AACountPerGene = codonAll$AACountPerGene,
                       relative_frequency = codonAll$relative_frequency)
  } else if (codSource == "riboSeq") {
    #
    stop('few things to correct')
    
    # Load dpn
    if (!is.null(dpn_path)) {
      checkDirectory(dpn_path)
      if (unlist(strsplit(dpn_path, ""))[length(unlist(strsplit(dpn_path, "")))] != "/") {
        dpn_path <- paste(dpn_path, "/", sep = "")
      }
    }
    dpn_files <- list.files(ifelse(is.null(dpn_path), ".", dpn_path), pattern = ".dpn$")
    #
    annotTmp <- annotBg
    annotTmp <- annotTmp[annotTmp$id %in% annotBgSel$id,]
    annotTmp$lenTmp <- as.numeric(sapply(annotTmp$CDS_seq, function(x) length(seqinr::s2c(x))))
    dataList <- lapply(dpn_files, readRiboDpn, dpn_path = dpn_path, annot = annotTmp, cds_filt = TRUE)
    names(dataList) <- dpn_files
    #
    dataListMerged <- list()
    for(i in unique(inputTable$condition)){
      dataListTmp <- dataList[inputTable$dpnFiles[which(inputTable$condition == i)]]
      #
      geneTmp <- as.character()
      for(s in 1:length(dataListTmp)){
        geneTmp <- append(geneTmp,names(dataListTmp[[s]])) 
      }
      #
      geneTmp<- unique(geneTmp)
      #
      cond_merged <- list()
      for(g in 1:length(geneTmp)){
        #Progress
        pb1 <- txtProgressBar(min=1, max=length(geneTmp), style=3)
        #Extract current transcirpt
        geneID <- geneTmp[g]
        #Extract all entries for that transcirpt in all samples in given condition
        conv <- lapply(dataListTmp, function(x) x[[geneID]])
        #Remove empty ones
        conv <- conv[lapply(conv,length)>0]
        #convert
        tmpConv <- t(plyr::ldply(conv, rbind,.id=NULL))
        #Select only reproducible peaks i.e. occur in selected number of samples in each condition
        vec <- apply(tmpConv,1,function(x) sum(x, na.rm=TRUE))
        #Add to output
        cond_merged[[geneID]] <- vec
        #Progress
        setTxtProgressBar(pb1, g)
      }
      dataListMerged[[i]] <- cond_merged
    }
    #calculate codons 
    #
    #
    seqs <- list()
    seqs[[1]] <- annotTmp$UTR5_seq
    seqs[[2]] <- annotTmp$CDS_seq
    #
    seqComb <- combSeq(seqIn = seqs)
    names(seqComb) <- annotTmp$geneID
    #
    dataCodon <- list()
    for(d in 1:length(dataListMerged)){
      dataMergedTmp <- dataListMerged[[d]]
      #
      outTmp <- list()
      for(j in 1:length(dataTmp)){
        #
        geneID <- names(dataTmp[j])
        print(geneID)
        
      
        #
        #cdsLengths_start <- as.numeric(sapply(annotTmp$UTR5_seq, function(x) length(seqinr::s2c(x))))+1
        #names(cdsLengths_start) <- annotTmp$geneID
        #
        #cdsLengths_end <- cdsLengths_start + annotTmp$lenTmp
        #names(cdsLengths_end) <- annotTmp$geneID
          
        #
        seqG <- seqComb[[geneID]]
        #seqCDS <- seqs[cdsLengths_start[geneID]:cdsLengths_end[geneID]]
        #
        codFreq <- seqinr::uco(seqG,index = "eff")
        names(codFreq) <- toupper(names(codFreq))
        if(!is.null(seqG)){
          #
          geneTmp <- dataTmp[[j]]
          codTmp <- sapply(as.numeric(names(geneTmp)), extract_seq, seqs=seqG)
          names(codTmp) <- as.numeric(geneTmp)
        }
        codObs <- as.numeric()
        for(cod in 1:length(codFreq)){
          codSum <- sum(as.numeric(names(codTmp[codTmp==names(codFreq[cod])])))
          codObs[cod] <- codSum
        }
        names(codObs) <- names(codFreq) 
        #rOut <- resid(lm(codObs ~ codFreq))
        #outTmp[[transID]] <- rOut
        outTmp[[geneID]] <- codObs
      }
      dataCodon[[d]] <- outTmp
    }
  }
  #
  #res <- list()
  #if (!is.null(ads) | !is.null(customBg)) {
  #  res[["background"]] <- annotBgSel$geneID
  #}

  #if (!is.null(ads)) {
  #  results <- anota2seqGetDirectedRegulations(ads)
  #  #
  #  resTmp <- vector("list", length = length(regulation))
  #  for (i in unique(contrast)) {
  #    resTmp[which(contrast == i)] <- results[[i]][regulation[contrast == i]]
  #  }
  #  names(resTmp) <- paste(regulation, paste("c", contrast, sep = ""), sep = "_")
  #  if (!is.null(geneList)) {
  #    resTmp <- append(resTmp, geneList)
  #  }
  #} else {
  #  resTmp <- geneList
  #}
  #res <- append(res, resTmp)
  #
  #if (length(res) < 2) {
  #  stop(
  #    "For comparison of codon composition between regulations, at least two regulations should be provided.",
  #    call. = FALSE
  #  )
  #}
  #
  #
  #indexes
  if(analysis == "codon" & codonN==1){
    if (species == "human") {
      codind <- read.delim(system.file("extdata/indexes/human/", "IndexesHuman.txt", package = "anota2seqUtils"))
    } else if (species == "mouse") {
      codind <- read.delim(system.file("extdata/indexes/mouse/", "IndexesMouse.txt", package = "anota2seqUtils"))
    } else {
      message('no available indexes for ', species)
    }
    
    #indexex
    indNames <- c('CAI', 'CBI', 'Fop','GC3s','tAI','L_aa')
    
    for(ind in indNames){
      #####
      index_sel <- codind[,which(colnames(codind)==ind)]
      names(index_sel) <- codind$external_gene_name
      
      resOutInd  <- resQuant(qvec = index_sel, a2sU = a2sU)
      colOutInd <- colPlot(a2sU)
      
      ##
      pdf(paste(nameOut, ind,'_index.pdf', sep = ""),width= 8,height=8, useDingbats = F)
      ylabel = paste(ind, 'index', sep = " ")
      plotUtils(resOutInd, colOutInd, comparisons, ylabel = ylabel, plotType = plotType_index)
      dev.off()
    }
  }
  #
  compTmpAll <- sort(unique(unlist(comparisons)),decreasing = F)
  compOut1 <- list()
  compOut2 <- list()
  compOut3 <- list()
  compOut4 <- list()
  for (i in compTmpAll) {
    if(i == 0){
      selTmp <- a2sU_bg(a2sU)
    } else {
      resTmp <- a2sU_dataIn(a2sU)
      selTmp <-  resTmp[[i]]
    }
    tmp <- codonAll[codonAll$geneID %in% selTmp, ]
    if (analysis == "codon") {
      #####
      tmp1 <- tmp %>% group_by(codon) %>% summarise(codonPerReg = sum(count))
      tmpSum <- tmp1$codonPerReg
      names(tmpSum) <- tmp1$codon
      
      if(i == 0){
        compOut1[["background"]] <- tmpSum
      } else {
        compOut1[[names(resTmp)[i]]] <- tmpSum
      }
      
      ####
      tmp2 <- tmp %>% group_by(codon) %>% summarise(freqPerReg = mean(frequency))
      tmpFreq <- tmp2$freqPerReg
      names(tmpFreq) <- tmp2$codon
      
      if(i == 0){
        compOut2[["background"]] <- tmpFreq
      } else {
        compOut2[[names(resTmp)[i]]] <- tmpFreq
      }
      
      #####
      tmp3 <- tmp %>% group_by(codon) %>% mutate(codonPerReg = sum(count))
      tmp3 <- tmp3 %>% group_by(AA) %>%  mutate(AAPerReg = sum(AACountPerGene))
      tmp3 <- subset(tmp3, !duplicated(codon))
      tmp3$codonNormAA <- tmp3$codonPerReg / tmp3$AAPerReg
      tmpcodonNormAA <- tmp3$codonNormAA
      names(tmpcodonNormAA) <- tmp3$codon
      
      if(i == 0){
        compOut3[["background"]] <- tmpcodonNormAA
      } else {
        compOut3[[names(resTmp)[i]]] <- tmpcodonNormAA
      }
      #
      tmp4 <- tmp %>%  group_by(codon) %>% summarise(freq = sum(frequency))
      tmpSumFreq <- tmp4$freq
      names(tmpSumFreq) <- tmp4$codon
      if(i == 0){
        compOut4[["background"]] <- tmpSumFreq
      } else {
        compOut4[[names(resTmp)[i]]] <- tmpSumFreq
      }
    } else if (analysis == "AA") {
      tmp1  <- tmp %>% group_by(AA) %>% summarise(AAPerReg = sum(count))
      #
      tmpSum <- tmp1$AAPerReg
      names(tmpSum) <- tmp1$AA
      
      if(i == 0){
        compOut1[["background"]] <- tmpSum
      } else {
        compOut1[[names(resTmp)[i]]] <- tmpSum
      }
      
      tmp4 <- tmp %>%  group_by(AA) %>% summarise(freq = sum(frequency))
      tmpSumFreq <- tmp4$freq
      names(tmpSumFreq) <- tmp4$AA
      if(i == 0){
        compOut4[["background"]] <- tmpSumFreq
      } else {
        compOut4[[names(resTmp)[i]]] <- tmpSumFreq
      }
    }
  }
  #
  codonsSel <- list()
  for (j in 1:length(comparisons)) {
    if (names(compOut1)[1] == 'background') {
      compTmp <- comparisons[[j]] + 1
    } else {
      compTmp <- comparisons[[j]]
    }
    resIn <- compOut1[compTmp]
    resIn <- do.call(rbind, resIn)
    
    #remove with 0 in both
    remInd <- as.numeric(which(apply(resIn, 2, sum)==0))
    resIn <- resIn[,-remInd]
    
    #
    if (sum(apply(resIn, 2, min) < 5) > 0) {
      if(isTRUE(rem5)){
        rem5Ind <- as.numeric(which(apply(resIn, 2, min) < 5))
        resIn <- resIn[,-rem5Ind]
      } else {
        warning("In the contingency table (counts of codons by geneset), some counts are lower than 5 which may invalidate the chi-squared test. It might be relevant to perform the analysis on a subset or groups of codons instead",
                call. = FALSE
        )
      }
    }
    
    chisqTest <- chisq.test(resIn)

    if (is.na(chisqTest$p.value)) {
      stop("The Chi-squared test could not be performed. Please check that all codons have at least one count for at least one geneset.")
    }


    StResiduals <- chisqTest$stdres
    signifLimit <- sqrt(qchisq(pAdj / (dim(resIn)[1] * dim(resIn)[2]), lower.tail = FALSE, df = 1))

    if (isTRUE(plotHeatmap)) {
      #
      residOut <- t(as.matrix(StResiduals))
      
      pdf(paste(nameOut, paste(colnames(residOut), collapse = "_"), "heatmap.pdf", sep = "_"), width = 20, height = 28, useDingbats = F)
      par(mar = c(10, 5, 5, 5), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 0.9, cex.main = 0.7, cex.lab = 0.7)
      col <- grDevices::colorRampPalette(c("blue", "white", "red"))(50)
      gplots::heatmap.2(residOut,
        col = col,
        breaks = c(seq(-ceiling(max(abs(residOut[, 1]))), ceiling(max(abs(residOut[, 1]))), length = 51)),
        margins = c(30, 50),
        key = TRUE,
        keysize = 0.5,
        dendrogram = "both",
        trace = "none",
        density.info = "none",
        labCol = colnames(residOut),
        key.par = list(cex = 0.90),
        cexCol = 2.1,
        cexRow = 1.4,
        key.xlab = "",
        lhei = c(3, 25),
        lwid = c(3, 13),
        na.rm = TRUE,
        main = ""
      )
      dev.off()
    }
    #
    signSel <- colnames(StResiduals)[which(abs(StResiduals[1, ]) > signifLimit)]
    
    ###
    if(length(signSel) == 0){
      message(paste("There are no significant codons for comparison ", j, sep = ""))
    } else {
      #
      if (analysis == "codon" & codonN==1) {
        resIn2 <- compOut2[compTmp]
        resIn2 <- as.data.frame(do.call(cbind, resIn2))
        resIn2$codon <- row.names(resIn2)
        resIn2$AA <- seqinr::aaa(seqinr::translate(seqinr::s2c(seqinr::c2s(rownames(resIn2)))))
        resIn2$AAcodon <- paste(seqinr::translate(seqinr::s2c(seqinr::c2s(rownames(resIn2)))), rownames(resIn2), sep = "_")

        #
        xylim <- roundNice(max(resIn2[, c(1, 2)]),direction = 'up')
        #
        pdf(paste(nameOut, paste(colnames(resIn2)[1:2], collapse = "_"), "averageFreq.pdf", sep = "_"), width = 8, height = 8, useDingbats = F)
        par(mar = c(10, 5, 5, 10), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 0.9, cex.main = 0.7, cex.lab = 0.7)
        pOut <- ggplot2::ggplot(resIn2, ggplot2::aes(!!sym(colnames(resIn2)[1]), !!sym(colnames(resIn2)[2]), col = AA)) +
          ggplot2::theme_bw() +
          ggplot2::xlim(0, xylim) +
          ggplot2::ylim(0, xylim) +
          ggplot2::geom_point(size = 0.5) +
          ggplot2::coord_fixed() +
          ggplot2::geom_line(ggplot2::aes(group = AA), col = "gray", linetype = 1, size = 0.2) +
          ggplot2::theme(legend.position = "none") +
          ggrepel::geom_text_repel(ggplot2::aes(label = AAcodon), size = 3, segment.size = 0.2) +
          ggplot2::labs(
            x = paste(colnames(resIn2)[1], "\n(average codon frequency - per thousand)", "\n is it really per thousand ??", sep = ""),
            y = paste(colnames(resIn2)[2], "\n(average codon frequency - per thousand)", "\n is it really per thousand ??", sep = "")
          )
        print(pOut)
        dev.off()

        #
        resIn3 <- compOut3[compTmp]
      
        resIn3 <- as.data.frame(do.call(cbind,  resIn3))
        resIn3$codon <- row.names( resIn3)
        resIn3$AA <- seqinr::aaa(seqinr::translate(seqinr::s2c(seqinr::c2s(rownames(resIn3)))))
        resIn3$AAcodon <- paste(seqinr::translate(seqinr::s2c(seqinr::c2s(rownames(resIn3)))), rownames(resIn3), sep = "_")

        #
        xylim <- roundNice(max(resIn3[, c(1, 2)]), direction='up')
        #
        pdf(paste(nameOut, paste(colnames(resIn3)[1:2], collapse = "_"), "codonUsage.pdf", sep = "_"), width = 8, height = 8, useDingbats = F)
        par(mar = c(10, 5, 5, 10), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 0.9, cex.main = 0.7, cex.lab = 0.7)
        pOut <- ggplot2::ggplot(resIn3, ggplot2::aes(!!sym(colnames(resIn3)[1]), !!sym(colnames(resIn3)[2]), col = AA)) +
          ggplot2::theme_bw() +
          ggplot2::xlim(0, xylim) +
          ggplot2::ylim(0, xylim) +
          ggplot2::geom_point(size = 0.5) +
          ggplot2::coord_fixed() +
          ggplot2::geom_line(ggplot2::aes(group = AA), col = "gray", linetype = 1, size = 0.2) +
          ggplot2::theme(legend.position = "none") +
          ggrepel::geom_text_repel(ggplot2::aes(label = AAcodon), size = 3, segment.size = 0.2) +
          ggplot2::labs(
            x = paste(colnames(resIn3)[1], "\n(Amino acid normalised codon usage)", sep = ""),
            y = paste(colnames(resIn3)[2], "\n(Amino acid normalised codon usage)", sep = "")
          )
          print(pOut)
          dev.off()
        #
     
        resIn4 <- data.frame(codon = colnames(resIn), t(resIn), row.names = NULL)
        resIn4$AA <- seqinr::aaa(seqinr::translate(seqinr::s2c(seqinr::c2s(resIn4$codon))))

        regComb <- combn(x = names(compOut1[compTmp]), m = 2)
        statOut <- log2(statOnDf(df = resIn4, regs = regComb, analysis = analysis)[signSel])
        #
        sumFreq <- compOut4[compTmp]
        sumFreq  <- rowSums(do.call(cbind, sumFreq))[signSel]

        #
        xlimT <- roundNice(max(abs(range(statOut))), direction='up')
        pdf(paste(nameOut, paste(regComb[,1], collapse = "_"), "codon_oddratio_vs_freq.pdf", sep = ""), width = 8, height = 8, useDingbats = F)
        m <- layout(mat = matrix(c(1, 2), nrow = 2, ncol = 1), heights = c(1, 5))
        #
        par(mar = c(0, 8, 0, 3),bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
        plot(0,1, xlim=c(-xlimT,xlimT), ylim=c(0,1),xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "", lwd = 1, bty = "n", font = 2,type="n", frame.plot = FALSE)
        legend(-xlimT, 0.5, paste(regComb[1],'vs',regComb[2],sep=" "), bty = "n", horiz = TRUE, xpd = T, y.intersp = 2.8, cex = 1.3)
        
        par(mar = c(5, 8, 0, 3), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.3)
        plot(statOut, sumFreq, col = "black", pch = 20, cex = 0.1, xlab = "", ylab = "", lwd = 1, bty = "n", font = 2, xlim = c(-xlimT, xlimT), ylim = range(sumFreq))#, xaxt = "n", yaxt = "n")
        text(statOut, sumFreq, names(sumFreq), col = "black", font = 2)

        mtext(side = 2, line = 3, "frequency ", col = "black", font = 2, cex = 1.7)
        mtext(side = 1, line = 3, "log2 odd ratio", col = "black", font = 2, cex = 1.7, at = 0)

        #axis(side = 2, seq(0, roundNice(range(sumFreq)[2], direction='up'), ifelse(roundNice(range(sumFreq)[2],direction='up') > 10, 10, 1)), font = 2, las = 2, lwd = 2)
        #axis(side = 1, seq(-xlimT, xlimT, 1), font = 2, lwd = 2)

        indUp <- as.numeric(which(statOut >= 0))
        statOut_up <- statOut[indUp]
        sumFreq_up <- sumFreq[indUp]
        codU <-  names(sumFreq_up)[sumFreq_up >= quantile(sumFreq_up, 1 - thresY1) & statOut_up >= quantile(statOut_up, 1 - thresX1)]
        if(length(codU) > 0){
          text(statOut_up[codU], sumFreq_up[codU], codU, col = "firebrick1", font = 2)
        } else {
          message('Non of the codons is labelled, please select more relaxed thresholds for thresX1 and thresY1')
        }
        #
        indDown <- as.numeric(which(statOut < 0))
        statOut_down <- statOut[indDown]
        sumFreq_down <- sumFreq[indDown]
        codD <- names(sumFreq_down)[sumFreq_down >= quantile(sumFreq_down, 1 - thresY2) & statOut_down <= quantile(statOut_down, thresX2)]
        if(length(codD) > 0){
          text(statOut_down[codD], sumFreq[codD], codD, col = "dodgerblue1", font = 2)
        } else {
          message('Non of the codons is labelled, please select more relaxed thresholds for thresX2 and thresY2')
        }
        dev.off()
      
      } else if (analysis == "AA") {

        resIn4 <- data.frame(AA = colnames(resIn), t(resIn), row.names = NULL)
        #
        regComb <- combn(x = names(compOut1[compTmp]), m = 2)
        statOut <- log2(statOnDf(df = resIn4, regs = regComb, analysis = analysis)[signSel])
      
        sumFreq <- compOut4[compTmp]
        sumFreq  <- rowSums(do.call(cbind, sumFreq))[signSel]
      
        xlimT <- roundNice(max(abs(range(statOut))), direction='up')
      
        pdf(paste(nameOut, paste(regComb[,1], collapse = "_"), "AA_oddratio_vs_freq.pdf", sep = ""), width = 8, height = 8, useDingbats = F)
        m <- layout(mat = matrix(c(1, 2), nrow = 2, ncol = 1), heights = c(1, 5))
        #
        par(mar = c(0, 8, 0, 3),bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
        plot(0,1, xlim=c(-xlimT,xlimT), ylim=c(0,1),xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "", lwd = 1, bty = "n", font = 2,type="n", frame.plot = FALSE)
        legend(-xlimT, 0.5, paste(regComb[1],'vs',regComb[2],sep=" "), bty = "n", horiz = TRUE, xpd = T, y.intersp = 2.8, cex = 1.3)
        
        par(mar = c(5, 8, 0, 3), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.3)
        plot(statOut, sumFreq, col = "black", pch = 20, cex = 0.1, xlab = "", ylab = "", lwd = 1, bty = "n", font = 2, xlim = c(-xlimT, xlimT), ylim = range(sumFreq))
        text(statOut, sumFreq, names(sumFreq), col = "black", font = 2)
        
        mtext(side = 2, line = 3, "frequency ", col = "black", font = 2, cex = 1.7)
        mtext(side = 1, line = 3, "log2 odd ratio", col = "black", font = 2, cex = 1.7, at = 0)
        
        #axis(side = 2, seq(0, roundNice(range(sumFreq)[2], direction='up'), ifelse(roundNice(range(sumFreq)[2],direction='up') > 10, 10, 1)), font = 2, las = 2, lwd = 2)
        #axis(side = 1, seq(-xlimT, xlimT, 1), font = 2, lwd = 2)
        
        indUp <- as.numeric(which(statOut >= 0))
        statOut_up <- statOut[indUp]
        sumFreq_up <- sumFreq[indUp]
        codU <-  names(sumFreq_up)[sumFreq_up >= quantile(sumFreq_up, 1 - thresY1) & statOut_up >= quantile(statOut_up, 1 - thresX1)]
        if(length(codU) > 0){
          text(statOut_up[codU], sumFreq_up[codU], codU, col = "firebrick1", font = 2)
        } else {
          message('Non of the AA is labelled, please select more relaxed thresholds for thresX1 and thresY1')
        }
        #
        indDown <- as.numeric(which(statOut < 0))
        statOut_down <- statOut[indDown]
        sumFreq_down <- sumFreq[indDown]
        codD <- names(sumFreq_down)[sumFreq_down >= quantile(sumFreq_down, 1 - thresY2) & statOut_down <= quantile(statOut_down, thresX2)]
        if(length(codD) > 0){
          text(statOut_down[codD], sumFreq[codD], codD, col = "dodgerblue1", font = 2)
        } else {
          message('Non of the AA is labelled, please select more relaxed thresholds for thresX2 and thresY2')
        }
        dev.off()
      }
    }
    codonsSel[[paste(regComb, collapse = "_vs_")]] <- list(Up = codU, Down = codD)
  }
  codonsOut <- new("anota2seqUtilsCodons",
                   codonsAll = codonsAllOut,
                   codonsSel  = codonsSel)

  
  a2sU@analysis@codons <- codonsOut
  return(a2sU)
}
