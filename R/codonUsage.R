codonUsage <- function(annotType = "ccds", # option: 'refseq' or 'ccds', 'custom', 'riboseq'
                       sourceSeq = "load", # option to 'load' available or create new 'create'
                       customFileCod = NULL,
                       analysis,
                       codonN = 1,
                       codSource = "sequence", # option: 'sequence' or if ribosomal profiling 'riboseq'
                       dpn_path = NULL, # path to dpn files if riboseq
                       cds_filt = TRUE,
                       pAdj = 0.01,
                       plotHeatmap = TRUE,
                       thresX1 = 0.3,
                       thresY1 = 0.3,
                       thresX2 = 0.3,
                       thresY2 = 0.3,
                       subregion = NULL, # number of nucleotides from start if positive or end if negative.
                       subregionSel= NULL, # select or exclude , required if subregion is not null.
                       inputTable=NULL,
                       comparisons,
                       pdfName = NULL) {
  #
  if (!checkUtils(a2sU)) {
    stop("a2sU is not a valid 'anota2seqUtilsData' object.")
  }
  #
  species <- anota2seqUtilsGetSpecies(a2sU)
  if(!is_valid_species(species)){
    stop("Please specify a species, at the moment only 'human' or 'mouse' are available).") 
  }
  #
  if(tolower(codSource)=='riboseq'){
    message("Please make sure the 'annot' was created using 'source == createFromFiles' option, and providing fasta file and corresponding postion file that was used for ribosomal profiling reads alignment. Otherwise positions in the dpn file will not match")
  }
  if(!is_annotType(annotType)){
    stop("Please provide 'annotType', i.e source of annotation'")
  }
  annotType <- tolower(annotType)
  
  if(annotType=="ccds"){
    if(tolower(codSource)=='riboseq'){
      stop("Please choose 'riboseq' if codSource is 'riboseq'")
    }
    sourceSeq <- tolower(sourceSeq)
    #if(is_valid_sourceSeq){
    #  stop("Please provide 'sourceSeq', i.e. 'load' or 'create'")
    #} else {
    #  if(!is_valid_species){
    #    stop("Please specify a species, at the moment only 'human' or 'mouse' are available")
    #  } else {
    #    species <- tolower(species)
    #  }
    #}
  } else if(annotType=="refseq" | annotType=='riboseq'){
    
    annot <- anota2seqUtilsGetCDS(a2sU)
    
  } else if(annotType=="custom"){
    if (is.null(customFileCod)) {
      stop("Please provide an customFileCod.")
    }
    annotTmp <- read.delim(customFileCod, stringsAsFactors = FALSE)
    checkAnnotCod (annotTmp)
    #
    lenTmp <- as.numeric(sapply(annotTmp$CDS_seq, function(x) length(seqinr::s2c(x))))
    annotTmp$lenTmp <- lenTmp
    #
    annotBgSelTmp <- isoSel(annot = annotTmp, method = anota2seqUtilsGetSelection(a2sU))
    #
    
    annot <- new("anota2seqUtilsRegion",
                       id = annotBgSelTmp$id,
                       geneID = annotBgSelTmp$geneID,
                       seq = annotBgSelTmp$CDS_seq)
  }
  if(!is_valid_analysis(analysis)){
    stop("Please provide an 'analysis'. It can only be 'codon' or 'AA'")
  } 
  checkSelection(selection)
  if(!is_number(codonN)){
    stop("Please provide numerical value for 'codonN'")
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
    if(length(which(unique(unlist(comparisons))==0))>0 && is.null(anota2seqUtilsGetBg(a2sU))){
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
  if (annotType == "ccds") {
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
        annotBgSelTmp <- isoSel(annot = annotTmp, method = anota2seqUtilsGetSelection(a2sU))
        #
        
        annot <- new("anota2seqUtilsRegion",
                           id = annotBgSelTmp$ccds_id_cl,
                           geneID = annotBgSelTmp$geneID,
                           seq = annotBgSelTmp$CDS_seq)
        
        #write.table(annot, file = "customDB_ccds.txt", col.names = T, row.names = F, sep = "\t", quote = F)
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
        annotBgSelTmp <- isoSel(annot = annotTmp, method = anota2seqUtilsGetSelection(a2sU))
        #
        
        annot <- new("anota2seqUtilsRegion",
                           id = annotBgSelTmp$ccds_id_cl,
                           geneID = annotBgSelTmp$geneID,
                           seq = annotBgSelTmp$CDS_seq)
        

        #write.table(annot, file = "customDB_ccds.txt", col.names = T, row.names = F, sep = "\t", quote = F)
      }
    } else if (sourceSeq == "load") {
      # list existing species
      currTmp <- list.files(system.file("extdata/annotation/ccds", package = "anota2seqUtils"))
      #
      species <- anota2seqUtilsGetSpecies(a2sU)
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
      annotBgSelTmp <- isoSel(annot = annotTmp, method = anota2seqUtilsGetSelection(a2sU))
      #
      
      annot <- new("anota2seqUtilsRegion",
                         id = annotBgSelTmp$ccds_id_cl,
                         geneID = annotBgSelTmp$geneID,
                         seq = annotBgSelTmp$CDS_seq)
    }
  }
  #
  if (!is.null(subregion)) {
    seqTmp <- anota2seqUtilsGetSeq(annot)
    #
    subSeq <- as.character(sapply(seqTmp, function(x) subset_seq(x, pos = subregion, subregionSel = subregionSel)))
    #
    seqTmp <- subSeq
  }

  if (codSource == "sequence") {
    geneID <- anota2seqUtilsGetGeneID(annot)
    #ptm <- proc.time()
    codonTmp <- list()
    for (i in 1:length(geneID)) {
      codonTmp[[i]] <- codonCount(gene = geneID[i], seq = seqTmp[i], codonN = codonN)
    }
    #proc.time() - ptm
    codonAll <- do.call(rbind, lapply(codonTmp, data.frame, stringsAsFactors = FALSE))
    if (codonN == 1) {
      codonAll <- codonAll[!codonAll$AA %in% c("Stp", "Trp", "Met"), ]
    }
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
    indNames <- c('CAI', 'CBI', 'Fop','GC3s','tai','L_aa')
    
    for(ind in indNames){
      #####
      index_sel <- codind[,which(colnames(codind)==indNames)]
      names(index_sel) <- codind$external_gene_name
      
      resOutInd  <- resQuant(qvec = index_sel, a2sU = a2sU)
      colOutInd <- colPlot(a2sU)
      
      ##
      pdf(paste(nameOut, indNames,'_index.pdf', sep = ""),width= 8,height=8, useDingbats = F)
      par(mar = c(12, 12, 5, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
      
      plotBoxplots(resOutInd, colOutInd, comparisons = comparisons, ylabel = paste(indNames,'index',sep='_')
                   
      dev.off()
    }
  }
  #
  compTmpAll <- sort(unique(unlist(comparisons)),decreasing = F)
  compOut1 <- list()
  compOut2 <- list()
  compOut3 <- list()
  
  for (i in compTmpAll) {
    if(i == 0){
      selTmp <- anota2seqUtilsGetBg(a2sU)
    } else {
      resTmp <- anota2seqUtilsGetDataIn(a2sU)
      selTmp <-  resTmp[[i]]
    }
    tmp <- codonAll[codonAll$geneID %in% selTmp, ]
    if (analysis == "codon") {
      
      #####
      tmp1 <- tmp %>% group_by(codon) %>% summarise(codonPerReg = sum(codonCount))
      tmpSum <- tmp1$codonPerReg
      names(tmpSum) <- tmp1$codon
      
      if(i == 0){
        compOut1[["background"]] <- tmpSum
      } else {
        compOut1[[names(resTmp)[i]]] <- tmpSum
      }
      
      ####
      tmp2 <- tmp %>% group_by(codon) %>% summarise(freqPerReg = mean(codonFreq))
      tmpFreq <- tmp2$freqPerReg
      names(tmpFreq) <- tmp2$codon
      
      if(i == 0){
        compOut2[["background"]] <- tmpFreq
      } else {
        compOut2[[names(resTmp)[i]]] <- tmpFreq
      }
      
      #####
      tmp3 <- tmp %>% group_by(codon) %>% mutate(codonPerReg = sum(codonCount))
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
    } else if (analysis == "AA") {
      tmp1  <- tmp %>% group_by(AA) %>% summarise(AAPerReg = sum(codonCount))
      #
      tmpSum <- tmp1$AAPerReg
      names(tmpSum) <- tmp1$AA
    }
    if(i == 0){
      compOut1[["background"]] <- tmpSum
    } else {
      compOut1[[names(resTmp)[i]]] <- tmpSum
    }
    
      resOutTmp <- data.frame(AA = colnames(resOut), t(resOut[compTmp, ]), row.names = NULL)
      #
      regComb <- combn(x = names(resSel), m = 2)
      statOut <- statOnDf(df = resOutTmp, regs = regComb, analysis = analysis)
    
      #
      freqTotal <- list()
        #
      gTmp <- union(resSelTmp[[1]], resSelTmp[[2]])
      tmp <- codonAll[codonAll$geneID %in% gTmp, ]
      #
      tmp2 <- tmp %>% group_by(AA) %>% summarise(freq = sum(codonFreq))
    
      #
     finalOut <- merge(tmp, data.frame(statOut), by.x = "AA", by.y = "row.names")
     finalOut <- finalOut[finalOut$AA %in% signSel, ]
    
    
  }

  
  #
  codonsOut <- list()
  for (j in 1:length(comparisons)) {
    if (!is.null(anota2seqUtilsGetBg(a2sU))) {
      compTmp <- comparisons[[j]] + 1
    } else {
      compTmp <- comparisons[[j]]
    }
    resIn <- resOut[compTmp, ]
    #
    chisqTest <- chisq.test(resIn)

    if (is.na(chisqTest$p.value)) {
      stop("The Chi-squared test could not be performed. Please check that all codons have at least one count for at least one geneset.")
    }

    if (sum(apply(resOut, 2, min) < 5) > 0) {
      warning("In the contingency table (counts of codons by geneset), some counts are lower than 5 which may invalidate the chi-squared test. It might be relevant to perform the analysis on a subset or groups of codons instead",
        call. = FALSE
      )
    }

    StResiduals <- chisqTest$stdres
    signifLimit <- sqrt(qchisq(pAdj / (dim(resOut)[1] * dim(resOut)[2]), lower.tail = FALSE, df = 1))

    if (isTRUE(plotHeatmap)) {
      #
      residOut <- t(as.matrix(StResiduals))
      
      pdf(paste(nameOut, paste(colnames(residOut), collapse = "_"), "_heatmap.pdf", sep = ""), width = 20, height = 28, useDingbats = F)
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
      if (analysis == "codon") {
        resSelTmp <- res[compTmp]

        compOut <- list()
        for (i in 1:2) {
          tmp <- codonAll[codonAll$geneID %in% resSelTmp[[i]], ]
          #
          tmp <- tmp %>%
            group_by(codon) %>%
            summarise(freqPerReg = mean(codonFreq))
          #
          tmpFreq <- tmp$freqPerReg
          names(tmpFreq) <- tmp$codon

          #
          compOut[[names(resSelTmp)[i]]] <- tmpFreq
        }
        compOut <- as.data.frame(do.call(cbind, compOut))
        compOut$codon <- row.names(compOut)
        compOut$AA <- seqinr::aaa(seqinr::translate(seqinr::s2c(seqinr::c2s(rownames(compOut)))))
        compOut$AAcodon <- paste(seqinr::translate(seqinr::s2c(seqinr::c2s(rownames(compOut)))), rownames(compOut), sep = "_")

      #
      xylim <- roundUpNice(max(compOut[, c(1, 2)]))
      #
      pdf(paste(nameOut, paste(colnames(compOut)[1:2], collapse = "_"), "_averageFreq.pdf", sep = ""), width = 8, height = 8, useDingbats = F)
      par(mar = c(10, 5, 5, 10), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 0.9, cex.main = 0.7, cex.lab = 0.7)
      pOut <- ggplot2::ggplot(compOut, ggplot2::aes(!!sym(colnames(compOut)[1]), !!sym(colnames(compOut)[2]), col = AA)) +
        ggplot2::theme_bw() +
        ggplot2::xlim(0, xylim) +
        ggplot2::ylim(0, xylim) +
        ggplot2::geom_point(size = 0.5) +
        ggplot2::coord_fixed() +
        ggplot2::geom_line(ggplot2::aes(group = AA), col = "gray", linetype = 1, size = 0.2) +
        ggplot2::theme(legend.position = "none") +
        ggrepel::geom_text_repel(ggplot2::aes(label = AAcodon), size = 3, segment.size = 0.2) +
        ggplot2::labs(
          x = paste(colnames(compOut)[1], "\n(average codon frequency - per thousand)", "\n is it really per thousand ??", sep = ""),
          y = paste(colnames(compOut)[2], "\n(average codon frequency - per thousand)", "\n is it really per thousand ??", sep = "")
        )
      print(pOut)
      dev.off()

      #
      compOut <- list()
      for (i in 1:2) {
        tmp <- codonAll[codonAll$geneID %in% resSelTmp[[i]], ]
        #
        tmp <- tmp %>%
          group_by(codon) %>%
          mutate(codonPerReg = sum(codonCount))
        tmp <- tmp %>%
          group_by(AA) %>%
          mutate(AAPerReg = sum(AACountPerGene))

        #
        tmp <- subset(tmp, !duplicated(codon))
        tmp$codonNormAA <- tmp$codonPerReg / tmp$AAPerReg

        #
        tmpcodonNormAA <- tmp$codonNormAA
        names(tmpcodonNormAA) <- tmp$codon

        #
        compOut[[names(resSelTmp)[i]]] <- tmpcodonNormAA
      }
      compOut <- as.data.frame(do.call(cbind, compOut))
      compOut$codon <- row.names(compOut)
      compOut$AA <- seqinr::aaa(seqinr::translate(seqinr::s2c(seqinr::c2s(rownames(compOut)))))
      compOut$AAcodon <- paste(seqinr::translate(seqinr::s2c(seqinr::c2s(rownames(compOut)))), rownames(compOut), sep = "_")

      #
      xylim <- roundUpNice(max(compOut[, c(1, 2)]))
      #
      pdf(paste(nameOut, paste(colnames(compOut)[1:2], collapse = "_"), "_codonUsage.pdf", sep = ""), width = 8, height = 8, useDingbats = F)
      par(mar = c(10, 5, 5, 10), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 0.9, cex.main = 0.7, cex.lab = 0.7)
      pOut <- ggplot2::ggplot(compOut, ggplot2::aes(!!sym(colnames(compOut)[1]), !!sym(colnames(compOut)[2]), col = AA)) +
        ggplot2::theme_bw() +
        ggplot2::xlim(0, xylim) +
        ggplot2::ylim(0, xylim) +
        ggplot2::geom_point(size = 0.5) +
        ggplot2::coord_fixed() +
        ggplot2::geom_line(ggplot2::aes(group = AA), col = "gray", linetype = 1, size = 0.2) +
        ggplot2::theme(legend.position = "none") +
        ggrepel::geom_text_repel(ggplot2::aes(label = AAcodon), size = 3, segment.size = 0.2) +
        ggplot2::labs(
          x = paste(colnames(compOut)[1], "\n(Amino acid normalised codon usage)", sep = ""),
          y = paste(colnames(compOut)[2], "\n(Amino acid normalised codon usage)", sep = "")
        )
      print(pOut)
      dev.off()

      resOutTmp <- data.frame(codon = colnames(resOut), t(resOut[compTmp, ]), row.names = NULL)
      resOutTmp$AA <- seqinr::aaa(seqinr::translate(seqinr::s2c(seqinr::c2s(resOutTmp$codon))))

      regComb <- combn(x = names(resSelTmp), m = 2)
      statOut <- statOnDf(df = resOutTmp, regs = regComb, analysis = analysis)

      #
      freqTotal <- list()
      gTmp <- union(resSelTmp[[1]], resSelTmp[[2]])
      tmp <- codonAll[codonAll$geneID %in% gTmp, ]
      #
      tmp <- tmp %>%
        group_by(codon) %>%
        summarise(freq = sum(codonFreq))

      #
      finalOut <- merge(tmp, data.frame(statOut), by.x = "codon", by.y = "row.names")
      finalOut <- finalOut[finalOut$codon %in% signSel, ]
      #
      xlimT <- roundUpNice(max(abs(range(log2(finalOut$statOut)))))
      pdf(paste(nameOut, paste(colnames(compOut)[1:2], collapse = "_"), "codon_oddratio_vs_freq.pdf", sep = ""), width = 8, height = 8, useDingbats = F)
      par(mar = c(5, 5, 5, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.3)
      plot(log2(finalOut$statOut), finalOut$freq, col = "black", pch = 20, cex = 0.1, xlab = "", ylab = "", lwd = 1, bty = "n", xaxt = "n", yaxt = "n", font = 2, xlim = c(-xlimT, xlimT), ylim = range(finalOut$freq))
      text(log2(finalOut$statOut), finalOut$freq, finalOut$codon, col = "black", font = 2)

      mtext(side = 2, line = 3, "frequency ", col = "black", font = 2, cex = 1.7)
      mtext(side = 1, line = 3, "log2 odd ratio", col = "black", font = 2, cex = 1.7, at = 0)

      axis(side = 2, seq(0, roundUpNice(range(finalOut$freq)[2]), ifelse(roundUpNice(range(finalOut$freq)[2]) > 10, 10, 1)), font = 2, las = 2, lwd = 2)
      axis(side = 1, seq(-xlimT, xlimT, 1), font = 2, lwd = 2)

      finalOut_u <- finalOut[finalOut$statOut >= 1, ]
      codU <- finalOut_u$codon[finalOut_u$freq >= quantile(finalOut_u$freq, 1 - thresY1) & finalOut_u$statOut >= quantile(finalOut_u$statOut, 1 - thresX1)]
      if(length(codU) > 0){
        finalOut_u <- finalOut_u[finalOut_u$codon %in% codU, ]
        text(log2(finalOut_u$statOut), finalOut_u$freq, finalOut_u$codon, col = "firebrick1", font = 2)
      } else {
        message('Non of the codons is labelled, please select more relaxed thresholds for thresX1 and thresY1')
      }
      #
      finalOut_d <- finalOut[finalOut$statOut < 1, ]
      codD <- finalOut_d$codon[finalOut_d$freq >= quantile(finalOut_d$freq, 1 - thresY2) & finalOut_d$statOut <= quantile(finalOut_d$statOut, thresX2)]
      if(length(codD) > 0){
        finalOut_d <- finalOut_d[finalOut_d$codon %in% codD, ]
        text(log2(finalOut_d$statOut), finalOut_d$freq, finalOut_d$codon, col = "dodgerblue1", font = 2)
      } else {
        message('Non of the codons is labelled, please select more relaxed thresholds for thresX2 and thresY2')
      }
      legend(-xlimT, roundUpNice(range(finalOut$freq)[2]), c(regComb[2], regComb[1]), fill = c("dodgerblue1", "firebrick1"), bty = "n", horiz = TRUE, xpd = T, y.intersp = 2.8, cex = 1.3)

      dev.off()
      
    } else if (analysis == "AA") {
      resOutTmp <- data.frame(AA = colnames(resOut), t(resOut[compTmp, ]), row.names = NULL)
      #
      regComb <- combn(x = names(resSel), m = 2)
      statOut <- statOnDf(df = resOutTmp, regs = regComb, analysis = analysis)

      #
      freqTotal <- list()
      #
      gTmp <- union(resSelTmp[[1]], resSelTmp[[2]])
      tmp <- codonAll[codonAll$geneID %in% gTmp, ]
      #
      tmp <- tmp %>%
        group_by(AA) %>%
        summarise(freq = sum(codonFreq))

      #
      finalOut <- merge(tmp, data.frame(statOut), by.x = "AA", by.y = "row.names")
      finalOut <- finalOut[finalOut$AA %in% signSel, ]
      #
      xlimT <- roundUpNice(max(abs(range(log2(finalOut$statOut)))))
      
      pdf(paste(nameOut, paste(colnames(compOut)[1:2], collapse = "_"), "AA_oddratio_vs_freq.pdf", sep = ""), width = 8, height = 8, useDingbats = F)
      par(mar = c(5, 8, 5, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.3)
      plot(finalOut$statOut, finalOut$freq, col = "black", pch = 20, cex = 0.1, xlab = "", ylab = "", lwd = 1, bty = "n", xaxt = "n", yaxt = "n", font = 2, xlim = c(-xlimT, xlimT), ylim = range(finalOut$freq))
      text(finalOut$statOut, finalOut$freq, finalOut$AA, col = "black")

      mtext(side = 2, line = 3, "frequency ", col = "black", font = 2, cex = 1.7)
      mtext(side = 1, line = 3, "log2 odd ratio", col = "black", font = 2, cex = 1.7, at = 0)

      axis(side = 2, seq(floor(range(finalOut$freq)[1]), ceiling(range(finalOut$freq)[2]), 10), font = 2, las = 2, lwd = 2)
      axis(side = 1, seq(-xlimT, xlimT, 21), font = 2, lwd = 2)

      finalOut_u <- finalOut[finalOut$statOut >= 1, ]
      codU <- finalOut_u$codon[finalOut_u$freq >= quantile(finalOut_u$freq, 1 - thresY1) & finalOut_u$statOut >= quantile(finalOut_u$statOut, 1 - thresX1)]
      finalOut_u <- finalOut_u[finalOut_u$codon %in% codU, ]
      text(log2(finalOut_u$statOut), finalOut_u$freq, finalOut_u$codon, col = "firebrick1", font = 2)

      finalOut_d <- finalOut[finalOut$statOut < 1, ]
      codD <- finalOut_d$codon[finalOut_d$freq >= quantile(finalOut_d$freq, 1 - thresY2) & finalOut_d$statOut <= quantile(finalOut_d$statOut, thresX2)]
      finalOut_d <- finalOut_d[finalOut_d$codon %in% codD, ]
      text(log2(finalOut_d$statOut), finalOut_d$freq, finalOut_d$codon, col = "dodgerblue1", font = 2)

      legend(-xlimT, roundUpNice(range(finalOut$freq)[2]), c(regComb[2], regComb[1]), fill = c("dodgerblue1", "firebrick1"), bty = "n", horiz = TRUE, xpd = T, y.intersp = 2.8, cex = 1.3)

      dev.off()
    }
    codonsOut[[paste(regComb, collapse = "_vs_")]] <- list(Up = codU, Down = codD)
  }
  codonsOut[["codonAll"]] <- codonAll
  return(codonsOut)
}
