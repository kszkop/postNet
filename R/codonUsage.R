codonUsage <- function(annot=NULL,
                       annotType = "ccds", # option: 'refseq' or 'ccds', 'custom
                       sourceCod = "load", # option to 'load' available or create new 'create'
                       customFileCod = NULL,
                       species = NULL, # it is required for ccds
                       analysis,
                       codonN = 1,
                       type = "sequence", # option: 'sequence' or if ribosomal profiling 'riboSeq'
                       dpn_path = NULL, # path to dpn files if riboSeq
                       cds_filt = TRUE,
                       pAdj = 0.01,
                       plotHeatmap = TRUE,
                       thresX1 = 0.3,
                       thresY1 = 0.3,
                       thresX2 = 0.3,
                       thresY2 = 0.3,
                       ads = NULL,
                       regulation = NULL,
                       contrast = NULL,
                       geneList = NULL,
                       geneListcolours = NULL,
                       customBg = NULL,
                       selection, # shortest, longest, random (default)
                       subregion = NULL, # number of nucleotides from start if positive or end if negative.
                       subregionSel, # select or exclude , required if subregion is not null.
                       comparisons,
                       outDir = NULL,
                       pdfName = NULL) {
  #
  if (!is.null(outDir)){
    dirTmp <- outDir
  } else {
    dirTmp <- paste('codonUsage', format(Sys.time(), "%Y%m%e_%X"),sep='_')
  }
  dir.create(dirTmp)
  #
  nameTmp <- ifelse(is.null(pdfName), analysis, paste(pdfName, analysis, sep = "_"))
  nameOut <- paste(dirTmp,nameTmp, sep='/')
  #
  if (is.null(comparisons)) {
    stop(
      "For further steps pairs of comparison of codon composition between regulations should be provided with direction for each.",
      call. = FALSE
    )
  }
  #
  if (annotType == "ccds") {
    if (sourceCod == "create") {
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
        gunzip("CCDS_nucleotide_human.fna.gz")
        ccdsSeq <- seqinr::read.fasta("CCDS_nucleotide_human.fna", seqtype = "AA")
        unlink("CCDS_nucleotide_human.fna")

        # Below, I just exclude these duplicated sequences
        ccdsSeq <- ccdsSeq[!duplicated(sapply(strsplit(names(ccdsSeq), "\\|"), function(x) x[1]))]

        ccdsSeq <- data.frame(id = sub("\\..*", "", names(ccdsSeq)), CDS_seq = t(as.data.frame(lapply(ccdsSeq, function(x) paste(x, collapse = "")))), row.names = NULL, stringsAsFactors = FALSE)

        annot <- merge(unique(ccds[, c("ccds_id_cl", "gene")]), ccdsSeq, by.x = "ccds_id_cl", by.y = "id", all.y = TRUE, all.x = FALSE)
        colnames(annot) <- c("id", "geneID", "CDS_seq")
        write.table(annot, file = "customDB_ccds.txt", col.names = T, row.names = F, sep = "\t", quote = F)
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
        gunzip("CCDS_nucleotide_mouse.fna.gz")
        ccdsSeq <- seqinr::read.fasta("CCDS_nucleotide_mouse.fna", seqtype = "AA")
        unlink("CCDS_nucleotide_mouse.fna")

        # Below, I just exclude these duplicated sequences
        ccdsSeq <- ccdsSeq[!duplicated(sapply(strsplit(names(ccdsSeq), "\\|"), function(x) x[1]))]

        ccdsSeq <- data.frame(id = sub("\\..*", "", names(ccdsSeq)), CDS_seq = t(as.data.frame(lapply(ccdsSeq, function(x) paste(x, collapse = "")))), row.names = NULL, stringsAsFactors = FALSE)

        annot <- merge(unique(ccds[, c("ccds_id_cl", "gene")]), ccdsSeq, by.x = "ccds_id_cl", by.y = "id", all.y = TRUE, all.x = FALSE)
        colnames(annot) <- c("id", "geneID", "CDS_seq")
        write.table(annot, file = "customDB_ccds.txt", col.names = T, row.names = F, sep = "\t", quote = F)
      }
    } else if (sourceCod == "load") {
      # list existing species
      currTmp <- list.files(system.file("extdata/annotation/ccds", package = "anota2seqUtils"))

      if (!species %in% currTmp) {
        stop("This option is only  available for species: human and mouse at the moment. Please use option createFromFile")
      }
      if (species == "human") {
        annot <- read.delim(system.file(paste("extdata/annotation/ccds/human", sep = "/"), "humanDB_ccds.txt.gz", package = "anota2seqUtils"), stringsAsFactors = FALSE)
      }
      if (species == "mouse") {
        annot <- read.delim(system.file(paste("extdata/annotation/ccds/mouse", sep = "/"), "mouseDB_ccds.txt.gz", package = "anota2seqUtils"), stringsAsFactors = FALSE) # }
      }
    }
  } else if (annotType == "refSeq") {
    annot <- annot[, c(1, 2, 4)]
    colnames(annot) <- c("id", "geneID", "CDS_seq")
  } else if (annotType == "custom") {
    annot <- read.delim(customFileCod, stringsAsFactors = FALSE)
    colnames(annot) <- c("id", "geneID", "CDS_seq")
  } else {
    stop("No correct option for annotation file provided")
  }
  if (type == "sequence") {
    # Subset annot for only expressed genes
    annotBg <- gSel(annot = annot, ads = ads, customBg = customBg, geneList = geneList)
    # Select region of interest
    annotTmp <- regSel(annot = annotBg, region = "CDS")
    # Per gene
    annotBgSel <- isoSel(annot = annotTmp, method = selection)
    #
    if (!is.null(subregion)) {
      #
      subSeq <- as.character(sapply(annotBgSel$seqTmp, function(x) subset_seq(x, pos = subregion, subregionSel = subregionSel)))
      #
      annotBgSel$seqTmp <- subSeq
    }
    annotBgSel <- annotBgSel[!is.na(annotBgSel$seqTmp), ]
    #
    codonTmp <- list()
    for (i in 1:nrow(annotBgSel)) {
      codonTmp[[i]] <- codonCount(gene = annotBgSel$geneID[i], seq = annotBgSel$seqTmp[i], codonN = codonN)
    }

    # ptm <- proc.time()
    # inTmp <- annotBgSel$seqTmp
    # names(inTmp) <- annotBgSel$geneID
    # test <- sapply(names(inTmp), codonCount, seq=inTmp,codonN=codonN )
    # Stop the clock
    # proc.time() - ptm

    codonAll <- do.call(rbind, lapply(codonTmp, data.frame, stringsAsFactors = FALSE))
    if (codonN == 1) {
      codonAll <- codonAll[!codonAll$AA %in% c("Stp", "Trp", "Met"), ]
    }
  } else if (type == "riboSeq") {
    # Remember to add information about annotation file so it is compatibile with the one used for alignement...
    # Load dpn
    if (!is.null(dpn_path)) {
      if (unlist(strsplit(dpn_path, ""))[length(unlist(strsplit(dpn_path, "")))] != "/") {
        dpn_path <- paste(dpn_path, "/", sep = "")
      }
    }
    #
    annot$lenTmp <- as.numeric(sapply(annot$CDS_seq, function(x) length(seqinr::s2c(x))))
    annotSel <- isoSel(annot = annot, method = selection)
    #
    dpn_files <- list.files(ifelse(is.null(dpn_path), ".", dpn_path), pattern = ".dpn$")
    dataList <- lapply(dpn_files, readRiboDpn, dpn_path = dpn_path, annot = annotSel, cds_filt = cds_filt)


    condition <- c(1, 1, 1, 1, 2, 2, 2, 2)
    InputTable <- data.frame(dpnFiles = dpn_files, condition = condition, stringsAsFactors = FALSE)
  }
  #
  res <- list()
  if (!is.null(ads) | !is.null(customBg)) {
    res[["background"]] <- annotBgSel$geneID
  }

  if (!is.null(ads)) {
    results <- anota2seqGetDirectedRegulations(ads)
    #
    resTmp <- vector("list", length = length(regulation))
    for (i in unique(contrast)) {
      resTmp[which(contrast == i)] <- results[[i]][regulation[contrast == i]]
    }
    names(resTmp) <- paste(regulation, paste("c", contrast, sep = ""), sep = "_")
    if (!is.null(geneList)) {
      resTmp <- append(resTmp, geneList)
    }
  } else {
    resTmp <- geneList
  }
  res <- append(res, resTmp)
  #
  if (length(res) < 2) {
    stop(
      "For comparison of codon composition between regulations, at least two regulations should be provided.",
      call. = FALSE
    )
  }
  #
  resOut <- list()
  for (i in 1:length(res)) {
    tmp <- codonAll[codonAll$geneID %in% res[[i]], ]
    if (analysis == "codon") {
      tmp <- tmp %>% group_by(codon) %>% summarise(codonPerReg = sum(codonCount))
      #
      tmpSum <- tmp$codonPerReg
      names(tmpSum) <- tmp$codon
    } else if (analysis == "AA") {
      tmp <- tmp %>% group_by(AA) %>% summarise(AAPerReg = sum(codonCount))
      #
      tmpSum <- tmp$AAPerReg
      names(tmpSum) <- tmp$AA
    }
    resOut[[names(res)[i]]] <- tmpSum
  }
  #
  resOut <- do.call(rbind, resOut)
  
  #indexes
  if(analysis == "codon"){
    if (species == "human") {
      codind <- data(codonIndexesHuman)
    }
    if (species == "mouse") {
      codind <- data(codonIndexesMouse)
    }
    #####CAI
    index_sel <- codind[,which(colnames(codind)=='CAI')]
    names(index_sel) <- codind$external_gene_name
  
    resOutInd <- resSel(vIn = index_sel, ads = ads, regulation = regulation, contrast = contrast, customBg = customBg, geneList = geneList)
    coloursOutInd <- coloursSel(ads = ads, regulation = regulation, geneList = geneList, geneListcolours = geneListcolours, customBg = customBg)
  
    ##
    pdf(paste(nameOut,'_CAI_index.pdf', sep = ""),width= 8,height=8, useDingbats = F)
    par(mar = c(12, 12, 5, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
    # 
    if (!is.null(regulation)) {
      xlimIn <- c(0.5, length(regulation) + ifelse(!is.null(geneList), length(geneList), 0) + 1.5)
    } else {
      xlimIn <- c(0.5, length(geneList) + 1.5)
    }
    plot(1, 1, xlim = xlimIn, ylim = c(0, range(as.numeric(unlist(resOutInd)))[2] + (1.25 * length(comparisons))), xaxt = "n", xlab = "", ylab = "", type = "n", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, frame.plot = FALSE)
  
    axis(side = 2, at=seq(0,1,0.2),labels = seq(0,1,0.2), font = 2, las = 2, lwd = 2)
    mtext(side = 2, line = 6, "CAI index", col = "black", font = 2, cex = 1.7, at = 0.5)
    text(1:length(resOut), par("usr")[3] - 0.45, labels = names(resOutInd), xpd = NA, cex = 0.9, srt = 45, adj = 1)
  
    if (!is.null(ads) | !is.null(customBg)) {
      abline(lty = 5, h = median(resOutInd[[1]]))
    }
    #
    for (i in 1:length(resOutInd)) {
      boxplot(resOutInd[[i]], add = TRUE, at = i, col = coloursOutInd[i], xaxt = "n", xlab = "", ylab = "", type = "n", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, frame.plot = FALSE, outcol = "grey65", whiskcol = "grey65", outline = FALSE, medcol = "black", staplelty = 0, whisklty = 1)
      text(i, 0, round(mean(antilog(resOutInd[[i]], 2), 0)), font = 2)
    }
    # Plot stats
    if (!is.null(comparisons)) {
      addStats(comparisons, ads, customBg, plotType, resOutInd, coloursOutInd)
    }
    dev.off()
  
    ######CBI
    index_sel <- codind[,which(colnames(codind)=='CBI')]
    names(index_sel) <- codind$external_gene_name
    
    resOutInd <- resSel(vIn = index_sel, ads = ads, regulation = regulation, contrast = contrast, customBg = customBg, geneList = geneList)
    coloursOutInd <- coloursSel(ads = ads, regulation = regulation, geneList = geneList, geneListcolours = geneListcolours, customBg = customBg)
  
    ##
    pdf(paste(nameOut,'_CBI_index.pdf', sep = ""),width= 8,height=8, useDingbats = F)
    par(mar = c(12, 12, 5, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
    # 
    if (!is.null(regulation)) {
      xlimIn <- c(0.5, length(regulation) + ifelse(!is.null(geneList), length(geneList), 0) + 1.5)
    } else {
      xlimIn <- c(0.5, length(geneList) + 1.5)
    }
    plot(1, 1, xlim = xlimIn, ylim = c(-1, range(as.numeric(unlist(resOutInd)))[2] + (1.25 * length(comparisons))), xaxt = "n", xlab = "", ylab = "", type = "n", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, frame.plot = FALSE)
  
    axis(side = 2, at=seq(-1,1,0.25),labels = seq(-1,1,0.25), font = 2, las = 2, lwd = 2)
    mtext(side = 2, line = 6, "CBI index", col = "black", font = 2, cex = 1.7, at = 0.5)
    text(1:length(resOut), par("usr")[3] - 0.45, labels = names(resOutInd), xpd = NA, cex = 0.9, srt = 45, adj = 1)
  
    if (!is.null(ads) | !is.null(customBg)) {
      abline(lty = 5, h = median(resOutInd[[1]]))
    }
    #
    for (i in 1:length(resOutInd)) {
      boxplot(resOutInd[[i]], add = TRUE, at = i, col = coloursOutInd[i], xaxt = "n", xlab = "", ylab = "", type = "n", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, frame.plot = FALSE, outcol = "grey65", whiskcol = "grey65", outline = FALSE, medcol = "black", staplelty = 0, whisklty = 1)
      text(i, 0, round(mean(antilog(resOutInd[[i]], 2), 0)), font = 2)
    }
    # Plot stats
    if (!is.null(comparisons)) {
      addStats(comparisons, ads, customBg, plotType, resOutInd, coloursOutInd)
    }
    dev.off()
  
    ######Fop
    index_sel <- codind[,which(colnames(codind)=='Fop')]
    names(index_sel) <- codind$external_gene_name
  
    resOutInd <- resSel(vIn = index_sel, ads = ads, regulation = regulation, contrast = contrast, customBg = customBg, geneList = geneList)
    coloursOutInd <- coloursSel(ads = ads, regulation = regulation, geneList = geneList, geneListcolours = geneListcolours, customBg = customBg)
  
    ##
    pdf(paste(nameOut,'_Fop_index.pdf', sep = ""),width= 8,height=8, useDingbats = F)
    par(mar = c(12, 12, 5, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
    # 
    if (!is.null(regulation)) {
      xlimIn <- c(0.5, length(regulation) + ifelse(!is.null(geneList), length(geneList), 0) + 1.5)
    } else {
      xlimIn <- c(0.5, length(geneList) + 1.5)
    }
    plot(1, 1, xlim = xlimIn, ylim = c(0, range(as.numeric(unlist(resOutInd)))[2] + (1.25 * length(comparisons))), xaxt = "n", xlab = "", ylab = "", type = "n", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, frame.plot = FALSE)
  
    axis(side = 2, at=seq(0,1,0.2),labels = seq(0,1,0.2), font = 2, las = 2, lwd = 2)
    mtext(side = 2, line = 6, "Fop index", col = "black", font = 2, cex = 1.7, at = 0.5)
    text(1:length(resOut), par("usr")[3] - 0.45, labels = names(resOutInd), xpd = NA, cex = 0.9, srt = 45, adj = 1)
  
    if (!is.null(ads) | !is.null(customBg)) {
      abline(lty = 5, h = median(resOutInd[[1]]))
    }
    #
    for (i in 1:length(resOutInd)) {
      boxplot(resOutInd[[i]], add = TRUE, at = i, col = coloursOutInd[i], xaxt = "n", xlab = "", ylab = "", type = "n", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, frame.plot = FALSE, outcol = "grey65", whiskcol = "grey65", outline = FALSE, medcol = "black", staplelty = 0, whisklty = 1)
      text(i, 0, round(mean(antilog(resOutInd[[i]], 2), 0)), font = 2)
    }
    # Plot stats
    if (!is.null(comparisons)) {
      addStats(comparisons, ads, customBg, plotType, resOutInd, coloursOutInd)
    }
    dev.off()
  
    #####GC3s
    index_sel <- codind[,which(colnames(codind)=='GC3s')]
    names(index_sel) <- codind$external_gene_name
  
    resOutInd <- resSel(vIn = index_sel, ads = ads, regulation = regulation, contrast = contrast, customBg = customBg, geneList = geneList)
    coloursOutInd <- coloursSel(ads = ads, regulation = regulation, geneList = geneList, geneListcolours = geneListcolours, customBg = customBg)
  
    ##
    pdf(paste(nameOut,'_GC3s_index.pdf', sep = ""),width= 8,height=8, useDingbats = F)
    par(mar = c(12, 12, 5, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
    # 
    if (!is.null(regulation)) {
      xlimIn <- c(0.5, length(regulation) + ifelse(!is.null(geneList), length(geneList), 0) + 1.5)
    } else {
      xlimIn <- c(0.5, length(geneList) + 1.5)
    }
    plot(1, 1, xlim = xlimIn, ylim = c(0, range(as.numeric(unlist(resOutInd)))[2] + (1.25 * length(comparisons))), xaxt = "n", xlab = "", ylab = "", type = "n", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, frame.plot = FALSE)
  
    axis(side = 2, at=seq(0,1,0.2),labels = seq(0,1,0.2), font = 2, las = 2, lwd = 2)
    mtext(side = 2, line = 6, "GC3s index", col = "black", font = 2, cex = 1.7, at = 0.5)
    text(1:length(resOut), par("usr")[3] - 0.45, labels = names(resOutInd), xpd = NA, cex = 0.9, srt = 45, adj = 1)
  
    if (!is.null(ads) | !is.null(customBg)) {
      abline(lty = 5, h = median(resOutInd[[1]]))
    }
    #
    for (i in 1:length(resOutInd)) {
      boxplot(resOutInd[[i]], add = TRUE, at = i, col = coloursOutInd[i], xaxt = "n", xlab = "", ylab = "", type = "n", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, frame.plot = FALSE, outcol = "grey65", whiskcol = "grey65", outline = FALSE, medcol = "black", staplelty = 0, whisklty = 1)
      text(i, 0, round(mean(antilog(resOutInd[[i]], 2), 0)), font = 2)
    }
    # Plot stats
    if (!is.null(comparisons)) {
      addStats(comparisons, ads, customBg, plotType, resOutInd, coloursOutInd)
    }
    dev.off()
  
    #######Nc
    tmp <- codind[!codind$Nc=='*****',]
    tmp <- tmp[!is.na(tmp$Nc),]
    index_sel <- as.numeric(levels(tmp$Nc))[tmp$Nc]
    names(index_sel) <- tmp$external_gene_name
  
    resOutInd <- resSel(vIn = index_sel, ads = ads, regulation = regulation, contrast = contrast, customBg = customBg, geneList = geneList)
    coloursOutInd <- coloursSel(ads = ads, regulation = regulation, geneList = geneList, geneListcolours = geneListcolours, customBg = customBg)
  
    pdf(paste(nameOut,'_Nc_index.pdf', sep = ""),width= 8,height=8, useDingbats = F)
    par(mar = c(12, 12, 5, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
    # 
    if (!is.null(regulation)) {
      xlimIn <- c(0.5, length(regulation) + ifelse(!is.null(geneList), length(geneList), 0) + 1.5)
    } else {
      xlimIn <- c(0.5, length(geneList) + 1.5)
    }
    plot(1, 1, xlim = xlimIn, ylim = c(0, range(as.numeric(unlist(resOutInd)))[2] + (1.25 * length(comparisons))), xaxt = "n", xlab = "", ylab = "", type = "n", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, frame.plot = FALSE)
  
    axis(side = 2, at=seq(0,100,20),labels = seq(0,100,20), font = 2, las = 2, lwd = 2)
    mtext(side = 2, line = 6, "Nc index", col = "black", font = 2, cex = 1.7, at = 0.5)
    text(1:length(resOut), par("usr")[3] - 0.45, labels = names(resOutInd), xpd = NA, cex = 0.9, srt = 45, adj = 1)
  
    if (!is.null(ads) | !is.null(customBg)) {
      abline(lty = 5, h = median(resOutInd[[1]]))
    }
    #
    for (i in 1:length(resOutInd)) {
      boxplot(resOutInd[[i]], add = TRUE, at = i, col = coloursOutInd[i], xaxt = "n", xlab = "", ylab = "", type = "n", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, frame.plot = FALSE, outcol = "grey65", whiskcol = "grey65", outline = FALSE, medcol = "black", staplelty = 0, whisklty = 1)
      text(i, 0, round(mean(antilog(resOutInd[[i]], 2), 0)), font = 2)
    }
    # Plot stats
    if (!is.null(comparisons)) {
      addStats(comparisons, ads, customBg, plotType, resOutInd, coloursOutInd)
    }
    dev.off()
  
    #####tAI
    index_sel <- codind[,which(colnames(codind)=='tai')]
    names(index_sel) <- codind$external_gene_name
  
    resOutInd <- resSel(vIn = index_sel, ads = ads, regulation = regulation, contrast = contrast, customBg = customBg, geneList = geneList)
    coloursOutInd <- coloursSel(ads = ads, regulation = regulation, geneList = geneList, geneListcolours = geneListcolours, customBg = customBg)
  
    ##
    pdf(paste(nameOut,'_tAI_index.pdf', sep = ""),width= 8,height=8, useDingbats = F)
    par(mar = c(12, 12, 5, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
    # 
    if (!is.null(regulation)) {
      xlimIn <- c(0.5, length(regulation) + ifelse(!is.null(geneList), length(geneList), 0) + 1.5)
    } else {
      xlimIn <- c(0.5, length(geneList) + 1.5)
    }
    plot(1, 1, xlim = xlimIn, ylim = c(0, range(as.numeric(unlist(resOutInd)))[2] + (1.25 * length(comparisons))), xaxt = "n", xlab = "", ylab = "", type = "n", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, frame.plot = FALSE)
  
    axis(side = 2, at=seq(0,1,0.2),labels = seq(0,1,0.2), font = 2, las = 2, lwd = 2)
    mtext(side = 2, line = 6, "tAI index", col = "black", font = 2, cex = 1.7, at = 0.5)
    text(1:length(resOut), par("usr")[3] - 0.45, labels = names(resOutInd), xpd = NA, cex = 0.9, srt = 45, adj = 1)
  
    if (!is.null(ads) | !is.null(customBg)) {
      abline(lty = 5, h = median(resOutInd[[1]]))
    }
    #
    for (i in 1:length(resOutInd)) {
      boxplot(resOutInd[[i]], add = TRUE, at = i, col = coloursOutInd[i], xaxt = "n", xlab = "", ylab = "", type = "n", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, frame.plot = FALSE, outcol = "grey65", whiskcol = "grey65", outline = FALSE, medcol = "black", staplelty = 0, whisklty = 1)
      text(i, 0, round(mean(antilog(resOutInd[[i]], 2), 0)), font = 2)
    }
    # Plot stats
    if (!is.null(comparisons)) {
      addStats(comparisons, ads, customBg, plotType, resOutInd, coloursOutInd)
    }
    dev.off()
  
    ###L_aa
    index_sel <- log2(codind[,which(colnames(codind)=='L_aa')])
    names(index_sel) <- codind$external_gene_name
  
    resOutInd <- resSel(vIn = index_sel, ads = ads, regulation = regulation, contrast = contrast, customBg = customBg, geneList = geneList)
    coloursOutInd <- coloursSel(ads = ads, regulation = regulation, geneList = geneList, geneListcolours = geneListcolours, customBg = customBg)
  
    ##
    pdf(paste(nameOut,'L_aa_index.pdf', sep = ""),width= 8,height=8, useDingbats = F)
    par(mar = c(12, 12, 5, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
    # 
    if (!is.null(regulation)) {
      xlimIn <- c(0.5, length(regulation) + ifelse(!is.null(geneList), length(geneList), 0) + 1.5)
    } else {
      xlimIn <- c(0.5, length(geneList) + 1.5)
    }
    plot(1, 1, xlim = xlimIn, ylim = c(0, range(as.numeric(unlist(resOutInd)))[2] + (1.25 * length(comparisons))), xaxt = "n", xlab = "", ylab = "", type = "n", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, frame.plot = FALSE)
  
    axis(side = 2, at=seq(0,10,2),labels = seq(0,10,2), font = 2, las = 2, lwd = 2)
    mtext(side = 2, line = 6, "log2 L_aa index", col = "black", font = 2, cex = 1.7, at = 0.5)
    text(1:length(resOut), par("usr")[3] - 0.45, labels = names(resOutInd), xpd = NA, cex = 0.9, srt = 45, adj = 1)
  
    if (!is.null(ads) | !is.null(customBg)) {
      abline(lty = 5, h = median(resOutInd[[1]]))
    }
    #
    for (i in 1:length(resOutInd)) {
      boxplot(resOutInd[[i]], add = TRUE, at = i, col = coloursOutInd[i], xaxt = "n", xlab = "", ylab = "", type = "n", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, frame.plot = FALSE, outcol = "grey65", whiskcol = "grey65", outline = FALSE, medcol = "black", staplelty = 0, whisklty = 1)
      text(i, 0, round(mean(antilog(resOutInd[[i]], 2), 0)), font = 2)
    }
    # Plot stats
    if (!is.null(comparisons)) {
      addStats(comparisons, ads, customBg, plotType, resOutInd, coloursOutInd)
    }
    dev.off()
  }
  #
  codonsOut <- list()
  for (j in 1:length(comparisons)) {
    if (!is.null(ads) | !is.null(customBg)) {
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
      finalOut_u <- finalOut_u[finalOut_u$codon %in% codU, ]
      text(log2(finalOut_u$statOut), finalOut_u$freq, finalOut_u$codon, col = "firebrick1", font = 2)

      finalOut_d <- finalOut[finalOut$statOut < 1, ]
      codD <- finalOut_d$codon[finalOut_d$freq >= quantile(finalOut_d$freq, 1 - thresY2) & finalOut_d$statOut <= quantile(finalOut_d$statOut, thresX2)]
      finalOut_d <- finalOut_d[finalOut_d$codon %in% codD, ]
      text(log2(finalOut_d$statOut), finalOut_d$freq, finalOut_d$codon, col = "dodgerblue1", font = 2)

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
