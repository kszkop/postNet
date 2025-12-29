codonUsage <- function(ptn,
                       annotType = "ptnCDS",
                       sourceSeq = "load",
                       analysis,
                       codonN = 1,
                       pAdj = 0.01,
                       rem5 = TRUE,
                       plotHeatmap = TRUE,
                       thresOddsUp = 0.3,
                       thresFreqUp = 0.3,
                       thresOddsDown = 0.3,
                       thresFreqDown = 0.3,
                       subregion = NULL,
                       subregionSel = NULL,
                       comparisons,
                       plotType_index = "boxplot",
                       setSeed = NULL,
                       pdfName = NULL) {
  #
  check_ptn(ptn)
  if (!is_annotType(annotType)) {
    stop("Please provide an input for 'annotType' specifying the source of annotation. The options are 'ccds', or 'ptnCDS'.")
  }
  if (annotType == "ccds") {
    species <- ptn_species(ptn)
    if (!is_valid_species(species)) {
      stop("Please specify a species. Currently, 'human' or 'mouse' are available.")
    }
    selectionTmp <- slot(ptn, "selection")
    check_selection(selectionTmp)
    if (!is_valid_sourceSeq(sourceSeq)) {
      stop("Please provide an input for 'sourceSeq'. The options are: 'load' or 'create'.")
    }
    if (sourceSeq == "create") {
      if (species == "human") {
        #
        download.file("ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS.current.txt", destfile = "CCDS_human.txt")
        ccds <- read.delim(file = "CCDS_human.txt", as.is = c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE))
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
        colnames(annotTmp)[2] <- "geneID"
        lenTmp <- as.numeric(sapply(annotTmp$CDS_seq, function(x) length(seqinr::s2c(x))))
        annotTmp$lenTmp <- lenTmp
        #
        annotSel <- isoSel(annot = annotTmp, method = ptn_selection(ptn), setSeed = setSeed )
        colnames(annotSel)[1:3] <- c("id", "geneID", "CDS_seq")
        #
        annot <- new("postNetRegion",
          id = annotSel$id,
          geneID = annotSel$geneID,
          sequences = annotSel$CDS_seq
        )

        ptn@annot@CCDS <- annot
      } else if (species == "mouse") {
        #
        download.file("ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_mouse/CCDS.current.txt", destfile = "CCDS_mouse.txt")
        ccds <- read.delim(file = "CCDS_mouse.txt", as.is = c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE))
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
        colnames(annotTmp)[2] <- "geneID"
        lenTmp <- as.numeric(sapply(annotTmp$CDS_seq, function(x) length(seqinr::s2c(x))))
        annotTmp$lenTmp <- lenTmp
        #
        annotSel <- isoSel(annot = annotTmp, method = ptn_selection(ptn), setSeed = setSeed)
        colnames(annotSel)[1:3] <- c("id", "geneID", "CDS_seq")
        #
        annot <- new("postNetRegion",
          id = annotSel$id,
          geneID = annotSel$geneID,
          sequences = annotSel$CDS_seq
        )

        ptn@annot@CCDS <- annot
      }
    } else if (sourceSeq == "load") {
      #
      currTmp <- list.files(system.file("extdata/annotation/ccds", package = "postNet"))
      #
      if (!species %in% currTmp) {
        stop("This option is currently only available for human and mouse. For analyses on other species or \
             annotations, please use the source = 'custom' option in the postNetStart() function to provide custom reference \
             sequences.")
      }
      if (species == "human") {
        annotTmp <- read.delim(system.file(paste("extdata/annotation/ccds/human", sep = "/"), "humanDB_ccds.txt.gz", package = "postNet"), stringsAsFactors = FALSE)
      }
      if (species == "mouse") {
        annotTmp <- read.delim(system.file(paste("extdata/annotation/ccds/mouse", sep = "/"), "mouseDB_ccds.txt.gz", package = "postNet"), stringsAsFactors = FALSE) # }
      }
      lenTmp <- as.numeric(sapply(annotTmp$CDS_seq, function(x) length(seqinr::s2c(x))))
      annotTmp$lenTmp <- lenTmp
      #
      annotSel <- isoSel(annot = annotTmp, method = ptn_selection(ptn), setSeed = setSeed)
      colnames(annotSel)[1:3] <- c("id", "geneID", "CDS_seq")
      #
      annot <- new("postNetRegion",
        id = annotSel$id,
        geneID = annotSel$geneID,
        sequences = annotSel$CDS_seq
      )

      ptn@annot@CCDS <- annot
    }
  }
  ####
  if (!is_valid_analysis(analysis)) {
    stop("Please provide an input for 'analysis'. The options are 'codon' or 'AA' (amino acid).")
  }
  if (!check_number(codonN)) {
    stop("Please provide an integer value for 'codonN' specifying the number of consecutive codons to consider in the analysis. \
    For example, 1 will consider each individual codon, while 2 will consider dicodons, etc.")
  }
  if (!check_logical(rem5)) {
    stop("The input for 'rem5' must be logical: TRUE or FALSE.")
  }
  #
  if (!is.null(subregion) && (!is.numeric(subregion) || !length(subregion) == 1)) {
    stop("The input for 'subregion' must be an integer.")
  }
  if (!is.null(subregionSel) && !subregionSel %in% c("select", "exclude")) {
    stop("The input for 'subregionSel' must be either 'select' or 'exclude'.")
  }

  if (!is.null(comparisons)) {
    if (!check_comparisons(comparisons)) {
      stop("The input for 'comparisons' must be a list of numeric vectors of paired comparisons. For example: list(c(0,2),c(0,1)). 0 always \ denotes the background gene set.")
    }
    #
    if (length(which(unique(unlist(comparisons)) == 0)) > 0 && is.null(ptn_background(ptn))) {
      stop("0 always denotes the background, but no background has been provided.")
    }
  } else {
    stop(
      "For further analysis to be performed, pairs of comparisons must be specified with the 'comparisons' parameter. \
      See the 'codonUsage()' help page for details: ?codonUsage.",
      call. = FALSE
    )
  }
  #
  nameTmp <- ifelse(is.null(pdfName), analysis, paste(pdfName, analysis, sep = "_"))
  nameOut <- nameTmp
  #
  if (annotType == "ccds") {
    seqTmp <- ptn_sequences(ptn, region = "CCDS")
    names(seqTmp) <- ptn_geneID(ptn, region = "CCDS")
  } else {
    seqTmp <- ptn_sequences(ptn, region = "CDS")
    names(seqTmp) <- ptn_geneID(ptn, region = "CDS")
  }
  # remove last codon (stop codon)
  seqTmp <- remove_last3(seqTmp)
  #
  if (!is.null(subregion)) {
    #
    subSeq <- sapply(seqTmp, function(x) subset_seq(x, pos = subregion, subregionSel = subregionSel))
    if (length(which(is.na(subSeq))) > 0) {
      message("For some sequences, the selected subregion is longer than the reference sequence region. These sequences will be removed from the analysis.")
    }
    seqTmp <- subSeq
  }
  seqTmp <- seqTmp[!is.na(seqTmp)]

  geneIDs <- names(seqTmp)

  codonTmp <- list()
  for (i in 1:length(geneIDs)) {
    codonTmp[[i]] <- codonCount(gene = geneIDs[i], seq = seqTmp[i], codonN = codonN)
  }
  codonAll <- data.table::rbindlist(codonTmp, use.names = TRUE, fill = TRUE)


  if (codonN == 1) {
    codonAll <- codonAll[!codonAll$AA %in% c("Stop", "W", "M", "O", "U"), ]
  }
  codonsAllOut <- new("postNetCodonsAll",
    geneID = codonAll$geneID,
    codon = codonAll$codon,
    AA = codonAll$AA,
    count = codonAll$count,
    frequency = codonAll$frequency,
    AACountPerGene = codonAll$AACountPerGene,
    relative_frequency = codonAll$relative_frequency
  )

  # indexes
  if (analysis == "codon" & codonN == 1) {
    species <- ptn_species(ptn)
    if (!is_valid_species(species)) {
      stop("Currently, codon usage indexes are only available for 'human' or 'mouse'.")
    }
    if (species == "human") {
      codind <- read.delim(system.file("extdata/indexes/human/", "IndexesHuman.txt", package = "postNetParcel"))
    } else if (species == "mouse") {
      codind <- read.delim(system.file("extdata/indexes/mouse/", "IndexesMouse.txt", package = "postNetParcel"))
    } else {
      message("No available codon indexes for ", species)
    }

    # indexes
    indNames <- c("CAI", "CBI", "Fop", "tAI", "L_aa")

    for (ind in indNames) {
      #####
      index_sel <- codind[, which(colnames(codind) == ind)]
      names(index_sel) <- codind$external_gene_name

      resOutInd <- resQuant(qvec = index_sel, ptn = ptn)
      colOutInd <- colPlot(ptn)

      ##
      pdf(paste(nameOut, ind, "_index.pdf", sep = ""), width = 8, height = 8, useDingbats = FALSE)
      ylabel <- paste(ind, "index", sep = " ")
      plotPostNet(resOutInd, colOutInd, comparisons, ylabel = ylabel, plotType = plotType_index)
      dev.off()
    }
  }
  #
  compTmpAll <- sort(unique(unlist(comparisons)), decreasing = F)
  compOut1 <- list()
  compOut2 <- list()
  compOut3 <- list()
  compOut4 <- list()
  for (i in compTmpAll) {
    if (i == 0) {
      selTmp <- ptn_background(ptn)
    } else {
      resTmp <- ptn_geneList(ptn)
      selTmp <- resTmp[[i]]
    }
    tmp <- codonAll[codonAll$geneID %in% selTmp, ]
    if (analysis == "codon") {
      #####
      tmp1 <- tmp %>%
        group_by(codon) %>%
        summarise(codonPerReg = sum(count))
      tmpSum <- tmp1$codonPerReg
      names(tmpSum) <- tmp1$codon

      if (i == 0) {
        compOut1[["background"]] <- tmpSum
      } else {
        compOut1[[names(resTmp)[i]]] <- tmpSum
      }

      ####
      tmp2 <- tmp %>%
        group_by(codon) %>%
        summarise(freqPerReg = mean(frequency))
      tmpFreq <- tmp2$freqPerReg
      names(tmpFreq) <- tmp2$codon

      if (i == 0) {
        compOut2[["background"]] <- tmpFreq
      } else {
        compOut2[[names(resTmp)[i]]] <- tmpFreq
      }

      #####
      tmp3 <- tmp %>%
        group_by(codon) %>%
        mutate(codonPerReg = sum(count))
      tmp3 <- tmp3 %>%
        group_by(AA) %>%
        mutate(AAPerReg = sum(AACountPerGene))
      tmp3 <- subset(tmp3, !duplicated(codon))
      tmp3$codonNormAA <- tmp3$codonPerReg / tmp3$AAPerReg
      tmpcodonNormAA <- tmp3$codonNormAA
      names(tmpcodonNormAA) <- tmp3$codon

      if (i == 0) {
        compOut3[["background"]] <- tmpcodonNormAA
      } else {
        compOut3[[names(resTmp)[i]]] <- tmpcodonNormAA
      }
      #
      tmp4 <- tmp %>%
        group_by(codon) %>%
        summarise(freq = sum(frequency))
      tmpSumFreq <- tmp4$freq
      names(tmpSumFreq) <- tmp4$codon
      if (i == 0) {
        compOut4[["background"]] <- tmpSumFreq
      } else {
        compOut4[[names(resTmp)[i]]] <- tmpSumFreq
      }
    } else if (analysis == "AA") {
      tmp1 <- tmp %>%
        group_by(AA) %>%
        summarise(AAPerReg = sum(count))
      #
      tmpSum <- tmp1$AAPerReg
      names(tmpSum) <- tmp1$AA

      if (i == 0) {
        compOut1[["background"]] <- tmpSum
      } else {
        compOut1[[names(resTmp)[i]]] <- tmpSum
      }

      tmp4 <- tmp %>%
        group_by(AA) %>%
        summarise(freq = sum(frequency))
      tmpSumFreq <- tmp4$freq
      names(tmpSumFreq) <- tmp4$AA
      if (i == 0) {
        compOut4[["background"]] <- tmpSumFreq
      } else {
        compOut4[[names(resTmp)[i]]] <- tmpSumFreq
      }
    }
  }
  #
  codonsSel <- list()
  for (j in 1:length(comparisons)) {
    if (names(compOut1)[1] == "background") {
      compTmp <- comparisons[[j]] + 1
    } else {
      compTmp <- comparisons[[j]]
    }
    resIn <- compOut1[compTmp]
    resIn <- do.call(rbind, resIn)

    # Remove columns with all zeros
    resIn <- resIn[, colSums(resIn) > 0]

    # Check if resIn has at least one positive entry
    if (all(resIn <= 0)) {
      stop("The contingency table (codon or AA counts by geneset) must have at least one positive entry \
           for the Chi-squared test to be performed.")
    }

    #
    if (sum(apply(resIn, 2, min) < 5) > 0) {
      if (isTRUE(rem5)) {
        rem5Ind <- as.numeric(which(apply(resIn, 2, min) < 5))
        resIn <- resIn[, -rem5Ind]
      } else {
        warning("In the contingency table (codon or AA counts by geneset), some counts are lower than 5 which may invalidate the Chi-squared test. \ You may wish to perform the analysis on a subset or groups of codons or AAs instead.",
          call. = FALSE
        )
      }
    }

    # Check if resIn has at least one positive entry
    if (all(resIn <= 0)) {
      stop("The contingency table (codon or AA counts by geneset) must have at least one positive entry \
           for the Chi-squared test to be performed.")
    }

    chisqTest <- chisq.test(resIn)

    if (is.na(chisqTest$p.value)) {
      stop("The Chi-squared test could not be performed. Please check that all codons or AAs have at least one count in at least one gene set.")
    }

    StResiduals <- chisqTest$stdres
    signifLimit <- sqrt(qchisq(pAdj / (dim(resIn)[1] * dim(resIn)[2]), lower.tail = FALSE, df = 1))

    if (isTRUE(plotHeatmap) & codonN == 1) {
      #
      residOut <- t(as.matrix(StResiduals))

      pdf(paste(nameOut, paste(colnames(residOut), collapse = "_"), "heatmap.pdf", sep = "_"), width = 20, height = 28, useDingbats = FALSE)
      par(mar = c(10, 5, 5, 5), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 0.9, cex.main = 0.7, cex.lab = 0.7)
      col <- grDevices::colorRampPalette(c("blue", "white", "red"))(50)

      max_abs <- ceiling(max(abs(residOut[, 1])))
      breaks <- sort(c(seq(-max_abs, max_abs, length = 51)))

      gplots::heatmap.2(residOut,
        col = col,
        breaks = breaks,
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
    if (length(signSel) == 0) {
      message(paste("There are no significant codons for the comparison ", j, sep = ""))
    } else {
      #
      if (analysis == "codon" & codonN == 1) {
        resIn2 <- compOut2[compTmp]
        resIn2 <- as.data.frame(do.call(cbind, resIn2))
        resIn2$codon <- row.names(resIn2)
        resIn2$AA <- seqinr::aaa(seqinr::translate(seqinr::s2c(seqinr::c2s(rownames(resIn2)))))
        resIn2$AAcodon <- paste(seqinr::translate(seqinr::s2c(seqinr::c2s(rownames(resIn2)))), rownames(resIn2), sep = "_")

        #
        xylim1 <- roundNice(max(resIn2[, c(1, 2)]), direction = "up")
        #
        pdf(paste(nameOut, paste(colnames(resIn2)[1:2], collapse = "_"), "averageFreq.pdf", sep = "_"), width = 8, height = 8, useDingbats = FALSE)
        par(mar = c(10, 5, 5, 10), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 0.9, cex.main = 0.7, cex.lab = 0.7)
        pOut1 <- ggplot2::ggplot(resIn2, ggplot2::aes(!!sym(colnames(resIn2)[1]), !!sym(colnames(resIn2)[2]), col = AA)) +
          ggplot2::theme_bw() +
          ggplot2::xlim(0, xylim1) +
          ggplot2::ylim(0, xylim1) +
          ggplot2::geom_abline(slope = 1, intercept = 0, col = "grey", linetype = "dashed", linewidth = 0.5) +
          ggplot2::geom_point() +
          ggplot2::coord_fixed() +
          ggplot2::geom_line(ggplot2::aes(group = AA), col = "gray", linetype = 1, linewidth = 0.2) +
          ggplot2::theme(legend.position = "none") +
          ggrepel::geom_text_repel(ggplot2::aes(label = AAcodon), size = 3, segment.size = 0.2, max.overlaps = 100) +
          ggplot2::labs(
            x = paste(colnames(resIn2)[1], "\n(Average codon frequency)", sep = ""),
            y = paste(colnames(resIn2)[2], "\n(Average codon frequency)", sep = "")
          )
        print(pOut1)
        dev.off()

        #
        resIn3 <- compOut3[compTmp]

        resIn3 <- as.data.frame(do.call(cbind, resIn3))
        resIn3$codon <- row.names(resIn3)
        resIn3$AA <- seqinr::aaa(seqinr::translate(seqinr::s2c(seqinr::c2s(rownames(resIn3)))))
        resIn3$AAcodon <- paste(seqinr::translate(seqinr::s2c(seqinr::c2s(rownames(resIn3)))), rownames(resIn3), sep = "_")

        #
        xylim2 <- roundNice(max(resIn3[, c(1, 2)]), direction = "up")
        #
        pdf(paste(nameOut, paste(colnames(resIn3)[1:2], collapse = "_"), "codonUsage.pdf", sep = "_"), width = 8, height = 8, useDingbats = FALSE)
        par(mar = c(10, 5, 5, 10), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 0.9, cex.main = 0.7, cex.lab = 0.7)
        pOut2 <- ggplot2::ggplot(resIn3, ggplot2::aes(!!sym(colnames(resIn3)[1]), !!sym(colnames(resIn3)[2]), col = AA)) +
          ggplot2::theme_bw() +
          ggplot2::xlim(0, xylim2) +
          ggplot2::ylim(0, xylim2) +
          ggplot2::geom_point() +
          ggplot2::geom_abline(slope = 1, intercept = 0, col = "grey", linetype = "dashed", linewidth = 0.5) +
          ggplot2::coord_fixed() +
          ggplot2::geom_line(ggplot2::aes(group = AA), col = "gray", linetype = 1, linewidth = 0.2) +
          ggplot2::theme(legend.position = "none") +
          ggrepel::geom_text_repel(ggplot2::aes(label = AAcodon), size = 3, segment.size = 0.2, max.overlaps = 100) +
          ggplot2::labs(
            x = paste(colnames(resIn3)[1], "\n(Amino acid normalised codon usage)", sep = ""),
            y = paste(colnames(resIn3)[2], "\n(Amino acid normalised codon usage)", sep = "")
          )
        print(pOut2)
        dev.off()
        #


        resIn4 <- data.frame(codon = colnames(resIn), t(resIn), row.names = NULL)
        resIn4$AA <- seqinr::aaa(seqinr::translate(seqinr::s2c(seqinr::c2s(resIn4$codon))))

        regComb <- combn(x = names(compOut1[compTmp]), m = 2)
        statOut <- log2(statOnDf(df = resIn4, regs = regComb, analysis = analysis)[signSel])
        #
        sumFreq <- compOut4[compTmp]
        sumFreq <- rowSums(do.call(cbind, sumFreq))[signSel]

        #
        xlimT <- roundNice(max(abs(range(statOut))), direction = "up")
        pdf(paste(nameOut, paste(regComb[, 1], collapse = "_"), "codon_oddratio_vs_freq.pdf", sep = ""), width = 8, height = 8, useDingbats = FALSE)
        m <- layout(mat = matrix(c(1, 2), nrow = 2, ncol = 1), heights = c(1, 5))
        #
        par(mar = c(0, 8, 0, 3), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
        plot(0, 1, xlim = c(-xlimT, xlimT), ylim = c(0, 1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "", lwd = 1, bty = "n", font = 2, type = "n", frame.plot = FALSE)
        legend(-xlimT, 0.5, paste(regComb[1], "vs", regComb[2], sep = " "), bty = "n", horiz = TRUE, xpd = TRUE, y.intersp = 2.8, cex = 1.3)

        par(mar = c(5, 8, 0, 3), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.3)
        plot(statOut, sumFreq, col = "black", pch = 20, cex = 0.1, xlab = "", ylab = "", lwd = 1, bty = "n", font = 2, xlim = c(-xlimT, xlimT), ylim = range(sumFreq)) # , xaxt = "n", yaxt = "n")
        text(statOut, sumFreq, names(sumFreq), col = "black", font = 2)

        mtext(side = 2, line = 3, "frequency ", col = "black", font = 2, cex = 1.7)
        mtext(side = 1, line = 3, "log2(Odds Ratio)", col = "black", font = 2, cex = 1.7, at = 0)

        # axis(side = 2, seq(0, roundNice(range(sumFreq)[2], direction='up'), ifelse(roundNice(range(sumFreq)[2],direction='up') > 10, 10, 1)), font = 2, las = 2, lwd = 2)
        # axis(side = 1, seq(-xlimT, xlimT, 1), font = 2, lwd = 2)

        indUp <- as.numeric(which(statOut >= 0))
        statOut_up <- statOut[indUp]
        sumFreq_up <- sumFreq[indUp]
        codU <- names(sumFreq_up)[sumFreq_up >= quantile(sumFreq_up, 1 - thresFreqUp) & statOut_up >= quantile(statOut_up, 1 - thresOddsUp)]
        if (length(codU) > 0) {
          text(statOut_up[codU], sumFreq_up[codU], codU, col = "firebrick1", font = 2)
        } else {
          message("None of the codons are labelled with the specified thresholds. \
                  Please select more relaxed thresholds for thresOddsUp and thresFreqUp.")
        }
        #
        indDown <- as.numeric(which(statOut < 0))
        statOut_down <- statOut[indDown]
        sumFreq_down <- sumFreq[indDown]
        codD <- names(sumFreq_down)[sumFreq_down >= quantile(sumFreq_down, 1 - thresFreqDown) & statOut_down <= quantile(statOut_down, thresOddsDown)]
        if (length(codD) > 0) {
          text(statOut_down[codD], sumFreq[codD], codD, col = "dodgerblue1", font = 2)
        } else {
          message("None of the codons are labelled with the specified thresholds. \
                  Please select more relaxed thresholds for thresOddsDown and thresFreqDown.")
        }
        dev.off()

        codonsSel[[paste(regComb, collapse = "_vs_")]] <- list(Up = codU, Down = codD)
      } else if (analysis == "AA") {
        resIn4 <- data.frame(AA = colnames(resIn), t(resIn), row.names = NULL)
        #
        regComb <- combn(x = names(compOut1[compTmp]), m = 2)
        statOut <- log2(statOnDf(df = resIn4, regs = regComb, analysis = analysis)[signSel])

        sumFreq <- compOut4[compTmp]
        sumFreq <- rowSums(do.call(cbind, sumFreq))[signSel]

        xlimT <- roundNice(max(abs(range(statOut))), direction = "up")

        pdf(paste(nameOut, paste(regComb[, 1], collapse = "_"), "AA_oddratio_vs_freq.pdf", sep = ""), width = 8, height = 8, useDingbats = FALSE)
        m <- layout(mat = matrix(c(1, 2), nrow = 2, ncol = 1), heights = c(1, 5))
        #
        par(mar = c(0, 8, 0, 3), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
        plot(0, 1, xlim = c(-xlimT, xlimT), ylim = c(0, 1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "", lwd = 1, bty = "n", font = 2, type = "n", frame.plot = FALSE)
        legend(-xlimT, 0.5, paste(regComb[1], "vs", regComb[2], sep = " "), bty = "n", horiz = TRUE, xpd = TRUE, y.intersp = 2.8, cex = 1.3)

        par(mar = c(5, 8, 0, 3), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.3)
        plot(statOut, sumFreq, col = "black", pch = 20, cex = 0.1, xlab = "", ylab = "", lwd = 1, bty = "n", font = 2, xlim = c(-xlimT, xlimT), ylim = range(sumFreq))
        text(statOut, sumFreq, names(sumFreq), col = "black", font = 2)

        mtext(side = 2, line = 3, "frequency ", col = "black", font = 2, cex = 1.7)
        mtext(side = 1, line = 3, "log2(Odds Ratio)", col = "black", font = 2, cex = 1.7, at = 0)

        # axis(side = 2, seq(0, roundNice(range(sumFreq)[2], direction='up'), ifelse(roundNice(range(sumFreq)[2],direction='up') > 10, 10, 1)), font = 2, las = 2, lwd = 2)
        # axis(side = 1, seq(-xlimT, xlimT, 1), font = 2, lwd = 2)

        indUp <- as.numeric(which(statOut >= 0))
        statOut_up <- statOut[indUp]
        sumFreq_up <- sumFreq[indUp]
        codU <- names(sumFreq_up)[sumFreq_up >= quantile(sumFreq_up, 1 - thresFreqUp) & statOut_up >= quantile(statOut_up, 1 - thresOddsUp)]
        if (length(codU) > 0) {
          text(statOut_up[codU], sumFreq_up[codU], codU, col = "firebrick1", font = 2)
        } else {
          message("None of the AAs are labelled with the specified thresholds. \
                  Please select more relaxed thresholds for thresOddsUp and thresFreqUp.")
        }
        #
        indDown <- as.numeric(which(statOut < 0))
        statOut_down <- statOut[indDown]
        sumFreq_down <- sumFreq[indDown]
        codD <- names(sumFreq_down)[sumFreq_down >= quantile(sumFreq_down, 1 - thresFreqDown) & statOut_down <= quantile(statOut_down, thresOddsDown)]
        if (length(codD) > 0) {
          text(statOut_down[codD], sumFreq[codD], codD, col = "dodgerblue1", font = 2)
        } else {
          message("None of the AAs are labelled with the specified thresholds. \
                  Please select more relaxed thresholds for thresOddsDown and thresFreqDown.")
        }
        dev.off()

        codonsSel[[paste(regComb, collapse = "_vs_")]] <- list(Up = codU, Down = codD)
      } else {
        codonsSel[[paste(row.names(resIn), collapse = "_vs_")]] <- signSel
      }
    }
  }
  codonsOut <- new("postNetCodons",
    codonAnalysis = codonsAllOut,
    codonSelection = codonsSel
  )

  ptn@analysis@codons <- codonsOut
  return(ptn)
}
