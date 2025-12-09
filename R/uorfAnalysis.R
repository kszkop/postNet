uorfAnalysis <- function(ptn,
                         startCodon = "ATG",
                         KozakContext = "strong",
                         onlyUTR5 = FALSE,
                         unitOut = "number",
                         comparisons = NULL,
                         plotOut = TRUE,
                         pdfName = NULL) {
  #
  check_ptn(ptn)
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
  if (!isStartCodon(startCodon)) {
    stop("The input for 'startCodon' must be a character vector of length one, and contain a 3-nucleotide sequence, ex. 'ATG'")
  }
  if (!isKozakContext(KozakContext)) {
    stop("The input for 'KozakContext' must be either: 'strong', 'adequate1', 'adequate2', 'weak', or 'any'.")
  }
  if (!check_logical(onlyUTR5)) {
    stop("The input for 'onlyUTR5' must be logical: TRUE or FALSE.")
  }
  if (!isUnitOut(unitOut)) {
    stop("The input for 'unitOut' must be either 'number' or 'position'.")
  }

  #
  KozakContext <- tolower(KozakContext)
  if (KozakContext == "strong") {
    context <- paste("[AG][ATGC][ATGC]", toupper(startCodon), "G", sep = "")
  } else if (KozakContext == "adequate1") {
    context <- paste("[AG][ATGC][ATGC]", toupper(startCodon), "[ATC]", sep = "")
  } else if (KozakContext == "adequate2") {
    context <- paste("[TC][ATGC][ATGC]", toupper(startCodon), "G", sep = "")
  } else if (KozakContext == "weak") {
    context <- paste("[TC][ATGC][ATGC]", toupper(startCodon), "[ATC]", sep = "")
  } else if (KozakContext == "any") {
    context <- paste("[ATGC][ATGC][ATGC]", toupper(startCodon), "[ATGC]", sep = "")
  } else {
    stop("Please provide a valid Kozak context.")
  }
  #
  uORFFinal <- list()

  seqTmp <- ptn_sequences(ptn, region = "UTR5")
  #
  if (!isTRUE(onlyUTR5)) {
    seq <- list()
    seq[[1]] <- ptn_sequences(ptn, region = "CDS")
    seq[[2]] <- ptn_sequences(ptn, region = "UTR3")
    extSeq <- combSeq(seqIn = seq)
    extSeq <- unlist(extSeq)
    #
    uorfOut <- mapply(calc_uORF, seqTmp = seqTmp, ext = extSeq, context = context, unit = tolower(unitOut), USE.NAMES = FALSE)
  } else {
    uorfOut <- sapply(seqTmp, function(x) calc_uORF(x, ext = NULL, context = context, unit = tolower(unitOut)), USE.NAMES = FALSE)
  }
  #
  names(uorfOut) <- ptn_geneID(ptn, region = "UTR5")
  #
  if (tolower(unitOut) == "number" & isTRUE(plotOut)) {
    #
    resOut <- resQuant(qvec = uorfOut, ptn = ptn)

    if (length(resOut) == 0) {
      stop("There are no regulated genes. Check the input or run without indicating 'regulation' and 'comparisons'.")
    }
    colOut <- colPlot(ptn)
    #
    resProp <- as.numeric()
    for (i in 1:length(resOut)) {
      resProp[i] <- length(resOut[[i]][resOut[[i]] > 0]) / length(resOut[[i]])
    }
    dataTmp <- as.numeric(unlist(resProp))
    ylimTmp2_1 <- 0
    ylimTmp2_2 <- roundNice(max(dataTmp), direction = "up")
    ylimTmp <- as.numeric(adjust_ylim(ylimTmp2_1, ylimTmp2_2))

    # Plot
    pdf(ifelse(is.null(pdfName), paste("uORFs_", KozakContext, ".pdf", sep = ""), paste(pdfName, "_uORFs_", KozakContext, ".pdf", sep = "")), width = 8, height = 8, useDingbats = F)
    m <- layout(mat = matrix(c(1, 2), nrow = 2, ncol = 1), heights = c(1, 5))
    xlimTmp <- c(0.5, length(resProp) + 0.5)
    #
    par(mar = c(0, 8, 3, 0), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
    ylimTmp1 <- ifelse(!is.null(comparisons), length(comparisons), 0)
    plot(1, ylimTmp1, xlim = xlimTmp, ylim = c(0, ylimTmp1), xaxt = "n", type = "n", yaxt = "n", xlab = "", ylab = "", main = "", lwd = 1, bty = "n", font = 2, frame.plot = FALSE)
    if (!is.null(comparisons)) {
      addStats(comparisons, plotType = "boxplot", resOut, colOut)
    }
    par(mar = c(8, 8, 0, 0), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
    plot(1, max(ylimTmp2_1, ylimTmp2_2), xlim = xlimTmp, ylim = ylimTmp, xaxt = "n", type = "n", yaxt = "n", xlab = "", ylab = "", main = "", lwd = 1, bty = "n", font = 2, frame.plot = FALSE)
    xpos <- seq_along(resProp)
    bar_width <- 1
    for (i in xpos) {
      rect(xleft = xpos[i] - bar_width / 2, xright = xpos[i] + bar_width / 2, ybottom = 0, ytop = resProp[i], col = colOut[i], border = NA)
      resPropTmp <- resProp[[i]] * 100
      resPropTmp <- ifelse(resPropTmp > -1 & resPropTmp < 1, round(resPropTmp, 2), round(resPropTmp, 0))
      text(i, 0.05, paste(resPropTmp, "%", sep = " "), font = 2, cex = 1.3)
    }
    axis(side = 2, font = 2, las = 2, lwd = 2) # ), at = seq(0, ylimTmp[2], 0.2), labels = seq(0, 1, 0.2))
    text(x = xpos, y = par("usr")[3] - 0.05 * diff(par("usr")[3:4]), labels = names(resOut), xpd = TRUE, cex = 0.9, srt = 45, adj = 1.1)
    mtext(side = 2, line = 5, paste("proportion of uORFs in \n", KozakContext, "Kozak context", sep = ""), col = "black", font = 2, cex = 1.7, at = mean(ylimTmp))

    dev.off()
  }
  uORFFinal[[paste("uORFs", startCodon, KozakContext, sep = "_")]] <- uorfOut
  #
  return(uORFFinal)
}
